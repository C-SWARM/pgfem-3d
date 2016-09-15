#include "vol_damage_int_alg.h"
#include <string.h>
#include "mkl_cblas.h"

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

#ifndef GET_DOF_IDS_ON_ELEM_H
#include "get_dof_ids_on_elem.h"
#endif

#ifndef GET_NDOF_ON_ELEM_H
#include "get_ndof_on_elem.h"
#endif

#ifndef CAST_MACROS_H
#include "cast_macros.h"
#endif

#ifndef TENSORS_H
#include "tensors.h"
#endif

#ifndef DEF_GRAD_H
#include "def_grad.h"
#endif

#ifndef NEW_POTENTIALS_H
#include "new_potentials.h"
#endif

#ifndef STABILIZED_H
#include "stabilized.h"
#endif

#ifndef DISP_BASED_ELEM_H
#include "displacement_based_element.h"
#endif

#ifndef VD_INT_ALG_DEBUG
#define VD_INT_ALG_DEBUG 0
#endif

#define MAX(a,b) (((a)<(b))? (b):(a))

static int integration_help(const int elem_id,
			    const int ip,
			    const int nne,
			    const int i,
			    const int j,
			    const int k,
			    const double *x,
			    const double *y,
			    const double *z,
			    const double *int_pt_ksi,
			    const double *int_pt_eta,
			    const double *int_pt_zet,
			    const double *disp,
			    const EPS *eps,
			    const SUPP sup,
			    const int stab_int,
			    double *Na,
			    double *N_x,
			    double *N_y,
			    double *N_z,
			    double *ST,
			    double *F,
			    double *C,
			    double *J,
			    double *wt,
			    const int analysis);

static int get_material_potential(double *Ybar,
				  const int elem_id,
				  const int ip,
				  const int nne,
				  const double kappa,
				  const double *Np,
				  const double *C,
				  const double J,
				  const double *dp,
				  const EPS *eps,
				  const SIG *sig,
				  const HOMMAT *mat,
				  const int stab_int,
				  const int analysis);

int vol_damage_int_alg(const int ne,
		       const int ndofn,
		       const double *d_r,
		       const double *r,
		       const ELEMENT *elem,
		       const NODE *node,
		       const HOMMAT *hommat,
		       const SUPP sup,
		       const double dt,
		       const int iter,
		       const MPI_Comm mpi_comm,
		       EPS *eps,
		       SIG *sig,
		       double *max_omega,
		       double *dissipation,
		       const int analysis,
		       const int mp_id)
{
  const int ndn = 3;
  int err = 0;
  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  double *F, *C;
  F = (double*) PGFEM_calloc(ndn*ndn,sizeof(double));
  C = (double*) PGFEM_calloc(ndn*ndn,sizeof(double));


  /* reset max_omega */
  *max_omega = 0.0;

  /* reset dissipation */
  *dissipation = 0.0;

  for(int elem_id=0; elem_id<ne; elem_id++){ /* for each element */

    /* if mu <= 0, element not allowed to damage, do not perform
       integration. Since the element is governed by one set of
       material properties, only need to check the parameter for the
       first integration point  */
    if(eps[elem_id].dam[0].params.mu <= 0.0){
      continue;
    }

    const int nne = elem[elem_id].toe;
    const int mat = elem[elem_id].mat[2];
    const double kappa = hommat[mat].E/(3*(1-2*hommat[mat].nu));
    const HOMMAT *ptrMat = &hommat[mat];

    long *nod = (long*) PGFEM_calloc(nne,sizeof(long));

    /* get the nodes on the element */
    elemnodes(elem_id,nne,nod,elem);

    /* get ndofs on element */
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node);

    /* allocate space */
    long *cn = (long*) PGFEM_calloc(ndofe,sizeof(double));
    double *disp = (double*) PGFEM_calloc(nne*ndn,sizeof(double));
    double *dp = (double*) PGFEM_calloc(nne,sizeof(double));
    double *r_e = (double*) PGFEM_calloc(ndofe,sizeof(double));
    double *r_en = (double*) PGFEM_calloc(ndofe,sizeof(double));
    double *x = (double*) PGFEM_calloc(nne,sizeof(double));
    double *y = (double*) PGFEM_calloc(nne,sizeof(double));
    double *z = (double*) PGFEM_calloc(nne,sizeof(double));

    /* get dof id numbers */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* get nodal coordinates */
    switch(analysis){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }

    /* get deformation on element */
    def_elem(cn,ndofe,d_r,elem,node,r_e,sup,0); /* increment of deformation */
    if(analysis == DISP){
      /*=== NOTE:
	disp-based element is TOTAL Lagrangian, thus we need
	the TOTAL deformation, not just the increment, in order to
	compute deformation gradients etc.
	=== */

      def_elem(cn,ndofe,r,elem,node,r_en,sup,1); /* deformation at tn */
      cblas_daxpy(ndofe,1.0,r_en,1,r_e,1);
    }

    /* filter displacement for mixed formulation elements */
    {
      int k=0;
      for(int n=0; n<nne; n++){
	for(int d=0; d<node[nod[n]].ndofn; d++){
	  if(d < ndn){
	    disp[n*ndn+d] = r_e[k+d];
	  } else if(d == ndn){
	    dp[n] = r_e[k+d];
	  }
	}
	k += node[nod[n]].ndofn;
      }
    }

    /*===================================================*/
    /*                  INTEGRATION                      */
    /*===================================================*/

    long npt_x, npt_y, npt_z;
    int_point(10,&npt_z);

    double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
    int_pt_ksi = aloc1(npt_z);
    int_pt_eta = aloc1(npt_z);
    int_pt_zet = aloc1(npt_z);
    weights = aloc1(npt_z);

    /* allocate space for the shape functions, derivatives etc */
    double *Na, *N_x, *N_y, *N_z, *ST;
    Na = aloc1 (nne);
    N_x = aloc1 (nne);
    N_y = aloc1 (nne);
    N_z = aloc1 (nne);
    ST = aloc1(3*3*ndn*nne);

    integrate(nne,&npt_x,&npt_y,&npt_z,
	      int_pt_ksi,int_pt_eta,int_pt_zet,
	      weights);

    double elem_diss = 0.0;
    int ip = 0;
    for(int ii=0; ii<npt_x; ii++){
      for(int j=0; j<npt_y; j++){
	for(int k=0; k<npt_z; k++){
	  double J;
	  double wt = weights[ip];

	  err = integration_help(elem_id,ip,nne,ii,j,k,x,y,z,
				 int_pt_ksi,int_pt_eta,int_pt_zet,
				 disp,eps,sup,0,Na,N_x,N_y,N_z,ST,F,C,&J,
				 &wt,analysis);

	  damage *ptrDam = &(eps[elem_id].dam[ip]);

	  double Ybar = 0.0;
	  double g = 0.0;
	  err += get_material_potential(&Ybar,elem_id,ip,nne,kappa,Na,C,J,
					dp,eps,sig,ptrMat,0,analysis);

	  g = damage_int_alg(ptrDam,Ybar,dt /* ,iter */);

	  if(VD_INT_ALG_DEBUG){
	    PGFEM_printerr("[%d] (elem,ip)::(%d,%d) wn+1 || Xn+1 || g\n"
		    "%1.12e || %1.12e || %1.12e\n",myrank,elem_id,ip,
		    ptrDam->w,ptrDam->X,g);
	  }

	  /* store max_omega */
	  *max_omega = MAX(ptrDam->w - ptrDam->wn,*max_omega);

	  /* integrate the element dissipation */
	  elem_diss += (ptrDam->w - ptrDam->wn) * Ybar *wt;

	  ip++;

	  /* inverted element detected, exit and return error */
	  if(err != 0) goto exit_function;
	}
      }
    }

    /* accumulate dissipation weighted by element volume */
    *dissipation += elem_diss;

    /*==========================================================*/
    /*             QUADRADIC INTEGRATION                        */
    /*==========================================================*/
    if(analysis == STABILIZED){ /* also require quadradic integration for
    			  pressure terms on stabilized element */
      integrate(10,&npt_x,&npt_y,&npt_z,
    		int_pt_ksi,int_pt_eta,int_pt_zet,
    		weights);
      ip = 0;
      for(int ii=0; ii<npt_x; ii++){
    	for(int j=0; j<npt_y; j++){
    	  for(int k=0; k<npt_z; k++){
    	    double J;
	    double wt = weights[ip];

    	    err = integration_help(elem_id,ip,nne,ii,j,k,x,y,z,
    				   int_pt_ksi,int_pt_eta,int_pt_zet,
    				   disp,eps,sup,1,Na,N_x,N_y,N_z,ST,F,C,&J,
				   &wt,analysis);

	    damage *ptrDam = &(eps[elem_id].st[ip].dam);

    	    double Ybar = 0.0;
    	    err += get_material_potential(&Ybar,elem_id,ip,nne,kappa,Na,C,J,
    					  dp,eps,sig,ptrMat,1,analysis);

    	    damage_int_alg(ptrDam,Ybar,dt /* ,iter */);

	    /* store max_omega */
	    *max_omega = MAX(ptrDam->w - ptrDam->wn,*max_omega);

    	    ip++;

    	    /* inverted element detected, exit and return error */
    	    if(err != 0) goto exit_function;
    	  }
    	}
      }
    } /* end analysis == STABILIZED */

    /*===================================================*/
    /*               END INTEGRATION                     */
    /*===================================================*/
  exit_function:
    /* free memory */
    free(nod);
    free(cn);
    free(r_e);
    free(r_en);
    free(disp);
    free(dp);
    free(x);
    free(y);
    free(z);

    free(int_pt_ksi);
    free(int_pt_eta);
    free(int_pt_zet);
    free(weights);
    free(Na);
    free(N_x);
    free(N_y);
    free(N_z);
    free(ST);

    if(err != 0) break;
  } /* each element */

  free(F);
  free(C);

  return err;
}

static int integration_help(const int elem_id,
			    const int ip,
			    const int nne,
			    const int i,
			    const int j,
			    const int k,
			    const double *x,
			    const double *y,
			    const double *z,
			    const double *int_pt_ksi,
			    const double *int_pt_eta,
			    const double *int_pt_zet,
			    const double *disp,
			    const EPS *eps,
			    const SUPP sup,
			    const int stab_int,
			    double *Na,
			    double *N_x,
			    double *N_y,
			    double *N_z,
			    double *ST,
			    double *F,
			    double *C,
			    double *J,
			    double *wt,
			    const int analysis)
{
  int err = 0;
  double ksi,eta,zet;
  const int ndn = 3;
  if(nne == 8){/* hexahedron */
    ksi = int_pt_ksi[i];
    eta = int_pt_eta[j];
    zet = int_pt_zet[k];
  } else { /* tetrahedron type */
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
  }

  double ****ST_tensor, **Fr_mat;
  ST_tensor = aloc4 (3,3,ndn,nne);
  Fr_mat = aloc2(3,3);

  shape_func(ksi,eta,zet,nne,Na);
  (*wt) *= deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
  shape_tensor (nne,ndn,N_x,N_y,N_z,ST_tensor);
  shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);

  def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
  mat2array(F,CONST_2(double) Fr_mat,3,3);

  /* add macro contribution to Fr.*/
  if(sup->multi_scale){
    cblas_daxpy(9,1.0,sup->F0,1,F,1);
  }

  *J = getJacobian(F,elem_id,&err);

  if(analysis != DISP){/* updated lagrangian formulation */
    /* get Fn */
    double *Fn;
    if(stab_int){ /* Fn from stab term */
      Fn = eps[elem_id].st[ip].Fpp;
    } else {
      Fn = eps[elem_id].il[ip].F;
    }
    int t_err = 0;

    /* J = Jr*Jn */
    (*J) *= getJacobian(Fn,elem_id,&t_err);
    err += t_err;

    /* F = FrFn */
    double *temp = aloc1(9);
    memcpy(temp,F,9*sizeof(double));
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,temp,3,Fn,3,0.0,F,3);
    free(temp);
  }

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
	      3,3,3,1.0,F,3,F,3,0.0,C,3);


  dealoc4(ST_tensor,3,3,ndn);
  dealoc2(Fr_mat,3);

  return err;
}/* integration_help */

static int get_material_potential(double *Ybar,
				  const int elem_id,
				  const int ip,
				  const int nne,
				  const double kappa,
				  const double *Np,
				  const double *C,
				  const double J,
				  const double *dp,
				  const EPS *eps,
				  const SIG *sig,
				  const HOMMAT *mat,
				  const int stab_int,
				  const int analysis)
{
  int err = 0;

  switch(analysis){
  case STABILIZED:
    {
      double Un_1, Jn_1;
      if(stab_int){/* get from eps.st */
	Un_1 = eps[elem_id].st[ip].Un_1;
	Jn_1 = eps[elem_id].st[ip].Jn_1;
      } else {
	Un_1 = eps[elem_id].il[ip].Un_1;
	Jn_1 = eps[elem_id].il[ip].Jn_1;
      }

      /* compute pressure terms */
      double Pn_1,Pn,P;
      Pn_1 = Pn = P = 0.0;
      for(int i=0; i<nne; i++){
	Pn_1 += Np[i]*sig[elem_id].pn_1[i];
	Pn += Np[i]*sig[elem_id].p[i];
	P += Np[i]*(sig[elem_id].p[i] + dp[i]);
      }

      err += stab_get_material_potential(Ybar,kappa,Un_1,Jn_1,
					 J,Pn_1,Pn,P,C,mat);
    }
    break;
  case DISP:
    err += DISP_get_material_potential(kappa,mat,C,J,Ybar);
    break;
  default:
    *Ybar = 0.0;
    break;
  }

  return err;
}/* get_material_potential() */
