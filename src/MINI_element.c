/* HEADER */
/**
 * AUTHORS:
 * Matthew Mosby
 */
#include "MINI_element.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "mkl_cblas.h"

#include "PGFEM_io.h"
#include "allocation.h"
#include "two_field_element.h"
#include "utils.h"
#include "incl.h"

/* need to update to new_potentials for consistencey */
#include "potential.h"

#include "tensors.h"
#include "def_grad.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "cast_macros.h"
#include "elem3d.h"
#include "index_macros.h"

#ifndef MINI_DEBUG
#define MINI_DEBUG 0
#endif

#ifndef MINI_P1_P1
#define MINI_P1_P1 0
#endif

void MINI_reset(ELEMENT *elem,
		const int nelem,
		const int npres,
		SIG *sig)
{
  for(int i=0; i<nelem; i++){
    for(int j=0; j<npres; j++){
      sig[i].d_p[i] = 0.0;
    }

    for(int j=0; j<elem[i].n_bub*elem[i].n_bub_dofs; j++){
      elem[i].bub_dofs[j] = 0.0;
    }
  }
}

int MINI_stiffmat_el(double *Ks,            /**< Element stiffmat */
		      const int ii,          /**< element id */
		      const int ndofn,
		      const int nne,
		      const double *x,
		      const double *y,
		      const double *z,
		      const ELEMENT *elem,
		      const HOMMAT *hommat,
		      const long *nod,
		      const NODE *node,
		      const EPS *eps,
		      const SIG *sig,
		      const double *r_e)    /**< dof values on elem */
{
  int count; /* ALWAYS reset before use */
  int err = 0;

  /* 3D */ 
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  
  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += node[nod[i]].ndofn;
  }
  
  const int n_bub = elem[ii].n_bub;
  const int n_bub_dofs = elem[ii].n_bub_dofs;
  const int total_nne = nne + n_bub;

  double *disp;
  disp = (double*) PGFEM_calloc (total_nne*ndn,sizeof(double));

  double *p;
  p = (double*) PGFEM_calloc (nne,sizeof(double));

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<node[nod[i]].ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
	disp[i*ndn + j] = r_e[count+j];
      } else if (j==ndn){                       /* pressure dof */
	p[i] = r_e[count+j];
      } else {                                         /* ERROR */
	PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
	PGFEM_Abort();
      }
    } /* ndofn */
    count += node[nod[i]].ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] = 
	elem[ii].bub_dofs[i*n_bub_dofs + j];
    }
  }

  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet;
  int_pt_ksi = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_eta = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_zet = (double*) PGFEM_calloc (npt_z,sizeof(double));

  double *weights;
  weights = (double*) PGFEM_calloc (npt_z,sizeof(double));

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
	     int_pt_ksi,int_pt_eta,int_pt_zet,weights);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

  /* allocate space for the stiffness tensors */
  double *Kuu, *Ktt, *Kut, *Kup, *Kpu, *Ktp, *Kpt, *Kpp;
  Kuu = aloc1(nne*ndn*nne*ndn);
  Ktt = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Kut = aloc1(nne*ndn*n_bub*n_bub_dofs);
  Kup = aloc1(nne*ndn*nne);
  Kpu = aloc1(nne*ndn*nne);
  Ktp = aloc1(n_bub*n_bub_dofs*nne);
  Kpt = aloc1(n_bub*n_bub_dofs*nne);
  Kpp = aloc1(nne*nne);

  /* allocate space for deformation gradients etc */
  double **Fn_mat, *Fn, **Fr_mat, *Fr, *Fr_I;
  double **C_mat, *C, *C_I, *AA;
  double ****dSdF_tensor, *dSdF, **S_mat, *S;
  Fn_mat = aloc2(3,3);
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_I = aloc1(9);
  C_mat = aloc2(3,3);
  C = aloc1(9);
  C_I = aloc1(9);
  AA = aloc1(9);
  dSdF_tensor = aloc4(3,3,3,3);
  dSdF = aloc1(81);
  S_mat = aloc2(3,3);
  S = aloc1(9);

  /* integration point in natural coords */
  double ksi, eta, zet, wt;

  /* pressure */
  double pressure, Upp;

  /* integration loop */
  /* NOTE: integration on tetrahedron has only one integration
     direction */
  for (int ip=0; ip<npt_z; ip++){
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
    wt = weights[ip];

    /* get shape functions and derivative of shape functions and
       Jacobian of element.  NOTE: x, y, z must also contain the
       element center */
    shape_func(ksi,eta,zet,total_nne,Na);
    J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
    get_bubble_grad(total_nne,ksi,eta,zet,x,y,z,N_x,N_y,N_z);
    shape_tensor (total_nne,ndn,N_x,N_y,N_z,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,total_nne);

    /* compute pressure */
    pressure = 0;
    for(int i=0; i<nne; i++){
      pressure += Na[i]*(sig[ii].p[i] + p[i]);
    }

    /* get the deformation gradients */
    memcpy(Fn,eps[ii].st[ip].Fpp,9*sizeof(double));
    array2mat(Fn,Fn_mat,3,3);
    double Jn = getJacobian(Fn,ii,&err);

    def_grad_get(total_nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
    mat2array(Fr,CONST_2(double) Fr_mat,3,3);
    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,Fr,3,Fn,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		3,3,3,1.0,AA,3,AA,3,0.0,C,3);
    array2mat(C,C_mat,3,3);

    // Get Deviatoric 2 P-K stress
    deviatoricStressFunctionPtr Stress;
    materialStiffnessFunctionPtr Stiffness;
    d2UdJ2FunctionPtr D2UDJ2;
    Stress = getDeviatoricStressFunction(1,&hommat[mat]);
    Stiffness = getMaterialStiffnessFunction(1,&hommat[mat]);
    D2UDJ2 = getd2UdJ2Function(1,&hommat[mat]);

    Stress(CCONST_2(double) C_mat,&hommat[mat],S_mat);
    Stiffness(CCONST_2(double) C_mat,&hommat[mat],
	      CCONST_2(double) Fn_mat,
	      CCONST_2(double) Fr_mat,
	      dSdF_tensor);
    D2UDJ2(Jn,Jr,&hommat[mat],&Upp);

    /* convert from mat/tensor to array */
    mat2array(S,CONST_2(double) S_mat,3,3);
    tensor4_2array(dSdF,CONST_4(double) dSdF_tensor,3,3,3,3);

    /* Begin computing matrices */
    /*** Kuu (12 x 12) ***/
    UL_Kuu_at_ip(Kuu,nne,total_nne,ST,Fn,Fr,
		 Fr_I,Jn,Jr,pressure,S,dSdF,J,wt,0);

    /*** Ktt (Kuu on bubble) (3 x 3) ***/
    UL_Kuu_at_ip(Ktt,nne,total_nne,ST,Fn,Fr,
		 Fr_I,Jn,Jr,pressure,S,dSdF,J,wt,1);

    /*** Kut (12 x 3) ***/
    UL_Kuu_at_ip(Kut,nne,total_nne,ST,Fn,Fr,
		 Fr_I,Jn,Jr,pressure,S,dSdF,J,wt,2);

    /*** Kup (12 x 4) ***/
    UL_Kup_at_ip(Kup,nne,total_nne,Na,ST,Fr_I,Jr,J,wt,0);

    /*** Kpu (4 x 12) ***/
    UL_Kpu_at_ip(Kpu,nne,total_nne,Na,ST,Fr_I,Jn,Jr,Upp,J,wt,0);

    /*** Ktp (3 x 4) ***/
    UL_Kup_at_ip(Ktp,nne,total_nne,Na,ST,Fr_I,Jr,J,wt,1);

    /*** Kpt (4 x 3) ***/
    UL_Kpu_at_ip(Kpt,nne,total_nne,Na,ST,Fr_I,Jn,Jr,Upp,J,wt,1);

    /*** Kpp (4 x 4) ***/
    UL_Kpp_at_ip(Kpp,nne,total_nne,Na,kappa,J,wt);

  } /* int z-dir */

  /* free un-needed */
  free(disp);
  free(p);
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  dealoc2(Fn_mat,2);
  free(Fn);
  dealoc2(Fr_mat,3);
  free(Fr);
  free(Fr_I);
  free(AA);
  dealoc2(C_mat,3);
  free(C);
  free(C_I);
  dealoc4(dSdF_tensor,3,3,3);
  free(dSdF);
  dealoc2(S_mat,3);
  free(S);

  /* stuff for debugging */
  char filename[50];
  FILE *debug_log;
  int err_rank = 0;
  if(MINI_DEBUG){
    PGFEM_Error_rank(&err_rank);
    sprintf(filename,"MINI_stiff_%d_%d.log",ii,err_rank);
    debug_log = fopen(filename,"a"); /* APPEND TO FILE */
    PGFEM_fprintf(debug_log,"========================================\n");
    PGFEM_fprintf(debug_log,"============= New Computation ==========\n");
    PGFEM_fprintf(debug_log,"========================================\n\n");
    PGFEM_fprintf(debug_log,"Kuu\n");
    print_array_d(debug_log,Kuu,nne*ndn*nne*ndn,nne*ndn,nne*ndn);
    PGFEM_fprintf(debug_log,"Ktt\n");
    print_array_d(debug_log,Ktt,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
		  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kut\n");
    print_array_d(debug_log,Kut,nne*ndn*n_bub*n_bub_dofs,nne*ndn,
		  n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kup\n");
    print_array_d(debug_log,Kup,nne*ndn*nne,nne*ndn,nne);
    PGFEM_fprintf(debug_log,"Kpu\n");
    print_array_d(debug_log,Kpu,nne*ndn*nne,nne,nne*ndn);
    PGFEM_fprintf(debug_log,"Ktp\n");
    print_array_d(debug_log,Ktp,nne*n_bub*n_bub_dofs,n_bub*n_bub_dofs,nne);
    PGFEM_fprintf(debug_log,"Kpt\n");
    print_array_d(debug_log,Kpt,nne*n_bub*n_bub_dofs,nne,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kpp\n");
    print_array_d(debug_log,Kpp,nne*nne,nne,nne);
  }

  /*** Construct Condensed Matrices ***/
  double *Ktt_I, *KutKttI, *KptKttI;
  Ktt_I = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  KutKttI = aloc1(nne*ndn*n_bub*n_bub_dofs);
  KptKttI = aloc1(nne*n_bub*n_bub_dofs);

  assert(n_bub*n_bub*n_bub_dofs*n_bub_dofs == 9);
  if(!MINI_P1_P1){
    inverse(Ktt,3,Ktt_I);
  }

  /*** OPTIMIZE Unroll these multiplications ***/

  /* KutKttI (12 x 3) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
	      nne*ndn,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
	      Kut,n_bub*n_bub_dofs,Ktt_I,n_bub*n_bub_dofs,
	      0.0,KutKttI,n_bub*n_bub_dofs);

  /* KptKttI (4 x 3) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
	      nne,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
	      Kpt,n_bub*n_bub_dofs,Ktt_I,n_bub*n_bub_dofs,
	      0.0,KptKttI,n_bub*n_bub_dofs);

  if(MINI_DEBUG){
    PGFEM_fprintf(debug_log,"Ktt_I\n");
    print_array_d(debug_log,Ktt_I,9,3,3);
    PGFEM_fprintf(debug_log,"Kut Ktt_I\n");
    print_array_d(debug_log,KutKttI,12*3,12,3);
    PGFEM_fprintf(debug_log,"Kpt Ktt_I\n");
    print_array_d(debug_log,KptKttI,4*3,4,3);
  }
 
  /*** Kuu  ***/
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
  	      nne*ndn,nne*ndn,n_bub*n_bub_dofs,-1.0,
  	      KutKttI,n_bub*n_bub_dofs,Kut,n_bub*n_bub_dofs,
  	      1.0,Kuu,nne*ndn);

  /*** Kup ***/
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
  	      nne*ndn,nne,n_bub*n_bub_dofs,-1.0,
  	      KutKttI,n_bub*n_bub_dofs,Ktp,nne,
  	      1.0,Kup,nne);

  /*** Kpu ***/
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
  	      nne,nne*ndn,n_bub*n_bub_dofs,-1.0,
  	      KptKttI,n_bub*n_bub_dofs,Kut,n_bub*n_bub_dofs,
  	      1.0,Kpu,nne*ndn);

  /*** Kpp ***/
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
  	      nne,nne,n_bub*n_bub_dofs,-1.0,
  	      KptKttI,n_bub*n_bub_dofs,Ktp,nne,
  	      1.0,Kpp,nne);

  if(MINI_DEBUG){
    PGFEM_fprintf(debug_log,"Kuu BAR\n");
    print_array_d(debug_log,Kuu,12*12,12,12);
    PGFEM_fprintf(debug_log,"Kup BAR\n");
    print_array_d(debug_log,Kup,12*4,12,4);
    PGFEM_fprintf(debug_log,"Kpp BAR\n");
    print_array_d(debug_log,Kpp,4*4,4,4);
  }

  /* free un-needed */
  free(Ktt_I);
  free(KutKttI);
  free(KptKttI);

  /*** Scatter the matrix into the proper form ***/
  for(int a=0; a<nne; a++){
    for(int b=0; b<ndofn; b++){
      for(int w=0; w<nne; w++){
	for(int g=0; g<ndofn; g++){
	  if ((b<ndn) && (g<ndn)){ /* from Kuu */
	    Ks[idx_K(a,b,w,g,nne,ndofn)] = Kuu[idx_K(a,b,w,g,nne,ndn)];
	  } else if ((b==ndn) && (g==ndn)){ /* from Kpp */
	    Ks[idx_K(a,b,w,g,nne,ndofn)] = Kpp[idx_K(a,0,w,0,nne,1)];
	  } else if ((b==ndn) || (g==ndn)){ /* from Kup */
	    /* Kpu = Trans(Kup) so get appropriate value */
	    if (g==ndn){ /* From Kup */
	      Ks[idx_K(a,b,w,g,nne,ndofn)] = 
		Kup[idx_K_gen(a,b,w,0,nne,ndn,nne,1)];
	    } else { /* From Kpu */
	      Ks[idx_K(a,b,w,g,nne,ndofn)] = 
		Kpu[idx_K_gen(a,0,w,g,nne,1,nne,ndn)];
	    }
	  } else {
	    PGFEM_printerr("ERROR constructing stiffness matrix!\n");
	    PGFEM_Abort();
	    abort();
	  }
	}
      }
    }
  }

  /* free remaining */
  free(Kuu);
  free(Ktt);
  free(Kut);
  free(Kup);
  free(Kpu);
  free(Ktp);
  free(Kpt);
  free(Kpp);

  /* print out the element stiffness */
  if (MINI_DEBUG){
    PGFEM_fprintf(debug_log,"Element Stiffness\n");
    print_array_d(debug_log,Ks,nne*ndofn*nne*ndofn,nne*ndofn,nne*ndofn);
    fclose(debug_log);
  }

  return err;
} /* stiffmat_el_lin_bub */

int MINI_resid_el(double *Res,         /**< Element residual */
		   const int ii,        /**< element id */
		   const int ndofn,
		   const int nne,
		   const double *x,
		   const double *y,
		   const double *z,
		   const ELEMENT *elem,
		   const long *nod,
		   const NODE *node,
		   const HOMMAT *hommat,
		   const EPS *eps,
		   const SIG *sig,
		   const double *r_e)    /**< dof values on elem */
{
  int count; /* ALWAYS reset before use */
  int err = 0;

  /* 3D */ 
  int ndn = 3;
  int mat = elem[ii].mat[2];
  
  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += node[nod[i]].ndofn;
  }
  
  int n_bub = elem[ii].n_bub;
  int n_bub_dofs = elem[ii].n_bub_dofs;
  int total_nne = nne + n_bub;

  double *disp;
  disp = (double*) PGFEM_calloc (total_nne*ndn,sizeof(double));

  double *p;
  p = (double*) PGFEM_calloc (nne,sizeof(double));

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<node[nod[i]].ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
	disp[i*ndn + j] = r_e[count+j];
      } else if (j==ndn){                       /* pressure dof */
	p[i] = r_e[count+j];
      } else {                                         /* ERROR */
	PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
	PGFEM_Abort();
      }
    } /* ndofn */
    count += node[nod[i]].ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] = 
	elem[ii].bub_dofs[i*n_bub_dofs + j];
    }
  }

  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet;
  int_pt_ksi = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_eta = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_zet = (double*) PGFEM_calloc (npt_z,sizeof(double));

  double *weights;
  weights = (double*) PGFEM_calloc (npt_z,sizeof(double));

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
	     int_pt_ksi,int_pt_eta,int_pt_zet,weights);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

  /* allocate space for deformation gradients etc */
  double **Fn_mat, *Fn, **Fr_mat, *Fr, *Fr_I;
  double **C_mat, *C, *AA;
  double ****dSdF_tensor, *dSdF, **S_mat, *S;
  Fn_mat = aloc2(3,3);
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_I = aloc1(9);
  AA = aloc1(9);
  C_mat = aloc2(3,3);
  C = aloc1(9);
  dSdF_tensor = aloc4(3,3,3,3);
  dSdF = aloc1(81);
  S_mat = aloc2(3,3);
  S = aloc1(9);

  /* Allocate space for the residuals */
  double *Ru, *Rt, *Rp;
  Ru = aloc1(nne*ndn);
  Rt = aloc1(n_bub*n_bub_dofs);
  Rp = aloc1(nne);

  /* allocate space for the stiffness tensors */
  double *Ktt, *Kut, *Kpt;
  Ktt = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Kut = aloc1(nne*ndn*n_bub*n_bub_dofs);
  Kpt = aloc1(n_bub*n_bub_dofs*nne);

  /* integration point in natural coords */
  double ksi, eta, zet, wt;

  /* pressure */
  double pressure, Up, Upp;

  /* integration loop */
  /* NOTE: integration on tetrahedron has only one integration
     direction */
  for (int ip=0; ip<npt_z; ip++){
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
    wt = weights[ip];

    /* get shape functions and derivative of shape functions and
       Jacobian of element.  NOTE: x, y, z must also contain the
       element center */
    shape_func(ksi,eta,zet,total_nne,Na);
    J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
    get_bubble_grad(total_nne,ksi,eta,zet,x,y,z,N_x,N_y,N_z);
    shape_tensor (total_nne,ndn,N_x,N_y,N_z,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,total_nne);

    /* compute total pressure */
    pressure = 0;
    for(int i=0; i<nne; i++){
      pressure += Na[i]*(sig[ii].p[i] + p[i]);
    }

    /* get the deformation gradients */
    memcpy(Fn,eps[ii].st[ip].Fpp,9*sizeof(double));
    array2mat(Fn,Fn_mat,3,3);
    double Jn = getJacobian(Fn,ii,&err);

    def_grad_get(total_nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
    mat2array(Fr,CONST_2(double) Fr_mat,3,3);
    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    /* Compute total F; F=FnFr */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,Fr,3,Fn,3,0.0,AA,3);
    /* C = F^t F */
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		3,3,3,1.0,AA,3,AA,3,0.0,C,3);
    array2mat(C,C_mat,3,3);

    // Get Deviatoric 2 P-K stress
    deviatoricStressFunctionPtr Stress;
    materialStiffnessFunctionPtr Stiffness;
    volumetricPressureFunctionPtr Pressure;
    d2UdJ2FunctionPtr D2UDJ2;
    Stress = getDeviatoricStressFunction(1,&hommat[mat]);
    Stiffness = getMaterialStiffnessFunction(1,&hommat[mat]);
    Pressure = getVolumetricPressureFunction(1,&hommat[mat]);
    D2UDJ2 = getd2UdJ2Function(1,&hommat[mat]);

    Stress(CCONST_2(double) C_mat,&hommat[mat],S_mat);
    Stiffness(CCONST_2(double) C_mat,
	      &hommat[mat],
	      CCONST_2(double) Fn_mat,
	      CCONST_2(double) Fr_mat,
	      dSdF_tensor);
    Pressure(Jn,Jr,&hommat[mat],&Up);
    D2UDJ2(Jn,Jr,&hommat[mat],&Upp);

    /* convert from mat/tensor to array */
    mat2array(S,CONST_2(double) S_mat,3,3);
    tensor4_2array(dSdF,CONST_4(double) dSdF_tensor,3,3,3,3);

    /* Compute residuals */
    /***  Ru (12 x 1) ***/
    UL_Ru_at_ip(Ru,nne,total_nne,ST,Fn,Fr,Fr_I,
		Jn,Jr,S,pressure,J,wt,0);

    /*** Rt (3 x 1) ***/
    UL_Ru_at_ip(Rt,nne,total_nne,ST,Fn,Fr,Fr_I,
		Jn,Jr,S,pressure,J,wt,1);

    /* Compute needed tangents */
    /*** Ktt (3 x 3) ***/
    UL_Kuu_at_ip(Ktt,nne,total_nne,ST,Fn,Fr,Fr_I,Jn,Jr,
		 pressure,S,dSdF,J,wt,1);

    /*** Kut (12 x 3) ***/
    UL_Kuu_at_ip(Kut,nne,total_nne,ST,Fn,Fr,Fr_I,Jn,Jr,
		 pressure,S,dSdF,J,wt,2);

    /* compute incremental pressure */
    pressure = 0;
    for(int i=0; i<nne; i++){
      pressure += Na[i]*(p[i]);
    }

    /*** Rp (4 x 1) ***/
    UL_Rp_at_ip(Rp,nne,Na,kappa,Up,pressure,J,wt);

    /*** Kpt (4 x 3) ***/
    UL_Kpu_at_ip(Kpt,nne,total_nne,Na,ST,Fr_I,Jn,Jr,Upp,J,wt,1);

  }/* end integration */

  /* free un-needed */
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  dealoc2(Fn_mat,2);
  free(Fn);
  dealoc2(Fr_mat,3);
  free(Fr);
  free(Fr_I);
  free(AA);
  dealoc2(C_mat,3);
  free(C);
  dealoc4(dSdF_tensor,3,3,3);
  free(dSdF);
  dealoc2(S_mat,3);
  free(S);

  /* stuff for debugging */
  char filename[50];
  FILE *debug_log;
  int err_rank = 0;
  if(MINI_DEBUG){
    PGFEM_Error_rank(&err_rank);
    sprintf(filename,"MINI_resid_%d_%d.log",ii,err_rank);
    debug_log = fopen(filename,"a"); /* APPEND TO FILE */
    PGFEM_fprintf(debug_log,"========================================\n");
    PGFEM_fprintf(debug_log,"============= New Computation ==========\n");
    PGFEM_fprintf(debug_log,"========================================\n\n");
    PGFEM_fprintf(debug_log,"Disp\n");
    print_array_d(debug_log,disp,total_nne*ndn,1,total_nne*ndn);
    PGFEM_fprintf(debug_log,"Pressure\n");
    print_array_d(debug_log,p,nne,1,nne);
  }
  free(disp);
  free(p);

  /*** Construct condensed residuals ***/
  double *Ktt_I, *KutKttI, *KptKttI;
  Ktt_I = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  KutKttI = aloc1(nne*ndn*n_bub*n_bub_dofs);
  KptKttI = aloc1(nne*n_bub*n_bub_dofs);

  assert(n_bub*n_bub*n_bub_dofs*n_bub_dofs == 9);
  if(!MINI_P1_P1){
    inverse(Ktt,3,Ktt_I);
  }

  /* KutKttI (12 x 3) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
	      nne*ndn,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
	      Kut,n_bub*n_bub_dofs,Ktt_I,n_bub*n_bub_dofs,
	      0.0,KutKttI,n_bub*n_bub_dofs);

  /* KptKttI (4 x 3) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
	      nne,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
	      Kpt,n_bub*n_bub_dofs,Ktt_I,n_bub*n_bub_dofs,
	      0.0,KptKttI,n_bub*n_bub_dofs);

  /*** Ru (12 x 1) ***/
  cblas_dgemv(CblasRowMajor,CblasNoTrans,
  	      nne*ndn,n_bub*n_bub_dofs,1.0,
  	      KutKttI,n_bub*n_bub_dofs,Rt,1,
  	      -1.0,Ru,1);

  /*** Rp (4 x 1) ***/
  cblas_dgemv(CblasRowMajor,CblasNoTrans,
  	      nne,n_bub*n_bub_dofs,1.0,
  	      KptKttI,n_bub*n_bub_dofs,Rt,1,
  	      -1.0,Rp,1);

  /* free un-needed */
  free(Ktt);
  free(Kpt);
  free(Kut);
  free(Ktt_I);
  free(KutKttI);
  free(KptKttI);
  free(Rt);

  int sign = -1;
  /*** Scatter the residuals (u,v,w,p) ***/
  /* NOTE residuals are negative because subtracted from RHS later */
  for (int a=0; a<nne; a++){
    for (int b=0; b<ndofn; b++){
      if (b<ndn){ /* from Ru */
	Res[a*ndofn+b] = sign*Ru[a*ndn+b];
      } else if (b == ndn){ /* from Rp */
	Res[a*ndofn+b] = sign*Rp[a];
      } else {
	PGFEM_printerr("ERROR constructing residual vector!\n");
	PGFEM_Abort();
	abort();
      }
    }
  }

  /* print out the element residual */
  if (MINI_DEBUG){
    PGFEM_fprintf(debug_log,"========== Element Residual ==========\n");
    print_array_d(debug_log,Res,nne*ndofn,1,nne*ndofn);
    fclose(debug_log);
  }

  /* free remaining */
  free(Ru);
  free(Rp);

  /* Temporary exit after first computation */
  /* PGFEM_Abort(); */
  /* abort(); */

  return err;

}

int MINI_update_bubble_el(ELEMENT *elem,
			   const int ii, /* id of element working on */
			   const int nne,
			   const NODE *node,
			   const int ndofn,
			   const double *x,
			   const double *y,
			   const double *z,
			   const long *nod, /* list of node ids on elem */
			   const EPS *eps,
			   const SIG *sig,
			   const HOMMAT *hommat,
			   const double *sol_e, /* accum. incr on elem */
			   const double *dsol_e) /* increment on element */
{
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  int err = 0;

  int count; /* ALWAYS reset before use */

  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += node[nod[i]].ndofn;
  }

  const int n_bub = elem[ii].n_bub;
  const int n_bub_dofs = elem[ii].n_bub_dofs;
  const int total_nne = nne + n_bub;

  /* stuff for debugging purposes */
  FILE * debug_log;
  int err_rank = 0;

  if(MINI_DEBUG){
    PGFEM_Error_rank(&err_rank);
    char fname[50];
    sprintf(fname,"MINI_up_%d_%d.log",err_rank,ii);
    debug_log = fopen(fname,"a");
    PGFEM_fprintf(debug_log,"***********************************\n");
    PGFEM_fprintf(debug_log,"Current bub dofs\n");
    print_array_d(debug_log,elem[ii].bub_dofs,
		  n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  /* disp contains displacements from nodes + bub */
  double *disp, *ddisp;
  disp = (double*) PGFEM_calloc (total_nne*ndn,sizeof(double));
  ddisp = (double*) PGFEM_calloc (nne*ndn,sizeof(double));

  double *p, *dp;
  p = (double*) PGFEM_calloc (nne,sizeof(double));
  dp = (double*) PGFEM_calloc (nne,sizeof(double));

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<node[nod[i]].ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
	disp[i*ndn + j] = sol_e[count+j];
	ddisp[i*ndn + j] = dsol_e[count+j];
      } else if (j==ndn){                       /* pressure dof */
	p[i] = sol_e[count+j];
	dp[i] = dsol_e[count+j];
      } else {                                         /* ERROR */
	PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
	PGFEM_Abort();
      }
    } /* ndofn */
    count += node[nod[i]].ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] = 
	elem[ii].bub_dofs[i*n_bub_dofs + j];
    }
  }

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet;
  int_pt_ksi = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_eta = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_zet = (double*) PGFEM_calloc (npt_z,sizeof(double));

  double *weights;
  weights = (double*) PGFEM_calloc (npt_z,sizeof(double));

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
	     int_pt_ksi,int_pt_eta,int_pt_zet,weights);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

  /* allocate space for the stiffness tensors */
  double *Ktt, *Kut, *Ktp, *Rt;
  Ktt = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Kut = aloc1(nne*ndn*n_bub*n_bub_dofs);
  Ktp = aloc1(n_bub*n_bub_dofs*nne);
  Rt = aloc1(n_bub*n_bub_dofs);

  /* allocate space for deformation gradients etc */
  double **Fn_mat, *Fn, **Fr_mat, *Fr, *Fr_I;
  double **C_mat, *C, *C_I, *AA;
  double ****dSdF_tensor, *dSdF, **S_mat, *S;
  Fn_mat = aloc2(3,3);
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_I = aloc1(9);
  C_mat = aloc2(3,3);
  C = aloc1(9);
  C_I = aloc1(9);
  AA = aloc1(9);
  dSdF_tensor = aloc4(3,3,3,3);
  dSdF = aloc1(81);
  S_mat = aloc2(3,3);
  S = aloc1(9);

  /* integration point in natural coords */
  double ksi, eta, zet, wt;

  /* pressure */
  double pressure, Upp;

  /* integration loop */
  /* NOTE: integration on tetrahedron has only one integration
     direction */
  for (int ip=0; ip<npt_z; ip++){
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
    wt = weights[ip];

    /* get shape functions and derivative of shape functions and
       Jacobian of element.  NOTE: x, y, z must also contain the
       element center */
    shape_func(ksi,eta,zet,total_nne,Na);
    J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
    get_bubble_grad(total_nne,ksi,eta,zet,x,y,z,N_x,N_y,N_z);
    shape_tensor (total_nne,ndn,N_x,N_y,N_z,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,total_nne);

    /* compute pressure */
    pressure = 0;
    for(int i=0; i<nne; i++){
      pressure += Na[i]*(sig[ii].p[i] + p[i]);
    }

    /* get the deformation gradients */
    memcpy(Fn,eps[ii].st[ip].Fpp,9*sizeof(double));
    array2mat(Fn,Fn_mat,3,3);
    double Jn = getJacobian(Fn,ii,&err);

    def_grad_get(total_nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
    mat2array(Fr,CONST_2(double) Fr_mat,3,3);
    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,Fr,3,Fn,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		3,3,3,1.0,AA,3,AA,3,0.0,C,3);
    array2mat(C,C_mat,3,3);

    // Get Deviatoric 2 P-K stress
    deviatoricStressFunctionPtr Stress;
    materialStiffnessFunctionPtr Stiffness;
    d2UdJ2FunctionPtr D2UDJ2;
    Stress = getDeviatoricStressFunction(1,&hommat[mat]);
    Stiffness = getMaterialStiffnessFunction(1,&hommat[mat]);
    D2UDJ2 = getd2UdJ2Function(1,&hommat[mat]);

    Stress(CCONST_2(double) C_mat,&hommat[mat],S_mat);
    Stiffness(CCONST_2(double) C_mat,
	      &hommat[mat],
	      CCONST_2(double) Fn_mat,
	      CCONST_2(double) Fr_mat,
	      dSdF_tensor);
    D2UDJ2(Jn,Jr,&hommat[mat],&Upp);

    /* convert from mat/tensor to array */
    mat2array(S,CONST_2(double) S_mat,3,3);
    tensor4_2array(dSdF,CONST_4(double) dSdF_tensor,3,3,3,3);

    /* Begin computing matrices */

    /*** Ktt (Kuu on bubble) (3 x 3) ***/
    UL_Kuu_at_ip(Ktt,nne,total_nne,ST,Fn,Fr,
		 Fr_I,Jn,Jr,pressure,S,dSdF,J,wt,1);

    /*** Kut (12 x 3) ***/
    UL_Kuu_at_ip(Kut,nne,total_nne,ST,Fn,Fr,
		 Fr_I,Jn,Jr,pressure,S,dSdF,J,wt,2);

    /*** Ktp (3 x 4) ***/
    UL_Kup_at_ip(Ktp,nne,total_nne,Na,ST,Fr_I,Jr,J,wt,1);

    /*** Rt (3 x 1) ***/
    UL_Ru_at_ip(Rt,nne,total_nne,ST,Fn,Fr,Fr_I,
		Jn,Jr,S,pressure,J,wt,1);

  } /* int z-dir */

  /* free un-needed */
  free(disp);
  free(p);
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  dealoc2(Fn_mat,2);
  free(Fn);
  dealoc2(Fr_mat,3);
  free(Fr);
  free(Fr_I);
  free(AA);
  dealoc2(C_mat,3);
  free(C);
  free(C_I);
  dealoc4(dSdF_tensor,3,3,3);
  free(dSdF);
  dealoc2(S_mat,3);
  free(S);

  /*** Compute the bubble update ***/
  double *Ktt_I, *ddt;
  ddt = aloc1(n_bub*n_bub_dofs);
  Ktt_I = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  if(!MINI_P1_P1){
    inverse(Ktt,3,Ktt_I);
  }

  if(MINI_DEBUG){
    PGFEM_fprintf(debug_log,"Ktt\n");
    print_array_d(debug_log,Ktt,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
		  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Ktt_I\n");
    print_array_d(debug_log,Ktt_I,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
		  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kut\n");
    print_array_d(debug_log,Kut,n_bub*n_bub_dofs*nne*ndn,
		  nne*ndn,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Ktp\n");
    print_array_d(debug_log,Ktp,n_bub*n_bub_dofs*nne,
		  n_bub*n_bub_dofs,nne);
    PGFEM_fprintf(debug_log,"ddisp\n");
    print_array_d(debug_log,ddisp,nne*ndn,1,nne*ndn);
    PGFEM_fprintf(debug_log,"dp\n");
    print_array_d(debug_log,dp,nne,1,nne);
    PGFEM_fprintf(debug_log,"Rt\n");
    print_array_d(debug_log,Rt,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  /* Rt += Kut'ddu */
  cblas_dgemv(CblasRowMajor,CblasTrans,12,3,1.0,Kut,3,ddisp,1,1.0,Rt,1);

  if(MINI_DEBUG){
    PGFEM_fprintf(debug_log,"Rt += Kut' ddisp\n");
    print_array_d(debug_log,Rt,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  /* Rt += Ktp ddp */
  cblas_dgemv(CblasRowMajor,CblasNoTrans,3,4,1.0,Ktp,4,dp,1,1.0,Rt,1);

  if(MINI_DEBUG){
    PGFEM_fprintf(debug_log,"Rt += Ktp dp\n");
    print_array_d(debug_log,Rt,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  /* ddt = -Ktt_I Rt */
  cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,-1.0,Ktt_I,3,Rt,1,0.0,ddt,1);

  if(MINI_DEBUG){
    PGFEM_fprintf(debug_log,"dt\n");
    print_array_d(debug_log,ddt,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  for (int i=0; i<n_bub*n_bub_dofs; i++){
    elem[ii].bub_dofs[i] += ddt[i];
  }

  if (MINI_DEBUG){
    PGFEM_fprintf(debug_log,"Updated bubble dofs\n");
    print_array_d(debug_log,elem[ii].bub_dofs,
  		  n_bub*n_bub_dofs,
  		  1,n_bub*n_bub_dofs);
    fclose(debug_log);
  }

  free(Ktt);
  free(Ktt_I);
  free(Kut);
  free(Ktp);
  free(Rt);
  free(ddt);
  free(ddisp);
  free(dp);

  return err;
}

int MINI_update_bubble(ELEMENT *elem,
		       const int nelem,
		       const NODE *node,
		       const int ndofn,
		       SUPP sup,
		       const EPS *eps,
		       const SIG *sig,
		       const HOMMAT *hommat,
		       const double *sol, /* accum. solution  on incr */
		       const double *dsol, /* sol from current iter */
		       const int iter)
{
  int err = 0;

  int nne, mat, nne_t;

  long *nod, *cn;
  double *x, *y, *z, *sup_def;
  double *sol_e, *dsol_e;

  sup_def = aloc1(sup->npd);

  /* for each element */
  for (int i=0; i<nelem; i++){
    nne = elem[i].toe;
    nne_t = nne + elem[i].n_bub;
    mat = elem[i].mat[2];

    /* allocate */
    nod = aloc1l(nne);
    cn = aloc1l(nne*ndofn);
    x = aloc1(nne_t);
    y = aloc1(nne_t);
    z = aloc1(nne_t);
    sol_e = aloc1(nne*ndofn);
    dsol_e = aloc1(nne*ndofn);

    /* get node ids on element */
    elemnodes(i,nne,nod,elem);

    /* get local dof ids on elemnt */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);

    /* Get the node and bubble coordinates */
    nodecoord_updated(nne,nod,node,x,y,z);
    element_center(nne,x,y,z);

    /*** compute the iteration and accumlated solution vectors on the
	 element ***/
    /* solution for this iteration */
    if (iter != 0){
      for (int j=0;j<sup->npd;j++){
	sup_def[j] = sup->defl_d[j];
	sup->defl_d[j] = 0.0;
      }
    }
    def_elem(cn,nne*ndofn,dsol,elem,node,dsol_e,sup,0);
    if (iter != 0){
      for (int j=0;j<sup->npd;j++)
	sup->defl_d[j] = sup_def[j];
    }
    /* accumulated solution for this load increment */
    if (iter == 0){
      for (int j=0;j<sup->npd;j++){
	sup_def[j] = sup->defl_d[j];
	sup->defl_d[j] = 0.0;
      }
    }
    def_elem(cn,nne*ndofn,sol,elem,node,sol_e,sup,0);
    if (iter == 0){
      for (int j=0;j<sup->npd;j++)
	sup->defl_d[j] = sup_def[j];
    }

    if(MINI_DEBUG){
      int err_rank = 0;
      PGFEM_Error_rank(&err_rank);
      if(err_rank==0){
	PGFEM_printf("Solution vector on element %d:\n",i);
	print_array_d(stdout,dsol_e,nne*ndofn,1,nne*ndofn);
      }
    }

    /* compute the bubble update on the element */
    err = MINI_update_bubble_el(elem,i,nne,node,ndofn,x,y,z,nod,
				eps,sig,hommat,sol_e,dsol_e);


    free(nod);
    free(cn);
    free(x);
    free(y);
    free(z);
    free(sol_e);
    free(dsol_e);
  }
  free(sup_def);

  return err;
}

void MINI_increment_el(ELEMENT *elem,
		       const int ii, /* id of element working on */
		       const int nne,
		       const NODE *node,
		       const long *nod,
		       const int ndofn,
		       const double *x,
		       const double *y,
		       const double *z,
		       EPS *eps,
		       SIG *sig,
		       const HOMMAT *hommat,
		       const double *sol_e)
{
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  int err = 0;

  int count; /* ALWAYS reset before use */

  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += node[nod[i]].ndofn;
  }

  const int n_bub = elem[ii].n_bub;
  const int n_bub_dofs = elem[ii].n_bub_dofs;
  const int total_nne = nne + n_bub;

  double *disp;
  disp = (double*) PGFEM_calloc (total_nne*ndn,sizeof(double));

  double *p;
  p = (double*) PGFEM_calloc (nne,sizeof(double));

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<node[nod[i]].ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
	disp[i*ndn + j] = sol_e[count+j];
      } else if (j==ndn){                       /* pressure dof */
	p[i] = sol_e[count+j];
      } else {                                         /* ERROR */
	PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
	PGFEM_Abort();
      }
    } /* ndofn */
    count += node[nod[i]].ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] = 
	elem[ii].bub_dofs[i*n_bub_dofs + j];
    }
  }

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet;
  int_pt_ksi = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_eta = (double*) PGFEM_calloc (npt_z,sizeof(double));
  int_pt_zet = (double*) PGFEM_calloc (npt_z,sizeof(double));

  double *weights;
  weights = (double*) PGFEM_calloc (npt_z,sizeof(double));

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
	     int_pt_ksi,int_pt_eta,int_pt_zet,weights);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

  /* allocate space for deformation gradients etc */
  double **Fn_mat, *Fn, **Fr_mat, *Fr, *Fr_I;
  double **C_mat, *C, *C_I, *E, *AA, *F_total, *F_total_I;
  double **S_mat, *S;
  Fn_mat = aloc2(3,3);
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_I = aloc1(9);
  C_mat = aloc2(3,3);
  C = aloc1(9);
  C_I = aloc1(9);
  E = aloc1(9);
  AA = aloc1(9);
  F_total = aloc1(9);
  F_total_I = aloc1(9);
  S_mat = aloc2(3,3);
  S = aloc1(9);

  /* identity */
  double *ident;
  ident = aloc1(9);
  ident[0] = ident[4] = ident[8] = 1.0;

  /* integration point in natural coords */
  double ksi, eta, zet, wt;

  /* pressure */
  double pressure;

  /* element volume */
  const double volume = Tetra_V(x,y,z);

  /* integration loop */
  /* NOTE: integration on tetrahedron has only one integration
     direction */
  for (int ip=0; ip<npt_z; ip++){
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
    wt = weights[ip];

    /* get shape functions and derivative of shape functions and
       Jacobian of element.  NOTE: x, y, z must also contain the
       element center */
    shape_func(ksi,eta,zet,total_nne,Na);
    J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
    get_bubble_grad(total_nne,ksi,eta,zet,x,y,z,N_x,N_y,N_z);
    shape_tensor (total_nne,ndofn-1,N_x,N_y,N_z,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,total_nne);

    /* compute pressure */
    pressure = 0;
    for(int i=0; i<nne; i++){
      pressure += Na[i]*(sig[ii].p[i] + p[i]);
    }

    /* get the deformation gradients */
    memcpy(Fn,eps[ii].st[ip].Fpp,9*sizeof(double));
    array2mat(Fn,Fn_mat,3,3);
    double Jn = getJacobian(Fn,ii,&err);

    def_grad_get(total_nne,ndofn-1,CONST_4(double) ST_tensor,disp,Fr_mat);
    mat2array(Fr,CONST_2(double) Fr_mat,3,3);
    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    /* total deformation gradient */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,Fr,3,Fn,3,0.0,F_total,3);

    /* Green-Lagrange deformation tensor */
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		3,3,3,1.0,F_total,3,F_total,3,0.0,C,3);
    array2mat(C,C_mat,3,3);
    inverse(C,3,C_I);

    // Get Deviatoric 2 P-K stress
    deviatoricStressFunctionPtr Stress;
    Stress = getDeviatoricStressFunction(1,&hommat[mat]);
    Stress(CCONST_2(double) C_mat,&hommat[mat],S_mat);
    mat2array(S,CONST_2(double) S_mat,3,3);

    /* Compute total stress S = dev(S) + p*Jn*Jr*C_I */
    cblas_daxpy(9,pressure*Jn*Jr,C_I,1,S,1);

    /* Compute G-L strain E = 0.5(C - 1) */
    for(int i=0; i<9; i++){
      E[i] = 0.5*(C[i] - ident[i]);
    }

    /* Compute Cauchy stress (Jn*Jr)^-1 F S F' */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,F_total,3,S,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
		3,3,3,1./(Jn*Jr),AA,3,F_total,3,0.0,S,3);

    sig[ii].el.o[0] += wt*J/volume*S[idx_2(0,0)];
    sig[ii].el.o[1] += wt*J/volume*S[idx_2(1,1)];
    sig[ii].el.o[2] += wt*J/volume*S[idx_2(2,2)];

    sig[ii].el.o[3] += wt*J/volume*S[idx_2(1,2)];
    sig[ii].el.o[4] += wt*J/volume*S[idx_2(0,2)];
    sig[ii].el.o[5] += wt*J/volume*S[idx_2(0,1)];

    /* Compute logarithmic strain */
    inverse(F_total,3,F_total_I);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		3,3,3,1.0,F_total_I,3,E,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,AA,3,F_total_I,3,0.0,E,3);

    eps[ii].el.o[0] += wt*J/volume*E[idx_2(0,0)];
    eps[ii].el.o[1] += wt*J/volume*E[idx_2(1,1)];
    eps[ii].el.o[2] += wt*J/volume*E[idx_2(2,2)];
	
    eps[ii].el.o[3] += 2.*wt*J/volume*E[idx_2(1,2)];
    eps[ii].el.o[4] += 2.*wt*J/volume*E[idx_2(0,2)];
    eps[ii].el.o[5] += 2.*wt*J/volume*E[idx_2(0,1)];

    /* update the deformation gradient at the integration point */
    memcpy(eps[ii].st[ip].Fpp,F_total,9*sizeof(double));
  } /* integration */

    /* update pressure in pressure */
  for (int i=0; i<nne; i++){
    sig[ii].p[i] += p[i];
  }

  /* null the bubble dofs */
  for (int i=0; i<n_bub*n_bub_dofs; i++){
    elem[ii].bub_dofs[i] = 0.0;
  }

  free(disp);
  free(p);
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  dealoc2(Fn_mat,2);
  free(Fn);
  dealoc2(Fr_mat,3);
  free(Fr);
  free(Fr_I);
  free(AA);
  free(F_total);
  free(F_total_I);
  dealoc2(C_mat,3);
  free(C);
  free(C_I);
  dealoc2(S_mat,3);
  free(S);
  free(E);
  free(ident);
}

void MINI_increment(ELEMENT *elem,
		    const int nelem,
		    NODE *node,
		    const int nnodes,
		    const int ndofn,
		    SUPP sup,
		    EPS *eps,
		    SIG *sig,
		    const HOMMAT *hommat,
		    const double *sol,
		    const MPI_Comm mpi_comm)
{
  const int ndn = 3;

  int nne, mat, nne_t, II;

  long *nod, *cn;
  double *x, *y, *z;
  double *sol_e;

  /* for each element */
  for (int i=0; i<nelem; i++){
    nne = elem[i].toe;
    nne_t = nne + elem[i].n_bub;
    mat = elem[i].mat[2];

    /* allocate */
    nod = aloc1l(nne);
    cn = aloc1l(nne*ndofn);
    x = aloc1(nne_t);
    y = aloc1(nne_t);
    z = aloc1(nne_t);
    sol_e = aloc1(nne*ndofn);

    /* get node ids on element */
    elemnodes(i,nne,nod,elem);

    /* get local dof ids on elemnt */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);

    /* Get the node and bubble coordinates */
    nodecoord_updated(nne,nod,node,x,y,z);
    element_center(nne,x,y,z);

    def_elem(cn,nne*ndofn,sol,elem,node,sol_e,sup,0);

    /* increment the element quantities */
    MINI_increment_el(elem,i,nne,node,nod,ndofn,
		      x,y,z,eps,sig,hommat,sol_e);

    free(nod);
    free(cn);
    free(x);
    free(y);
    free(z);
    free(sol_e);
  } /* For each element */

  /*** Update Coordinates ***/
  for (int ii=0;ii<nnodes;ii++){
    for (int i=0;i<ndn;i++){
      II = node[ii].id[i];
      if (II > 0){
	if (i == 0) node[ii].x1 += sol[II-1];
	if (i == 1) node[ii].x2 += sol[II-1];
	if (i == 2) node[ii].x3 += sol[II-1];
      }
      if (II < 0){
	if (i == 0) node[ii].x1 += sup->defl_d[abs(II)-1];
	if (i == 1) node[ii].x2 += sup->defl_d[abs(II)-1];
	if (i == 2) node[ii].x3 += sup->defl_d[abs(II)-1];
      }
    }
  }/* end ii < nn */

  int myrank = 0;
  double PL, GPL;
  MPI_Comm_rank(mpi_comm,&myrank);
  PL = T_VOLUME (nelem,ndofn-1,elem,node);
  /* Gather Volume from all domains */
  MPI_Allreduce(&PL,&GPL,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  if (myrank == 0) {
    PGFEM_printf ("AFTER DEF - VOLUME = %12.12f\n",GPL);
  }
}

void MINI_check_resid(const int ndofn,
		      const int ne,
		      const ELEMENT *elem,
		      const NODE *node,
		      const HOMMAT *hommat,
		      const EPS *eps,
		      const SIG *sig,
		      const double *d_r,
		      SUPP sup,
		      const double *RR,
		      const long *DomDof,
		      const int ndofd,
		      const int GDof,
		      const COMMUN comm,
		      const MPI_Comm mpi_comm)
{

  /* compute the norm of the residauls for each variable */
  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);
  const int ndn = 3;

  int count; /* ALWAYS reset before use */
  int err = 0;
  long *nod, *cn;
  double *r_e, *x, *y, *z, *fe;
  double nrT, nrB;
  nrT = nrB = 0.0;

  double *f, *f_u, *f_p;
  f = aloc1(ndofd);
  f_u = aloc1(ndofd);
  f_p = aloc1(ndofd);

  for(int ii=0; ii<ne; ii++){

    /* Number of element nodes */
    const int nne = elem[ii].toe;
    const int nPres = nne;
    const int total_nne = nne + elem[ii].n_bub;
    /* Nodes on element */
    nod = aloc1l (nne);
    elemnodes (ii,nne,nod,elem);
    /* Element Dof */
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node);

    const int n_bub = elem[ii].n_bub;
    const int n_bub_dofs = elem[ii].n_bub_dofs;
    const int mat = elem[ii].mat[2];
    const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

    /* allocation */
    r_e = aloc1 (ndofe); /* container for accum solution */
    x = aloc1 (total_nne);
    y = aloc1 (total_nne);
    z = aloc1 (total_nne);
    fe = aloc1 (ndofe); /* container for element global res */
    cn = aloc1l (ndofe);

    /* coordinates */
    nodecoord_updated (nne,nod,node,x,y,z);
    element_center(nne,x,y,z);

    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    /* deformation on element */
    def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);

    /*============== COMPUTE ELEMENT RESIDUALS ================*/
    double *disp, *p;
    disp = aloc1(total_nne*ndn);
    p = aloc1(nPres);

    count = 0;
    for (int i=0; i<nne; i++){
      for (int j=0; j<node[nod[i]].ndofn; j++){

	/* filter displacement and pressure from element unknowns */
	if (j<ndn){                           /* displacement dof */
	  disp[i*ndn + j] = r_e[count+j];
	} else if (j==ndn){                       /* pressure dof */
	  p[i] = r_e[count+j];
	} else {                                         /* ERROR */
	  PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
	  PGFEM_Abort();
	}
      } /* ndofn */
      count += node[nod[i]].ndofn;
    } /* nne */

    /* displacements from bubble */
    for (int i=0; i<n_bub; i++){
      for (int j=0; j<n_bub_dofs; j++){
	disp[nne*ndn + i*n_bub_dofs + j] = 
	  elem[ii].bub_dofs[i*n_bub_dofs + j];
      }
    }

    /*=== BEGIN INTEGRATION ===*/
    long npt_x, npt_y, npt_z;
    int_point(total_nne,&npt_z);

    double *int_pt_ksi, *int_pt_eta, *int_pt_zet;
    int_pt_ksi = (double*) PGFEM_calloc (npt_z,sizeof(double));
    int_pt_eta = (double*) PGFEM_calloc (npt_z,sizeof(double));
    int_pt_zet = (double*) PGFEM_calloc (npt_z,sizeof(double));

    double *weights;
    weights = (double*) PGFEM_calloc (npt_z,sizeof(double));

    /* get the integration points and weights */
    integrate (total_nne,&npt_x,&npt_y,&npt_z,
	       int_pt_ksi,int_pt_eta,int_pt_zet,weights);

    /* allocate space for the shape functions, derivatives etc */
    double *Na, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
    Na = aloc1 (total_nne);
    N_x = aloc1 (total_nne);
    N_y = aloc1 (total_nne);
    N_z = aloc1 (total_nne);
    ST_tensor = aloc4 (3,3,ndn,total_nne);
    ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

    /* allocate space for deformation gradients etc */
    double **Fn_mat, *Fn, **Fr_mat, *Fr, *Fr_I;
    double **C_mat, *C, *AA;
    double **S_mat, *S;
    Fn_mat = aloc2(3,3);
    Fn = aloc1(9);
    Fr_mat = aloc2(3,3);
    Fr = aloc1(9);
    Fr_I = aloc1(9);
    AA = aloc1(9);
    C_mat = aloc2(3,3);
    C = aloc1(9);
    S_mat = aloc2(3,3);
    S = aloc1(9);

    /* Allocate space for the residuals */
    double *Ru, *Rt, *Rp;
    Ru = aloc1(nne*ndn);
    Rt = aloc1(n_bub*n_bub_dofs);
    Rp = aloc1(nne);

    /* integration point in natural coords */
    double ksi, eta, zet, wt;

    /* pressure */
    double pressure, Up, Upp;

    /* integration loop */
    /* NOTE: integration on tetrahedron has only one integration
       direction */
    for (int ip=0; ip<npt_z; ip++){
      ksi = int_pt_ksi[ip];
      eta = int_pt_eta[ip];
      zet = int_pt_zet[ip];
      wt = weights[ip];

      /* get shape functions and derivative of shape functions and
	 Jacobian of element.  NOTE: x, y, z must also contain the
	 element center */
      shape_func(ksi,eta,zet,total_nne,Na);
      J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
      get_bubble_grad(total_nne,ksi,eta,zet,x,y,z,N_x,N_y,N_z);
      shape_tensor (total_nne,ndn,N_x,N_y,N_z,ST_tensor);
      shapeTensor2array(ST,CONST_4(double) ST_tensor,total_nne);

      /* compute total pressure */
      pressure = 0;
      for(int i=0; i<nPres; i++){
	pressure += Na[i]*(sig[ii].p[i] + p[i]);
      }


    /* get the deformation gradients */
    memcpy(Fn,eps[ii].st[ip].Fpp,9*sizeof(double));
    array2mat(Fn,Fn_mat,3,3);
    double Jn = getJacobian(Fn,ii,&err);

    def_grad_get(total_nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
    mat2array(Fr,CONST_2(double) Fr_mat,3,3);
    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    /* Compute total F; F=FnFr */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		3,3,3,1.0,Fr,3,Fn,3,0.0,AA,3);
    /* C = F^t F */
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		3,3,3,1.0,AA,3,AA,3,0.0,C,3);
    array2mat(C,C_mat,3,3);

    // Get Deviatoric 2 P-K stress
    deviatoricStressFunctionPtr Stress;
    volumetricPressureFunctionPtr Pressure;
    d2UdJ2FunctionPtr D2UDJ2;
    Stress = getDeviatoricStressFunction(1,&hommat[mat]);
    Pressure = getVolumetricPressureFunction(1,&hommat[mat]);
    D2UDJ2 = getd2UdJ2Function(1,&hommat[mat]);

    Stress(CCONST_2(double) C_mat,&hommat[mat],S_mat);
    Pressure(Jn,Jr,&hommat[mat],&Up);
    D2UDJ2(Jn,Jr,&hommat[mat],&Upp);

    /* convert from mat/tensor to array */
    mat2array(S,CONST_2(double) S_mat,3,3);

    /* Compute residuals */
    /***  Ru (12 x 1) ***/
    UL_Ru_at_ip(Ru,nne,total_nne,ST,Fn,Fr,Fr_I,
		Jn,Jr,S,pressure,J,wt,0);

    /*** Rt (3 x 1) ***/
    UL_Ru_at_ip(Rt,nne,total_nne,ST,Fn,Fr,Fr_I,
		Jn,Jr,S,pressure,J,wt,1);

    /* compute incremental pressure */
    pressure = 0;
    for(int i=0; i<nne; i++){
      pressure += Na[i]*(p[i]);
    }

    /*** Rp (4 x 1) ***/
    UL_Rp_at_ip(Rp,nne,Na,kappa,Up,pressure,J,wt);

    }/* end integration */

    /* free un-needed */
    free(int_pt_ksi);
    free(int_pt_eta);
    free(int_pt_zet);
    free(weights);
    free(Na);
    free(N_x);
    free(N_y);
    free(N_z);
    dealoc4(ST_tensor,3,3,ndn);
    free(ST);
    dealoc2(Fn_mat,2);
    free(Fn);
    dealoc2(Fr_mat,3);
    free(Fr);
    free(Fr_I);
    free(AA);
    dealoc2(C_mat,3);
    free(C);
    dealoc2(S_mat,3);
    free(S);
    free(disp);
    free(p);
    /*=== END INTEGRATION ===*/

    int sign = 1;
    /*** Scatter the residuals (u,v,w,p) ***/
    /* NOTE residuals are negative because subtracted from RHS later */
    for (int a=0; a<nne; a++){
      for (int b=0; b<ndofn; b++){
	if (b<ndn){ /* from Ru */
	  fe[a*ndofn+b] = sign*Ru[a*ndn+b];
	} else if (b == ndn){ /* from Rp */
	  fe[a*ndofn+b] = sign*Rp[a];
	} else {
	  PGFEM_printerr("ERROR constructing residual vector!\n");
	  PGFEM_Abort();
	  abort();
	}
      }
    }

    /* Compute norm of LOCAL residuals */
    nrB += cblas_ddot(n_bub*n_bub_dofs,Rt,1,Rt,1);
    /*============== END COMPUTE ELEMENT RESIDUALS ============*/

   /* Assemble global residual for U and P, but separate into
       separate arrays */
    int j = 0;
    for (int k=0;k<nne;k++){
      for (int kk=0;kk<node[nod[k]].ndofn;kk++){
	int II = node[nod[k]].id[kk]-1;
	if (II < 0) continue;
	if (kk<node[nod[k]].ndofn-1){
	  f_u[II] += fe[j+kk];
	} else {
	  f_p[II] += fe[j+kk];
	}
      }/*end kk*/
      j += node[nod[k]].ndofn;
    }/*end k*/

    free(r_e);
    free(x);
    free(y);
    free(z);
    free(fe);
    free(cn);

  }/* end ii < ne*/

  /* Compute norms */
  double normal, tmp;
  double *BS_f;
  BS_f = aloc1(DomDof[myrank]);

  /*=== Norm of Ru ===*/
  for (int i=0;i<ndofd;i++) f[i] = RR[i] - f_u[i];
  LToG(f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
  normal = cblas_ddot(DomDof[myrank],BS_f,1,BS_f,1);
  MPI_Allreduce(&normal,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Ru = %12.12e\n",normal);
  memset(BS_f,0,DomDof[myrank]*sizeof(double));

  /*=== Norm of Rp ===*/
  LToG(f_p,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
  normal = cblas_ddot(DomDof[myrank],BS_f,1,BS_f,1);
  MPI_Allreduce(&normal,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Rp = %12.12e\n",normal);

  /*=== Norm of Rt ===*/
  MPI_Allreduce(&nrB,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Rt = %12.12e\n",normal);

  free(f_u);
  free(f_p);
  free(f);

}/* MINI_chec_resid */
