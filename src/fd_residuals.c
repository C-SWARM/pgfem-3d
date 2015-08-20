#include "fd_residuals.h"
#include "assert.h"
#include "enumerations.h"
#include "utils.h"
#include "allocation.h"
#include "resid_on_elem.h"
#include "stabilized.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "matice.h"
#include "elem3d.h"
#include <math.h>
#include <string.h>
#include "three_field_element.h"
#include "condense.h"
#include "new_potentials.h"
#include "tensors.h"
#include "cast_macros.h"
#include "index_macros.h"
#include "def_grad.h"
#include "mkl_cblas.h"
#include "dynamics.h"

#include "constitutive_model.h"
#include "plasticity_model.h"

#define ndn 3
			   
int update_residuals_from_constitutive_model(double *f,
        const int ii,
        const int ndofn,
        const int nne,
        const int nsd,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        MATGEOM matgeom,
        const long *nod,
        const NODE *node,
        double dt,
        SIG *sig,
        EPS *eps,
        const SUPP sup,
        double *r_e)
{
  int total_Lagrangian = 0;
  
  int err = 0;

  if (PGFEM3D_DEV_TEST) {
    //assert(0);
  }
    
  double *u;
  u = aloc1(nne*nsd);

  for(int a=0;a<nne;a++)
  {
  	for(int b=0; b<nsd;b++)
  		u[a*nsd+b] = r_e[a*ndofn+b];	
  }

  Matrix(double) Fn, Fr, Fnp1, pFn;
  Matrix_construct_redim(double,Fn ,3,3);
  Matrix_construct_redim(double,Fr ,3,3);
  Matrix_construct_redim(double,Fnp1 ,3,3); 
  
  Matrix_construct_redim(double,pFn ,3,3);     
       
  Matrix(double) eFnp1,pFnp1,pFnp1_I,L,S;  
  Matrix_construct_redim(double,eFnp1,3,3);  
  Matrix_construct_redim(double,pFnp1,3,3);
  Matrix_construct_redim(double,pFnp1_I,3,3);      
  Matrix_construct_redim(double,L ,81,1);  
  Matrix_construct_redim(double,S    ,3,3);
  
  Matrix(double) pFnI, eFn, M, eFnM;  
  Matrix_construct_redim(double,pFnI,3,3); 
  Matrix_construct_redim(double, eFn,3,3); 
  Matrix_construct_redim(double,   M,3,3);
  Matrix_construct_redim(double,eFnM,3,3);
            
  Matrix(double) ST_ab, AA, sAA;
  Matrix_construct(double,ST_ab);
  Matrix_construct_redim(double,AA,3,3);
  Matrix_construct_redim(double,sAA, 3,3);

  Matrix(double) MTeFnT_sAA, MTeFnT_sAA_eFnM;
  Matrix_construct_redim(double,MTeFnT_sAA        ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFnM    ,3,3);

  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);      
  int compute_stiffness = 0;
   
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,Fr);
    
    Constitutive_model *m = &(eps[ii].model[ip-1]);
    Matrix(double) *Fs = (m->vars).Fs;    


//    if(updated_Lagrangian)
//    {
      Matrix_AeqB(Fn,1.0,Fs[TENSOR_Fn]);
      Matrix_AeqB(pFn,1.0,Fs[TENSOR_pFn]);      
//    }   
//    else
//    {  
//      Matrix_eye(Fn,3);
//      Matrix_eye(pFn,3);
//    } 
    

    Matrix_inv(pFn, pFnI);
    Matrix_AxB(eFn,1.0,0.0,Fn,0,pFnI,0); 
   
    // --> update plasticity part
    if(total_Lagrangian)
    {
      Matrix(double) FnI;
      Matrix_construct_redim(double, FnI,3,3);
      Matrix_inv(Fn,FnI);
      Matrix_AeqB(Fnp1,1.0,Fr);  // Fn+1 
      Matrix_AxB(Fr,1.0,0.0,Fnp1,0,FnI,0);  // Fn+1          
      Matrix_cleanup(FnI);          
    }
    else
    {
      Matrix_AxB(Fnp1,1.0,0.0,Fr,0,Fn,0);  // Fn+1    
    }      
    constitutive_model_update_plasticity(&pFnp1,&Fnp1,&eFn,m,dt);
    Matrix_AxB(M,1.0,0.0,pFnI,0,pFnp1,0);    
    // <-- update plasticity part

    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(S,0.0);    
    
    Matrix_inv(pFnp1, pFnp1_I);
    Matrix_AxB(eFnp1,1.0,0.0,Fnp1,0,pFnp1_I,0);
    constitutive_model_update_elasticity(m,&eFnp1,dt,&L,&S,compute_stiffness);
    // <-- update elasticity part
            
    Matrix_AxB(eFnM,1.0,0.0,eFn,0,M,0);
    double Jn; Matrix_det(Fn, Jn);
    
    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];
        Matrix_init_w_array(ST_ab,3,3,ptrST_ab);
        Matrix_AxB(AA,1.0,0.0,Fr,1,ST_ab,0); 
        Matrix_symmetric(AA,sAA);

        Matrix_AxB(MTeFnT_sAA,1.0,0.0,eFnM,1,sAA,0);
        Matrix_AxB(MTeFnT_sAA_eFnM,1.0,0.0,MTeFnT_sAA,0,eFnM,0); 
        double MTeFnT_sAA_eFnM_S = 0.0; 
        Matrix_ddot(MTeFnT_sAA_eFnM,S,MTeFnT_sAA_eFnM_S);          
        
        int fe_id = a*ndofn + b;              
        f[fe_id] += 1.0/Jn*fe.detJxW*MTeFnT_sAA_eFnM_S;
      }
    }       
  }
  
  free(u);
  Matrix_cleanup(Fn);
  Matrix_cleanup(Fr);
  Matrix_cleanup(Fnp1);
  Matrix_cleanup(pFn);

  Matrix_cleanup(eFnp1);    
  Matrix_cleanup(pFnp1);
  Matrix_cleanup(pFnp1_I);        
  Matrix_cleanup(L);  
  Matrix_cleanup(S);

  Matrix_cleanup(pFnI);
  Matrix_cleanup(eFn);  
  Matrix_cleanup(M);
  Matrix_cleanup(eFnM);  
  
  Matrix_cleanup(ST_ab);
  Matrix_cleanup(AA);
  Matrix_cleanup(sAA);
  Matrix_cleanup(MTeFnT_sAA);
  Matrix_cleanup(MTeFnT_sAA_eFnM);
    
  FEMLIB_destruct(&fe);
  return err;
}        
        

int fd_residuals (double *f_u,
		  long ne,
		  int n_be,
		  long ndofn,
		  long npres,
		  double *d_r,
		  double *r,
		  NODE *node,
		  ELEMENT *elem,
		  BOUNDING_ELEMENT *b_elems,
		  MATGEOM matgeom,
		  HOMMAT *hommat,
		  SUPP sup,
		  EPS *eps,
		  SIG *sig,
		  double nor_min,
		  CRPL *crpl,
		  double dt,
		  double t,
		  double stab,
		  long nce,
		  COEL *coel,
		  MPI_Comm mpi_comm,
		  const PGFem3D_opt *opts,
		  double alpha, double *r_n, double *r_n_1)
{
/* make decision to include ineria*/
  const int mat = elem[0].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;
  
  if(fabs(rho)<MIN_DENSITY)
    include_inertia = 0;
   
/* decision end*/   
  int err = 0;
  /* long i,j,nne,ndofe,ndofc,k,kk,II,*nod,P,R,*cn; */
  /* double *r_e,*x,*y,*z,*fe,*X,*Y; */

  /* int nne_t; */

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);
  
  for (int i=0;i<ne;i++){
    /* Number of element nodes */
    int nne = elem[i].toe;
    int nne_t = nne + elem[i].n_bub;
    /* Nodes on element */
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node);

    /* allocation */
    double *r_e = aloc1 (ndofe);

    double *x,*y,*z;
    if(opts->analysis_type == MINI
       || opts->analysis_type == MINI_3F){/* P1+B/P1 */
      x = aloc1 (nne_t);
      y = aloc1 (nne_t);
      z = aloc1 (nne_t);
    } else {
      x = aloc1 (nne);
      y = aloc1 (nne);
      z = aloc1 (nne);
    }

    double *fe = aloc1 (ndofe);
    long *cn = aloc1l (ndofe);
    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    /* coordinates */
    if(sup->multi_scale)
    {
	    nodecoord_total(nne,nod,node,x,y,z);
	    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
    } 
    else
    {
      switch(opts->analysis_type)
      {
      case DISP:
        if(PGFEM3D_DEV_TEST)
        {  
	        nodecoord_updated(nne,nod,node,x,y,z);
	        def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);
	      } 
        else
        {  
          nodecoord_total(nne,nod,node,x,y,z);
	        def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
	      }  
        break; 
      case TF:
	      nodecoord_total(nne,nod,node,x,y,z);
	      def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
	      break;
      default:
	      nodecoord_updated(nne,nod,node,x,y,z);
	      /* deformation on element */
	      def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);    
	      break;
      }
    }

    if(opts->analysis_type == MINI
       || opts->analysis_type == MINI_3F){/* P1+B/P1 */
      element_center(nne,x,y,z);
    }

		int nVol = N_VOL_TREE_FIELD;
    int nsd = 3;

    if(include_inertia)
    {
      residuals_w_inertia_el(fe,i,nne,ndofn,npres,nVol,ndofe,r_e,node,elem,hommat,sup,eps,sig,
		                         nod,cn,x,y,z,dt,t,opts,alpha,r_n,r_n_1);		
    }
    else
    {    
    /* Residuals on element */
    switch(opts->analysis_type){
    case STABILIZED:
      err = resid_st_elem (i,ndofn,nne,elem,nod,node,hommat,
			   x,y,z,eps,sig,sup,r_e,nor_min,fe,dt,stab);
      break;
    case MINI:
      err = MINI_resid_el(fe,i,ndofn,nne,x,y,z,elem,
			  nod,node,hommat,eps,sig,r_e);
      break;
    case MINI_3F:
      err = MINI_3f_resid_el(fe,i,ndofn,nne,x,y,z,elem,
			     nod,node,hommat,eps,sig,r_e);
      break;
    case DISP:
    {
      double *bf = aloc1(ndofe);
	    memset(bf, 0, sizeof(double)*ndofe);
      DISP_resid_body_force_el(bf,i,ndofn,nne,x,y,z,elem,hommat,node,dt,t);   

      if (PGFEM3D_DEV_TEST) {
        update_residuals_from_constitutive_model(fe,i,ndofn,nne,nsd,elem,hommat,matgeom,nod,node,
                                                 dt,sig,eps,sup,r_e);
      } else {
        err =  DISP_resid_el(fe,i,ndofn,nne,x,y,z,elem,
                             hommat,nod,node,eps,sig,sup,r_e);
      }
        			   
	    for(long a = 0; a<ndofe; a++)
	      fe[a] += -bf[a];			        

      dealoc1(bf);			   
      break;
    }  
    case TF:
    { 
	    double *bf = aloc1(ndofe);
	    memset(bf, 0, sizeof(double)*ndofe);
      DISP_resid_body_force_el(bf,i,ndofn,nne,x,y,z,elem,hommat,node,dt,t); 

//      residuals_3f_w_inertia_el(fe,i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,
//	                                dt,sig,eps,-1.0,r_e,r_e);


      residuals_3f_el(fe,i,ndofn,nne,npres,nVol,nsd,
			       x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,r_e);
	    for(long a = 0; a<ndofe; a++)
	      fe[a] += -bf[a];			        

      dealoc1(bf);
			break;
		}	
    default:
      err = resid_on_elem (i,ndofn,nne,nod,elem,node,matgeom,
			   hommat,x,y,z,eps,sig,r_e,npres,
			   nor_min,fe,crpl,dt,opts->analysis_type);

      break;
    }
  }

    /* Assembly */
    {
      int j = 0;
      int II = 0;
      for (int k=0;k<nne;k++){
	for (int kk=0;kk<node[nod[k]].ndofn;kk++){
	  II = node[nod[k]].id[kk]-1;
	  if (II < 0) continue;
	  f_u[II] += fe[j+kk];
	}/*end kk*/
	j += node[nod[k]].ndofn;
      }/*end k*/
    }

    dealoc1l (nod);
    dealoc1 (r_e);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1 (fe);
    dealoc1l (cn);

    /*** RETURN on error ***/
    if(err != 0) return err;

  }/* end i < ne*/

  /**** COHESIVE ELEMENT RESIDUALS ****/
  if (opts->cohesive == 1){
    int ndofc = 3;

    for (int i=0;i<nce;i++){
      
      int nne = coel[i].toe/2;
      int ndofe = coel[i].toe*ndofc;
      
      long *nod = aloc1l (coel[i].toe);
      double *r_e = aloc1 (ndofe);
      double *x = aloc1 (coel[i].toe);
      double *y = aloc1 (coel[i].toe);
      double *z = aloc1 (coel[i].toe);
      double *fe = aloc1 (ndofe);
      long *cn = aloc1l (ndofe);
      
      for (int j=0;j<coel[i].toe;j++)
	nod[j] = coel[i].nod[j];

      nodecoord_updated (coel[i].toe,nod,node,x,y,z);
      
      /* code numbers on element */
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nod,node,cn);
      
      /* deformation on element */
      def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);

      /* Residuals on element */
      resid_co_elem (i,ndofc,nne,nod,x,y,z,coel,r_e,fe,nor_min,myrank);
      
      /* Localization */
      {
	int II = 0;
	for (int k=0;k<coel[i].toe;k++){
	  for (int kk=0;kk<ndofc;kk++){
	    II = node[nod[k]].id[kk]-1;
	    if (II < 0)
	      continue;
	    f_u[II] += fe[k*ndofc+kk];
	  }/*end kk*/
	}/*end k*/
      }
      
      dealoc1l (nod);
      dealoc1 (r_e);
      dealoc1 (x);
      dealoc1 (y);
      dealoc1 (z);
      dealoc1 (fe);
      dealoc1l (cn);
    }/* end i < nce */
  }/* end coh == 1 */

  /*===============================================
    |             BOUNDARY ELEMENTS               |
    ===============================================*/
    
  for (int i=0; i<n_be; i++){

    /* get the coordinates and dof id's on the element nodes */
    const long *ptr_vnodes = elem[b_elems[i].vol_elem_id].nod; /* --"-- */
    const ELEMENT *ptr_elem = &elem[b_elems[i].vol_elem_id]; /* --"-- */
    const BOUNDING_ELEMENT *ptr_be = &b_elems[i];
    const int nne = ptr_elem->toe;

    double *x = aloc1(nne);
    double *y = aloc1(nne);
    double *z = aloc1(nne);

    switch(opts->analysis_type){
    case DISP:
      nodecoord_total(nne,ptr_vnodes,node,x,y,z);
      break;
    default:
      nodecoord_updated(nne,ptr_vnodes,node,x,y,z);
      break;
    }

    int ndofe = get_ndof_on_bnd_elem(node,ptr_be,elem);
    double *r_e = aloc1(ndofe);
    double *fe = aloc1(ndofe);
    long *cn = aloc1l(ndofe);
    long *Gcn = aloc1l(ndofe);

    get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn);
    get_dof_ids_on_bnd_elem(1,ndofn,node,ptr_be,elem,Gcn);

    /* compute the deformation on the element */
    def_elem(cn,ndofe,d_r,NULL,NULL,r_e,sup,0);

    /* TOTAL LAGRANGIAN formulation */
    if(opts->analysis_type == DISP){
      double *r_en = aloc1(ndofe);
      def_elem(cn,ndofe,r,NULL,NULL,r_en,sup,1);
      vvplus(r_e,r_en,ndofe);
      free(r_en);
    }

    /* for debugging */
    double *RR = aloc1(ndofe);
    def_elem(cn,ndofe,f_u,NULL,NULL,RR,sup,2);

    if(opts->analysis_type == DISP){
      err += DISP_resid_bnd_el(fe,i,ndofn,ndofe,x,y,z,b_elems,elem,
			       hommat,node,eps,sig,sup,r_e);
    } else {
      /* not implemented/needed so do nothing */
    }

    /* Assembly to local part of the residual vector */
    {
      for(int j=0; j<ndofe; j++){
	int II = cn[j] - 1;
	if(II >= 0){
	  f_u[II] += fe[j];
	}
      }
    }

    def_elem(cn,ndofe,f_u,NULL,NULL,RR,sup,2);

    free(x);
    free(y);
    free(z);

    free(r_e);
    free(fe);
    free(cn);
    free(Gcn);
    free(RR);

    /*** RETURN on error ***/
    if(err != 0) return err;
  } /* for each bounding element */
  return err;
}
