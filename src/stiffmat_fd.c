/* HEADER */
/**
 * AUTHORS:
 * Matthew Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * Karel Matous, University of Notre Dame, kmatous [at] nd.edu
 */
#include "stiffmat_fd.h"
#include "assert.h"
#include "enumerations.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "elem3d.h"
#include "allocation.h"
#include "PLoc_Sparse.h"
#include "stabilized.h"
#include "stiffmatel_fd.h"
#include "utils.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "matice.h"

#include "three_field_element.h"
#include "condense.h"
#include "new_potentials.h"
#include "tensors.h"
#include "cast_macros.h"
#include "index_macros.h"
#include "def_grad.h"

#include "mkl_cblas.h"
#include "femlib.h"
#include "dynamics.h"

#include "constitutive_model.h"
#include "plasticity_model.h"

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#define ndn 3

#ifndef PGFEM3D_DEV_TEST
#define PGFEM3D_DEV_TEST 0
#endif

static const int periodic = 0;

int update_stiffness_from_constitutive_model(double *lk,
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

  
  Matrix(double) Fr, Fnp1;
  Matrix(double) FreFn,eFnp1,pFnp1,pFnp1_I,L,dMdu,S;  
  Matrix(double) pFnI, eFn, M, eFnM;  
  Matrix(double) ST_ab, ST_wg, AA, BB, CC;
  Matrix(double) sAA, sBB, sCC;
  Matrix(double) MTeFnT_sAA, MTeFnT_sAA_eFn,MTeFnT_sAA_eFnM,FrTFr,MTeFnT_FrTFr,MTeFnT_FrTFreFn,MTeFnT_FrTFreFndMdu,dCdu,MTeFnT_sBB;
  Matrix(double) L_dCdu,MTeFnT_sCC,MTeFnT_sCC_eFnM,MTeFnT_sAA_eFndMdu,sMTeFnT_sAA_eFndMdu;
      
  Matrix_construct_redim(double,Fr ,3,3);
  Matrix_construct_redim(double,Fnp1 ,3,3);      
  Matrix_construct_redim(double,L ,81,1);  
  Matrix_construct_redim(double,dMdu ,3,3);
  Matrix_construct_redim(double,FreFn,3,3);      
  Matrix_construct_redim(double,eFnp1,3,3);
  Matrix_construct_redim(double,pFnp1,3,3);
  Matrix_construct_redim(double,pFnp1_I,3,3);    
  Matrix_construct_redim(double,S    ,3,3);     
  Matrix_construct_redim(double,pFnI,3,3); 
  Matrix_construct_redim(double, eFn,3,3); 
  Matrix_construct_redim(double,   M,3,3);
  Matrix_construct_redim(double,eFnM,3,3); 
  Matrix_construct(double,ST_ab);
  Matrix_construct(double,ST_wg);
  Matrix_construct_redim(double,AA,3,3);
  Matrix_construct_redim(double,BB,3,3);
  Matrix_construct_redim(double,CC,3,3);  
  Matrix_construct_redim(double,sAA, 3,3);
  Matrix_construct_redim(double,sBB, 3,3);  
  Matrix_construct_redim(double,sCC, 3,3);    
  Matrix_construct_redim(double,MTeFnT_sAA        ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFn     ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFnM    ,3,3);
  Matrix_construct_redim(double,FrTFr            ,3,3);
  Matrix_construct_redim(double,MTeFnT_FrTFr      ,3,3);
  Matrix_construct_redim(double,MTeFnT_FrTFreFn    ,3,3);
  Matrix_construct_redim(double,MTeFnT_FrTFreFndMdu,3,3);
  Matrix_construct_redim(double,dCdu             ,3,3);
  Matrix_construct_redim(double,MTeFnT_sBB        ,3,3);
  Matrix_construct_redim(double,L_dCdu           ,3,3);
  Matrix_construct_redim(double,MTeFnT_sCC        ,3,3);
  Matrix_construct_redim(double,MTeFnT_sCC_eFnM    ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFndMdu ,3,3);
  Matrix_construct_redim(double,sMTeFnT_sAA_eFndMdu,3,3);  

  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0);
  int compute_stiffness = 1;      

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,Fr);

    Constitutive_model *m = &(eps[ii].model[ip-1]);
    Matrix(double) *Fs = (m->vars).Fs;
    Matrix_AxB(FrTFr,1.0,0.0,Fr,1,Fr,0);
    
    Matrix_inv(Fs[TENSOR_pFn], pFnI);
    Matrix_AxB(eFn,1.0,0.0,Fs[TENSOR_Fn],0,pFnI,0); 
    
    // --> update plasticity part
    Matrix_AxB(Fnp1,1.0,0.0,Fs[TENSOR_Fn],0,Fr,0);  // Fn+1    
    
   constitutive_model_update_plasticity(&pFnp1,&Fnp1,&eFn,m,dt);
/////////////////////////////////////////////////////////////////////// 
    Matrix_AxB(M,1.0,0.0,pFnI,0,pFnp1,0);   
    // <-- update plasticity part
/////////////////////////////////////////////////////////////////////////////////            
/////////////////////////////////////////////////////////////////////////////////
            
//            Matrix_eye(M,3);
//            Matrix_eye(pFnp1,3);            
            
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////  

    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(S,0.0);    
    
    Matrix_inv(pFnp1, pFnp1_I);
    Matrix_AxB(eFnp1,1.0,0.0,Fnp1,0,pFnp1_I,0);
    constitutive_model_update_elasticity(m,&eFnp1,dt,&L,&S,compute_stiffness);
    // <-- update elasticity part

    // --> start computing tagent
    Matrix_AxB(FreFn,1.0,0.0,Fr,0,eFn,0);
    Matrix_AxB(eFnp1,1.0,0.0,FreFn,0,M,0);    
    Matrix_AxB(eFnM,1.0,0.0,eFn,0,M,0);
    
    double Jn; Matrix_det(Fs[TENSOR_Fn], Jn);
    
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
        Matrix_AxB(MTeFnT_sAA_eFn,1.0,0.0,MTeFnT_sAA,0,eFn,0);
        Matrix_AxB(MTeFnT_sAA_eFnM,1.0,0.0,MTeFnT_sAA_eFn,0,M,0);  

        for(int w=0; w<nne; w++)
        {
          for(int g=0; g<nsd; g++)
          { 
            const double* const ptrST_wg = &(fe.ST)[idx_4_gen(w,g,0,0,
                                                      nne,nsd,nsd,nsd)]; 
            Matrix_init_w_array(ST_wg,3,3,ptrST_wg); 
            
            // --> update stiffness w.r.t plasticity
            constitutive_model_update_dMdu(m,&dMdu,&eFnp1,&S,&L,&ST_wg,dt);
/////////////////////////////////////////////////////////////////////////////////            
/////////////////////////////////////////////////////////////////////////////////
//            Matrix_print(dMdu);
//            Matrix_init(dMdu,0.0);
            
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////            
            // <-- update stiffness w.r.t plasticity
            Matrix_AxB(BB,1.0,0.0,Fr,1,ST_wg,0); 
            Matrix_symmetric(BB,sBB);
            Matrix_AxB(CC, 1.0,0.0,ST_ab,1,ST_wg,0);
            Matrix_symmetric(CC,sCC);
            
            // compute dCdu
            Matrix_AxB(MTeFnT_FrTFr,1.0,0.0,eFnM,1,FrTFr,0);
            Matrix_AxB(MTeFnT_FrTFreFn,1.0,0.0,MTeFnT_FrTFr,0,eFn,0);
            Matrix_AxB(MTeFnT_FrTFreFndMdu,1.0,0.0,MTeFnT_FrTFreFn,0,dMdu,0);            
            Matrix_symmetric(MTeFnT_FrTFreFndMdu,dCdu);
            
            Matrix_AxB(MTeFnT_sBB,1.0,0.0,eFnM,1,sBB,0);
            Matrix_AxB(dCdu,1.0,1.0,MTeFnT_sBB,0,eFnM,0);
            
            // compute MTeFnT_sAA_eFnM:L:dCdu
            Matrix_Tns4_dd_Tns2(L_dCdu,L,dCdu);
            double MTeFnT_sAA_eFnM_L_dCdu = 0.0;
            Matrix_ddot(MTeFnT_sAA_eFnM,L_dCdu,MTeFnT_sAA_eFnM_L_dCdu);
            
            // compute MTeFnT_sCC_eFnM
            Matrix_AxB(MTeFnT_sCC,1.0,0.0,eFnM,1,sCC,0);
            Matrix_AxB(MTeFnT_sCC_eFnM,1.0,0.0,MTeFnT_sCC,0,eFnM,0);
            
            // compute MTeFnT_sCC_eFnM:S
            double MTeFnT_sCC_eFnM_S = 0.0;
            Matrix_ddot(MTeFnT_sCC_eFnM,S,MTeFnT_sCC_eFnM_S);
            
            // compute MTeFnT_sAA_eFndMdu
            Matrix_AxB(MTeFnT_sAA_eFndMdu,1.0,0.0,MTeFnT_sAA_eFn,0,dMdu,0);    
            Matrix_symmetric(MTeFnT_sAA_eFndMdu, sMTeFnT_sAA_eFndMdu);        

            // compute MTeFnT_sAA_eFndMdu:S
            double sMTeFnT_sAA_eFndMdu_S = 0.0;            
            Matrix_ddot(sMTeFnT_sAA_eFndMdu,S,sMTeFnT_sAA_eFndMdu_S);
            
            const int lk_idx = idx_K(a,b,w,g,nne,nsd);  
                      
            lk[lk_idx] += 1.0/Jn*fe.detJxW*(MTeFnT_sAA_eFnM_L_dCdu + 2.0*sMTeFnT_sAA_eFndMdu_S + MTeFnT_sCC_eFnM_S);            
          }
        }
      }
    }
  }
  free(u);
  
  Matrix_cleanup(Fr);
  Matrix_cleanup(Fnp1);      
  Matrix_cleanup(L);  
  Matrix_cleanup(dMdu);
  Matrix_cleanup(FreFn);      
  Matrix_cleanup(eFnp1);
  Matrix_cleanup(pFnp1);
  Matrix_cleanup(pFnp1_I);
  Matrix_cleanup(S);    
  Matrix_cleanup(pFnI); 
  Matrix_cleanup(eFn); 
  Matrix_cleanup(M);
  Matrix_cleanup(eFnM); 
  Matrix_cleanup(ST_ab);
  Matrix_cleanup(ST_wg);
  Matrix_cleanup(AA);
  Matrix_cleanup(BB);
  Matrix_cleanup(CC);  
  Matrix_cleanup(sAA);
  Matrix_cleanup(sBB);  
  Matrix_cleanup(sCC);    
  Matrix_cleanup(MTeFnT_sAA);
  Matrix_cleanup(MTeFnT_sAA_eFn);
  Matrix_cleanup(MTeFnT_sAA_eFnM);
  Matrix_cleanup(FrTFr);
  Matrix_cleanup(MTeFnT_FrTFr);
  Matrix_cleanup(MTeFnT_FrTFreFn);
  Matrix_cleanup(MTeFnT_FrTFreFndMdu);
  Matrix_cleanup(dCdu);
  Matrix_cleanup(MTeFnT_sBB);
  Matrix_cleanup(L_dCdu);
  Matrix_cleanup(MTeFnT_sCC);
  Matrix_cleanup(MTeFnT_sCC_eFnM);
  Matrix_cleanup(MTeFnT_sAA_eFndMdu);
  Matrix_cleanup(sMTeFnT_sAA_eFndMdu);
    
  FEMLIB_destruct(&fe);
  
  return err;
}        

int el_compute_stiffmat(int i, double *lk, long ndofn, long nne, long npres, int nVol, int nsd,
                        ELEMENT *elem, NODE *node, HOMMAT *hommat, MATGEOM matgeom, SIG *sig, EPS *eps, SUPP sup,
                        double dt, double nor_min, double stab, CRPL *crpl, long FNR, double lm,
	                      double *x, double *y, double *z, double *fe, long *nod, double *r_n, double *r_e, 
	                      double alpha, int include_inertia, const int analysis)
{
  int err = 0;
  if(include_inertia)
  {
    stiffmat_disp_w_inertia_el(lk,i,ndofn,nne,npres,nVol,nsd,x, y, z,	
                               elem,hommat,nod,node,dt,
                               sig,eps,sup,analysis,alpha,r_n,r_e);          
  }
  else
  {
    switch(analysis){
    case STABILIZED:
      err += stiffmatel_st(i,ndofn,nne,x,y,z,elem,hommat,nod,node,sig,eps,
			  sup,r_e,npres,nor_min,lk,dt,stab,FNR,lm,fe);
      break;
    case MINI:
      err += MINI_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			    hommat,nod,node,eps,sig,r_e);
      break;
    case MINI_3F:
      err += MINI_3f_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			      hommat,nod,node,eps,sig,r_e);
      break;
    case DISP:
    {      
      if (PGFEM3D_DEV_TEST){
        err += update_stiffness_from_constitutive_model(lk,i,ndofn,nne,nsd,elem,hommat,matgeom,nod,node,
                                                 dt,sig,eps,sup,r_e);
      } else {
        err += DISP_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
                               hommat,nod,node,eps,sig,sup,r_e);
      }
      break;
    }
    case TF:
      stiffmat_3f_el(lk,i,ndofn,nne,npres,nVol,nsd,
                  x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,-1.0,r_e);
      break;
    default:
      err += stiffmatel_fd (i,ndofn,nne,nod,x,y,z,elem,matgeom,
			  hommat,node,sig,eps,r_e,npres,
			  nor_min,lk,dt,crpl,FNR,lm,fe,analysis);
      break;
    } // switch (analysis)
  } // if(include_inertia)
  
  return err;
}        
/* This function may not be used outside this file */
static int el_stiffmat(int i, /* Element ID */
			double **Lk,
			int *Ap,
			int *Ai,
			long ndofn,
			ELEMENT *elem,
			NODE *node,
			HOMMAT *hommat,
			MATGEOM matgeom,
			SIG *sig,
			EPS *eps,
			double *d_r,
			double *r,
			long npres,
			SUPP sup,
			long iter,
			double nor_min,
			double dt,
			CRPL *crpl,
			double stab,
			long FNR,
			double lm,
			double *f_u,
			int myrank,
			int nproc,
			long GDof,
			COMMUN comm,
			int *Ddof,
			int interior,
			const int analysis,
			PGFEM_HYPRE_solve_info *PGFEM_hypre,
			double alpha, double *r_n, double *r_n_1)
{
  if(-1==i && 0==iter && PGFEM3D_DEV_TEST)
    constitutive_model_test(hommat);

/* make a decision to include ineria*/
  const int mat = elem[i].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;
  int nsd = 3;
  
  if(fabs(rho)<MIN_DENSITY)
  {
    include_inertia = 0;
  }
/* decision end*/ 
  
  
  int err = 0;
  long j,l,nne,ndofe,*cnL,*cnG,*nod,II;
  double *x,*y,*z,*r_e,*sup_def,*fe;
  long kk;
  
  /* Number of element nodes */
  nne = elem[i].toe;
    
  /* Nodes on element */
  nod = aloc1l (nne);
  elemnodes (i,nne,nod,elem);
    
  /* Element Dof */
  ndofe = get_ndof_on_elem_nodes(nne,nod,node);

  /* allocation */
  cnL = aloc1l (ndofe);
  cnG = aloc1l (ndofe);

  const int nne_t = nne + elem[i].n_bub;

  if (analysis == MINI
      || analysis == MINI_3F){ /* P1+B/P1 element */
    x = aloc1 (nne_t);
    y = aloc1 (nne_t);
    z = aloc1 (nne_t);
  } else {
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);
  }
  r_e = aloc1 (ndofe);
  fe = aloc1 (ndofe);
  if(sup->npd>0){
    sup_def = aloc1(sup->npd);
  } else {
    sup_def = NULL;
  }
    
  /* Coordinates of nodes */
  /*
   * The coordinates define the configuration for computing gradients. 
   */
  if(sup->multi_scale){
    /* multi-scale analysis is TOTAL LAGRANGIAN for all element
       formulations */
    nodecoord_total (nne,nod,node,x,y,z);
  } else {
    switch(analysis){
    case DISP: case TF:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
  }

  /* if P1+B/P1, get element centroid coords */
  if (analysis == MINI || analysis == MINI_3F){
    element_center(nne,x,y,z);
  }
    
  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cnL);
  get_dof_ids_on_elem_nodes(1,nne,ndofn,nod,node,cnG);
    
  /*=== deformation on element ===*/

  /* set the increment of applied def=0 on first iter */
  if (iter == 0) {
    for (j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }
  }

  /* get the deformation on the element */
  if(sup->multi_scale){
    /* multi-scale analysis is TOTAL LAGRANGIAN for all element
       formulations */
    def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
  } else {
    switch(analysis){
    case DISP: case TF: /* TOTAL LAGRANGIAN */
      def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
      break;
    default:
      def_elem (cnL,ndofe,d_r,elem,node,r_e,sup,0);
      break;
    }
  }

  /* recover thei increment of applied def on first iter */
  if (iter == 0){
    for (j=0;j<sup->npd;j++)
      sup->defl_d[j] = sup_def[j];
   }
    
  int nVol = N_VOL_TREE_FIELD;
  Matrix(double) lk;
  Matrix_construct_redim(double,lk,ndofe,ndofe);   
  Matrix_init(lk, 0.0);
  
  err += el_compute_stiffmat(i,lk.m_pdata,ndofn,nne,npres,nVol,nsd,
                        elem,node,hommat,matgeom,sig,eps,sup,
                        dt,nor_min,stab,crpl,FNR,lm,
	                      x,y,z,fe,nod,r_n,r_e,
	                      alpha,include_inertia,analysis);  

/*
  if(include_inertia)
  {
    stiffmat_disp_w_inertia_el(lk.m_pdata,i,ndofn,nne,npres,nVol,nsd,x, y, z,	
                               elem,hommat,nod,node,dt,
                               sig,eps,sup,analysis,alpha,r_n,r_e);          
  }
  else
  { 
           
  switch(analysis){
  case STABILIZED:
    err = stiffmatel_st(i,ndofn,nne,x,y,z,elem,hommat,nod,node,sig,eps,
			 sup,r_e,npres,nor_min,lk.m_pdata,dt,stab,FNR,lm,fe);
    break;
  case MINI:
    err = MINI_stiffmat_el(lk.m_pdata,i,ndofn,nne,x,y,z,elem,
			   hommat,nod,node,eps,sig,r_e);
    break;
  case MINI_3F:
    err = MINI_3f_stiffmat_el(lk.m_pdata,i,ndofn,nne,x,y,z,elem,
			      hommat,nod,node,eps,sig,r_e);
    break;
  case DISP:
    {
      
      if (PGFEM3D_DEV_TEST) {
        update_stiffness_from_constitutive_model(lk.m_pdata,i,ndofn,nne,nsd,elem,hommat,matgeom,nod,node,
                                                 dt,sig,eps,sup,r_e);

      } else {
        err = DISP_stiffmat_el(lk.m_pdata,i,ndofn,nne,x,y,z,elem,
                               hommat,nod,node,eps,sig,sup,r_e);
      }
    }
    break;
  case TF:
    stiffmat_3f_el(lk.m_pdata,i,ndofn,nne,npres,nVol,nsd,
                  x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,-1.0,r_e);
        //stiffmat_3f_el(lk.m_pdata,i,ndofn,nne,npres,nVol,nsd,
        //          x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,r_e);
        break;
  default:
    err = stiffmatel_fd (i,ndofn,nne,nod,x,y,z,elem,matgeom,
			 hommat,node,sig,eps,r_e,npres,
			 nor_min,lk.m_pdata,dt,crpl,FNR,lm,fe,analysis);
    break;
  } // switch (analysis)
  } // if(include_inertia)
    */
  if (PFEM_DEBUG){
    char filename[50];
    switch(analysis){
    case STABILIZED:
      sprintf(filename,"stab_stiff_%d.log",myrank);
      break;
    case MINI:
      sprintf(filename,"MINI_stiff_%d.log",myrank);
      break;
    case MINI_3F:
      sprintf(filename,"MINI_3f_stiff_%d.log",myrank);
      break;
    default:
      sprintf(filename,"stiff_%d.log",myrank);
      break;
    }

    FILE *output;
    output = fopen(filename,"a");
    print_array_d(output,lk.m_pdata,ndofe*ndofe,ndofe,ndofe);
    fclose(output);
  }

  /* Localization of TANGENTIAL LOAD VECTOR */
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    for (l=0;l<nne;l++){
      for (kk=0;kk<node[nod[l]].ndofn;kk++){
	II = node[nod[l]].id[kk]-1;
	if (II < 0)  continue;
	f_u[II] += fe[l*node[nod[l]].ndofn+kk];
      }/*end l */
    }/*end kk */
  }/* end periodic */

  /* Assembly */
  PLoc_Sparse (Lk,lk.m_pdata,Ai,Ap,cnL,cnG,ndofe,Ddof,GDof,
	       myrank,nproc,comm,interior,PGFEM_hypre,analysis);

  /*  dealocation  */
  free (cnL);
  free (cnG);
  free (nod);
  Matrix_cleanup(lk);
  free (x);
  free (y);
  free (z);
  free (fe);
  free (r_e);
  free (sup_def);

  return err;

} /* ELEMENT STIFFNESS */

/* This function may not be used outside of this file */
static void coel_stiffmat(int i, /* coel ID */
			  double **Lk,
			  int *Ap,
			  int *Ai,
			  long ndofc,
			  ELEMENT *elem,
			  NODE *node,
			  EPS *eps,
			  double *d_r,
			  double *r,
			  long npres,
			  SUPP sup,
			  long iter,
			  double nor_min,
			  double dt,
			  CRPL *crpl,
			  double stab,
			  COEL *coel,
			  long FNR,
			  double lm,
			  double *f_u,
			  int myrank,
			  int nproc,
			  long *DomDof,
			  long GDof,
			  COMMUN comm,
			  int *Ddof,
			  int interior,
			  const int analysis,
			  PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  long j,l,nne,ndofe,*cnL,*cnG,*nod,P,R,II;
  double *lk,*x,*y,*z,*r_e,*sup_def,*fe, *X, *Y;
  long kk;

  nne = coel[i].toe/2;
  ndofe = coel[i].toe*ndofc;
      
  /* Alocation */
  nod = aloc1l (coel[i].toe);
  lk = aloc1 (ndofe*ndofe);
  x = aloc1 (coel[i].toe); 
  y = aloc1 (coel[i].toe);
  z = aloc1 (coel[i].toe);
  r_e = aloc1 (ndofe);
  fe = aloc1 (ndofe);
  cnL = aloc1l (ndofe);
  cnG = aloc1l (ndofe);

  if(sup->npd > 0){
    sup_def = aloc1(sup->npd);
  } else {
    sup_def = NULL;
  }

  /* Element Node */
  for (j=0;j<coel[i].toe;j++)
    nod[j] = coel[i].nod[j];

  /* Coordinates */
  /* 
   * NOTE: get updated coordinates for all formulations, computes jump
   * based on both coordinates and deformation.
   */
  nodecoord_updated (coel[i].toe,nod,node,x,y,z);
      
  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nod,node,cnL);
  get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nod,node,cnG);

  /* deformation on element */
  if (iter == 0){
    for (j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }
  }

  def_elem (cnL,ndofe,d_r,elem,node,r_e,sup,0);
  if (iter == 0)
    for (j=0;j<sup->npd;j++)
      sup->defl_d[j] = sup_def[j];
      
  if (periodic == 1){/* Periodic */

    X = aloc1 (3);
    Y = aloc1 (3);
	
    for (j=0;j<coel[i].toe;j++){
      X[0] = x[j];
      X[1] = y[j];
      X[2] = z[j];
      for (P=0;P<3;P++){
	Y[P] = 0.0;
	for (R=0;R<3;R++){
	  Y[P] += eps[0].F[P][R]*X[R];
	}
      }
      for (P=0;P<3;P++){
	II = node[nod[j]].id[P];
	    
	if (P == 0)
	  x[j] = Y[P];
	if (P == 1)
	  y[j] = Y[P];
	if (P == 2)
	  z[j] = Y[P];
	    
	if (II > 0){
	  if (P == 0)
	    x[j] += r[II-1];
	  if (P == 1)
	    y[j] += r[II-1];
	  if (P == 2)
	    z[j] += r[II-1];
	}
      }/* P < 3 */
    }/* j < coel[i].toe */

    dealoc1 (X);
    dealoc1 (Y);
	
  }/* end periodic */

  /**************************/
  /*     FOR PERIODIC       */
  /*                        */
  /* xn+1 = Fn+1*X + un + u */
  /*    x = Fn+1*X + un     */
  /**************************/

  nulld (lk,ndofe*ndofe); 
  stiff_mat_coh (i,ndofc,nne,nod,x,y,z,coel,r_e,lk,
		 nor_min,eps,FNR,lm,fe,myrank);
      
  /* Assembly */
  PLoc_Sparse (Lk,lk,Ai,Ap,cnL,cnG,ndofe,Ddof,
	       GDof,myrank,nproc,comm,interior,PGFEM_hypre,analysis);

  /* Localization of TANGENTIAL LOAD VECTOR */
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    for (l=0;l<coel[i].toe;l++){
      for (kk=0;kk<ndofc;kk++){
	II = node[nod[l]].id[kk]-1;
	if (II < 0)  continue;
	f_u[II] += fe[l*ndofc+kk];
      }/*end l */
    }/*end kk */
  }/* end periodic */

  /*  dealocation  */
  dealoc1l (cnL);
  dealoc1l (cnG);
  dealoc1l (nod);
  dealoc1 (lk); 
  dealoc1 (x);
  dealoc1(y); 
  dealoc1 (z);
  dealoc1 (r_e);
  dealoc1 (fe);
  free(sup_def);

} /* COHESIVE ELEMENT STIFFNESS */

static int bnd_el_stiffmat(int belem_id,
			   double **Lk,
			   int *Ap,
			   int *Ai,
			   long ndofn,
			   ELEMENT *elem,
			   BOUNDING_ELEMENT *b_elems,
			   NODE *node,
			   HOMMAT *hommat,
			   MATGEOM matgeom,
			   SIG *sig,
			   EPS *eps,
			   double *d_r,
			   double *r,
			   long npres,
			   SUPP sup,
			   long iter,
			   double nor_min,
			   double dt,
			   CRPL *crpl,
			   double stab,
			   long FNR,
			   double lm,
			   double *f_u,
			   int myrank,
			   int nproc,
			   long GDof,
			   COMMUN comm,
			   int *Ddof,
			   int interior,
			   const int analysis,
			   PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  int err = 0;
  const BOUNDING_ELEMENT *ptr_be = &b_elems[belem_id];
  const ELEMENT *ptr_ve = &elem[ptr_be->vol_elem_id];
  const long *ptr_vnodes = ptr_ve->nod;
  const int nn_ve = ptr_ve->toe;

  /* get coordinated for BOUNDING ELEMENT */
  double *x = aloc1(nn_ve);
  double *y = aloc1(nn_ve);
  double *z = aloc1(nn_ve);
  switch(analysis){
  case DISP:
    nodecoord_total(nn_ve,ptr_vnodes,node,x,y,z);
    break;
  default:
    nodecoord_updated(nn_ve,ptr_vnodes,node,x,y,z);
    break;
  }

  /* get the local and global dof id's */
  int ndof_ve = get_ndof_on_bnd_elem(node,ptr_be,elem);

  long *cn_ve = aloc1l(ndof_ve);
  long *Gcn_ve = aloc1l(ndof_ve);

  get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn_ve);
  get_dof_ids_on_bnd_elem(1,ndofn,node,ptr_be,elem,Gcn_ve);

  /* compute the deformation on the element */
  double *v_disp = aloc1(ndof_ve);

  if(iter == 0){
    /* on iter == 0, null increment of deflection */
    double *sup_def = aloc1(sup->npd);
    for (int j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }

    def_elem(cn_ve,ndof_ve,d_r,NULL,NULL,v_disp,sup,0);

    for (int j=0;j<sup->npd;j++){
      sup->defl_d[j] = sup_def[j];
    }
    free(sup_def);
  } else {
    def_elem(cn_ve,ndof_ve,d_r,NULL,NULL,v_disp,sup,0);
  }

  if(analysis == DISP){ /* TOTAL LAGRANGIAN formulation */
    double *ve_n = aloc1(ndof_ve);
    def_elem(cn_ve,ndof_ve,r,NULL,NULL,ve_n,sup,1);
    vvplus(v_disp,ve_n,ndof_ve);
    free(ve_n);
  }

  /* compute the local stiffness matrix */
  double *lk = aloc1(ndof_ve*ndof_ve);
  if(analysis == DISP){
    err += DISP_stiffmat_bnd_el(lk,belem_id,ndofn,ndof_ve,
				x,y,z,b_elems,elem,hommat,node,eps,
				sig,sup,v_disp);
  } else {
    /* Not implemented, do nothing */
  }

  /* only assemble to global stiffness if no error */
  if(err == 0){
    /* PLoc_Sparse_rec(Lk,lk,Ai,Ap,Gcn_be,Gcn_ve,ndof_be,ndof_ve,Ddof, */
    /* 		   GDof,myrank,nproc,comm,interior); */
    PLoc_Sparse(Lk,lk,Ai,Ap,cn_ve,Gcn_ve,ndof_ve,Ddof,
		GDof,myrank,nproc,comm,interior,PGFEM_hypre,analysis);
  }


  free(x);
  free(y);
  free(z);

  free(cn_ve);
  free(Gcn_ve);

  free(v_disp);
  free(lk);

  return err;
} /* Bounding element stiffnes matrix */

/* This is the re-written function which computes elem stiffness on
   boundaries first, then interior elem stiffnesses before
   assembly. */
int stiffmat_fd(int *Ap,
		 int *Ai,
		 long ne,
		 int n_be,
		 long ndofn,
		 ELEMENT *elem,
		 BOUNDING_ELEMENT *b_elems,
		 long nbndel,
		 long *bndel,
		 NODE *node,
		 HOMMAT *hommat,
		 MATGEOM matgeom,
		 SIG *sig,
		 EPS *eps,
		 double *d_r,
		 double *r,
		 long npres,
		 SUPP sup,
		 long iter,
		 double nor_min,
		 double dt,
		 CRPL *crpl,
		 double stab,
		 long nce,
		 COEL *coel,
		 long FNR,
		 double lm,
		 double *f_u,
		 int myrank,
		 int nproc,
		 long *DomDof,
		 long GDof,
		 COMMUN comm,
		 MPI_Comm mpi_comm,
		 PGFEM_HYPRE_solve_info *PGFEM_hypre,
		 const PGFem3D_opt *opts,double alpha, double *r_n, double *r_n_1)
{
  int err = 0;
  long i,ndofc;
  int *Ddof;
  double **Lk,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;

  /* interior element counters */
  int idx = 0;
  int skip = 0;

  if(opts->solverpackage == 0){
    PGFEM_printf("BlockSolve no longer supported\n");
    abort();
  }

  err += init_and_post_stiffmat_comm(&Lk,&recieve,&req_r,&sta_r,
				     mpi_comm,comm);

  /* Lk = (double**) PGFEM_calloc (nproc,sizeof(double*)); */
  /* for (i=0;i<nproc;i++) { */
  /*   if (myrank == i || comm->S[i] == 0) */
  /*     k = 1; */
  /*   else  */
  /*     k = comm->AS[i]; */
  /*   Lk[i] = (double*) PGFEM_calloc (k,sizeof(double)); */
  /* } */
  /* if (Lk == NULL){ */
  /*   PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__); */
  /*   fflush(stdout);  */
  /*   PGFEM_Comm_code_abort (mpi_comm,i); */
  /* } */
  
  /* /\* Allocate recieve *\/ */
  /* recieve = (double**) PGFEM_calloc (nproc,sizeof(double*)); */
  /* for (i=0;i<nproc;i++) { */
  /*   if (comm->AR[i] == 0) */
  /*     KK = 1; */
  /*   else */
  /*     KK = comm->AR[i]; */
  /*   recieve[i] = (double*) PGFEM_calloc (KK,sizeof(double)); */
  /* } */
  
  /* /\* Allocate request fields *\/ */
  /* if (comm->Nr == 0) */
  /*   KK = 1; */
  /* else */
  /*   KK = comm->Nr; */
  /* sta_r = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status)); */
  /* req_r = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request)); */
  
  /* /\* Receive data *\/ */
  /* for (i=0;i<comm->Nr;i++){ */
  /*   KK = comm->Nrr[i]; */
  /*   MPI_Irecv (recieve[KK],comm->AR[KK],MPI_DOUBLE,KK, */
  /* 	       MPI_ANY_TAG,mpi_comm,&req_r[i]); */
  /* }/\* end i < nproc *\/ */
  
  /* Allocate */
  Ddof = aloc1i (nproc);
  
  /* Set Ddof */
  Ddof[0] = DomDof[0];
  for (i=1;i<nproc;i++)
    Ddof[i] = Ddof[i-1] + DomDof[i];
  
  /***** COMM BOUNDARY ELEMENTS *****/
  for(i=0; i<nbndel; i++){
    err += el_stiffmat(bndel[i],Lk,Ap,Ai,ndofn,elem,node,hommat,
		       matgeom,sig,eps,d_r,r,npres,sup,iter,nor_min,
		       dt,crpl,stab,FNR,lm,f_u,myrank,nproc,GDof,comm,
		       Ddof,0,opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);

    /* If there is an error, complete communication and exit */
    if(err != 0) goto send;
  }

  /***** COHESIVE ELEMENTS *****/

  /* Need to split into boundary and interior parts as with the
     regular elements */
  if (opts->cohesive == 1){/* COHESIVE STIFFNESS */
    ndofc = 3; 
    
    /* WHY IS THIS HERE */
    if (iter == 0) FNR = 0;
    
    if (nor_min < 1.e-10)
      nor_min = 1.e-10;
    
    for (i=0;i<nce;i++){
      coel_stiffmat(i,Lk,Ap,Ai,ndofc,elem,node,eps,
		    d_r,r,npres,sup,iter,nor_min,dt,crpl,
		    stab,coel,FNR,lm,f_u,myrank,nproc,DomDof,
		    GDof,comm,Ddof,0,opts->analysis_type,PGFEM_hypre);
    }
  }


  /**** BOUNDING ELEMENTS ******/
  /* In the future, these elements will be listed as with the
     volumetric elements to properly overlay computation and
     communication. For now, this is the best place for them as they
     are a proportinally smaller group than the interior volume
     elements for typical problems and are guarenteed to be on the
     communication boundary for periodic domains. */

  /* temporary for compile testing */
  /* int n_be = 0; */
  /* BOUNDING_ELEMENT *b_elems = NULL; */
  for(i=0; i<n_be; i++){
    err += bnd_el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,b_elems,node,hommat,
			   matgeom,sig,eps,d_r,r,npres,sup,iter,nor_min,
			   dt,crpl,stab,FNR,lm,f_u,myrank,nproc,GDof,
			   comm,Ddof,0,opts->analysis_type,PGFEM_hypre);

    /* If there is an error, complete communication and exit */
    if(err != 0) goto send;
  }

  if (PFEM_DEBUG){
    char ofile[50];
    switch(opts->analysis_type){
    case STABILIZED:
      sprintf(ofile,"stab_el_stiff_send_%d.log",myrank);
      break;
    case MINI:
      sprintf(ofile,"MINI_el_stiff_send_%d.log",myrank);
      break;
    case MINI_3F:
      sprintf(ofile,"MINI_3f_el_stiff_send_%d.log",myrank);
      break;
    default:
      sprintf(ofile,"el_stiff_send_%d.log",myrank);
      break;
    }

    FILE *out;
    out = fopen(ofile,"a");
    for(int send_proc=0; send_proc<nproc; send_proc++){
      if(send_proc==myrank || comm->S[send_proc] ==0) continue;
      for(int n_data=0; n_data<comm->AS[send_proc]; n_data++){
	PGFEM_fprintf(out,"%12.12e ",Lk[send_proc][n_data]);
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n");
    fclose(out);
  }

  /**********************************/
  /***** SEND BOUNDARY AND COEL *****/
  /**********************************/
 send:
  err += send_stiffmat_comm(&sta_s,&req_s,Lk,mpi_comm,comm);

  /* If error, complete communication and exit */
  if(err != 0) goto wait;

   /**********************************/
  /*****    CONTINUE WORKING    *****/
  /**********************************/

  /***** INTERIOR ELEMENTS *****/

  if(nbndel > 0){/* this is 99% of the time */
    for(i=0; i<ne; i++){
      if(idx < nbndel-1){
	if(i == 0 && idx == 0 && bndel[idx] == 0){
	  idx++;
	  skip++;
	  continue;
	} else if(i == bndel[idx]){
	  idx++;
	  skip++;
	  continue;
	} else if (idx == 0 && i < bndel[idx]){
	  err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,
			    sig,eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,
			    stab,FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
	} else if (idx > 0 && bndel[idx-1] < i && i < bndel[idx]){
	  err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			    eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,
			    FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
	} else {
	  PGFEM_printf("[%d]ERROR: problem in determining if element %ld"
		 " is on interior.\n", myrank, i);
	}
      } else {
	if(i != bndel[nbndel-1]){
	  err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			    eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,
			    FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
	}
      }

      /* If there is an error, complete communication and exit */
      if(err != 0) goto wait;
    }

    /* Check to make sure I got all of them */
    if(skip != nbndel - 1){
      PGFEM_printf("[%d]WARNING: number skipped elem != nbndel, check code\n",myrank);
    }
  } else { /* communication by coheisve elements only, nbndel = 0 */
    for(i=0; i<ne; i++){
      err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,FNR,
			lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
      /* If there is an error, complete communication and exit */
      if(err != 0) goto wait;
    }
  }

  /**********************************/
  /*****    RECEIVE AND ADD     *****/
  /**********************************/
  
 wait:
  err += assemble_nonlocal_stiffmat(comm,sta_r,req_r,PGFEM_hypre,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,comm);
    
 /* exit_function: */
  /* Deallocate recieve */
  for (i=0;i<nproc;i++)
    free (recieve[i]);
  free (recieve);
  
  /*  dealocation  */
  for (i=0;i<nproc;i++)
    free(Lk[i]);
  
  free (Lk);
  free (Ddof);
  free (sta_s);
  free (sta_r);
  free (req_s);
  free (req_r);

  return err;
}


/** Assemble non-local parts as they arrive */
int assemble_nonlocal_stiffmat(const COMMUN pgfem_comm,
			       MPI_Status *sta_r,
			       MPI_Request *req_r,
			       PGFEM_HYPRE_solve_info *PGFEM_hypre,
			       double **recv)
{
  int err = 0;
  int comm_idx = 0;
  int n_received = 0;
  while (n_received < pgfem_comm->Nr){
    /* get the communication index */
    err += MPI_Waitany(pgfem_comm->Nr,req_r,&comm_idx,sta_r);

    /* convert communication index to proc_id */
    const int proc = pgfem_comm->Nrr[comm_idx];

    /* get number of rows */
    const int nrows = pgfem_comm->R[proc];

    /* allocate rows and cols to receive */
    int *row_idx = PGFEM_calloc(nrows,sizeof(int));
    int *ncols = PGFEM_calloc(nrows,sizeof(int));
    int *col_idx = PGFEM_calloc(pgfem_comm->AR[proc],sizeof(int));

    /* get row and column ids */
    int idx = 0;
    for(int j=0; j<pgfem_comm->R[proc]; j++){
      row_idx[j] = pgfem_comm->RGID[proc][j];
      ncols[j] = pgfem_comm->RAp[proc][j];
      for(int k=0; k<ncols[j]; k++){
	col_idx[idx] = pgfem_comm->RGRId[proc][idx];
	++idx;
      }
    }

    /* assemble to local part of global stiffness */
    err += HYPRE_IJMatrixAddToValues(PGFEM_hypre->hypre_k,
				     nrows,ncols,row_idx,col_idx,
				     recv[proc]);

    /* free memory */
    free(row_idx);
    free(ncols);
    free(col_idx);

    /* increment counter */
    n_received++;
  }
  return err;
}
