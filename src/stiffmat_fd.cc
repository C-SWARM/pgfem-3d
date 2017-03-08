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
#include "PGFem3D_data_structure.h"
#include "macro_micro_functions.h"

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#define ndn 3

static const int periodic = 0;

int el_compute_stiffmat(int i,
                        double *lk,
                        long ndofn,
                        long nne,
                        long npres,
                        int nVol,
                        int nsd,
                        ELEMENT *elem,
                        NODE *node,
                        HOMMAT *hommat,
                        MATGEOM matgeom,
                        SIG *sig,
                        EPS *eps,
                        SUPP sup,
                        double dt,
                        double nor_min,
                        double stab,
                        CRPL *crpl,
                        long FNR,
                        double lm,
                        double *x,
                        double *y,
                        double *z,
                        double *fe,
                        long *nod,
                        double *r_n,
                        double *r_e,
                        double alpha,
                        int include_inertia,
                        const int analysis,
                        int cm)
{
  int err = 0;
  if(include_inertia)
  {
    stiffmat_disp_w_inertia_el(lk,i,ndofn,nne,npres,nVol,nsd,x, y, z,	
                               elem,hommat,nod,node,dt,
                               sig,eps,sup,analysis,cm,alpha,r_n,r_e);          
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
        err += DISP_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
                                hommat,nod,node,eps,sig,sup,r_e,dt);
      break;
    case TF:
      stiffmat_3f_el(lk,i,ndofn,nne,npres,nVol,nsd,
                  x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,-1.0,r_e);
      break;
    case CM:
      {
        switch(cm) {
        case UPDATED_LAGRANGIAN: // intentionally left to flow
        case TOTAL_LAGRANGIAN:          
          err += stiffness_el_crystal_plasticity(lk,i,ndofn,nne,nsd,elem,
                                                 nod,node,dt,eps,sup,r_e,
                                                 cm);
          break;

        case MIXED_ANALYSIS_MODE:
          err += stiffness_el_crystal_plasticity(lk,i,ndofn,nne,nsd,elem,
                                                 nod,node,dt,eps,sup,r_e,
                                                 1 /* TL */);
          break;

        default: assert(0 && "should never get here"); break;
        }
        break;
      }
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
			const int cm,
			PGFEM_HYPRE_solve_info *PGFEM_hypre,
			double alpha, double *r_n, double *r_n_1,
			const int mp_id)
{
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
  ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

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
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    case TF:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    case CM:
    {  
      switch(cm) {
      case UPDATED_LAGRANGIAN:
        nodecoord_updated (nne,nod,node,x,y,z);        
        break;        
      case TOTAL_LAGRANGIAN:
        nodecoord_total (nne,nod,node,x,y,z);
        break;
      case MIXED_ANALYSIS_MODE: // Total_lagrangian
        nodecoord_total (nne,nod,node,x,y,z);
        break;        
      default: assert(0 && "undefined CM type"); break;
      }
      break;      
    }    
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
  get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cnL,mp_id);
  get_dof_ids_on_elem_nodes(1,nne,ndofn,nod,node,cnG,mp_id);
    
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
    case DISP:
      def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
    case TF: /* TOTAL LAGRANGIAN */
      def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
      break;
    case CM:
    {  
      switch(cm) {
      case UPDATED_LAGRANGIAN:
        def_elem (cnL,ndofe,d_r,elem,node,r_e,sup,0);
        break;                
      case TOTAL_LAGRANGIAN:
        def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
        break;
      case MIXED_ANALYSIS_MODE:
        def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
        break;        

      default: assert(0 && "undefined CM type"); break;
      }
      break;      
    }
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
                             alpha,include_inertia,analysis,cm);

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
      for (kk=0;kk<ndofn;kk++){
	II = node[nod[l]].id_map[mp_id].id[kk]-1;
	if (II < 0)  continue;
	f_u[II] += fe[l*ndofn+kk];
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
			  PGFEM_HYPRE_solve_info *PGFEM_hypre,
			  const int mp_id)
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
  get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nod,node,cnL,mp_id);
  get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nod,node,cnG,mp_id);

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
	II = node[nod[j]].id_map[mp_id].id[P];
	    
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
	II = node[nod[l]].id_map[mp_id].id[kk]-1;
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
			   PGFEM_HYPRE_solve_info *PGFEM_hypre,
			   const int mp_id)
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
  int ndof_ve = get_ndof_on_bnd_elem(node,ptr_be,elem,ndofn);

  long *cn_ve = aloc1l(ndof_ve);
  long *Gcn_ve = aloc1l(ndof_ve);

  get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn_ve, mp_id);
  get_dof_ids_on_bnd_elem(1,ndofn,node,ptr_be,elem,Gcn_ve,mp_id);

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
   
   
int el_compute_stiffmat_MP(FEMLIB *fe,
                           double *lk,
                           GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *fv,
                           SOLVER_OPTIONS *sol,
                           LOADING_STEPS *load,
                           CRPL *crpl,
                           const PGFem3D_opt *opts,
                           MULTIPHYSICS *mp,
                           int mp_id,
                           double dt,
                           double lm,
                           double *be,
                           double *r_e)
{
  int err = 0;  
  int eid = fe->curt_elem_id;

  const int mat_id = grid->element[eid].mat[2];
  double rho = mat->hommat[mat_id].density;

  long *nod = (fe->node_id).m_pdata; // list of node ids in this element    

  // make a decision to include ineria
  long include_inertia = 1;
  if(fabs(rho)<MIN_DENSITY)
  {
    include_inertia = 0;
  }
  
  int nVol = N_VOL_TREE_FIELD;
  SUPP sup = load->sups[mp_id];  
  
  double *x = (fe->temp_v).x.m_pdata;
  double *y = (fe->temp_v).y.m_pdata;
  double *z = (fe->temp_v).z.m_pdata;    
  
  if(include_inertia)
  {
    err += stiffness_with_inertia(fe,lk,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);
    /*
    stiffmat_disp_w_inertia_el(lk,eid,fv->ndofn,fe->nne,fv->npres,nVol,grid->nsd,x, y, z,
            grid->element,mat->hommat,nod,grid->node,dt,
            fv->sig,fv->eps,sup,opts->analysis_type,opts->cm,sol->alpha,fv->u_n,r_e);*/
  }
  else
  {
    switch(opts->analysis_type){
      case STABILIZED:
        err += stiffmatel_st(eid,fv->ndofn,fe->nne,x,y,z,grid->element,mat->hommat,nod,grid->node,fv->sig,fv->eps,
                sup,r_e,fv->npres,sol->nor_min,lk,dt,opts->stab,sol->FNR,lm,be);
        break;
      case MINI:
        err += MINI_stiffmat_el(lk,eid,fv->ndofn,fe->nne,x,y,z,grid->element,
                mat->hommat,nod,grid->node,fv->eps,fv->sig,r_e);
        break;
      case MINI_3F:
        err += MINI_3f_stiffmat_el(lk,eid,fv->ndofn,fe->nne,x,y,z,grid->element,
                mat->hommat,nod,grid->node,fv->eps,fv->sig,r_e);
        break;
      case DISP:
        err += DISP_stiffmat_el(lk,eid,fv->ndofn,fe->nne,x,y,z,grid->element,
                mat->hommat,nod,grid->node,fv->eps,fv->sig,sup,r_e,dt);
        break;
      case TF:
        stiffmat_3f_el(lk,eid,fv->ndofn,fe->nne,fv->npres,nVol,grid->nsd,
                x,y,z,grid->element,mat->hommat,nod,grid->node,dt,fv->sig,fv->eps,sup,-1.0,r_e);
        break;
      case CM:
        err += stiffness_el_constitutive_model(fe,lk,r_e,grid,mat,fv,sol,load,crpl,
                                               opts,mp,mp_id,dt);
            
        break;
      default:
        err += stiffmatel_fd (eid,fv->ndofn,fe->nne,nod,x,y,z,grid->element,mat->matgeom,
                mat->hommat,grid->node,fv->sig,fv->eps,r_e,fv->npres,
                sol->nor_min,lk,dt,crpl,sol->FNR,lm,be,opts->analysis_type);
        break;
    } // switch (analysis)
  } // if(include_inertia)
  
  if(err>0)
  {
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    
    printf("Error is detected: rank= %d, eid = %d\n",myrank, eid);
  }
  return err;
}



/* This function may not be used outside this file */
static int el_stiffmat_MP(int eid,
        double **Lk,
        int *Ddof,
        int interior,
        GRID *grid,
        MATERIAL_PROPERTY *mat,
        FIELD_VARIABLES *fv,
        SOLVER_OPTIONS *sol,
        LOADING_STEPS *load,
        COMMUNICATION_STRUCTURE *com,
        CRPL *crpl,
        MPI_Comm mpi_comm,
        const PGFem3D_opt *opts,
        MULTIPHYSICS *mp,
        int mp_id,
        double dt,
        long iter,
        int myrank)
{
  int err = 0;
  double lm = 0.0;
  int intg_order = 0;
  
  SUPP sup = load->sups[mp_id];

  // multi-scale simulation, displacement base and three-fied-mixed (elastic only)
  // are total Lagrangian based  
  int total_Lagrangian = 0;  
  switch(opts->analysis_type)
  {
    case DISP: // intented to flow
    case TF:
      total_Lagrangian = 1;
      break;
    case CM:
      if(opts->cm != UPDATED_LAGRANGIAN)
        total_Lagrangian = 1;
      
      break;
  }
  
  if(sup->multi_scale)
    total_Lagrangian = 1;
    
  // set FEMLIB  
  FEMLIB fe;
  if (opts->analysis_type == MINI || opts->analysis_type == MINI_3F)
    FEMLIB_initialization_by_elem_w_bubble(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);    
  else
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);  
  
  long *nod = (fe.node_id).m_pdata; // list of node ids in this element  
  
  /* Element Dof */
  long ndofe = get_ndof_on_elem_nodes(fe.nne,nod,grid->node,fv->ndofn);
  
  
  long*cnL,*cnG;
  double *r_e,*sup_def,*be;
  
  /* allocation */
  cnL = aloc1l (ndofe);
  cnG = aloc1l (ndofe);
  r_e = aloc1 (ndofe);
  be = aloc1 (ndofe);

  if(sup->npd>0)
    sup_def = aloc1(sup->npd);
  else
    sup_def = NULL;
    
  
  //global local ids on element
  get_dof_ids_on_elem_nodes(0,fe.nne,fv->ndofn,nod,grid->node,cnL,mp_id);
  get_dof_ids_on_elem_nodes(1,fe.nne,fv->ndofn,nod,grid->node,cnG,mp_id);
  
  // deformation on element  
  // set the increment of applied def=0 on first iter
  if (iter == 0) 
  {
    for(int j=0;j<sup->npd;j++)
    {
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }
  }
  
  // get the deformation on the element
  if(total_Lagrangian)
    def_elem_total(cnL,ndofe,fv->u_np1,fv->d_u,grid->element,grid->node,sup,r_e);
  else 
    def_elem (cnL,ndofe,fv->d_u,grid->element,grid->node,r_e,sup,0);
  
  /* recover thei increment of applied def on first iter */
  if (iter == 0)
  {
    for(int j=0; j<sup->npd; j++)
      sup->defl_d[j] = sup_def[j];
  }
  
  int nVol = N_VOL_TREE_FIELD;
  Matrix(double) lk;
  Matrix_construct_redim(double,lk,ndofe,ndofe);
  Matrix_init(lk, 0.0);
  
  err += el_compute_stiffmat_MP(&fe,lk.m_pdata,grid,mat,fv,sol,load,
                                crpl,opts,mp,mp_id,dt,lm,be,r_e);

  FEMLIB_destruct(&fe);          
  
  if (PFEM_DEBUG){
    char filename[50];
    switch(opts->analysis_type){
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
  if (periodic == 1 && (sol->FNR == 2 || sol->FNR == 3))
  {
    for(int l=0; l<fe.nne; l++)
    {
      for(int kk=0; kk<fv->ndofn; kk++)
      {
        long II = grid->node[nod[l]].id_map[mp_id].id[kk]-1;
        if(II < 0)
          continue;
          
        fv->f_u[II] += be[l*fv->ndofn+kk];
      }
    }
  }// end periodic
  
  // Assembly
  PLoc_Sparse (Lk,lk.m_pdata,com->Ai,com->Ap,cnL,cnG,ndofe,Ddof,com->GDof,
          myrank,com->nproc,com->comm,interior,sol->PGFEM_hypre,opts->analysis_type);
  
  //  dealocation
  free (cnL);
  free (cnG);
  Matrix_cleanup(lk);
  free (be);
  free (r_e);
  free (sup_def);
  
  return err;
  
} // ELEMENT STIFFNESS

//////////////////////////////////////////////////////////////////////////////////////////////



/// Compute stiffnes
///
/// Computes element stiffness matrices and assembles local
/// part. Off-process portions of the matrix are communicated via
/// non-blocking point-to-point send/receives using information in
/// COMMUN. Elements with global DOFs are computed first to overlap
/// communication with computation of fully local elements.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com communication object
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dt time step
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int stiffmat_fd_MP(GRID *grid,
        MATERIAL_PROPERTY *mat,
        FIELD_VARIABLES *fv,
        SOLVER_OPTIONS *sol,
        LOADING_STEPS *load,
        COMMUNICATION_STRUCTURE *com,
        CRPL *crpl,
        MPI_Comm mpi_comm,
        const PGFem3D_opt *opts,
        MULTIPHYSICS *mp,
        int mp_id,
        double dt,
        long iter,
        int myrank)
{
  int err = 0;
  int total_Lagrangian  = 0;
  int intg_order        = 0;
  int do_assemble       = 1; // udate stiffness matrix
  int compute_load4pBCs = 0; // if 1, compute load due to Dirichlet BCs.
  // if 0, update only stiffness
  double lm = 0.0;
  
  double **Lk,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;
  
  err += init_and_post_stiffmat_comm(&Lk,&recieve,&req_r,&sta_r,
          mpi_comm,com->comm);
  
  Matrix(int) Ddof;
  Matrix_construct_redim(int, Ddof,com->nproc,1);
  
  Ddof.m_pdata[0] = com->DomDof[0];
  for (int ia=1; ia<com->nproc; ia++)
    Ddof.m_pdata[ia] = Ddof.m_pdata[ia-1] + com->DomDof[ia];
  
  for(int eid=0; eid<com->nbndel; eid++)
  {
    err += el_stiffmat_MP(com->bndel[eid],Lk,Ddof.m_pdata,0,grid,mat,fv,sol,load,com,crpl,
            mpi_comm,opts,mp,mp_id,dt,iter,myrank);

    if(err != 0)
      break;
  }
  
  if(err==0)
  {
    // COHESIVE ELEMENTS
    // Need to split into boundary and interior parts as with the
    //  regular elements
    if (opts->cohesive == 1)
    {
      long ndofc = 3;
      if (sol->nor_min < 1.e-10)
        sol->nor_min = 1.e-10;
      
      for(int eid=0; eid<grid->nce; eid++)
      {
        coel_stiffmat(eid,Lk,com->Ap,com->Ai,ndofc,grid->element,grid->node,fv->eps,
                fv->d_u,fv->u_np1,fv->npres,load->sups[mp_id],iter,sol->nor_min,dt,crpl,
                opts->stab,grid->coel,0,0,fv->f_u,myrank,com->nproc,com->DomDof,
                com->GDof,com->comm,Ddof.m_pdata,0,opts->analysis_type,sol->PGFEM_hypre, mp_id);
      }
    }
  }
  
  if(err==0)
  {
    // BOUNDING ELEMENTS
    // In the future, these elements will be listed as with the
    // volumetric elements to properly overlay computation and
    // communication. For now, this is the best place for them as they
    // are a proportinally smaller group than the interior volume
    // elements for typical problems and are guarenteed to be on the
    // communication boundary for periodic domains.
    
    // temporary for compile testing
    // int n_be = 0;
    // BOUNDING_ELEMENT *b_elems = NULL;
    for(int eid=0; eid<grid->n_be; eid++){
      err += bnd_el_stiffmat(eid,Lk,com->Ap,com->Ai,fv->ndofn,grid->element,grid->b_elems,grid->node,mat->hommat,
              mat->matgeom,fv->sig,fv->eps,fv->d_u,fv->u_np1,fv->npres,load->sups[mp_id],iter,sol->nor_min,
              dt,crpl,opts->stab,sol->FNR,lm,fv->f_u,myrank,com->nproc,com->GDof,
              com->comm,Ddof.m_pdata,0,opts->analysis_type,sol->PGFEM_hypre,mp_id);
      
      // If there is an error, complete exit the loop
      if(err != 0)
        break;
    }
  }
  
  if (PFEM_DEBUG)
  {
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
    for(int send_proc=0; send_proc<com->nproc; send_proc++)
    {
      if(send_proc==myrank || (com->comm)->S[send_proc] ==0)
        continue;
      
      for(int n_data=0; n_data<(com->comm)->AS[send_proc]; n_data++)
      {
        PGFEM_fprintf(out,"%12.12e ",Lk[send_proc][n_data]);
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n");
    fclose(out);
  }
  
  /// SEND BOUNDARY AND COEL
  err += send_stiffmat_comm(&sta_s,&req_s,Lk,mpi_comm,com->comm);
  
  int skip = 0;
  int idx  = 0;
  
  for(int eid=0; eid<grid->ne; eid++)
  {
    int is_it_in = is_element_interior(eid,&idx,&skip,com->nbndel,com->bndel,myrank);
    
    if(is_it_in==-1)
    {
      err = 1;
      break;
    }
    
    if(is_it_in==0)
      continue;
    
    
    // do volume integration at an element
    int interior = 1;
    
    err += el_stiffmat_MP(eid,Lk,Ddof.m_pdata,1,grid,mat,fv,sol,load,com,crpl,
            mpi_comm,opts,mp,mp_id,dt,iter,myrank);
    
    if(err != 0)
      break;
  }
  
  err += assemble_nonlocal_stiffmat(com->comm,sta_r,req_r,sol->PGFEM_hypre,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,com->comm);
  
  // stiffnes build is completed
  // deallocate memory
  for(int ia=0; ia<com->nproc; ia++)
  {
    free (recieve[ia]);
    free(Lk[ia]);
  }
  free (recieve);
  free (Lk);
  free (sta_s);
  free (sta_r);
  free (req_s);
  free (req_r);
  
  Matrix_cleanup(Ddof);
  return err;
}

/// Multiscale simulation interface to compute stiffness matrix
///
/// In computing stiffness matrix, first perform data mapping from multiscale to multiphysics,
/// and call function computing stiffness matrix using mulphysics data structures. 
/// 
/// \param[in] c structure of macroscale information
/// \param[in,out] s contains the information for the history-dependent solution
/// \param[in] opts structure PGFem3D option
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] nor_min nonlinear convergence tolerance
/// \param[in] FNR if 1: Full Newton-Raphson
///                   0: only compute stiffnes at the 1st iteration
/// \param[in] myrank current process rank
/// \param[in] nproc   number of total process
/// \return non-zero on internal error
int stiffmat_fd_multiscale(COMMON_MACROSCALE *c,
                           MACROSCALE_SOLUTION *s,
                           const PGFem3D_opt *opts,
                           long iter,
                           double nor_min,
                           long FNR,
                           int myrank,
                           int nproc)
{
  int err = 0;
  int mp_id = 0;
  
  // initialize and define multiphysics
  MULTIPHYSICS mp;  
  int id = MULTIPHYSICS_MECHANICAL;
  int ndim = c->ndofn;
  int write_no = 0;

  int *coupled_ids = (int *) malloc(sizeof(int));
  char *physicsname = (char *) malloc(sizeof(char)*1024);
  {
    coupled_ids[0] = 0;
    sprintf(physicsname, "Mechanical");
 
    mp.physicsno      = 1;
    mp.physicsname    = &physicsname;
    mp.physics_ids    = &id;
    mp.ndim           = &ndim;
    mp.write_no       = &write_no;
    mp.write_ids      = NULL;
    mp.coupled_ids    = &coupled_ids;
    mp.total_write_no = 0;
  }
  
  // initialize and define mesh object  
  GRID grid;
  grid_initialization(&grid);
  {
    grid.ne          = c->ne;
    grid.nn          = c->nn;
    grid.element     = c->elem;
    grid.node        = c->node;
    grid.nce         = c->nce;
    grid.coel        = c->coel;
  }
  
  // initialize and define field variables
  FIELD_VARIABLES fv;
  {
    field_varialbe_initialization(&fv);
    fv.ndofn  = c->ndofn;
    fv.ndofd  = c->ndofd;
    fv.npres  = c->npres;
    fv.sig    = s->sig_e;
    fv.eps    = s->eps;
    fv.u_np1  = s->r;
    fv.f      = s->f;
    fv.d_u    = s->d_r;
    fv.dd_u   = s->rr;
    fv.R      = s->R;
    fv.f_defl = s->f_defl;
    fv.RR     = s->RR;
    fv.f_u    = s->f_u;
    fv.RRn    = s->RRn;
    fv.BS_x   = s->BS_x;
    fv.BS_f   = s->BS_f;
    fv.BS_RR  = s->BS_RR;
    fv.BS_f_u = s->BS_f_u;
  }                  

  /// initialize and define iterative solver object
  SOLVER_OPTIONS sol;
  {
    solution_scheme_initialization(&sol);
    sol.FNR          = FNR;
    sol.PGFEM_hypre  = c->SOLVER;
    sol.err          = c->lin_err;
    sol.alpha        = 0.0;
    sol.nor_min      = nor_min;              
  }
  
  // initialize and define loading steps object
  LOADING_STEPS load;    
  {
    loading_steps_initialization(&load);
    load.sups     = &(c->supports);
  }
  
  // initialize and define material properties
  MATERIAL_PROPERTY mat;
  {
    material_initialization(&mat);
    mat.hommat  = c->hommat;
    mat.matgeom = c->matgeom;
  }
  
  /// initialize and define communication structures
  COMMUNICATION_STRUCTURE com;
  {
    communication_structure_initialization(&com); 
    com.nproc  = nproc;
    com.Ap     = c->Ap;
    com.Ai     = c->Ai;
    com.DomDof = c->DomDof;
    com.comm   = c->pgfem_comm;
    com.GDof   = c->GDof;
    com.nbndel = c->nbndel;
    com.bndel  = c->bndel;
  }
  
  err += stiffmat_fd_MP(&grid,&mat,&fv,&sol,&load,&com,s->crpl,
                        c->mpi_comm,opts,&mp,mp_id,s->dt,iter,myrank); 
  
  free(coupled_ids);
  free(physicsname);
  
  return err;
}                                        
