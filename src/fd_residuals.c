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
#include "PGFem3D_data_structure.h"

#define ndn 3

/* assemble the element residual to the local portion of the global
   residual vector */
static void fd_res_assemble(double *f_u,
                            const double *fe,
                            const NODE *node,
                            const int nne,
                            const int ndofn,
                            const long *nod,
                            const int mp_id)
{
  for (int k = 0; k < nne; k++) {
    for (int kk = 0; kk < ndofn; kk++){
      int II = node[nod[k]].id_map[mp_id].id[kk] - 1;
      if (II < 0) continue;
      f_u[II] += fe[k * ndofn + kk];
    }
  }
}

/* compute the local residual on a single element */
static int fd_res_elem(double *fe,
                       const int i,
                       ELEMENT *elem,
                       const long ndofn,
                       const long npres,
                       const double *d_r,
                       const double *r,
                       NODE *node,
                       MATGEOM matgeom,
                       HOMMAT *hommat,
                       SUPP sup,
                       EPS *eps,
                       SIG *sig,
                       const double nor_min,
                       CRPL *crpl,
                       const double *dts,
                       const double t,
                       const double stab,
                       MPI_Comm mpi_comm,
                       const PGFem3D_opt *opts,
                       const double alpha,
                       double *r_n,
                       double *r_n_1,
                       const int include_inertia,
                       const int mp_id)
{
  int err = 0;
  double dt = dts[DT_NP1];

  /* Number of element nodes */
  int nne = elem[i].toe;
  int nne_t = nne + elem[i].n_bub;
  /* Nodes on element */
  long *nod = aloc1l (nne);
  elemnodes (i,nne,nod,elem);
  /* Element Dof */
  int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

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

  long *cn = aloc1l (ndofe);
  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

  /* coordinates */
  if(sup->multi_scale) {
    nodecoord_total(nne,nod,node,x,y,z);
    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
  } else {
    switch(opts->analysis_type) {
    case DISP: case TF:
      /* total Lagrangian */
      nodecoord_total(nne,nod,node,x,y,z);
      def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
      break;

    case CM:
      switch(opts->cm) {
      case UPDATED_LAGRANGIAN:
        nodecoord_updated (nne,nod,node,x,y,z);
        def_elem(cn,ndofe,d_r,elem,node,r_e,sup,0);
        break;
      case TOTAL_LAGRANGIAN:        
        nodecoord_total (nne,nod,node,x,y,z);
        def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
        break;
      case MIXED_ANALYSIS_MODE:
        nodecoord_total (nne,nod,node,x,y,z);
        def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
        break;
      default: assert(0 && "undefined CM type"); break;
      }
      break;

    default: /* updated Lagrangian */
      nodecoord_updated(nne,nod,node,x,y,z);
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

  if(include_inertia) {
    err += residuals_w_inertia_el(fe,i,nne,ndofn,npres,nVol,
                           ndofe,r_e,node,elem,hommat,sup,eps,sig,
                           nod,cn,x,y,z,dts,t,opts,alpha,r_n,r_n_1);
  } else {
    /* Residuals on element */
    switch(opts->analysis_type) {
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

        err =  DISP_resid_el(fe,i,ndofn,nne,x,y,z,elem,
                             hommat,nod,node,eps,sig,sup,r_e,dt);
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
    case CM:
      {
        switch(opts->cm)
          {
          case UPDATED_LAGRANGIAN: // intentionally left to flow
          case TOTAL_LAGRANGIAN:
            err += residuals_el_crystal_plasticity(fe,i,ndofn,nne,nsd,elem,nod,node,
                                                   dt,eps,sup,r_e, opts->cm);
            break;
          case MIXED_ANALYSIS_MODE:
            err += residuals_el_crystal_plasticity(fe,i,ndofn,nne,nsd,elem,nod,node,
                                                   dt,eps,sup,r_e, 1 /* TL */);
            break;
          default: assert(0 && "undefined CM type"); break;
          }
        break;
      }
    default:
      err = resid_on_elem (i,ndofn,nne,nod,elem,node,matgeom,
                           hommat,x,y,z,eps,sig,r_e,npres,
                           nor_min,fe,crpl,dt,opts->analysis_type);

      break;
    }
  }


  dealoc1l (nod);
  dealoc1 (r_e);
  dealoc1 (x);
  dealoc1 (y);
  dealoc1 (z);
  dealoc1l (cn);

  return err;
}

/* compute the residual for a single cohesive element */
static int fd_res_coel(double *fe,
                       const int i,
                       NODE *node,
                       COEL *coel,
                       SUPP sup,
                       const int ndofc,
                       const double *d_r,
                       const double nor_min,
                       const int myrank,
                       const int mp_id)
{
  int err = 0;
  const int nne = coel[i].toe/2;
  const int nnet = 2 * nne;
  const int ndofe = coel[i].toe*ndofc;

  long *nod = aloc1l (nnet);
  double *r_e = aloc1 (ndofe);
  double *x = aloc1 (nnet);
  double *y = aloc1 (nnet);
  double *z = aloc1 (nnet);
  long *cn = aloc1l (ndofe);

  for (int j = 0; j < nnet; j++)
    nod[j] = coel[i].nod[j];

  nodecoord_updated (nnet,nod,node,x,y,z);

  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,nnet,ndofc,nod,node,cn,mp_id);

  /* deformation on element */
  def_elem (cn,ndofe,d_r,NULL,node,r_e,sup,0);

  /* Residuals on element */
  resid_co_elem (i,ndofc,nne,nod,x,y,z,coel,r_e,fe,nor_min,myrank);

  dealoc1l (nod);
  dealoc1 (r_e);
  dealoc1 (x);
  dealoc1 (y);
  dealoc1 (z);
  dealoc1l (cn);

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
                  double *dts,
                  double t,
                  double stab,
                  long nce,
                  COEL *coel,
                  MPI_Comm mpi_comm,
                  const PGFem3D_opt *opts,
                  double alpha,
                  double *r_n,
                  double *r_n_1,
                  const int mp_id)
{
  /* make decision to include ineria*/
  const int mat = elem[0].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;

  if(fabs(rho)<MIN_DENSITY)
    include_inertia = 0;

  /* decision end*/
  int err = 0;

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  for (int i=0;i<ne;i++) {
    const int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    const int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);
    double *fe = aloc1 (ndofe);

    err += fd_res_elem(fe, i, elem, ndofn, npres, d_r, r, node,
                       matgeom, hommat, sup, eps, sig, nor_min,
                       crpl, dts, t, stab, mpi_comm, opts, alpha,
                       r_n, r_n_1, include_inertia,mp_id);

    fd_res_assemble(f_u, fe, node, nne, ndofn, nod, mp_id);

    dealoc1l (nod);
    dealoc1 (fe);

    /*** RETURN on error ***/
    if(err != 0) return err;

  }/* end i < ne*/

  /**** COHESIVE ELEMENT RESIDUALS ****/
  if (opts->cohesive == 1){
    const int ndofc = 3;

    for (int i=0;i<nce;i++){

      int ndofe = coel[i].toe*ndofc;
      long *nod = aloc1l (coel[i].toe);
      double *fe = aloc1 (ndofe);
      for (int j=0;j<coel[i].toe;j++)
        nod[j] = coel[i].nod[j];

      err += fd_res_coel(fe, i, node, coel, sup, ndofc, d_r, nor_min, myrank,mp_id);
      fd_res_assemble(f_u, fe, node, coel[i].toe, ndofc, nod, mp_id);

      dealoc1l (nod);
      dealoc1 (fe);

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

    int ndofe = get_ndof_on_bnd_elem(node,ptr_be,elem,ndofn);
    double *r_e = aloc1(ndofe);
    double *fe = aloc1(ndofe);
    long *cn = aloc1l(ndofe);
    long *Gcn = aloc1l(ndofe);

    get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn ,mp_id);
    get_dof_ids_on_bnd_elem(1,ndofn,node,ptr_be,elem,Gcn,mp_id);

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


/// Compute residuals on a single element
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \return non-zero on internal error

static int fd_res_elem_MP(double *be,
        const int eid,
        GRID *grid,
        MATERIAL_PROPERTY *mat,
        FIELD_VARIABLES *fv,
        SOLVER_OPTIONS *sol,
        LOADING_STEPS *load,
        CRPL *crpl,
        MPI_Comm mpi_comm,
        const PGFem3D_opt *opts,
        MULTIPHYSICS *mp,
        int mp_id,
        double t,
        double *dts,
        int include_inertia,
        int updated_deformation)
{
  int err = 0;
  int intg_order = 0;
  double dt = dts[DT_NP1];
  
  ELEMENT *elem = grid->element;
  SUPP sup = load->sups[mp_id];
  
  double *f;
  if(updated_deformation)
    f = fv->f;
  else
    f = fv->d_u;
  
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
  
  long *cn;
  double *r_e;
  
  /* allocation */
  cn  = aloc1l (ndofe);
  r_e = aloc1 (ndofe);
  
  //global local ids on element
  get_dof_ids_on_elem_nodes(0,fe.nne,fv->ndofn,nod,grid->node,cn,mp_id);
  
  // get the deformation on the element
  if(total_Lagrangian)
    def_elem_total(cn,ndofe,fv->u_np1,f,grid->element,grid->node,sup,r_e);
  else
    def_elem (cn,ndofe,f,grid->element,grid->node,r_e,sup,0);
  
  int nVol = N_VOL_TREE_FIELD;
  int nsd = 3;
  
  double *x = (fe.temp_v).x.m_pdata;
  double *y = (fe.temp_v).y.m_pdata;
  double *z = (fe.temp_v).z.m_pdata;
  
  if(include_inertia) {
    err += residuals_w_inertia_el(be,eid,fe.nne,fv->ndofn,fv->npres,nVol,
                                  ndofe,r_e,grid->node,elem,mat->hommat,sup,fv->eps,fv->sig,
                                  nod,cn,x,y,z,dts,t,opts,sol->alpha,fv->u_n,fv->u_nm1);
  } else {
    /* Residuals on element */
    switch(opts->analysis_type) {
      case STABILIZED:
        err = resid_st_elem (eid,fv->ndofn,fe.nne,elem,nod,grid->node,mat->hommat,
                             x,y,z,fv->eps,fv->sig,sup,r_e,sol->nor_min,be,dt,opts->stab);
        break;
      case MINI:
        err = MINI_resid_el(be,eid,fv->ndofn,fe.nne,x,y,z,elem,
                            nod,grid->node,mat->hommat,fv->eps,fv->sig,r_e);
        break;
      case MINI_3F:
        err = MINI_3f_resid_el(be,eid,fv->ndofn,fe.nne,x,y,z,elem,
                               nod,grid->node,mat->hommat,fv->eps,fv->sig,r_e);
        break;
      case DISP:
      {
        double *bf = aloc1(ndofe);
        memset(bf, 0, sizeof(double)*ndofe);
        DISP_resid_body_force_el(bf,eid,fv->ndofn,fe.nne,x,y,z,elem,mat->hommat,grid->node,dt,t);
        
        err =  DISP_resid_el(be,eid,fv->ndofn,fe.nne,x,y,z,elem,
                             mat->hommat,nod,grid->node,fv->eps,fv->sig,sup,r_e,dt);
        for(long a = 0; a<ndofe; a++)
          be[a] += -bf[a];
        
        dealoc1(bf);
        break;
      }
      case TF:
      {
        double *bf = aloc1(ndofe);
        memset(bf, 0, sizeof(double)*ndofe);
        DISP_resid_body_force_el(bf,eid,fv->ndofn,fe.nne,x,y,z,elem,mat->hommat,grid->node,dt,t);
                
        residuals_3f_el(be,eid,fv->ndofn,fe.nne,fv->npres,nVol,nsd,
                        x,y,z,elem,mat->hommat,nod,grid->node,dt,fv->sig,fv->eps,sup,r_e);
        for(long a = 0; a<ndofe; a++)
          be[a] += -bf[a];
        
        dealoc1(bf);
        break;
      }
      case CM:
        err += residuals_el_crystal_plasticity(be,eid,fv->ndofn,fe.nne,nsd,elem,nod,grid->node,
                                               dt,fv->eps,sup,r_e, total_Lagrangian);
        break;      
      default:
        err = resid_on_elem (eid,fv->ndofn,fe.nne,nod,elem,grid->node,mat->matgeom,
                             mat->hommat,x,y,z,fv->eps,fv->sig,r_e,fv->npres,
                             sol->nor_min,be,crpl,dt,opts->analysis_type);
        
        break;
    }
  }
  
  
  dealoc1 (r_e);
  dealoc1l (cn);
  
  return err;
}

/// Compute residuals
///
/// Compute redidual vector for mechanical problem. 
/// Integration algorithm will be perfromed based on constitutive model.
/// If either integration algorithm if faild to converge or jacobian of the 
/// deformation gradient is small than zero, element loop will be stopped and
/// return non-zero value.    
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \return non-zero on internal error
long fd_residuals_MP(GRID *grid,
                     MATERIAL_PROPERTY *mat,
                     FIELD_VARIABLES *fv,
                     SOLVER_OPTIONS *sol,
                     LOADING_STEPS *load,
                     CRPL *crpl,
                     MPI_Comm mpi_comm,
                     const PGFem3D_opt *opts,
                     MULTIPHYSICS *mp,
                     int mp_id,
                     double t,
                     double *dts,
                     int updated_deformation)
{
  int err = 0;
  ELEMENT *elem = grid->element;
  BOUNDING_ELEMENT *b_elems = grid->b_elems;
  SUPP sup = load->sups[mp_id];
   
  double *f;
  if(updated_deformation)
    f = fv->f;
  else
    f = fv->d_u;  
  
  /* make decision to include ineria*/
  const int mat_id = grid->element[0].mat[2];
  double rho = mat->hommat[mat_id].density;
  long include_inertia = 1;

  if(fabs(rho)<MIN_DENSITY)
    include_inertia = 0;

  /* decision end*/

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  for (int i=0;i<grid->ne;i++) {
    const int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    const int ndofe = get_ndof_on_elem_nodes(nne,nod,grid->node,fv->ndofn);
    double *fe = aloc1 (ndofe);

    err += fd_res_elem_MP(fe, i, grid, mat, fv, sol, load, crpl,
                          mpi_comm, opts, mp, mp_id, t, dts, 
                          include_inertia, updated_deformation);

    fd_res_assemble(fv->f_u, fe, grid->node, nne, fv->ndofn, nod, mp_id);

    dealoc1l (nod);
    dealoc1 (fe);

    /*** RETURN on error ***/
    if(err != 0) return err;

  }/* end i < grid->ne*/

  /**** COHESIVE ELEMENT RESIDUALS ****/
  if (opts->cohesive == 1){
    COEL *coel = grid->coel;
    const int ndofc = 3;

    for (int i=0;i<grid->nce;i++){

      int ndofe = coel[i].toe*ndofc;
      long *nod = aloc1l (coel[i].toe);
      double *fe = aloc1 (ndofe);
      for (int j=0;j<coel[i].toe;j++)
        nod[j] = coel[i].nod[j];

      err += fd_res_coel(fe, i, grid->node, coel, sup, ndofc, f, sol->nor_min, myrank,mp_id);
      fd_res_assemble(fv->f_u, fe, grid->node, coel[i].toe, ndofc, nod, mp_id);

      dealoc1l (nod);
      dealoc1 (fe);

    }/* end i < grid->nce */
  }/* end coh == 1 */

  /*===============================================
    |             BOUNDARY ELEMENTS               |
    ===============================================*/

  for (int i=0; i<grid->n_be; i++){

    /* get the coordinates and dof id's on the element nodes */
    const long *ptr_vnodes = elem[b_elems[i].vol_elem_id].nod; /* --"-- */
    const ELEMENT *ptr_elem = elem + b_elems[i].vol_elem_id; /* --"-- */
    const BOUNDING_ELEMENT *ptr_be = &b_elems[i];
    const int nne = ptr_elem->toe;

    double *x = aloc1(nne);
    double *y = aloc1(nne);
    double *z = aloc1(nne);

    switch(opts->analysis_type){
    case DISP:
      nodecoord_total(nne,ptr_vnodes,grid->node,x,y,z);
      break;
    default:
      nodecoord_updated(nne,ptr_vnodes,grid->node,x,y,z);
      break;
    }

    int ndofe = get_ndof_on_bnd_elem(grid->node,ptr_be,elem,fv->ndofn);
    double *r_e = aloc1(ndofe);
    double *fe = aloc1(ndofe);
    long *cn = aloc1l(ndofe);
    long *Gcn = aloc1l(ndofe);

    get_dof_ids_on_bnd_elem(0,fv->ndofn,grid->node,ptr_be,elem,cn ,mp_id);
    get_dof_ids_on_bnd_elem(1,fv->ndofn,grid->node,ptr_be,elem,Gcn,mp_id);

    /* compute the deformation on the element */
    def_elem(cn,ndofe,f,NULL,NULL,r_e,sup,0);

    /* TOTAL LAGRANGIAN formulation */
    if(opts->analysis_type == DISP){
      double *r_en = aloc1(ndofe);
      def_elem(cn,ndofe,fv->u_np1,NULL,NULL,r_en,sup,1);
      vvplus(r_e,r_en,ndofe);
      free(r_en);
    }

    /* for debugging */
    double *RR = aloc1(ndofe);
    def_elem(cn,ndofe,fv->f_u,NULL,NULL,RR,sup,2);

    if(opts->analysis_type == DISP){
      err += DISP_resid_bnd_el(fe,i,fv->ndofn,ndofe,x,y,z,b_elems,elem,
			       mat->hommat,grid->node,fv->eps,fv->sig,sup,r_e);
    } else {
      /* not implemented/needed so do nothing */
    }

    /* Assembly to local part of the residual vector */
    {
      for(int j=0; j<ndofe; j++){
	int II = cn[j] - 1;
	if(II >= 0){
	  fv->f_u[II] += fe[j];
	}
      }
    }

    def_elem(cn,ndofe,fv->f_u,NULL,NULL,RR,sup,2);

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

/* compute the reaction force for each magnitude of prescribed
   deflection. CAVEATS: Does not include contributions from cohesive
   or boundary elements. */
int fd_res_compute_reactions(const long ndofn,
                             const long npres,
                             const double *d_r,
                             const double *r,
                             ELEMENT *elem,
                             NODE *node,
                             MATGEOM matgeom,
                             HOMMAT *hommat,
                             SUPP sup,
                             EPS *eps,
                             SIG *sig,
                             const double nor_min,
                             CRPL *crpl,
                             const double *dts,
                             const double t,
                             const double stab,
                             MPI_Comm mpi_comm,
                             const PGFem3D_opt *opts,
                             const double alpha,
                             double *r_n,
                             double *r_n_1,
                             const int mp_id)
{
  int err = 0;

  /* make decision to include ineria*/
  const long include_inertia = 0; //computing reaction in quasi-steady state gives correct results

  const int ne = sup->nde;
  const long *el_id = sup->lepd;
  const int n_rxn = sup->npd + 1;
  double *rxn = PGFEM_calloc(n_rxn, sizeof(*rxn));
  double *RXN = PGFEM_calloc(n_rxn, sizeof(*RXN));
  for (int i = 0; i < ne; i++) {
    const int nne = elem[el_id[i]].toe;
    long *nod = aloc1l (nne);
    elemnodes (el_id[i],nne,nod,elem);
    const int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);
    double *fe = aloc1 (ndofe);

    err += fd_res_elem(fe, el_id[i], elem, ndofn, npres, d_r, r, node,
                       matgeom, hommat, sup, eps, sig, nor_min, crpl,
                       dts, t, stab, mpi_comm, opts, alpha, r_n, r_n_1,
                       include_inertia,mp_id);

    /* Previous may have called integration algorithm. Need to reset
       state variables to retain consistent tangent and to ensure we
       didn't play with any rate sensitive behavior */
    if (opts->analysis_type == CM) constitutive_model_reset_state(eps, ne, elem);

    long *cn = aloc1l (ndofe);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    for (int j = 0; j < ndofe; j++) {
      if (cn[j] <= 0) {
        rxn[labs(cn[j])] += fe[j];
      }
    }
    free(nod);
    free(cn);
    free(fe);
  }

  /* communicate reactions on all domains */
  int myrank = -1;
  err += MPI_Comm_rank(mpi_comm, &myrank);
  err += MPI_Reduce(rxn, RXN, n_rxn, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
  if (myrank == 0) {
    PGFEM_printf("Reactions: (fixed 1 ... n)\n");
    print_array_d(PGFEM_stdout, RXN, n_rxn, 1, n_rxn);
  }

  free(rxn);
  free(RXN);
  return err;
}
