#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "PGFem3D_data_structure.h"
#include "allocation.h"
#include "cast_macros.h"
#include "condense.h"
#include "constitutive_model.h"
#include "def_grad.h"
#include "displacement_based_element.h"
#include "dynamics.h"
#include "elem3d.h"
#include "enumerations.h"
#include "fd_residuals.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "index_macros.h"
#include "macro_micro_functions.h"
#include "matice.h"
#include "mkl_cblas.h"
#include "new_potentials.h"
#include "resid_on_elem.h"
#include "solver_file.h"
#include "stabilized.h"
#include "tensors.h"
#include "three_field_element.h"
#include "utils.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <ttl/ttl.h>

using namespace multiscale::net;

namespace {
using pgfem3d::Solver;
using pgfem3d::CommunicationStructure;
using pgfem3d::MultiscaleCommon;
using pgfem3d::MULTISCALE_SOLUTION;
}


template <int R, class S = double>
using Tensor = ttl::Tensor<R, 3, S>;

template <int R, class S = double *>
using TensorA = ttl::Tensor<R, 3, S>;
  
static constexpr ttl::Index<'i'> i;
static constexpr ttl::Index<'j'> j;

  
/* assemble the element residual to the local portion of the global
   residual vector */
static void fd_res_assemble(double *f_u,
                            const double *fe,
                            const Node *node,
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

/* compute the residual for a single cohesive element */
static int fd_res_coel(double *fe,
                       const int i,
                       Node *node,
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

void compute_Neumann_boundary_conditions(FEMLIB *fe,
                                         double *be,
                                         const double *r_e,
                                         Grid *grid,
                                         FieldVariables *fv,
                                         int mp_id,
                                         double t,
                                         double *dts,
                                         const int include_inertia,
                                         const double alpha){
  // no NBE available
  if(grid->NBE.m_row*grid->NBE.m_col == 0)
    return;
    
  if(grid->NBE(mp_id).nbe_no == 0)
    return;
  
  // check validation of this integration  
  const int eid = fe->curt_elem_id;
  if(grid->NBE(mp_id).element_check(eid, 0)==0)
    return;
  
  int this_bnd_elem_id = grid->NBE(mp_id).element_check(eid, 1);
    
  // number of boundary element to be integrated
  const int bnd_elem_no = grid->NBE(mp_id).element_ids(this_bnd_elem_id, 1);

  // number of spatial dimensions  
  const int nsd   = fe->nsd;
  const int ndofn = fv->ndofn;
  
  // mid point rule
  double dt_1_minus_alpha =  -1.0;
  double dt_alpha         =  0.0;
    
  // displacement at t[n] and t[n-1]
  Matrix<double> Un, Unm1;  
  
  if(include_inertia){ 
    dt_1_minus_alpha = dts[DT_NP1]*(1.0-alpha);
    dt_alpha = dts[DT_N]*alpha;    
    
      Un.initialization(fe->nne*nsd, 1, 0.0);
    Unm1.initialization(fe->nne*nsd, 1, 0.0);
    
    for(int ia=0; ia<fe->nne; ++ia){
      for(int ib=0; ib<nsd; ++ib){
        Un(  ia*nsd+ib) = fv->u_n[  fe->node_id(ia)*ndofn+ib];
        Unm1(ia*nsd+ib) = fv->u_nm1[fe->node_id(ia)*ndofn+ib];
      }
    }
  }

  Constitutive_model *m = fv->eps[eid].model;

  for(int iA = 0; iA<bnd_elem_no; ++iA){
    const int face_id    = grid->NBE(mp_id).bnd_elements(this_bnd_elem_id)(iA, 0);
    const int feature_id = grid->NBE(mp_id).bnd_elements(this_bnd_elem_id)(iA, 1);
    const int load_type  = grid->NBE(mp_id).load_type(feature_id); // same as number of boundary loads

    double Pnp1[3] = {}, Pn[3] = {}, Pnm1[3] = {};

    // compute pressure at t_np1
    for(int ia=0; ia<load_type; ++ia){
      Pnp1[ia] = grid->NBE(mp_id).load_values_np1(feature_id, ia);
    }
        
    // compute pressure at t[n] and t[n-1] if transient

    if(include_inertia){                
      for(int ia=0; ia<load_type; ++ia){
        Pn[  ia] = grid->NBE(mp_id).load_values_n(  feature_id, ia);
        Pnm1[ia] = grid->NBE(mp_id).load_values_nm1(feature_id, ia);
      }
    }
    
    if(load_type == 1){ // pressure bounary condition
      for(int ia=2; ia<nsd; ++ia){
        Pnp1[ia] = Pnp1[0];
        Pn[  ia] = Pn[  0];
        Pnm1[ia] = Pnm1[0];        
      }
    }    

    // start boundary element FEM integration for
    FemLibBoundary fes(fe, face_id, fe->intg_order);    

    // apply quadrature rule 
    for(int ip=0; ip<fes.nint; ++ip){      
      fes.elem_basis_S(ip);
      fes.update_shape_tensor();
      
      // compute pressure terms at t(n+1), t(n), and t(n-1)
      double PNnp1[3], PNn[3], PNnm1[3];      
      if(load_type == 1){ // pressure boundary condition
        TensorA<1> N(fes.normal);
        for(int ia=0; ia<nsd; ++ia){
          PNnp1[ia] = Pnp1[0]*fes.normal[ia];
          PNn[  ia] = Pn[  0]*fes.normal[ia];
          PNnm1[ia] = Pnm1[0]*fes.normal[ia];
        }
      } else{             // traction boundary condition
        for(int ia=0; ia<nsd; ++ia){
          PNnp1[ia] = Pnp1[ia];
          PNn[  ia] = Pn[  ia];
          PNnm1[ia] = Pnm1[ia];
        }
      }
      
      // compute deformation gradient * pressure at t(n+a) and t(n-1-a)
      Tensor<1> FITNnpa = {}, FITNnma = {};
      double Jnpa = {}, Jnma = {};
      
      if(include_inertia){
        Tensor<1> PNnpa = {}, PNnma = {};
        mid_point_rule(PNnpa.data, PNn,   Pnp1, alpha, nsd);
        mid_point_rule(PNnma.data, PNnm1, Pn,   alpha, nsd);
        Matrix<double> Unpa(fe->nne*nsd, 1, 0.0), Unma(fe->nne*nsd, 1, 0.0);
        mid_point_rule(Unpa.m_pdata,   Un.m_pdata, r_e,        alpha, fe->nne*nsd);
        mid_point_rule(Unma.m_pdata, Unm1.m_pdata, Un.m_pdata, alpha, fe->nne*nsd);
        Tensor<2> Fnpa = {}, Fnma = {}, FnpaI = {}, FnmaI = {};
        fes.update_deformation_gradient(nsd, Unpa.m_pdata, Fnpa.data, (m->param)->pFI);
        fes.update_deformation_gradient(nsd, Unma.m_pdata, Fnma.data, (m->param)->pFI);

        Jnpa = ttl::det(Fnpa);
        FnpaI = ttl::inverse(Fnpa);
        Jnma = ttl::det(Fnma);
        FnmaI = ttl::inverse(Fnma);
        FITNnpa = FnpaI(j, i)*PNnpa(j);
        FITNnma = FnmaI(j, i)*PNnma(j);
         
      } else{
        Tensor<2> Fnpa = {}, FnpaI = {};
        TensorA<1> PNnpa(PNnp1);      
        fes.update_deformation_gradient(nsd, r_e, Fnpa.data, (m->param)->pFI);
        Jnpa = ttl::det(Fnpa);
        FnpaI = ttl::inverse(Fnpa);
        FITNnpa = FnpaI(j, i)*PNnpa(j);
      }

      for(int ia = 0; ia<fes.nne_bnd; ++ia){
        for(int ib = 0; ib<nsd; ++ib){
          int nid = fes.Volume2Boundary[ia];
          be[nid*ndofn+ib] += fes.N(nid)*fes.detJxW*(
                              dt_1_minus_alpha*Jnpa*FITNnpa.data[ib] + 
                                      dt_alpha*Jnma*FITNnma.data[ib]);
        } 
      }
    }
  }  
}

/// Compute residuals on a single element
///
/// \param[out]    be              computed element residual
/// \param[in]     eid             element id
/// \param[in]     grid            a mesh object
/// \param[in]     mat             a material object
/// \param[in,out] fv              variables object for field variables
/// \param[in]     sol             object for solution scheme
/// \param[in]     load            object for loading
/// \param[in]     crpl            object for lagcy crystal plasticity
/// \param[in]     com             communicatio structure
/// \param[in]     opts            structure PGFem3D option
/// \param[in]     mp_id           mutiphysics id
/// \param[in]     t               time
/// \param[in]     dts             time step sizes a n, and n+1
/// \return non-zero on internal error
static int fd_res_elem_MP(double *be,
                          const int eid,
                          Grid *grid,
                          MaterialProperty *mat,
                          FieldVariables *fv,
                          Solver *sol,
                          LoadingSteps *load,
                          CRPL *crpl,
                          const CommunicationStructure *com,
                          const PGFem3D_opt *opts,
                          const Multiphysics& mp,
                          int mp_id,
                          double t,
                          double *dts,
                          int include_inertia,
                          int updated_deformation)
{
  int err = 0;
  int intg_order = 0;
  double dt = dts[DT_NP1];

  Element *elem = grid->element;
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
   case CM:   // intented to flow
   case CM3F:
    if(opts->cm != UPDATED_LAGRANGIAN)
      total_Lagrangian = 1;

    break;
  }

  if(sup->multi_scale)
    total_Lagrangian = 1;

  // set FEMLIB
  FEMLIB fe;

  if(opts->analysis_type == CM || opts->analysis_type == CM3F)
  {
    Constitutive_model *m = (fv->eps[eid]).model;
    double *pFI = m->param->pFI;
    fe.initialization(eid,grid->element,grid->node,intg_order,total_Lagrangian,pFI);
  }
  else
  {
    if (opts->analysis_type == MINI || opts->analysis_type == MINI_3F)
      fe.initialization(eid,grid->element,grid->node,intg_order,total_Lagrangian,NULL,true);
    else
      fe.initialization(eid,grid->element,grid->node,intg_order,total_Lagrangian,NULL);
  }

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

  double *x = (fe.temp_v).x.m_pdata;
  double *y = (fe.temp_v).y.m_pdata;
  double *z = (fe.temp_v).z.m_pdata;

  if(include_inertia) {
    err += residual_with_inertia(&fe,be,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dts,t,com->rank);
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

         residuals_3f_el(&fe, be, r_e, grid, mat, fv);

         for(long a = 0; a<ndofe; a++)
           be[a] += -bf[a];

         dealoc1(bf);
         break;
       }
     case CM:  // intented to flow
     case CM3F:
      err += residuals_el_constitutive_model(&fe,be,r_e,grid,mat,fv,sol,load,crpl,
                                             opts,mp,mp_id,dts,t);
      break;
     default:
      err = resid_on_elem (eid,fv->ndofn,fe.nne,nod,elem,grid->node,mat->matgeom,
                           mat->hommat,x,y,z,fv->eps,fv->sig,r_e,fv->npres,
                           sol->nor_min,be,crpl,dt,opts->analysis_type);

      break;
    }
  }
  
  // applying NB conditions  
  compute_Neumann_boundary_conditions(&fe, be, r_e,
                                      grid, fv, mp_id, t, dts, include_inertia,sol->alpha);

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
/// \param[in] com communication structure
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \return non-zero on internal error
long fd_residuals_MP(Grid *grid,
                     MaterialProperty *mat,
                     FieldVariables *fv,
                     Solver *sol,
                     LoadingSteps *load,
                     CRPL *crpl,
                     const CommunicationStructure *com,
                     const PGFem3D_opt *opts,
                     const Multiphysics& mp,
                     int mp_id,
                     double t,
                     double *dts,
                     int updated_deformation)
{
  int err = 0;
  Element *elem = grid->element;
  BoundingElement *b_elems = grid->b_elems;
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

  int myrank = com->rank;
  
  // Neumann boundary condition loading values
  if(grid->NBE.m_row*grid->NBE.m_col > 0)
    grid->NBE(mp_id).compute_load_values(t, t-dts[DT_NP1], t-dts[DT_NP1]-dts[DT_N]);

  for (int i=0;i<grid->ne;i++) {
    const int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    const int ndofe = get_ndof_on_elem_nodes(nne,nod,grid->node,fv->ndofn);
    double *fe = aloc1 (ndofe);

    err += fd_res_elem_MP(fe, i, grid, mat, fv, sol, load, crpl,
                          com, opts, mp, mp_id, t, dts,
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
    const Element *ptr_elem = elem + b_elems[i].vol_elem_id; /* --"-- */
    const BoundingElement *ptr_be = &b_elems[i];
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
      for(int j=0; j<ndofe; j++)
      {
        int II = cn[j] - 1;
        if(II >= 0)
          fv->f_u[II] += fe[j];
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

/// compute the reaction force for multiphysics mode (for each magnitude of prescribed
/// deflection). CAVEATS: Does not include contributions from cohesive
/// or boundary elements.
///
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] com object for communication
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes at n, and n+1
/// \return non-zero on internal error
int fd_res_compute_reactions_MP(Grid *grid,
                                MaterialProperty *mat,
                                FieldVariables *fv,
                                Solver *sol,
                                LoadingSteps *load,
                                CRPL *crpl,
                                const CommunicationStructure *com,
                                const PGFem3D_opt *opts,
                                const Multiphysics& mp,
                                int mp_id,
                                double t,
                                double *dts)
{
  int err = 0;

  //computing reaction in quasi-steady state gives correct results
  const long include_inertia = 0;

  const int updated_deformation = 0;

  Element *elem = grid->element;
  SUPP sup = load->sups[mp_id];

  const int ne = sup->nde;
  const long *el_id = sup->lepd;
  const int n_rxn = sup->npd + 1;
  double *rxn = PGFEM_calloc(double, n_rxn);
  double *RXN = PGFEM_calloc(double, n_rxn);

  for(int i = 0; i < ne; i++)
  {
    const int nne = elem[el_id[i]].toe;
    long *nod = aloc1l(nne);
    elemnodes(el_id[i],nne,nod,elem);
    const int ndofe = get_ndof_on_elem_nodes(nne,nod,grid->node,fv->ndofn);
    double *fe = aloc1(ndofe);

    err += fd_res_elem_MP(fe, el_id[i], grid, mat, fv, sol, load, crpl,
                          com, opts, mp, mp_id, t, dts,
                          include_inertia, updated_deformation);

    // Previous may have called integration algorithm. Need to reset
    // state variables to retain consistent tangent and to ensure we
    // didn't play with any rate sensitive behavior
    if (opts->analysis_type == CM || opts->analysis_type == CM3F)
      constitutive_model_reset_state(fv->eps, grid->ne, grid->element);

    long *cn = aloc1l (ndofe);
    get_dof_ids_on_elem_nodes(0,nne,fv->ndofn,nod,grid->node,cn,mp_id);

    for (int j = 0; j < ndofe; j++)
    {
      if (cn[j] <= 0)
        rxn[labs(cn[j])] += fe[j];
    }
    free(nod);
    free(cn);
    free(fe);
  }

  /* communicate reactions on all domains */
  int myrank = com->rank;
  try {
    com->net->reduce(rxn, RXN, n_rxn, NET_DT_DOUBLE, NET_OP_SUM, 0, com->comm);
  } catch(...) {
    err += 1;
  }
  if (myrank == 0)
  {
    PGFEM_printf("Reactions: (fixed 1 ... n)\n");
    print_array_d(PGFEM_stdout, RXN, n_rxn, 1, n_rxn);
  }

  free(rxn);
  free(RXN);
  return err;
}

/// Multiscale simulation interface computing reaction forces
///
/// \param[in] c structure of macroscale information
/// \param[in] s contains the information for the history-dependent solution
/// \param[in] solver_file structure for storing/updating the data
/// \param[in] opts structure PGFem3D option
/// \param[in] dts time step sizes at n, and n+1
/// \return non-zero on internal error
int fd_res_compute_reactions_multiscale(MultiscaleCommon *c,
                                        MULTISCALE_SOLUTION *s,
                                        SOLVER_FILE *solver_file,
                                        const PGFem3D_opt *opts,
                                        double *dts)
{
  int err = 0;

  // initialize and define multiphysics
  Multiphysics mp;
  int id = MULTIPHYSICS_MECHANICAL;
  int ndim = c->ndofn;
  int write_no = 0;

  vector<int> coupled_ids;
  char *physicsname = (char *) malloc(sizeof(char)*1024);
  {
    coupled_ids.push_back(0);
    sprintf(physicsname, "Mechanical");

    mp.physicsno      = 1;
    mp.physicsname    = &physicsname;
    mp.physics_ids    = &id;
    mp.ndim           = &ndim;
    mp.write_no       = &write_no;
    mp.coupled_ids.push_back(coupled_ids);
    mp.total_write_no = 0;
  }

  Grid grid;
  grid_initialization(&grid);
  {
    grid.ne          = c->ne;
    grid.nn          = c->nn;
    grid.element     = c->elem;
    grid.b_elems     = NULL;
    grid.node        = c->node;
    grid.nce         = c->nce;
    grid.coel        = c->coel;
    grid.n_be        = 0;
  }

  // initialize and define field variables
  FieldVariables fv;
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
    fv.BS_f   = s->BS_f;
    fv.BS_RR  = s->BS_RR;
    fv.BS_f_u = s->BS_f_u;
    fv.sig_n  = s->sig_n;
    fv.NORM   = s->NORM;
  }

  /// initialize and define iterative solver object
  Solver sol{};
  {
    sol.nor_min  = solver_file->nonlin_tol;
    sol.alpha   = 0.0;
  }

  // initialize and define loading steps object
  LoadingSteps load;
  {
    loading_steps_initialization(&load);
    load.sups     = &(c->supports);
  }

  // initialize and define material properties
  MaterialProperty mat;
  {
    material_initialization(&mat);
    mat.hommat  = c->hommat;
    mat.matgeom = c->matgeom;
  }

  err += fd_res_compute_reactions_MP(&grid,&mat,&fv,&sol,&load,s->crpl,c,opts,mp,
                                     0,s->times[s->tim+1],dts);

  free(physicsname);

  return err;
}
