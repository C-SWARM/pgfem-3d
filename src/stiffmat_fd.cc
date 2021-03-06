/**
 * AUTHORS:
 * Matthew Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * Karel Matous, University of Notre Dame, kmatous [at] nd.edu
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "PGFem3D_data_structure.h"
#include "PLoc_Sparse.h"
#include "allocation.h"
#include "cast_macros.h"
#include "condense.h"
#include "constitutive_model.h"
#include "def_grad.h"
#include "displacement_based_element.h"
#include "dynamics.h"
#include "elem3d.h"
#include "enumerations.h"
#include "femlib.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "incl.h"
#include "index_macros.h"
#include "macro_micro_functions.h"
#include "matice.h"
#include "new_potentials.h"
#include "stabilized.h"
#include "stiffmat_fd.h"
#include "stiffmatel_fd.h"
#include "tensors.h"
#include "three_field_element.h"
#include "utils.h"
#include <mkl_cblas.h>
#include <cassert>
#include <limits>
#include <ttl/ttl.h>

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

using namespace multiscale::net;

namespace {
using pgfem3d::SparseComm;
using pgfem3d::Solver;
using pgfem3d::solvers::SparseSystem;
using pgfem3d::CommunicationStructure;
using pgfem3d::MultiscaleCommon;
using pgfem3d::MULTISCALE_SOLUTION;

const constexpr int periodic = 0;
}

template <int R, class S = double>
using Tensor = ttl::Tensor<R, 3, S>;

template <int R, class S = double *>
using TensorA = ttl::Tensor<R, 3, S>;
  
static constexpr ttl::Index<'i'> i;
static constexpr ttl::Index<'j'> j;
static constexpr ttl::Index<'k'> k;
static constexpr ttl::Index<'l'> l;

/* This function may not be used outside of this file */
static void coel_stiffmat(int i, /* coel ID */
                          double **Lk,
                          int *Ap,
                          Ai_t *Ai,
                          long ndofc,
                          Element *elem,
                          Node *node,
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
                          SparseComm *comm,
                          long *Ddof,
                          int interior,
                          const int analysis,
                          SparseSystem *system,
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
    for (j=0;j<sup->npd;j++){
      assert(sup_def != NULL && "sup_def can't be null");
      sup->defl_d[j] = sup_def[j];
    }

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
               GDof,myrank,nproc,comm,interior,system,analysis);

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
  PGFEM_free(sup_def);

} /* COHESIVE ELEMENT STIFFNESS */

static int bnd_el_stiffmat(int belem_id,
                           double **Lk,
                           int *Ap,
                           Ai_t *Ai,
                           long ndofn,
                           Element *elem,
                           BoundingElement *b_elems,
                           Node *node,
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
                           SparseComm *comm,
                           long *Ddof,
                           int interior,
                           const int analysis,
                           SparseSystem *system,
                           const int mp_id)
{
  int err = 0;
  const BoundingElement *ptr_be = &b_elems[belem_id];
  const Element *ptr_ve = &elem[ptr_be->vol_elem_id];
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
  assert(0 < ndof_ve and ndof_ve <= std::numeric_limits<int>::max());

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
    PGFEM_free(sup_def);
  } else {
    def_elem(cn_ve,ndof_ve,d_r,NULL,NULL,v_disp,sup,0);
  }

  if(analysis == DISP){ /* TOTAL LAGRANGIAN formulation */
    double *ve_n = aloc1(ndof_ve);
    def_elem(cn_ve,ndof_ve,r,NULL,NULL,ve_n,sup,1);
    vvplus(v_disp,ve_n,ndof_ve);
    PGFEM_free(ve_n);
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
    PLoc_Sparse(Lk,lk,Ai,Ap,cn_ve,Gcn_ve,ndof_ve,Ddof,
                GDof,myrank,nproc,comm,interior,system,analysis);
  }


  PGFEM_free(x);
  PGFEM_free(y);
  PGFEM_free(z);

  PGFEM_free(cn_ve);
  PGFEM_free(Gcn_ve);

  PGFEM_free(v_disp);
  PGFEM_free(lk);

  return err;
} /* Bounding element stiffnes matrix */

void compute_Neumann_boundary_conditions_stiffnes_mat(FEMLIB *fe,
                                                      double *lk,
                                                      const double *r_e,
                                                      Grid *grid,
                                                      FieldVariables *fv,
                                                      int mp_id,
                                                      double dt,
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
  double dt_1_minus_alpha = -1.0;
    
  // displacement at t[n] and t[n-1]
  Matrix<double> Un, Unm1;  
  
  if(include_inertia){ 
    dt_1_minus_alpha = dt*alpha*(1.0-alpha);
    
    Un.initialization(fe->nne*nsd, 1, 0.0);
    
    for(int ia=0; ia<fe->nne; ++ia){
      for(int ib=0; ib<nsd; ++ib){
        Un(  ia*nsd+ib) = fv->u_n[  fe->node_id(ia)*ndofn+ib];
      }
    }
  }

  Constitutive_model *m = fv->eps[eid].model;

  for(int iA = 0; iA<bnd_elem_no; ++iA){
    const int face_id    = grid->NBE(mp_id).bnd_elements(this_bnd_elem_id)(iA, 0);
    const int feature_id = grid->NBE(mp_id).bnd_elements(this_bnd_elem_id)(iA, 1);
    const int load_type  = grid->NBE(mp_id).load_type(feature_id); // same as number of boundary loads

    double Pnp1[3] = {}, Pn[3] = {};

    // compute pressure at t_np1
    for(int ia=0; ia<load_type; ++ia){
      Pnp1[ia] = grid->NBE(mp_id).load_values_np1(feature_id, ia);
    }
        
    // compute pressure at t[n] if transient

    if(include_inertia){                
      for(int ia=0; ia<load_type; ++ia)
        Pn[  ia] = grid->NBE(mp_id).load_values_n(feature_id, ia);
    }
    
    if(load_type == 1){ // pressure bounary condition
      for(int ia=2; ia<nsd; ++ia){
        Pnp1[ia] = Pnp1[0];
        Pn[  ia] = Pn[  0];
      }
    }    

    // start boundary element FEM integration for
    FemLibBoundary fes(fe, face_id, fe->intg_order);    

    // apply quadrature rule 
    for(int ip=0; ip<fes.nint; ++ip){      
      fes.elem_basis_S(ip);
      fes.update_shape_tensor();
      
      // compute pressure terms at t(n+1), t(n)
      double PNnp1[3] = {}, PNn[3] = {};
      if(load_type == 1){ // pressure boundary condition
        TensorA<1> N(fes.normal);
        for(int ia=0; ia<nsd; ++ia){
          PNnp1[ia] = Pnp1[ia]*fes.normal[ia];
          PNn[  ia] = Pn[  ia]*fes.normal[ia];
        }
      } else{             // traction boundary condition
        for(int ia=0; ia<nsd; ++ia){
          PNnp1[ia] = Pnp1[ia];
          PNn[  ia] = Pn[  ia];
        }
      }
      
      // compute deformation gradient * pressure at t(n+a)
      Tensor<2> FnpaI = {};
      Tensor<1> PNnpa = {};
      double Jnpa = {};
      
      if(include_inertia){

        mid_point_rule(PNnpa.data, PNn,   Pnp1, alpha, nsd);
        Matrix<double> Unpa(fe->nne*nsd, 1, 0.0);
        
        mid_point_rule(Unpa.m_pdata,   Un.m_pdata, r_e, alpha, fe->nne*nsd);
        Tensor<2> Fnpa = {}, FnpaI = {};
        fes.update_deformation_gradient(nsd, Unpa.m_pdata, Fnpa.data, (m->param)->pFI);

        Jnpa = ttl::det(Fnpa);
        FnpaI = ttl::inverse(Fnpa);         
      } else{       
        Tensor<2> Fnpa = {}, FnpaI = {};
        for(int ia=0; ia<nsd; ++ia)
          PNnpa.data[ia] = PNnp1[ia];
          
        fes.update_deformation_gradient(nsd, r_e, Fnpa.data, (m->param)->pFI);
        Jnpa = ttl::det(Fnpa);
        FnpaI = ttl::inverse(Fnpa);
      }

      for(int ia = 0; ia<fes.nne_bnd; ++ia){
        int a_nid = fes.Volume2Boundary[ia];
        for(int ic=0; ic<fes.nne_bnd; ++ic){
          int c_nid = fes.Volume2Boundary[ic];
          for(int id=0; id<nsd; ++id){
            const int id_cd = idx_4_gen(ic,id,0,0,fe->nne,nsd,nsd,nsd);
            TensorA<2> ST_cd((fes.ST)+id_cd);
            Tensor<4> dFIdF = {};
            dFIdF(i,j,k,l) = -FnpaI(j,k)*FnpaI(l,i);
            Tensor<1> FIT = (FnpaI(l,k)*ST_cd(k,l)*FnpaI(j,i) - dFIdF(i,j,k,l)*ST_cd(k,l))*PNnpa(j);

            const int lk_idx = idx_K(a_nid,id,c_nid,id,fe->nne,nsd);              
            lk[lk_idx] += fes.N(a_nid)*fes.detJxW*dt_1_minus_alpha*Jnpa*FIT.data[id];
          }
        }
      }
    }
  }  
}

/// compute element stiffness matrix
///
/// Computes an element stiffness matrices. Each type of analysis will
/// be excuted in this function. Contribution to the residual
/// calculated using this function when prescribed displacement is
/// applied.
///
/// \param[in] fe finite element object
/// \param[out] lk computed local(element level) stiffness matrix
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step
/// \param[in] lm Load multiplier level in Arc Length scheme
/// \param[in] be tangential load vector when periodic and solution scheme is Arc Length
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \return non-zero on internal error
int el_compute_stiffmat_MP(FEMLIB *fe,
                           double *lk,
                           Grid *grid,
                           MaterialProperty *mat,
                           FieldVariables *fv,
                           Solver *sol,
                           LoadingSteps *load,
                           CRPL *crpl,
                           const PGFem3D_opt *opts,
                           const Multiphysics& mp,
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

  SUPP sup = load->sups[mp_id];

  double *x = (fe->temp_v).x.m_pdata;
  double *y = (fe->temp_v).y.m_pdata;
  double *z = (fe->temp_v).z.m_pdata;

  if(include_inertia)
    err += stiffness_with_inertia(fe,lk,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);
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
      stiffmat_3f_el(fe,lk,r_e,grid,mat,fv,-1,dt);
      break;
     case CM:  // intened to flow
     case CM3F:
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

  compute_Neumann_boundary_conditions_stiffnes_mat(fe, lk, r_e,
                                                   grid, fv, mp_id, dt, 
                                                   include_inertia,sol->alpha);
  if(err>0)
  {
    PGFEM_printf("Error is detected: eid = %d\n", eid);
  }
  return err;
}


/// Prepare finite element objects for computing stiffness matrix
///
/// Based on the element type, this function allocates a FEMLIB object,
/// and call actual function to compute stiffness matrix.
///
/// \param[in] eid element id
/// \param[out] Lk computed stiffness matrix
/// \param[in] Ddof shifted degree of freedom (global end dof id in each domain)
/// \param[in] interior if 1, local element
///                        0, element on the communication boundary
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com communication object
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] iter number of Newton Raphson interataions
/// \return non-zero on internal error
static int el_stiffmat_MP(int eid,
                          double **Lk,
                          long *Ddof,
                          int interior,
                          Grid *grid,
                          MaterialProperty *mat,
                          FieldVariables *fv,
                          Solver *sol,
                          LoadingSteps *load,
                          const CommunicationStructure *com,
                          CRPL *crpl,
                          const PGFem3D_opt *opts,
                          const Multiphysics& mp,
                          int mp_id,
                          double dt,
                          long iter)
{
  int err = 0;
  double lm = 0.0;
  int intg_order = 0;
  int myrank = com->rank;

  SUPP sup = load->sups[mp_id];

  // multi-scale simulation, displacement base and three-fied-mixed (elastic only)
  // are total Lagrangian based
  int total_Lagrangian = 0;
  switch(opts->analysis_type)
  {
   case DISP: // intended to flow
   case TF:
    total_Lagrangian = 1;
    break;
   case CM:   // intended to flow
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
    for(int j=0; j<sup->npd; j++){
      assert(sup_def != NULL && "sup_def can't be null");
      sup->defl_d[j] = sup_def[j];
    }
  }

  Matrix<double> lk(ndofe,ndofe,0.0);

  err += el_compute_stiffmat_MP(&fe,lk.m_pdata,grid,mat,fv,sol,load,
                                crpl,opts,mp,mp_id,dt,lm,be,r_e);

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
               myrank,com->nproc,com->spc,interior,sol->system,opts->analysis_type);

  //  dealocation
  PGFEM_free (cnL);
  PGFEM_free (cnG);
  PGFEM_free (be);
  PGFEM_free (r_e);
  PGFEM_free (sup_def);

  return err;

}


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
/// \param[in,out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com communication object
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dt time step
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int stiffmat_fd_MP(Grid *grid,
                   MaterialProperty *mat,
                   FieldVariables *fv,
                   Solver *sol,
                   LoadingSteps *load,
                   const CommunicationStructure *com,
                   CRPL *crpl,
                   const PGFem3D_opt *opts,
                   const Multiphysics& mp,
                   int mp_id,
                   double dt,
                   long iter)
{
  int myrank = com->rank;
  int err = 0;
  // if 0, update only stiffness
  double lm = 0.0;

  double **Lk,**receive;

  try {
    com->spc->post_stiffmat(&Lk,&receive);
  } catch(...) {
    err++;
  }

  Matrix<long> Ddof(com->nproc,1);

  Ddof.m_pdata[0] = com->DomDof[0];
  for (int ia=1; ia<com->nproc; ia++)
    Ddof.m_pdata[ia] = Ddof.m_pdata[ia-1] + com->DomDof[ia];

  for(int eid=0; eid<com->nbndel; eid++)
  {
    err += el_stiffmat_MP(com->bndel[eid],Lk,Ddof.m_pdata,0,grid,mat,fv,sol,load,com,crpl,
              opts,mp,mp_id,dt,iter);

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
                      opts->stab,grid->coel,0,0,fv->f_u,com->rank,com->nproc,com->DomDof,
                      com->GDof,com->spc,Ddof.m_pdata,0,opts->analysis_type,sol->system, mp_id);
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
    // BoundingElement *b_elems = NULL;
    for(int eid=0; eid<grid->n_be; eid++){
      err += bnd_el_stiffmat(eid,Lk,com->Ap,com->Ai,fv->ndofn,grid->element,grid->b_elems,grid->node,mat->hommat,
                             mat->matgeom,fv->sig,fv->eps,fv->d_u,fv->u_np1,fv->npres,load->sups[mp_id],iter,sol->nor_min,
                             dt,crpl,opts->stab,sol->FNR,lm,fv->f_u,com->rank,com->nproc,com->GDof,
                             com->spc,Ddof.m_pdata,0,opts->analysis_type,sol->system,mp_id);

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
      if(send_proc==myrank || (com->spc)->S[send_proc] ==0)
        continue;

      for(int n_data=0; n_data<(com->spc)->AS[send_proc]; n_data++)
      {
        PGFEM_fprintf(out,"%12.12e ",Lk[send_proc][n_data]);
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n");
    fclose(out);
  }

  /// SEND BOUNDARY AND COEL
  try {
    com->spc->send_stiffmat();
  } catch(...) {
    err++;
  }

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
    err += el_stiffmat_MP(eid,Lk,Ddof.m_pdata,1,grid,mat,fv,sol,load,com,crpl,
                          opts,mp,mp_id,dt,iter);

    if(err != 0)
      break;
  }

  try {
    com->spc->assemble_nonlocal_stiffmat(sol->system);
    com->spc->finalize_stiffmat();
  } catch(...) {
    err++;
  }

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
int stiffmat_fd_multiscale(MultiscaleCommon *c,
                           MULTISCALE_SOLUTION *s,
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

  // initialize and define mesh object
  Grid grid;
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
    fv.BS_x   = s->BS_x;
    fv.BS_f   = s->BS_f;
    fv.BS_RR  = s->BS_RR;
    fv.BS_f_u = s->BS_f_u;
  }

  /// initialize and define iterative solver object
  Solver sol{};
  {
    sol.FNR     = FNR;
    sol.system  = c->SOLVER;
    sol.err     = c->lin_err;
    sol.alpha   = 0.0;
    sol.nor_min = nor_min;
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

  /// initialize and define communication structures
  CommunicationStructure com;
  {
    communication_structure_initialization(&com);
    com.Ap     = c->Ap;
    com.Ai     = c->Ai;
    com.DomDof = c->DomDof;
    com.GDof   = c->GDof;
    com.nbndel = c->nbndel;
    com.bndel  = c->bndel;
    com.boot   = c->boot;
    com.net    = c->net;
    com.comm   = c->comm;
    com.rank   = c->rank;
    com.nproc  = c->nproc;
    com.spc    = c->spc;
  }

  err += stiffmat_fd_MP(&grid,&mat,&fv,&sol,&load,&com,s->crpl,
                        opts,mp,mp_id,s->dt,iter);

  free(physicsname);

  return err;
}
