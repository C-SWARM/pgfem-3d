#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "dynamics.h"
#include "allocation.h"
#include "constitutive_model.h"
#include "displacement_based_element.h"
#include "enumerations.h"
#include "femlib.h"
#include "three_field_element.h"
#include "utils.h"

using pgfem3d::Solver;

#ifndef VERIFICATION_USING_MMS
#define INTG_ORDER 0

void MMS_body_force(double *b, const HOMMAT  *hommat, HyperElasticity *elast, double t, double X, double Y, double Z, const bool is4cm)
{
  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
}

#else
#define INTG_ORDER 0
#include "../verification_MMS/MMS.h"
#endif

void DISP_resid_body_force_el(double *f,
                              const int ii,
                              const int ndofn,
                              const int nne,
                              const Element *elem,
                              const HOMMAT *hommat,
                              const Node *node,
                              double dt, 
                              double t,
                              HyperElasticity *elast,
                              bool is4cm)
{
  const int mat = elem[ii].mat[2];
  // double rho = hommat[mat].density;
  int ndofe = nne*ndofn;


  /* make sure the f vector contains all zeros */
  memset(f,0,ndofe*sizeof(double));

  FEMLIB fe(ii, elem, node, INTG_ORDER,1);

  Matrix<double> bf(ndofn, 1, 0.0);
  
  for(int ip = 0; ip<fe.nint; ip++)
  {
    fe.elem_basis_V(ip);
    bf.set_values(0.0);
    MMS_body_force(bf.m_pdata, &hommat[mat], elast, t,  fe.x_ip(0), fe.x_ip(1), fe.x_ip(2), is4cm);

    for(long a = 0; a<nne; a++)
    {
      for(long b=0; b<3; b++)
      {
        long id = a*ndofn + b;
        f[id] += bf(b)*fe.N(a)*fe.detJxW;
      }
    }
  }
}

void DISP_resid_w_inertia_el(FEMLIB *fe,
                             double *f,
                             const int ii,
                             const int ndofn,
                             const int nne,
                             const Element *elem,
                             const HOMMAT *hommat,
                             const Node *node, const double *dts, double t,
                             double *r_2, double* r_1, double *r_0, double alpha,
                             HyperElasticity *elast,
                             const bool is4cm)
{
  const int mat = elem[ii].mat[2];
  double rho = hommat[mat].density;
  int ndofe = nne*ndofn;
  
  #ifdef VERIFICATION_USING_MMS
    if(!is4cm){
      MATERIAL_ELASTICITY mat_e;
      elast = new HyperElasticity;
      set_properties_using_E_and_nu(&mat_e,hommat[mat].E,hommat[mat].nu);
      mat_e.m01 = hommat[mat].m01;
      mat_e.m10 = hommat[mat].m10;
      mat_e.G   = hommat[mat].G;
      mat_e.kappa      = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));
      mat_e.devPotFlag = hommat[mat].devPotFlag;
      mat_e.volPotFlag = hommat[mat].volPotFlag;
      elast->construct_elasticity(&mat_e, true);
    }
  #endif    

  /* make sure the f vector contains all zeros */
  memset(f,0,ndofe*sizeof(double));

  Matrix<double> du(3,1);

  Matrix<double> bf0(ndofn, 1, 0.0), bf1(ndofn, 1, 0.0), bf2(ndofn, 1, 0.0);
  Matrix<double> bf_npa(ndofn, 1, 0.0), bf_nma(ndofn, 1, 0.0); 
    

  for(int ip = 0; ip<fe->nint; ++ip)
  {
    fe->elem_basis_V(ip);

    du.set_values(0.0);       
    for(int ia=0; ia<ndofn; ++ia)
      bf0(ia) = bf1(ia) = bf2(ia) = bf_npa(ia) = bf_nma(ia) = 0.0;   

    for(long a = 0; a<nne; a++){
      for(long b = 0; b<3; b++){
        long id = a*ndofn + b;
        du(b) += fe->N(a)*(dts[DT_N]*r_2[id]-(dts[DT_NP1]+dts[DT_N])*r_1[id]+dts[DT_NP1]*r_0[id]);
      }
    }

    double t1 = t-dts[DT_NP1];
    double t0 = t - dts[DT_NP1] - dts[DT_N];

    //    if(t0>=0)
    MMS_body_force(bf0.m_pdata, &hommat[mat], elast, t0, fe->x_ip(0), fe->x_ip(1), fe->x_ip(2), is4cm);

    //    if(t1>=0)
    MMS_body_force(bf1.m_pdata, &hommat[mat], elast, t1, fe->x_ip(0), fe->x_ip(1), fe->x_ip(2), is4cm);

    MMS_body_force(bf2.m_pdata, &hommat[mat], elast, t,  fe->x_ip(0), fe->x_ip(1), fe->x_ip(2), is4cm);

    mid_point_rule(bf_npa.m_pdata, bf1.m_pdata, bf2.m_pdata, alpha, ndofn);
    mid_point_rule(bf_nma.m_pdata, bf0.m_pdata, bf1.m_pdata, alpha, ndofn);    

    for(long a = 0; a<nne; a++)
    {
      for(long b=0; b<3; b++)
      {
        long id = a*ndofn + b;
        f[id] += rho/dts[DT_NP1]/dts[DT_N]*fe->N(a)*du(b)*fe->detJxW;
        f[id] -= (1.0-alpha)*dts[DT_NP1]*bf_npa(b)*fe->N(a)*fe->detJxW;
        f[id] -=       alpha*dts[DT_N  ]*bf_nma(b)*fe->N(a)*fe->detJxW;
      }
    }
  }
  
  #ifdef VERIFICATION_USING_MMS
    if(!is4cm){
      delete elast;
    }
  #endif   
}

/// compute element residual vector in transient
///
/// Displacement based function is used to compute acceleration term and each type of analysis will
/// be excuted and added their constribution to the residual.
/// CAVEAT Acceleration term is based on the total Lagrangian, but other part support
/// total and updated Lagrangian.
///
/// \param[in] fe finite element helper object
/// \param[out] be computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1)
///                                                dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] t current time
/// \return non-zero on internal error
int residual_with_inertia(FEMLIB *fe,
                          double *be,
                          double *r_e,
                          Grid *grid,
                          MaterialProperty *mat,
                          FieldVariables *fv,
                          Solver *sol,
                          LoadingSteps *load,
                          CRPL *crpl,
                          const PGFem3D_opt *opts,
                          const Multiphysics& mp,
                          int mp_id,
                          double *dts,
                          double t,
                          int myrank)
{
  int err = 0;

  if(!(opts->analysis_type == DISP
       || opts->analysis_type == TF
       || opts->analysis_type == CM
       || opts->analysis_type == CM3F))
  {
    if(myrank==0)
      PGFEM_printf("With inertia is supported only for analysis type = DISP or TF or CM\n");

    return err;
  }

  int nsd = 3;
  int eid = fe->curt_elem_id;

  // const int mat_id = grid->element[eid].mat[2];
  // double rho = mat->hommat[mat_id].density;

  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int ndofe = nne*ndofn;

  double *r_n_a   = aloc1(ndofe);
  double *r_n_1_a = aloc1(ndofe);
  double *r0      = aloc1(ndofe);
  double *r0_     = aloc1(ndofe);
  double *f_i     = aloc1(ndofe);
  memset(be,  0, sizeof(double)*ndofe);
  memset(f_i, 0, sizeof(double)*ndofe);

  const long *nod = (fe->node_id).m_pdata;

  for(int I=0;I<nne;I++)
  {
    for(int J=0; J<nsd; J++)
    {
      r0[I*ndofn + J] =   fv->u_n[nod[I]*ndofn + J];
      r0_[I*ndofn + J] = fv->u_nm1[nod[I]*ndofn + J];
    }
  }

  mid_point_rule(r_n_1_a,r0_,r0,  sol->alpha, ndofe);
  mid_point_rule(r_n_a,  r0, r_e, sol->alpha, ndofe);

  double *x = (fe->temp_v).x.m_pdata;
  double *y = (fe->temp_v).y.m_pdata;
  double *z = (fe->temp_v).z.m_pdata;
  SUPP sup = load->sups[mp_id];

  bool is4cm = false; 
  HyperElasticity *elast = NULL;
  
  if(opts->analysis_type == CM || opts->analysis_type == CM3F){
    elast = (fv->eps[eid].model[0].param)->cm_elast;
    if(fv->eps[eid].model[0].param->type == MANUFACTURED_SOLUTIONS)
      is4cm = true;
  }
  DISP_resid_w_inertia_el(fe, f_i,eid,ndofn,nne,grid->element,mat->hommat,grid->node,dts,t,r_e, r0, r0_, sol->alpha, elast, is4cm);

  switch(opts->analysis_type)
  {
   case DISP:
     {
       double *f_n_a   = aloc1(ndofe);
       double *f_n_1_a = aloc1(ndofe);

       err +=  DISP_resid_el(f_n_1_a,eid,ndofn,nne,x,y,z,grid->element,
                             mat->hommat,nod,grid->node,fv->eps,fv->sig,sup,r_n_1_a, dts[DT_N]);

       err +=  DISP_resid_el(f_n_a,eid,ndofn,nne,x,y,z,grid->element,
                             mat->hommat,nod,grid->node,fv->eps,fv->sig,sup,r_n_a, dts[DT_NP1]);

       for(long a = 0; a<ndofe; a++)
         be[a] = -f_i[a] - (1.0-(sol->alpha))*dts[DT_NP1]*f_n_a[a] - (sol->alpha)*dts[DT_N]*f_n_1_a[a];

       PGFEM_free(f_n_a);
       PGFEM_free(f_n_1_a);

       break;
     }
   case TF:
     {
       if(0<(sol->alpha) && (sol->alpha)<1.0)
         residuals_3f_w_inertia_el(fe,be,r_e,grid,mat,fv,sol,dts,r_n_1_a,r_n_a);

       for(long a = 0; a<ndofe; a++)
         be[a] -= f_i[a];

       break;

     }
   case CM:   //intended to flow
   case CM3F:
     {
       err += residuals_el_constitutive_model_w_inertia(fe,be,r_e,r_n_a,r_n_1_a,grid,mat,fv,sol,load,crpl,opts,mp,dts,mp_id,t);

       for(long a = 0; a<ndofe; a++)
         be[a] -= f_i[a]; // - (1.0-alpha)*dt and - alpha*dt are included in be[a]
         
       break;
     }
   default:
    PGFEM_printf("Only displacement based element and three field element are supported\n");
    break;
  }

  PGFEM_free(r_n_a);
  PGFEM_free(r_n_1_a);
  PGFEM_free(r0);
  PGFEM_free(r0_);
  PGFEM_free(f_i);
  return err;
}

/// compute element stiffness matrix in transient
///
/// Displacement based function is used to compute acceleration term and each type of analysis will
/// be excuted and added their constribution to the residual.
/// CAVEAT Acceleration term is based on the total Lagrangian, but other part support
/// total and updated Lagrangian.
///
/// \param[in] fe finite element helper object
/// \param[out] Ke computed computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_with_inertia(FEMLIB *fe,
                           double *Ks,
                           double *r_e,
                           Grid *grid,
                           MaterialProperty *mat,
                           FieldVariables *fv,
                           Solver *sol,
                           LoadingSteps *load,
                           CRPL *crpl,
                           const PGFem3D_opt *opts,
                           const Multiphysics& mp,
                           int mp_id,
                           double dt)
{

  int err = 0;
  int eid = fe->curt_elem_id;

  const int mat_id = grid->element[eid].mat[2];
  double rho = mat->hommat[mat_id].density;

  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int ndofe = nne*ndofn;

  Matrix<double> Kuu_K,Kuu_I, u_npa, u_n;
  Kuu_I.initialization(ndofe,ndofe,0.0);
  Kuu_K.initialization(ndofe,ndofe,0.0);
  u_npa.initialization(ndofe,1,0.0);
  u_n.initialization(ndofe,1,0.0);

  /* make sure the stiffenss matrix contains all zeros */
  memset(Ks,0,ndofe*ndofe*sizeof(double));

  const long *nod = (fe->node_id).m_pdata;
  for (long I=0;I<nne;I++)
  {
    for(long J=0; J<ndofn; J++)
      u_n(I*ndofn+J,0) = fv->u_n[nod[I]*ndofn + J];
  }

  mid_point_rule(u_npa.m_pdata, u_n.m_pdata, r_e, sol->alpha, ndofe);

  if(opts->analysis_type == DISP ||
     opts->analysis_type == TF   ||
     opts->analysis_type == CM   ||
     opts->analysis_type == CM3F)
  {
    for(int ip = 0; ip<fe->nint; ip++)
    {
      fe->elem_basis_V(ip);

      for(long a = 0; a<nne; a++)
      {
        for(long c=0; c<nne; c++)
        {
          for(long b=0; b<ndofn; b++)
            Kuu_I(a*ndofn+b,c*ndofn+b) += rho/dt*fe->N(a)*fe->N(c)*fe->detJxW;
        }
      }
    }
  }

  double *x = (fe->temp_v).x.m_pdata;
  double *y = (fe->temp_v).y.m_pdata;
  double *z = (fe->temp_v).z.m_pdata;

  SUPP sup = load->sups[mp_id];

  switch(opts->analysis_type)
  {
   case DISP:
    err = DISP_stiffmat_el(Kuu_K.m_pdata,eid,ndofn,nne,x,y,z,grid->element,
                           mat->hommat,nod,grid->node,fv->eps,fv->sig,sup,u_npa.m_pdata,dt);

    for(long a = 0; a<ndofe*ndofe; a++)
      Ks[a] = -Kuu_I.m_pdata[a]-(sol->alpha)*(1.0-(sol->alpha))*dt*Kuu_K.m_pdata[a];
    
    break;

   case TF:
    if(0<sol->alpha && sol->alpha<1.0)
    {
      stiffmat_3f_el(fe,Kuu_K.m_pdata,u_npa.m_pdata,grid,mat,fv,sol->alpha,dt);
    }
    for(long a = 0; a<ndofe*ndofe; a++)
      Ks[a] = -Kuu_I.m_pdata[a] + Kuu_K.m_pdata[a];

    break;
   case CM:  // intended to flow
     err += stiffness_el_constitutive_model_w_inertia(fe,Kuu_K.m_pdata,r_e,u_npa.m_pdata,
                                                      grid,mat,fv,sol,load,crpl,
                                                      opts,mp,mp_id,dt);

     for(long a = 0; a<ndofe*ndofe; a++)
       Ks[a] = -Kuu_I.m_pdata[a]-(sol->alpha)*(1.0-(sol->alpha))*dt*Kuu_K.m_pdata[a];

     break;
   case CM3F:
     err += stiffness_el_constitutive_model_w_inertia(fe,Kuu_K.m_pdata,r_e,u_npa.m_pdata,
                                                      grid,mat,fv,sol,load,crpl,
                                                      opts,mp,mp_id,dt);

     for(long a = 0; a<ndofe*ndofe; a++)
       Ks[a] = -Kuu_I.m_pdata[a] + Kuu_K.m_pdata[a];

     break;     

   default:
    PGFEM_printf("Only displacement based element and three field element are supported\n");
    break;
  }
  return err;
}
