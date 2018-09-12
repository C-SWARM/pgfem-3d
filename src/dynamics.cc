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

void MMS_body_force(double *b, const HOMMAT  *hommat, ELASTICITY *elast, double t, double X, double Y, double Z, const bool is4cm)
{
  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
}

#else
#define INTG_ORDER 1
#include "../verification_MMS/MMS.h"
#endif

void DISP_resid_body_force_el(double *f,
                              const int ii,
                              const int ndofn,
                              const int nne,
                              const double *x,
                              const double *y,
                              const double *z,
                              const Element *elem,
                              const HOMMAT *hommat,
                              const Node *node,
                              double dt, 
                              double t,
                              ELASTICITY *elast,
                              bool is4cm)
{
  const int mat = elem[ii].mat[2];
  // double rho = hommat[mat].density;
  int ndofe = nne*ndofn;

  

  /* make sure the f vector contains all zeros */
  memset(f,0,ndofe*sizeof(double));

  FEMLIB fe(ii, elem, node, INTG_ORDER,1);

  double *bf = aloc1(ndofn);

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);
    double X[3];
    X[0] = X[1] = X[2] = 0.0;

    memset(bf,0,ndofn*sizeof(double));

    for(long a = 0; a<nne; a++)
    {
      X[0] += fe.N(a+1)*x[a];
      X[1] += fe.N(a+1)*y[a];
      X[2] += fe.N(a+1)*z[a];
    }

    MMS_body_force(bf, &hommat[mat], elast, t,  X[0], X[1], X[2], is4cm);

    for(long a = 0; a<nne; a++)
    {
      for(long b=0; b<3; b++)
      {
        long id = a*ndofn + b;
        f[id] += bf[b]*fe.N(a+1)*fe.detJxW;
      }
    }
  }
  dealoc1(bf);
}

void DISP_resid_w_inertia_el(FEMLIB *fe,
                             double *f,
                             const int ii,
                             const int ndofn,
                             const int nne,
                             const double *x,
                             const double *y,
                             const double *z,
                             const Element *elem,
                             const HOMMAT *hommat,
                             const Node *node, const double *dts, double t,
                             double *r_2, double* r_1, double *r_0, double alpha,
                             ELASTICITY *elast,
                             const bool is4cm)
{
  const int mat = elem[ii].mat[2];
  double rho = hommat[mat].density;
  int ndofe = nne*ndofn;

  /* make sure the f vector contains all zeros */
  memset(f,0,ndofe*sizeof(double));

  Matrix<double> du(3,1);

  double *bf0, *bf1, *bf2, *bf_n1a, *bf;
  bf0 = aloc1(ndofn);
  bf1 = aloc1(ndofn);
  bf2 = aloc1(ndofn);
  bf_n1a = aloc1(ndofn);
  bf     = aloc1(ndofn);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    fe->elem_basis_V(ip);

    du.set_values(0.0);
    double X[3];
    X[0] = X[1] = X[2] = 0.0;

    memset(bf0,   0,ndofn*sizeof(double));
    memset(bf1,   0,ndofn*sizeof(double));
    memset(bf2,   0,ndofn*sizeof(double));
    memset(bf_n1a,0,ndofn*sizeof(double));
    memset(bf,    0,ndofn*sizeof(double));

    for(long a = 0; a<nne; a++)
    {
      X[0] += fe->N(a+1)*x[a];
      X[1] += fe->N(a+1)*y[a];
      X[2] += fe->N(a+1)*z[a];
      for(long b = 0; b<3; b++)
      {
        long id = a*ndofn + b;
        du(b+1) += fe->N(a+1)*(dts[DT_N]*r_2[id]-(dts[DT_NP1]+dts[DT_N])*r_1[id]+dts[DT_NP1]*r_0[id]);
      }
    }

    double t1 = t-dts[DT_NP1];
    double t0 = t - dts[DT_NP1] - dts[DT_N];

    //    if(t0>=0)
    MMS_body_force(bf0, &hommat[mat], elast, t0, X[0], X[1], X[2], is4cm);

    //    if(t1>=0)
    MMS_body_force(bf1, &hommat[mat], elast, t1, X[0], X[1], X[2], is4cm);

    MMS_body_force(bf2, &hommat[mat], elast, t,  X[0], X[1], X[2], is4cm);

    mid_point_rule(bf_n1a, bf0, bf1, alpha, ndofn);
    mid_point_rule(bf, bf1, bf2, alpha, ndofn);

    for(long a = 0; a<nne; a++)
    {
      for(long b=0; b<3; b++)
      {
        long id = a*ndofn + b;
        f[id] += rho/dts[DT_NP1]/dts[DT_N]*fe->N(a+1)*du(b+1)*fe->detJxW;
        f[id] -= (1.0-alpha)*dts[DT_NP1]*bf[b]*fe->N(a+1)*fe->detJxW;
        f[id] -= alpha*dts[DT_N]*bf_n1a[b]*fe->N(a+1)*fe->detJxW;
      }
    }
  }

  dealoc1(bf0);
  dealoc1(bf1);
  dealoc1(bf2);
  dealoc1(bf_n1a);
  dealoc1(bf);
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
  ELASTICITY *elast = NULL;
  
  if(opts->analysis_type == CM || opts->analysis_type == CM3F){
    elast = (fv->eps[eid].model[0].param)->cm_elast;
    if(fv->eps[eid].model[0].param->type == MANUFACTURED_SOLUTIONS)
      is4cm = true;
  }
  DISP_resid_w_inertia_el(fe,f_i,eid,ndofn,nne,x,y,z,grid->element,mat->hommat,grid->node,dts,t,r_e, r0, r0_, sol->alpha, elast, is4cm);

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

       free(f_n_a);
       free(f_n_1_a);

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
       err += residuals_el_constitutive_model_w_inertia(fe,be,r_e,r_n_a,r_n_1_a,grid,mat,fv,sol,load,crpl,opts,mp,dts,mp_id,dts[DT_NP1],t);

       for(long a = 0; a<ndofe; a++)
         be[a] -= f_i[a]; // - (1.0-alpha)*dt and - alpha*dt are included in be[a]
         
       break;
     }
   default:
    PGFEM_printf("Only displacement based element and three field element are supported\n");
    break;
  }

  free(r_n_a);
  free(r_n_1_a);
  free(r0);
  free(r0_);
  free(f_i);
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
      u_n(I*ndofn+J+1,1) = fv->u_n[nod[I]*ndofn + J];
  }

  mid_point_rule(u_npa.m_pdata, u_n.m_pdata, r_e, sol->alpha, ndofe);

  if(opts->analysis_type == DISP ||
     opts->analysis_type == TF   ||
     opts->analysis_type == CM   ||
     opts->analysis_type == CM3F)
  {
    for(int ip = 1; ip<=fe->nint; ip++)
    {
      fe->elem_basis_V(ip);

      for(long a = 0; a<nne; a++)
      {
        for(long c=0; c<nne; c++)
        {
          for(long b=1; b<=ndofn; b++)
            Kuu_I(a*ndofn+b,c*ndofn+b) += rho/dt*fe->N(a+1)*fe->N(c+1)*fe->detJxW;
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
