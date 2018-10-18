/// This file defines the implementation for finite strain J2
/// plasticity with damage employing the stress split.
///
/// Authors:
///  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
///  Sangmin Lee, University of Notre Dame, <slee43@nd.edu>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "allocation.h"
#include "cm_j2_plasticity.h"
#include "cm_iso_viscous_damage.h"
#include "constitutive_model.h"
#include "index_macros.h"
#include "utils.h"
#include <math.h>
#include <string.h> 
#include <assert.h>
#include "J2_plasticity.h"
#include "continuum_damage_model.h"
#include <ttl/ttl.h>

// Define constant dimensions. Note cannot use `static const` with
// initialization list
namespace {
  template <int R, class S = double>
  using Tensor = ttl::Tensor<R, 3, S>;
    
  template <int R, class S = double *>
  using TensorA = ttl::Tensor<R, 3, S>;  
    
  template <int R>
  using Delta = ttl::Tensor<R, 3, double>;
  
  template<class T1, class T2> int inv(T1 &A, T2 &AI)
  {
    int err = inv3x3(A.data, AI.data);
    return err;
  }
  
  // ttl indexes
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'m'> m;
  static constexpr ttl::Index<'n'> n;
    
  constexpr int            DIM_3 = 3;
  constexpr int          DIM_3x3 = 9;
  constexpr int      DIM_3x3x3x3 = 81;
  constexpr double DAMAGE_THRESH = 0.9999;
  constexpr double   DELTA_W_MAX = 0.05;    
      
  enum {TENSOR_Fnm1,
        TENSOR_Fn, 
        TENSOR_Fnp1,
        TENSOR_pSnm1,
        TENSOR_pSn, 
        TENSOR_pSnp1,         
        TENSOR_end};
        
  enum {VAR_w_np1,
        VAR_H_np1,
        VAR_X_np1,
        VAR_ep_np1,
        VAR_gam_np1,
        VAR_w_n,
        VAR_X_n,
        VAR_H_n,
        VAR_ep_n,
        VAR_gam_n,        
        VAR_w_nm1,
        VAR_X_nm1,
        VAR_H_nm1,
        VAR_ep_nm1,
        VAR_gam_nm1,        
        VAR_end};
              
  enum {FLAG_damaged_n,
        FLAG_damaged,
        FLAG_end};
                
  enum {PARAM_G, 
        PARAM_nu, 
        PARAM_beta, 
        PARAM_hp, 
        PARAM_k0,
        PARAM_mu, 
        PARAM_w_max, 
        PARAM_P1, 
        PARAM_P2, 
        PARAM_Yin,
        PARAM_da,
        PARAM_db,
        PARAM_va,
        PARAM_vb,
        PARAM_NO};      
}

double j2d_compute_w_npa(const Constitutive_model *m,
                         const int npa,
                         const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_w_nm1, VAR_w_n, VAR_w_np1, npa, alpha);
}

double j2d_compute_H_npa(const Constitutive_model *m,
                         const int npa,
                         const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_H_nm1, VAR_H_n, VAR_H_np1, npa, alpha);
}

double j2d_compute_X_npa(const Constitutive_model *m,
                         const int npa,
                         const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_X_nm1, VAR_X_n, VAR_X_np1, npa, alpha);
}

double j2d_compute_ep_npa(const Constitutive_model *m,
                          const int npa,
                          const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_ep_nm1, VAR_ep_n, VAR_ep_np1, npa, alpha);
}

double j2d_compute_gam_npa(const Constitutive_model *m,
                           const int npa,
                           const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_gam_nm1, VAR_gam_n, VAR_gam_np1, npa, alpha);
}

int J2d_determine_damaged_npa(const Constitutive_model *m,
                              const int npa)
{
  int *flags = m->vars_list[0][m->model_id].flags;
  int flag_npa = flags[FLAG_damaged];
  switch(npa)
  {
    case 0:
      flag_npa = flags[FLAG_damaged_n];
      break;
    case 1:
      flag_npa = flags[FLAG_damaged];
      break;
  }
  return flag_npa;
}

// Y = G/2 *(tr(bbar) - 3) + kU 
double j2d_compute_Y0(ELASTICITY *elast,
                      double *F_in)
{
  Tensor<2> bbar;
  TensorA<2> F(F_in);
  double J = 0.0;
  J = det3x3(F_in);
  
  double J23 = pow(J,-2.0/3.0);
  bbar(i,j) = J23*F(i,k)*F(j,k);
  
  const double G = elast->mat->G;
  const double kappa = elast->mat->kappa;
  
  double U = 0.0;
  elast->compute_u(&U, J);

  return 0.5*G*(bbar.data[0] + bbar.data[4] + bbar.data[8] - 3.0) + kappa*U;
}

// see box 4 of Simo and Ju (1989)
int CM_J2P_PARAM::integration_algorithm(Constitutive_model *m,
                                        CM_Ctx &cm_ctx)
const
{
  int err = 0;

  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  int    *flags      = m->vars_list[0][m->model_id].flags;
      
  // store the deformation gradient
  memcpy(Fs[TENSOR_Fnp1].m_pdata, cm_ctx.F, DIM_3x3 * sizeof(double));
  
  err += J2_plasticity_integration_alg(Fs[TENSOR_pSnp1].m_pdata,
                                       vars + VAR_ep_np1,
                                       vars + VAR_gam_np1,
                                       Fs[TENSOR_Fnp1].m_pdata,
                                       Fs[TENSOR_Fn].m_pdata,
                                       Fs[TENSOR_pSn].m_pdata,
                                       vars[VAR_ep_n],
                                       this->cm_mat->mat_J2p,
                                       this->cm_mat->mat_e); 
                                       

  double Y0 = j2d_compute_Y0(this->cm_elast, Fs[TENSOR_Fnp1].m_pdata);
  err += continuum_damage_integration_alg_public((this->cm_mat)->mat_d,this->cm_elast,
                                                  vars+VAR_w_np1, vars+VAR_X_np1, vars+VAR_H_np1, flags + FLAG_damaged,                                          
                                                  vars[VAR_w_n], vars[VAR_X_n], cm_ctx.dt, Y0);  
  
  return err;
}

int CM_J2P_PARAM::compute_dev_stress(const Constitutive_model *m,
                                     CM_Ctx &cm_ctx,
                                     double *stress)
const
{
  return 0;
}

int CM_J2P_PARAM::compute_dudj(const Constitutive_model *m,
                               CM_Ctx &cm_ctx,
                               double *dudj)
const
{
  int err = 0;
  double eJ = det3x3(cm_ctx.F);
  
  double *vars  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;

  ELASTICITY *elasticity = this->cm_elast;
  *dudj = damage_compute_dudj(elasticity, eJ, vars[VAR_w_np1]);
  
  return err;
}

double CM_J2P_PARAM::compute_dudj(const Constitutive_model *m,
                                  double theta_e,
                                  const int npa,
                                  const double alpha)
const 
{
  double w_npa = j2d_compute_w_npa(m,npa,alpha);
  ELASTICITY *elasticity = this->cm_elast;
  return damage_compute_dudj(elasticity, theta_e, w_npa);          
}

int CM_J2P_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                      CM_Ctx &cm_ctx,
                                      double *L)
const
{
  return 0;
}

int CM_J2P_PARAM::compute_d2udj2(const Constitutive_model *m,
                                 CM_Ctx &cm_ctx,
                                 double *d2udj2)
const
{
  int err = 0;
  double eJ = det3x3(cm_ctx.F);
  
  double *vars  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;

  // scale by damage variable 
  ELASTICITY *elasticity = this->cm_elast;
  *d2udj2 = damage_compute_d2udj2(elasticity, eJ, vars[VAR_w_np1]);
  
  return err;
}

double CM_J2P_PARAM::compute_d2udj2(const Constitutive_model *m,
                                    double theta_e,
                                    const int npa,
                                    const double alpha)
const 
{
  double w_npa = j2d_compute_w_npa(m,npa,alpha);
  ELASTICITY *elasticity = this->cm_elast;
  return damage_compute_d2udj2(elasticity, theta_e, w_npa);           
}

int CM_J2P_PARAM::update_elasticity(const Constitutive_model *m,
                                    CM_Ctx &cm_ctx,
                                    double *L_in,
                                    double *S,
                                    const int compute_stiffness)
const
{
  int err = 0;
  
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;

  ELASTICITY *elasticity = this->cm_elast;
  
  double *S_tmp = elasticity->S;
  double *L_tmp = elasticity->L;
  
  elasticity->S = S;
  elasticity->L = L_in;
  
  double w_npa = j2d_compute_w_npa(m, cm_ctx.npa, cm_ctx.alpha);
  double H_npa = j2d_compute_H_npa(m, cm_ctx.npa, cm_ctx.alpha);
  double g_npa = j2d_compute_gam_npa(m, cm_ctx.npa, cm_ctx.alpha);
  
  int is_damaged  = J2d_determine_damaged_npa(m, cm_ctx.npa);
  
  Tensor<2> sp_npa;
  State_variables *sv = m->vars_list[0] + m->model_id;
  err += sv->compute_Fs_npa(sp_npa.data, TENSOR_pSnm1, TENSOR_pSnm1, TENSOR_pSnm1,  cm_ctx.npa, cm_ctx.alpha);
  
  double dam = 1.0 - w_npa;

  Tensor<2> S0, Sbar;
      
  if(cm_ctx.eFnpa){
    err += compute_S0_Sbar_public(S0.data, Sbar.data, cm_ctx.eFnpa, sp_npa.data, elasticity);    

    if(compute_stiffness){
      err += compute_Lbar_public(L_in, cm_ctx.eFnpa, 
                                 Fs[TENSOR_Fn].m_pdata,
                                 Fs[TENSOR_pSn].m_pdata,
                                 g_npa,
                                 this->cm_mat->mat_J2p,
                                 elasticity);
    }        
  } else {
    Tensor<2> eF = {};

    if(cm_ctx.is_coulpled_with_thermal){
      TensorA<2> Fnp1(Fs[TENSOR_Fnp1].m_pdata);
      TensorA<2> hFnp1(cm_ctx.hFnp1);
      Tensor<2> hFnp1_I;

      err += inv(hFnp1, hFnp1_I);
      eF = Fnp1(i,k)*hFnp1_I(k,j);  
      err += compute_S0_Sbar_public(S0.data, Sbar.data, eF.data, sp_npa.data, elasticity);
      if(compute_stiffness){
        err += compute_Lbar_public(L_in, eF.data, 
                                   Fs[TENSOR_Fn].m_pdata,
                                   Fs[TENSOR_pSn].m_pdata,
                                   g_npa,
                                   this->cm_mat->mat_J2p,
                                   elasticity);
      }
    }
    else{
      err += compute_S0_Sbar_public(S0.data, Sbar.data, Fs[TENSOR_Fnp1].m_pdata, sp_npa.data, elasticity);
      if(compute_stiffness){
        err += compute_Lbar_public(L_in, Fs[TENSOR_Fnp1].m_pdata, 
                                   Fs[TENSOR_Fn].m_pdata,
                                   Fs[TENSOR_pSn].m_pdata,
                                   g_npa,
                                   this->cm_mat->mat_J2p,
                                   elasticity);
      }
    }      
  }
  
  err += apply_damage_on_stress(S, Sbar.data, w_npa);
  
  if(compute_stiffness)
  {
    for(int ia=0; ia<DIM_3x3x3x3; ++ia)
      L_in[ia] *= dam;
    
    if(is_damaged){
      double dt_mu = cm_ctx.dt*(this->cm_mat->mat_d->mu);
      double evo = dt_mu*H_npa/(1.0+dt_mu);        
      TensorA<4> L(L_in);
      L(i,j,k,l) = L(i,j,k,l) - 0.5*evo*(Sbar(i,j)*S0(k,l) + S0(i,j)*Sbar(k,l));
    }
  }  
  
  elasticity->S = S_tmp;
  elasticity->L = L_tmp;   
  return err;
}

int CM_J2P_PARAM::update_elasticity_dev(const Constitutive_model *m,
                                        double *eFnpa,
                                        double *L_in,
                                        double *S,
                                        const int npa,
                                        const double alpha,
                                        const double dt,
                                        const int compute_stiffness)  
const
{
  int err = 0;  
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;

  ELASTICITY *elasticity = this->cm_elast;
    
  double *S_tmp = elasticity->S;
  double *L_tmp = elasticity->L;
  
  elasticity->S = S;
  elasticity->L = L_in;
  
  double w_npa = j2d_compute_w_npa(m, npa, alpha);
  double H_npa = j2d_compute_H_npa(m, npa, alpha);
  double g_npa = j2d_compute_gam_npa(m, npa, alpha);
  
  int is_damaged  = J2d_determine_damaged_npa(m, npa);
  
  Tensor<2> sp_npa;
  State_variables *sv = m->vars_list[0] + m->model_id;
  err += sv->compute_Fs_npa(sp_npa.data, TENSOR_pSnm1, TENSOR_pSnm1, TENSOR_pSnm1,  npa, alpha);
  
  double dam = 1.0 - w_npa;

  Tensor<2> S0, Sbar;
      
  if(eFnpa){
    err += compute_S0_Sbar_dev_public(S0.data, Sbar.data, eFnpa, sp_npa.data, elasticity);    

    if(compute_stiffness){
      err += compute_Lbar_dev_public(L_in, eFnpa, 
                                     Fs[TENSOR_Fn].m_pdata,
                                     Fs[TENSOR_pSn].m_pdata,
                                     g_npa,
                                     this->cm_mat->mat_J2p,
                                     elasticity);
    }        
  } else {

    err += compute_S0_Sbar_dev_public(S0.data, Sbar.data, Fs[TENSOR_Fnp1].m_pdata, sp_npa.data, elasticity);
    if(compute_stiffness){
      err += compute_Lbar_dev_public(L_in, Fs[TENSOR_Fnp1].m_pdata, 
                                     Fs[TENSOR_Fn].m_pdata,
                                     Fs[TENSOR_pSn].m_pdata,
                                     g_npa,
                                     this->cm_mat->mat_J2p,
                                     elasticity);
    }      
  }
  
  err += apply_damage_on_stress(S, Sbar.data, w_npa);
  
  if(compute_stiffness)
  {
    for(int ia=0; ia<DIM_3x3x3x3; ++ia)
      L_in[ia] *= dam;
    
    if(is_damaged){
      double dt_mu = dt*(this->cm_mat->mat_d->mu);
      double evo = dt_mu*H_npa/(1.0+dt_mu);        
      TensorA<4> L(L_in);
      L(i,j,k,l) = L(i,j,k,l) - 0.5*evo*(Sbar(i,j)*S0(k,l) + S0(i,j)*Sbar(k,l));
    }
  }  
  
  elasticity->S = S_tmp;
  elasticity->L = L_tmp;   
  return err; 
}

int CM_J2P_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  
  // deformation gradients
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;  
  Fs[TENSOR_Fnm1] = Fs[TENSOR_Fn];
  Fs[TENSOR_Fn]   = Fs[TENSOR_Fnp1];
  Fs[TENSOR_pSnm1] = Fs[TENSOR_pSn];
  Fs[TENSOR_pSn]   = Fs[TENSOR_pSnp1];  

  // state variables 
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;  
  vars[VAR_w_nm1] = vars[VAR_w_n];
  vars[VAR_H_nm1] = vars[VAR_H_n];
  vars[VAR_X_nm1] = vars[VAR_X_n];
  
  vars[VAR_ep_nm1]  = vars[VAR_ep_n];
  vars[VAR_gam_nm1] = vars[VAR_gam_n];
    
  vars[VAR_w_n] = vars[VAR_w_np1];
  vars[VAR_H_n] = vars[VAR_H_np1];
  vars[VAR_X_n] = vars[VAR_X_np1];
  
  vars[VAR_ep_n]  = vars[VAR_ep_np1];
  vars[VAR_gam_n] = vars[VAR_gam_np1];

  // flags 
  int *flags = m->vars_list[0][m->model_id].flags;  
  flags[FLAG_damaged_n] = flags[FLAG_damaged];

  return err;
}

int CM_J2P_PARAM::reset_state_vars(Constitutive_model *m)
const 
{
  int err = 0;

  // deformation gradients 
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[TENSOR_Fnp1]  = Fs[TENSOR_Fn];
  Fs[TENSOR_pSnp1]  = Fs[TENSOR_pSn];

  // state variables
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;  
  vars[VAR_w_np1]   = vars[VAR_w_n];
  vars[VAR_H_np1]   = vars[VAR_H_n];
  vars[VAR_X_np1]   = vars[VAR_X_n];
  vars[VAR_ep_np1]  = vars[VAR_ep_n];
  vars[VAR_gam_np1] = vars[VAR_gam_n];  

  // flags
  int *flags = m->vars_list[0][m->model_id].flags;  
  flags[FLAG_damaged] = flags[FLAG_damaged_n];

  return err;
}

int CM_J2P_PARAM::get_subdiv_param(const Constitutive_model *m,
                                   double *subdiv_param,
                                   double dt)
const
{
  int err = 0;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  *subdiv_param = compute_subdiv_param(vars[VAR_w_n], vars[VAR_w_np1]);
  return err;
}

int CM_J2P_PARAM::get_var_info(Model_var_info &info)
const
{
  const CMVariableNames variable_names[] = {{VAR_w_np1,  "w"},
                                            {VAR_H_np1,  "H"},
                                            {VAR_X_np1,  "X" },
                                            {VAR_ep_np1, "ep" },
                                            {VAR_gam_np1,"gam"},
                                            {VAR_w_n,    "w_n"},
                                            {VAR_X_n,    "X_n"},
                                            {VAR_H_n,    "H_n"},
                                            {VAR_ep_n,   "ep_n"},
                                            {VAR_gam_n,  "gam_n"},                                            
                                            {VAR_w_nm1,  "w_nm1"},
                                            {VAR_X_nm1,  "X_nm1"},
                                            {VAR_H_nm1,  "H_nm1"},
                                            {VAR_ep_nm1, "ep_nm1"},
                                            {VAR_gam_nm1,"gam_nm1"}                                            
                                           };
  
  
  const CMVariableNames tensor_names[] = {{TENSOR_Fnm1, "Fnm1"},
                                          {TENSOR_Fn,   "Fn"},
                                          {TENSOR_Fnp1, "F"},
                                          {TENSOR_pSnm1,"pSnm1"},
                                          {TENSOR_pSn,  "pSn"},
                                          {TENSOR_pSnp1,"pS"} 
                                         };
  const CMVariableNames flag_names[] = {{FLAG_damaged_n, "damaged_n"},
                                        {FLAG_damaged,   "damaged"  }
                                       };
  
  int err = constitutive_model_info(info, VAR_end,    variable_names,
                                          TENSOR_end, tensor_names,
                                          FLAG_end,   flag_names);

  return err;
}

int J2D_get_F(const Constitutive_model *m,
              double *F,
              const int stepno)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_J2P_PARAM::get_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const
{
  return J2D_get_F(m,F,stepno);
}

int CM_J2P_PARAM::set_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_J2P_PARAM::get_eF(const Constitutive_model *m,
                         double *F,
                         const int stepno)
const
{
  return J2D_get_F(m,F,stepno);
}

int CM_J2P_PARAM::get_pF(const Constitutive_model *m,
                         double *F,
                         const int stepno)
const 
{
  F[0] = F[4] = F[8] = 1.0;
  F[1] = F[2] = F[3] = F[5] = F[6] = F[7] = 0.0;
  
  return 0;
}

int CM_J2P_PARAM::get_eF_of_hF(const Constitutive_model *m,
                               double *eF_in,
                               double *hFI_in,
                               const int stepno)
const
{
  int err = 0;
  TensorA<2> eF(eF_in), hFI(hFI_in);
  
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
    {  
      TensorA<2> F(Fs[TENSOR_Fnm1].m_pdata);
      eF = F(i,k)*hFI(k,j);
      break;
    }
    case 1: // n
    {  
      TensorA<2> F(Fs[TENSOR_Fn].m_pdata);
      eF = F(i,k)*hFI(k,j);
      break;
    }
    case 2: // n+1
    {  
      TensorA<2> F(Fs[TENSOR_Fnp1].m_pdata);
      eF = F(i,k)*hFI(k,j);
      break;
    }
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);

  return err;
}

int CM_J2P_PARAM::reset_state_vars_using_temporal(const Constitutive_model *m, 
                                                  State_variables *var)
const
{
  int err = 0;
  Matrix<double> *Fs    = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs_in = var->Fs;
  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fn    ].m_pdata[ia] = Fs_in[TENSOR_Fn    ].m_pdata[ia];
    Fs[TENSOR_Fnm1  ].m_pdata[ia] = Fs_in[TENSOR_Fnm1  ].m_pdata[ia];
    Fs[TENSOR_pSn   ].m_pdata[ia] = Fs_in[TENSOR_pSn   ].m_pdata[ia];
    Fs[TENSOR_pSnm1 ].m_pdata[ia] = Fs_in[TENSOR_pSnm1 ].m_pdata[ia];    
  }

  double *vars     = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *vars_in  = var->state_vars[0].m_pdata;
  vars[VAR_w_n]   = vars_in[VAR_w_n];
  vars[VAR_H_n]   = vars_in[VAR_H_n];
  vars[VAR_X_n]   = vars_in[VAR_X_n];
  vars[VAR_ep_n]  = vars_in[VAR_ep_n];
  vars[VAR_gam_n] = vars_in[VAR_gam_n];  
  
  vars[VAR_w_nm1]   = vars_in[VAR_w_nm1];
  vars[VAR_H_nm1]   = vars_in[VAR_H_nm1];
  vars[VAR_X_nm1]   = vars_in[VAR_X_nm1];
  vars[VAR_ep_nm1]  = vars_in[VAR_ep_nm1];
  vars[VAR_gam_nm1] = vars_in[VAR_gam_nm1];  
  
  
  int *flags    = m->vars_list[0][m->model_id].flags;
  int *flags_in = var->flags;
  
  flags[FLAG_damaged_n] = flags_in[FLAG_damaged_n];
  
  return err;
}

int CM_J2P_PARAM::update_np1_state_vars_to_temporal(const Constitutive_model *m, 
                                                    State_variables *var)
const
{
  int err = 0;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs    = var->Fs;
  for(int ia=0; ia<DIM_3x3; ia++){
    Fs[TENSOR_Fnp1 ].m_pdata[ia] = Fs_in[TENSOR_Fnp1 ].m_pdata[ia];
    Fs[TENSOR_pSnp1].m_pdata[ia] = Fs_in[TENSOR_pSnp1].m_pdata[ia];
  }    

  double *vars_in = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *vars    = var->state_vars[0].m_pdata;
  vars[VAR_w_np1]   = vars_in[VAR_w_np1];
  vars[VAR_H_np1]   = vars_in[VAR_H_np1];
  vars[VAR_X_np1]   = vars_in[VAR_X_np1];
  vars[VAR_ep_np1]  = vars_in[VAR_ep_np1];
  vars[VAR_gam_np1] = vars_in[VAR_gam_np1];  
  
  int *flags_in = m->vars_list[0][m->model_id].flags;
  int *flags    = var->flags;
  
  flags[FLAG_damaged] = flags_in[FLAG_damaged];
  
  return err;
}

int CM_J2P_PARAM::save_state_vars_to_temporal(const Constitutive_model *m, 
                                              State_variables *var)
const
{
  int err = 0;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs    = var->Fs;
  for(int ia=0; ia<TENSOR_end; ia++)
  {
    for(int ib=0; ib<DIM_3x3; ib++)    
      Fs[ia].m_pdata[ib] = Fs_in[ia].m_pdata[ib];
  }

  double *vars_in = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *vars    = var->state_vars[0].m_pdata;
  for(int ia=0; ia<VAR_end; ia++)
    vars[ia] = vars_in[ia];  
    
  int *flags_in = m->vars_list[0][m->model_id].flags;
  int *flags    = var->flags;
  
  for(int ia=0; ia<FLAG_end; ia++)
    flags[ia] = flags_in[ia];  
  
  return err;
}

int CM_J2P_PARAM::get_hardening(const Constitutive_model *m,
                                double *var,
                                const int stepno)
const
{
  int err = 0;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  switch(stepno)
  {
    case 0: // n-1
      *var = vars[VAR_w_nm1];
      break;
    case 1: // n
      *var = vars[VAR_w_n];
      break;
    case 2: // n+1
      *var = vars[VAR_w_np1];
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;  
}

int CM_J2P_PARAM::get_plast_strain_var(const Constitutive_model *m,
                                       double *chi)
const                                   
{
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  *chi = vars[VAR_ep_n];
  return 0;
}

int CM_J2P_PARAM::write_restart(FILE *out, 
                                const Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *vars  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  const int *flags = (m->vars_list[0][m->model_id]).flags;
  
  err += cm_write_tensor_restart(out, Fs[TENSOR_Fn   ].m_pdata);
  err += cm_write_tensor_restart(out, Fs[TENSOR_pSn  ].m_pdata);
  err += cm_write_tensor_restart(out, Fs[TENSOR_Fnm1 ].m_pdata);
  err += cm_write_tensor_restart(out, Fs[TENSOR_pSnm1].m_pdata);
    
  fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", 
                             vars[VAR_w_n  ], vars[VAR_w_nm1 ], 
                             vars[VAR_X_n  ], vars[VAR_X_nm1 ], 
                             vars[VAR_H_n  ], vars[VAR_H_nm1 ],
                             vars[VAR_ep_n ], vars[VAR_ep_nm1], 
                             vars[VAR_gam_n], vars[VAR_gam_nm1]);
                                                          
  fprintf(out, "%d\n", flags[FLAG_damaged_n]);

  return err;
}

int CM_J2P_PARAM::read_restart(FILE *in, 
                                Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *vars  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  int *flags = (m->vars_list[0][m->model_id]).flags;
  
  err += cm_read_tensor_restart(in, Fs[TENSOR_Fn   ].m_pdata);
  err += cm_read_tensor_restart(in, Fs[TENSOR_pSn  ].m_pdata);
  err += cm_read_tensor_restart(in, Fs[TENSOR_Fnm1 ].m_pdata);
  err += cm_read_tensor_restart(in, Fs[TENSOR_pSnm1].m_pdata);  
    
  fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                                        vars+VAR_w_n,   vars+VAR_w_nm1, 
                                        vars+VAR_X_n,   vars+VAR_X_nm1, 
                                        vars+VAR_H_n,   vars+VAR_H_nm1,
                                        vars+VAR_ep_n,  vars+VAR_ep_nm1,
                                        vars+VAR_gam_n, vars+VAR_gam_nm1);                                        
                                        
  fscanf(in, "%d", flags + FLAG_damaged_n);

  err += this->reset_state_vars(m);
  return err;
}

int CM_J2P_PARAM::set_init_vals(Constitutive_model *m)
const
{
  int err = 0;
  double I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
  double zero[9] = {};
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  memcpy(Fs[TENSOR_Fnm1].m_pdata, I, DIM_3x3*sizeof(double));
  memcpy(Fs[TENSOR_Fn  ].m_pdata, I, DIM_3x3*sizeof(double));
  memcpy(Fs[TENSOR_Fnp1].m_pdata, I, DIM_3x3*sizeof(double));
  
  memcpy(Fs[TENSOR_pSnm1].m_pdata, zero, DIM_3x3*sizeof(double));
  memcpy(Fs[TENSOR_pSn  ].m_pdata, zero, DIM_3x3*sizeof(double));
  memcpy(Fs[TENSOR_pSnp1].m_pdata, zero, DIM_3x3*sizeof(double));    
  

  double *vars = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  for(int ia=0; ia<VAR_end; ia++)
    vars[ia] = 0.0;
    
  int *flags = m->vars_list[0][m->model_id].flags;
  
  for(int ia=0; ia<FLAG_end; ia++)
    flags[ia] = 0;
    
  return err;
}

int CM_J2P_PARAM::read_param(FILE *in)
const
{
  int err = 0;
  // get pointer to parameter data
  double *param = this->model_param;
  assert(param != NULL); // check the pointer

  //  There are three (3) sets of parameters:
  // 
  //  1. elastic properties (G, nu) **NOTE**: These are currently read
  //     from the HOMMAT object.
  //  2. plastic properties (hp, beta, k0)
  //  3. damage properties (mu, p1, p2, Yin)
  // 
  // 
  // Get the elastic properties from the HOMMAT object.
  // 
  // This is to maintain consistency with the other implementations and
  // to ensure that the correct bulk modulous is computed elsewhere.
  // Note that we **DO NOT** use the deviatoric constitutive law
  // prescribed in the HOMMAT object, but use a spatial formulation of
  // the Neo-Hookean model to maintain consistency with the J2+damage
  // formulation in Simo and Ju.

  int match = 2;
  param[PARAM_G]  = this->p_hmat->G;
  param[PARAM_nu] = this->p_hmat->nu;

  // PLASTIC PROPERTIES 
  err += scan_for_valid_line(in);
  match += fscanf(in, "%lf %lf %lf",
                  param + PARAM_beta, param + PARAM_hp, param + PARAM_k0);

  err += scan_for_valid_line(in);

  // DAMAGE PROPERTIES 
  match += fscanf(in, "%lf %lf %lf %lf %lf",
                  param + PARAM_mu, param + PARAM_w_max, param + PARAM_P1,
                  param + PARAM_P2, param + PARAM_Yin);

  if (match != PARAM_NO) err++;
  assert(match == PARAM_NO && "Did not read expected number of parameters");

  // PARAM_w_max in [0, 1) 
  if (param[PARAM_w_max] >= 1) param[PARAM_w_max] = DAMAGE_THRESH;

  // scan past any other comment/blank lines in the block 
  err += scan_for_valid_line(in);

  // not expecting EOF, check and return error if encountered 
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");
  
  this->cm_mat->mat_d   = new MATERIAL_CONTINUUM_DAMAGE;
  set_damage_parameters(this->cm_mat->mat_d, param[PARAM_P1],  param[PARAM_P2],
                              param[PARAM_Yin], param[PARAM_mu], param[PARAM_w_max]);  

  this->cm_mat->mat_J2p = new MATERIAL_J2_PLASTICITY;
  set_J2_plasticity_parameters(this->cm_mat->mat_J2p, 
                               param[PARAM_hp], param[PARAM_beta], param[PARAM_k0]);  
  return err;
}


int CM_J2P_PARAM::model_dependent_initialization(void)
{
  int err = 0;
  this->type              = J2_PLASTICITY_DAMAGE;
  this->n_param           = PARAM_NO;
  this->model_param       = new double[PARAM_NO]();
  this->n_param_index     = 0;
  this->model_param_index = NULL;
  this->gcm_solver_info   = NULL;

  return err;
}

int CM_J2P_PARAM::model_dependent_finalization(void)
{
  int err = 0;
  delete (this->cm_mat)->mat_d;
  return err;
}