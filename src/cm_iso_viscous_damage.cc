/// This file defines the implementation for the isotropic viscous
/// damage model.
/// 
/// REFERENCES:
/// 
/// Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
/// finite strains." Computational Mechanics 5.5 (1989): 375-400.
/// 
/// Mosby, Matthew and K. Matous. "On mechanics and material length scales of failure
/// in heterogeneous interfaces using a finite strain high performance
/// solver." Modelling and Simulation in Materials Science and
/// Engineering 23.8 (2015): 085014.
/// 
/// Authors:
///  Sangmin Lee, University of Notre Dame, <slee43@nd.edu>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "allocation.h"
#include "cm_iso_viscous_damage.h"
#include "constitutive_model.h"
#include "index_macros.h"
#include "utils.h"
#include <math.h>
#include <string.h>
#include <assert.h>
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
        TENSOR_end};
        
  enum {VAR_w_np1,
        VAR_H_np1,
        VAR_X_np1,
        VAR_w_n,
        VAR_X_n,
        VAR_H_n,
        VAR_w_nm1,
        VAR_X_nm1,
        VAR_H_nm1,
        VAR_end};
              
  enum {FLAG_damaged_n,
        FLAG_damaged,
        FLAG_end};
                
  enum {PARAM_mu, 
        PARAM_w_max, 
        PARAM_P1, 
        PARAM_P2, 
        PARAM_Yin,
        PARAM_NO};      
}

double ivd_compute_w_npa(const Constitutive_model *m,
                         const int npa,
                         const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_w_nm1, VAR_w_n, VAR_w_np1, npa, alpha);
}

double ivd_compute_H_npa(const Constitutive_model *m,
                         const int npa,
                         const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_H_nm1, VAR_H_n, VAR_H_np1, npa, alpha);
}

double ivd_compute_X_npa(const Constitutive_model *m,
                         const int npa,
                         const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_X_nm1, VAR_X_n, VAR_X_np1, npa, alpha);
}

int determine_damaged_npa(const Constitutive_model *m,
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

int CM_IVD_PARAM::integration_algorithm(Constitutive_model *m,
                                        CM_Ctx &cm_ctx)
const
{
  
  int err = 0;
  
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  int    *flags      = m->vars_list[0][m->model_id].flags;

  // integration algorithm
  ELASTICITY *elasticity = this->cm_elast;
  MATERIAL_CONTINUUM_DAMAGE *mat_d = this->cm_mat->mat_d;
      
  // store the deformation gradient
  memcpy(Fs[TENSOR_Fnp1].m_pdata, cm_ctx.F, DIM_3x3 * sizeof(double));
  
  err += continuum_damage_integration_alg(mat_d,elasticity,
                                          vars+VAR_w_np1,
                                          vars+VAR_X_np1,
                                          vars+VAR_H_np1,
                                          flags+FLAG_damaged,                                          
                                          vars[VAR_w_n],
                                          vars[VAR_X_n],
                                          cm_ctx.dt, Fs[TENSOR_Fnp1].m_pdata); 

  return err;
}

int CM_IVD_PARAM::compute_dev_stress(const Constitutive_model *m,
                                     CM_Ctx &cm_ctx,
                                     double *stress)
const
{
  return 0;
}

int CM_IVD_PARAM::compute_dudj(const Constitutive_model *m,
                               CM_Ctx &cm_ctx,
                               double *dudj)
const
{
  int err = 0;
  double eJ = det3x3(cm_ctx.F);
  
  double *vars  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;

  // scale by damage variable 
  ELASTICITY *elasticity = this->cm_elast;
  *dudj = damage_compute_dudj(elasticity, eJ, vars[VAR_w_np1]);
  
  return err;
}

double CM_IVD_PARAM::compute_dudj(const Constitutive_model *m,
                                  double theta_e,
                                  const int npa,
                                  const double alpha)
const 
{
  double w_npa = ivd_compute_w_npa(m,npa,alpha);
  ELASTICITY *elasticity = this->cm_elast;
  return damage_compute_dudj(elasticity, theta_e, w_npa);          
}

int CM_IVD_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                      CM_Ctx &cm_ctx,
                                      double *L)
const
{
  return 0;
}

int CM_IVD_PARAM::compute_d2udj2(const Constitutive_model *m,
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

double CM_IVD_PARAM::compute_d2udj2(const Constitutive_model *m,
                                    double theta_e,
                                    const int npa,
                                    const double alpha)
const 
{
  double w_npa = ivd_compute_w_npa(m,npa,alpha);
  ELASTICITY *elasticity = this->cm_elast;
  return damage_compute_d2udj2(elasticity, theta_e, w_npa);           
}

int CM_IVD_PARAM::update_elasticity(const Constitutive_model *m,
                                    CM_Ctx &cm_ctx,
                                    double *L,
                                    double *S,
                                    const int compute_stiffness)
const
{
  int err = 0;  
  //double *vars  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  //int *flags = m->vars_list[0][m->model_id].flags;

  ELASTICITY *elasticity = this->cm_elast;
  MATERIAL_CONTINUUM_DAMAGE *mat_d = this->cm_mat->mat_d;
  
  double *S_tmp = elasticity->S;
  double *L_tmp = elasticity->L;
  
  elasticity->S = S;
  elasticity->L = L;
  
  double w_npa = ivd_compute_w_npa(m, cm_ctx.npa, cm_ctx.alpha);
  double H_npa = ivd_compute_H_npa(m, cm_ctx.npa, cm_ctx.alpha);
  int is_damaged  = determine_damaged_npa(m, cm_ctx.npa);
    
  if(cm_ctx.eFnpa)
  {
    err += update_damage_elasticity(mat_d, elasticity, 
                                    w_npa,
                                    H_npa,
                                    is_damaged,
                                    cm_ctx.dt, 
                                    cm_ctx.eFnpa, 
                                    compute_stiffness);
  }    
  else
  {
    // shorthand of deformation gradients
    Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;  
    Tensor<2> eF = {};  
    
    if(cm_ctx.is_coulpled_with_thermal){
      TensorA<2> Fnp1(Fs[TENSOR_Fnp1].m_pdata);
      TensorA<2> hFnp1(cm_ctx.hFnp1);
      Tensor<2> hFnp1_I;

      err += inv(hFnp1, hFnp1_I);
      eF = Fnp1(i,k)*hFnp1_I(k,j);
  
      err += update_damage_elasticity(mat_d, elasticity, 
                                      w_npa,
                                      H_npa,
                                      is_damaged,
                                      cm_ctx.dt, 
                                      eF.data, 
                                      compute_stiffness);
    }
    else{
      err += update_damage_elasticity(mat_d, elasticity, 
                                      w_npa,
                                      H_npa,
                                      is_damaged,
                                      cm_ctx.dt, 
                                      Fs[TENSOR_Fnp1].m_pdata, 
                                      compute_stiffness);
    }      
  }
  
  elasticity->S = S_tmp;
  elasticity->L = L_tmp;   
  return err;
}

int CM_IVD_PARAM::update_elasticity_dev(const Constitutive_model *m,
                                        double *eFnpa,
                                        double *L,
                                        double *S,
                                        const int npa,
                                        const double alpha,
                                        const double dt,
                                        const int compute_stiffness)  
const
{
  int err = 0;
  
  //double *vars = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  //int *flags = m->vars_list[0][m->model_id].flags;

  ELASTICITY *elasticity = this->cm_elast;
  MATERIAL_CONTINUUM_DAMAGE *mat_d = this->cm_mat->mat_d;
  
  double *S_tmp = elasticity->S;
  double *L_tmp = elasticity->L;
  
  elasticity->S = S;
  elasticity->L = L;
  
  double w_npa = ivd_compute_w_npa(m, npa, alpha);
  double H_npa = ivd_compute_H_npa(m, npa, alpha);
  int is_damaged  = determine_damaged_npa(m, npa);
    
  err += update_damage_elasticity_dev(mat_d, elasticity, 
                                      w_npa,
                                      H_npa,
                                      is_damaged,
                                      dt, 
                                      eFnpa, 
                                      compute_stiffness);
  elasticity->S = S_tmp;
  elasticity->L = L_tmp;   
  return err; 
}

int CM_IVD_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  
  // deformation gradients
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;  
  Fs[TENSOR_Fnm1] = Fs[TENSOR_Fn];
  Fs[TENSOR_Fn]   = Fs[TENSOR_Fnp1];

  // state variables 
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;  
  vars[VAR_w_nm1] = vars[VAR_w_n];
  vars[VAR_H_nm1] = vars[VAR_H_n];
  vars[VAR_X_nm1] = vars[VAR_X_n];
    
  vars[VAR_w_n] = vars[VAR_w_np1];
  vars[VAR_H_n] = vars[VAR_H_np1];
  vars[VAR_X_n] = vars[VAR_X_np1];

  // flags 
  int *flags = m->vars_list[0][m->model_id].flags;  
  flags[FLAG_damaged_n] = flags[FLAG_damaged];

  return err;
}

int CM_IVD_PARAM::reset_state_vars(Constitutive_model *m)
const 
{
  int err = 0;

  // deformation gradients 
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[TENSOR_Fnp1]  = Fs[TENSOR_Fn];

  // state variables
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;  
  vars[VAR_w_np1] = vars[VAR_w_n];
  vars[VAR_H_np1] = vars[VAR_H_n];
  vars[VAR_X_np1] = vars[VAR_X_n];

  // flags
  int *flags = m->vars_list[0][m->model_id].flags;  
  flags[FLAG_damaged] = flags[FLAG_damaged_n];

  return err;
}

int CM_IVD_PARAM::get_subdiv_param(const Constitutive_model *m,
                                   double *subdiv_param,
                                   double dt)
const
{
  int err = 0;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  *subdiv_param = compute_subdiv_param(vars[VAR_w_n], vars[VAR_w_np1]);
  return err;
}

int CM_IVD_PARAM::get_var_info(Model_var_info &info)
const
{
  const CMVariableNames variable_names[] = {{VAR_w_np1,"w"  },
                                            {VAR_H_np1,"H"  },
                                            {VAR_X_np1,"X"  },
                                            {VAR_w_n,  "w_n"},
                                            {VAR_X_n,  "X_n"},
                                            {VAR_H_n,  "H_n"},
                                            {VAR_w_nm1,"w_nm1"},
                                            {VAR_X_nm1,"X_nm1"},
                                            {VAR_H_nm1,"H_nm1"},
                                           };
  
  
  const CMVariableNames tensor_names[] = {{TENSOR_Fnm1, "Fnm1"},
                                          {TENSOR_Fn,   "Fn"},
                                          {TENSOR_Fnp1, "F"} 
                                         };
  const CMVariableNames flag_names[] = {{FLAG_damaged_n, "damaged_n"},
                                        {FLAG_damaged,   "damaged"  },
                                       };
  
  int err = constitutive_model_info(info, VAR_end,    variable_names,
                                          TENSOR_end, tensor_names,
                                          FLAG_end,   flag_names);

  return err;
}

int CM_IVD_PARAM::get_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_IVD_PARAM::set_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->set_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_IVD_PARAM::get_eF(const Constitutive_model *m,
                         double *F,
                         const int stepno)
const
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_IVD_PARAM::get_pF(const Constitutive_model *m,
                         double *F,
                         const int stepno)
const 
{
  F[0] = F[4] = F[8] = 1.0;
  F[1] = F[2] = F[3] = F[5] = F[6] = F[7] = 0.0;
  
  return 0;
}

int CM_IVD_PARAM::get_eF_of_hF(const Constitutive_model *m,
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

int CM_IVD_PARAM::reset_state_vars_using_temporal(const Constitutive_model *m, 
                                                  State_variables *var)
const
{
  int err = 0;
  Matrix<double> *Fs    = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs_in = var->Fs;
  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fn   ].m_pdata[ia] = Fs_in[TENSOR_Fn   ].m_pdata[ia];
    Fs[TENSOR_Fnm1 ].m_pdata[ia] = Fs_in[TENSOR_Fnm1 ].m_pdata[ia];
  }

  double *vars     = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *vars_in  = var->state_vars[0].m_pdata;
  vars[VAR_w_n] = vars_in[VAR_w_n];
  vars[VAR_H_n] = vars_in[VAR_H_n];
  vars[VAR_X_n] = vars_in[VAR_X_n];
  
  vars[VAR_w_nm1] = vars_in[VAR_w_nm1];
  vars[VAR_H_nm1] = vars_in[VAR_H_nm1];
  vars[VAR_X_nm1] = vars_in[VAR_X_nm1];
  
  
  int *flags    = m->vars_list[0][m->model_id].flags;
  int *flags_in = var->flags;
  
  flags[FLAG_damaged_n] = flags_in[FLAG_damaged_n];
  
  return err;
}


int CM_IVD_PARAM::update_np1_state_vars_to_temporal(const Constitutive_model *m, 
                                                    State_variables *var)
const
{
  int err = 0;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs    = var->Fs;
  for(int ia=0; ia<DIM_3x3; ia++)
    Fs[TENSOR_Fnp1].m_pdata[ia] = Fs_in[TENSOR_Fnp1].m_pdata[ia];

  double *vars_in = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *vars    = var->state_vars[0].m_pdata;
  vars[VAR_w_np1] = vars_in[VAR_w_np1];
  vars[VAR_H_np1] = vars_in[VAR_H_np1];
  vars[VAR_X_np1] = vars_in[VAR_X_np1];
  
  int *flags_in = m->vars_list[0][m->model_id].flags;
  int *flags    = var->flags;
  
  flags[FLAG_damaged] = flags_in[FLAG_damaged];
  
  return err;
}

int CM_IVD_PARAM::save_state_vars_to_temporal(const Constitutive_model *m, 
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


int CM_IVD_PARAM::get_hardening(const Constitutive_model *m,
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

int CM_IVD_PARAM::get_plast_strain_var(const Constitutive_model *m,
                                       double *chi)
const                                   
{
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  *chi = vars[VAR_X_n];
  return 0;
}

int CM_IVD_PARAM::write_restart(FILE *out, 
                                const Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *vars  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  const int *flags = (m->vars_list[0][m->model_id]).flags;
  
  err += cm_write_tensor_restart(out, Fs[TENSOR_Fn].m_pdata);
  err += cm_write_tensor_restart(out, Fs[TENSOR_Fnm1].m_pdata);
    
  fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e\n", 
                             vars[VAR_w_n], vars[VAR_w_nm1], 
                             vars[VAR_X_n], vars[VAR_X_nm1], 
                             vars[VAR_H_n], vars[VAR_H_nm1]);
                                                          
  fprintf(out, "%d\n", flags[FLAG_damaged_n]);

  return err;
}

int CM_IVD_PARAM::read_restart(FILE *in, 
                                Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *vars  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  int *flags = (m->vars_list[0][m->model_id]).flags;
  
  cm_read_tensor_restart(in, Fs[TENSOR_Fn].m_pdata);
  cm_read_tensor_restart(in, Fs[TENSOR_Fnm1].m_pdata);
    
  fscanf(in, "%lf %lf %lf %lf %lf %lf", vars+VAR_w_n, vars+VAR_w_nm1, 
                                        vars+VAR_X_n, vars+VAR_X_nm1, 
                                        vars+VAR_H_n, vars+VAR_H_nm1);                                        
                                        
  fscanf(in, "%d", flags + FLAG_damaged_n);

  err += this->reset_state_vars(m);
  return err;
}

int CM_IVD_PARAM::set_init_vals(Constitutive_model *m)
const
{
  int err = 0;
  double I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};

  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  for(int ia=0; ia<TENSOR_end; ia++)
  {
    for(int ib=0; ib<DIM_3x3; ib++)
      Fs[ia].m_pdata[ib] = I[ib];
  }

  double *vars = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  for(int ia=0; ia<VAR_end; ia++)
    vars[ia] = 0.0;
    
  int *flags = m->vars_list[0][m->model_id].flags;
  
  for(int ia=0; ia<FLAG_end; ia++)
    flags[ia] = 0;
    
  return err;
}

int CM_IVD_PARAM::read_param(FILE *in)
const
{
  int err = 0;
  // get pointer to parameter data
  double *param     = this->model_param;
  //int    *param_idx = this->model_param_index;
  assert(param     != NULL); // check the pointer
  //assert(param_idx != NULL); // check the pointer

  // scan to non-blank/comment line 
  err += scan_for_valid_line(in);

  // READ PROPERTIES IN ALPHABETICAL ORDER 
  int match = fscanf(in, "%lf %lf %lf %lf %lf",
                     param + PARAM_mu, param + PARAM_w_max, param + PARAM_P1,
                     param + PARAM_P2, param + PARAM_Yin);                     

  if (match != PARAM_NO) err++;
    assert(match == PARAM_NO && "Did not read expected number of parameters");

  // PARAM_w_max in [0, 1) 
  if (param[PARAM_w_max] >= 1.0) param[PARAM_w_max] = DAMAGE_THRESH;

  // scan past any other comment/blank lines in the block 
  err += scan_for_valid_line(in);

  // not expecting EOF, check and return error if encountered 
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");
  
  this->cm_mat->mat_d = new MATERIAL_CONTINUUM_DAMAGE;
  set_damage_parameters(this->cm_mat->mat_d, param[PARAM_P1],  param[PARAM_P2],
                              param[PARAM_Yin], param[PARAM_mu], param[PARAM_w_max]);
                                 
  return err;
}

int CM_IVD_PARAM::model_dependent_initialization(void)
{
  int err = 0;
  this->type              = ISO_VISCOUS_DAMAGE;
  this->n_param           = PARAM_NO;
  this->model_param       = new double[PARAM_NO]();
  this->n_param_index     = 0;
  this->model_param_index = NULL;
  this->gcm_solver_info   = NULL;

  return err;
}

int CM_IVD_PARAM::model_dependent_finalization(void)
{
  int err = 0;
  delete (this->cm_mat)->mat_d;
  return err;
}
