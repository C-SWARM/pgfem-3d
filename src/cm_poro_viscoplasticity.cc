/// Define poro-visco plasticity model functions for the constitutive model interface
/// 
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#include "cm_poro_viscoplasticity.h"

#include "plasticity_model.h"
#include "constitutive_model.h"
#include "cm_placeholder_functions.h"
#include "new_potentials.h"
#include "elem3d.h"
#include "gen_path.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "material_properties.h"
#include "hyperelasticity.h"

#include <iostream>
#include <iomanip>
#include <time.h>
#include <ttl/ttl.h>
#include"poro_visco_plasticity.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX_D_ALPHA 0.0002

static const int FLAG_end = 0;

// ttl declarations
namespace {
  /// Only dealing with 3 dimensional double data, but it is sometimes const.
  template <int R, class S = double>
  using Tensor = ttl::Tensor<R, 3, S>;
    
  template <int R, class S = double *>
  using TensorA = ttl::Tensor<R, 3, S>;  
  
  template <int R, class S = double>
  using Tensor9x9 = ttl::Tensor<R, 9, S>;
  
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
}

enum variable_names {
  VAR_pc_n,
  VAR_pc_np1,
  VAR_pc_nm1,
  VAR_end
};

enum tensor_names {
  TENSOR_Fn,
  TENSOR_pFn,
  TENSOR_Fnp1,
  TENSOR_pFnp1,
  TENSOR_Fnm1,
  TENSOR_pFnm1,
  TENSOR_end
};

enum param_names {
  PARAM_yf_M,       // Yield function parameters
  PARAM_yf_alpha,   //   :
  PARAM_flr_m,      // Flow rule parameters
  PARAM_flr_gamma0, //   :
  PARAM_hr_a1,      // Hardening rule parameters
  PARAM_hr_a2,      //   :
  PARAM_hr_Lambda1, //   :
  PARAM_hr_Lambda2, //   :
  PARAM_c_inf,      // Cohesion rule parameters
  PARAM_c_Gamma,    //   :
  PARAM_d_B,        // Transition rule parameters
  PARAM_d_pcb,      //   :
  PARAM_mu_0,       // Shear modulus parameters
  PARAM_mu_1,       //   :
  PARAM_K_p0,       // Bulk modulus parameters
  PARAM_K_kappa,    //   :
  PARAM_pl_n,       // Power law exponent
  PARAM_cf_g0,      // Compaction function parameters
  PARAM_cf_pcinf,   //   :
  PARAM_pc_0,       // initial pc
  PARAM_pJ,         // initial plastic deformation
  PARAM_tol,        // convergence tolerance
  PARAM_computer_zero, // computer zero
  PARAM_NO
};

enum param_index_names {
  PARAM_max_itr_stag,
  PARAM_max_itr_M,
  PARAM_max_subdivision,
  PARAM_INX_NO
};

int test_cm_poro_viscoplasticity_model(Constitutive_model *m);

double pc_npa(const Constitutive_model *m,
              const int npa,
              const double alpha)
{

  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_pc_nm1, VAR_pc_n, VAR_pc_np1, npa, alpha);
}              

int CM_PVP_PARAM::integration_algorithm(Constitutive_model *m,
                                        CM_Ctx &cm_ctx)
const
{
  int err = 0;

  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  
  memcpy(Fs[TENSOR_Fnp1].m_pdata, cm_ctx.F, DIM_3x3*sizeof(double));

  GcmPvpIntegrator integrator;
  integrator.mat = (m->param)->cm_mat->mat_pvp;
  integrator.solver_info = (m->param)->gcm_solver_info;
  integrator.set_tensors(Fs[TENSOR_Fnp1].m_pdata,
                         Fs[TENSOR_Fn].m_pdata,
                         Fs[TENSOR_Fnm1].m_pdata,                         
                         Fs[TENSOR_pFnp1].m_pdata,
                         Fs[TENSOR_pFn].m_pdata,
                         cm_ctx.hFnp1,
                         cm_ctx.hFn);
  integrator.pc_np1 = vars+VAR_pc_np1;
  integrator.pc_n   = vars[VAR_pc_n];
  err += integrator.run_integration_algorithm(cm_ctx.dtn, cm_ctx.dt);

  return err;
}


int CM_PVP_PARAM::compute_dev_stress(const Constitutive_model *m,
                                     CM_Ctx &cm_ctx,
                                     double *stress)
const
{
  int err = 0;
//  double eC[DIM_3x3] = {};
      
//  err += cp_compute_eC(cm_ctx.F,eC);
//  err += cp_dev_stress(eC,m->param->p_hmat,stress->m_pdata);

  return err;
}

int CM_PVP_PARAM::compute_dudj(const Constitutive_model *m,
                               CM_Ctx &cm_ctx,
                               double *dudj)
const
{
  int err = 0;
//  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
//  double J = det3x3(cm_ctx.F);
//  Pressure(J,m->param->p_hmat,dudj);
  return err;
}

double CM_PVP_PARAM::compute_dudj(const Constitutive_model *m,
                                  double theta_e,
                                  const int npa,
                                  const double alpha)
const 
{
  double pc = pc_npa(m, npa, alpha);
  MaterialPoroViscoPlasticity *mat_pvp = (m->param)->cm_mat->mat_pvp;      
  return poro_visco_plasticity_intf_compute_dudj(theta_e, pc, mat_pvp);           
}

int CM_PVP_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                      CM_Ctx &cm_ctx,
                                      double *L)
const
{
  int err = 0;
//  double eC[DIM_3x3] = {};
  
//  err += cp_compute_eC(cm_ctx.F,eC);
//  err += cp_dev_tangent(eC,m->param->p_hmat,tangent->m_pdata);
  return err;
}

int CM_PVP_PARAM::compute_d2udj2(const Constitutive_model *m,
                                 CM_Ctx &cm_ctx,
                                 double *d2udj2)
const
{
  int err = 0;
//  double J = det3x3(cm_ctx.F);

//  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
//  D_Pressure(J,m->param->p_hmat,d2udj2);
  return err;
}

double CM_PVP_PARAM::compute_d2udj2(const Constitutive_model *m,
                                    double theta_e,
                                    const int npa,
                                    const double alpha)
const 
{
  double pc = pc_npa(m, npa, alpha);  
  MaterialPoroViscoPlasticity *mat_pvp = (m->param)->cm_mat->mat_pvp;
  return poro_visco_plasticity_intf_compute_d2udj2(theta_e, pc, mat_pvp);           
}

int CM_PVP_PARAM::update_elasticity(const Constitutive_model *m,
                                    CM_Ctx &cm_ctx,
                                    double *L_in,
                                    double *S,
                                    const int compute_stiffness)
const
{
  int err = 0;
  double pc = pc_npa(m, cm_ctx.npa, cm_ctx.alpha);      
  MaterialPoroViscoPlasticity *mat_pvp = (m->param)->cm_mat->mat_pvp;    
  
  double *L = NULL;
  if(compute_stiffness)
    L = L_in;
  
  if(cm_ctx.eFnpa)
  {
    poro_visco_plasticity_update_elasticity(S, L, mat_pvp, 
                                            cm_ctx.eFnpa, pc, compute_stiffness);
  }    
  else
  {
    // shorthand of deformation gradients
    Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;  
    Tensor<2> eF = {}, pFnp1_I;  
    TensorA<2> pFnp1(Fs[TENSOR_pFnp1].m_pdata), Fnp1(Fs[TENSOR_Fnp1].m_pdata);
    
    err += inv(pFnp1, pFnp1_I);    
    if(cm_ctx.is_coulpled_with_thermal)
    {
      TensorA<2> hFnp1(cm_ctx.hFnp1);
      Tensor<2> hFnp1_I;

      err += inv(hFnp1, hFnp1_I);
      eF = Fnp1(i,k)*hFnp1_I(k,l)*pFnp1_I(l,j);    
    }
    else
      eF = Fnp1(i,k)*pFnp1_I(k,j);
    
    poro_visco_plasticity_update_elasticity(S, L, mat_pvp,
                                            eF.data, pc, compute_stiffness);                                           
  }     
  return err;
}

int CM_PVP_PARAM::update_elasticity_dev(const Constitutive_model *m,
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
  double pc = pc_npa(m, npa, alpha);
  MaterialPoroViscoPlasticity *mat_pvp = (m->param)->cm_mat->mat_pvp;    
  
  double *L = NULL;
  if(compute_stiffness)
    L = L_in;
  
  poro_visco_plasticity_update_elasticity_dev(S, L, mat_pvp, eFnpa, pc, compute_stiffness);           
  return err; 
}

int CM_PVP_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Fs[TENSOR_Fnm1] = Fs[TENSOR_Fn];    
  Fs[TENSOR_Fn]   = Fs[TENSOR_Fnp1];
  Fs[TENSOR_pFnm1]= Fs[TENSOR_pFn];
  Fs[TENSOR_pFn]  = Fs[TENSOR_pFnp1];
  
  state_var[VAR_pc_nm1] = state_var[VAR_pc_n];
  state_var[VAR_pc_n]   = state_var[VAR_pc_np1];
  return err;
}

int CM_PVP_PARAM::reset_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Fs[TENSOR_Fnp1]  = Fs[TENSOR_Fn];
  Fs[TENSOR_pFnp1] = Fs[TENSOR_pFn];  
  state_var[VAR_pc_np1] = state_var[VAR_pc_n];
  return err;
}

int CM_PVP_PARAM::get_subdiv_param(const Constitutive_model *m,
                                   double *subdiv_param,
                                   double dt)
const
{
  int err = 0;  

  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  GcmSolverInfo *solver_info = m->param->gcm_solver_info;
  
  MaterialPoroViscoPlasticity *mat_pvp = (m->param)->cm_mat->mat_pvp;
  
  double gamma_d = 0.0;
  double gamma_v = 0.0;  
  poro_visco_plasticity_intf_compute_gammas(gamma_d,
                                            gamma_v,
                                            Fs[TENSOR_Fnp1 ].m_pdata,
                                            Fs[TENSOR_Fn   ].m_pdata,
                                            Fs[TENSOR_pFnp1].m_pdata,
                                            Fs[TENSOR_pFn  ].m_pdata,
                                            state_var[VAR_pc_np1],
                                            state_var[VAR_pc_n],
                                            dt,
                                            mat_pvp,
                                            solver_info);

  double alpha = (fabs(gamma_d)>fabs(gamma_d))?fabs(gamma_d):fabs(gamma_v);
  *subdiv_param = dt*alpha/MAX_D_ALPHA;
  return err;
}

int CM_PVP_PARAM::get_var_info(Model_var_info &info)
const 
{
  const CMVariableNames variable_names[] = {{VAR_pc_n  , "pc_n"  },
                                            {VAR_pc_np1, "pc_np1"},
                                            {VAR_pc_nm1, "pc_nm1"}
                                           };

  const CMVariableNames tensor_names[] = {{TENSOR_Fn   , "Fn"   }, 
                                          {TENSOR_pFn  , "pFn"  },  
                                          {TENSOR_Fnp1 , "Fnp1" }, 
                                          {TENSOR_pFnp1, "pFnp1"},
                                          {TENSOR_Fnm1 , "Fnm1" }, 
                                          {TENSOR_pFnm1, "pFnm1"},
                                         };

  const CMVariableNames flag_names[] = {};                                       
  int err = constitutive_model_info(info, VAR_end,    variable_names,
                                          TENSOR_end, tensor_names,
                                          FLAG_end,   flag_names);
  return err;
}

int CM_PVP_PARAM::get_pF(const Constitutive_model *m,
                         double *F,
                         const int stepno)
const 
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
      memcpy(F,Fs[TENSOR_pFnm1].m_pdata,DIM_3x3*sizeof(double));
      break;
    case 1: // n
      memcpy(F,Fs[TENSOR_pFn].m_pdata,  DIM_3x3*sizeof(double));
      break;      
    case 2: // n+1
      memcpy(F,Fs[TENSOR_pFnp1].m_pdata,DIM_3x3*sizeof(double));
      break;      
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;
}

int CM_PVP_PARAM::get_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const 
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_PVP_PARAM::set_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_PVP_PARAM::get_eF(const Constitutive_model *m,
                         double *eF_in,
                         const int stepno)
const
{
  int err = 0;
  Tensor<2> pFI;
  TensorA<2> eF(eF_in);
  
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
    {  
      TensorA<2> pF(Fs[TENSOR_pFnm1].m_pdata), F(Fs[TENSOR_Fnm1].m_pdata);
      inv(pF, pFI);
      eF = F(i,k)*pFI(k,j);
      break;
    }
    case 1: // n
    {  
      TensorA<2> pF(Fs[TENSOR_pFn].m_pdata), F(Fs[TENSOR_Fn].m_pdata);
      inv(pF, pFI);
      eF = F(i,k)*pFI(k,j);
      break;
    }
    case 2: // n+1
    {  
      TensorA<2> pF(Fs[TENSOR_pFnp1].m_pdata), F(Fs[TENSOR_Fnp1].m_pdata);
      inv(pF, pFI);
      eF = F(i,k)*pFI(k,j);
      break;
    }
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);

  return err;
}

int CM_PVP_PARAM::get_eF_of_hF(const Constitutive_model *m,
                               double *eF_in,
                               double *hFI_in,
                               const int stepno)
const
{
  int err = 0;
  Tensor<2> pFI;
  TensorA<2> eF(eF_in), hFI(hFI_in);
  
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
    {  
      TensorA<2> pF(Fs[TENSOR_pFnm1].m_pdata), F(Fs[TENSOR_Fnm1].m_pdata);
      inv(pF, pFI);
      eF = F(i,k)*hFI(k,l)*pFI(l,j);
      break;
    }
    case 1: // n
    {  
      TensorA<2> pF(Fs[TENSOR_pFn].m_pdata), F(Fs[TENSOR_Fn].m_pdata);
      inv(pF, pFI);
      eF = F(i,k)*hFI(k,l)*pFI(l,j);      
      break;
    }
    case 2: // n+1
    {  
      TensorA<2> pF(Fs[TENSOR_pFnp1].m_pdata), F(Fs[TENSOR_Fnp1].m_pdata);
      inv(pF, pFI);
      eF = F(i,k)*hFI(k,l)*pFI(l,j);
      break;
    }
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);

  return err;      
}

int CM_PVP_PARAM::reset_state_vars_using_temporal(const Constitutive_model *m, 
                                                  State_variables *var)
const {
  int err = 0;
  Matrix<double> *Fs    = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs_in = var->Fs;
  double *state_var     = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *state_var_in  = var->state_vars[0].m_pdata;

  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fn   ].m_pdata[ia] = Fs_in[TENSOR_Fn   ].m_pdata[ia];
    Fs[TENSOR_pFn  ].m_pdata[ia] = Fs_in[TENSOR_pFn  ].m_pdata[ia];
    Fs[TENSOR_Fnm1 ].m_pdata[ia] = Fs_in[TENSOR_Fnm1 ].m_pdata[ia];
    Fs[TENSOR_pFnm1].m_pdata[ia] = Fs_in[TENSOR_pFnm1].m_pdata[ia];
  }
  
  state_var[VAR_pc_n]   = state_var_in[VAR_pc_n];
  state_var[VAR_pc_nm1] = state_var_in[VAR_pc_nm1];  
    
  return err;
}

int CM_PVP_PARAM::update_np1_state_vars_to_temporal(const Constitutive_model *m, 
                                                    State_variables *var)
const{
  int err = 0;
  Matrix<double> *Fs    = var->Fs;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
  double *state_var     = var->state_vars[0].m_pdata;
  double *state_var_in  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
    
  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fnp1 ].m_pdata[ia] = Fs_in[TENSOR_Fnp1 ].m_pdata[ia];
    Fs[TENSOR_pFnp1].m_pdata[ia] = Fs_in[TENSOR_pFnp1].m_pdata[ia];
  }
  
  state_var[VAR_pc_np1]   = state_var_in[VAR_pc_np1];
  

  return err;
}

int CM_PVP_PARAM::save_state_vars_to_temporal(const Constitutive_model *m, 
                                              State_variables *var)
const{
  int err = 0;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs    = var->Fs;
  double *state_var     = var->state_vars[0].m_pdata;
  double *state_var_in  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;  
  
  for(int ia=0; ia<TENSOR_end; ia++)
  {
    for(int ib=0; ib<DIM_3x3; ib++)
      Fs[ia].m_pdata[ib] = Fs_in[ia].m_pdata[ib];
  }
    
  for(int ia=0; ia<VAR_end; ia++)
    state_var[ia] = state_var_in[ia];
    
  return err;
}

int CM_PVP_PARAM::get_hardening(const Constitutive_model *m,
                                double *var,
                                const int stepno)
const 
{
  int err = 0;
  double *s_var = m->vars_list[0][m->model_id].state_vars->m_pdata;
  switch(stepno)
  {
    case 0: // n-1
      *var = s_var[VAR_pc_nm1];
      break;
    case 1: // n
      *var = s_var[VAR_pc_n];
      break;
    case 2: // n+1
      *var = s_var[VAR_pc_np1];
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;  
}

int CM_PVP_PARAM::get_plast_strain_var(const Constitutive_model *m,
                                       double *var)
const 
{
  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  *var = det3x3(Fs[TENSOR_pFn].m_pdata);
  return err;  
}


int CM_PVP_PARAM::write_restart(FILE *fp, const Constitutive_model *m)
const
{

  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);
                                                   
  fprintf(fp, "%.17e %.17e\n", state_var[VAR_pc_n], state_var[VAR_pc_nm1]);
  
  return err;
}

int CM_PVP_PARAM::read_restart(FILE *fp, Constitutive_model *m)
const
{
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  cm_read_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);
                                                        
  fscanf(fp, "%lf %lf\n", state_var+VAR_pc_n, state_var+VAR_pc_nm1);

  this->reset_state_vars(m);
  return 0;  
}

int CM_PVP_PARAM::compute_dMdu(const Constitutive_model *m,
                               CM_Ctx &cm_ctx,
                               double *Grad_op,
                               const int nne,
                               const int nsd,
                               double *dM_du)
const
{
  for(int ia=0; ia<nne*nsd*DIM_3x3; ia++)
    dM_du[ia] = 0.0;
  
  int err = 0;
  const double dt = cm_ctx.dt;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  MaterialPoroViscoPlasticity *mat_pvp = (m->param)->cm_mat->mat_pvp;
  GcmSolverInfo *solver_info = (m->param)->gcm_solver_info;

  if(cm_ctx.alpha<0){
    poro_visco_plasticity_compute_dMdu(dM_du, Grad_op, mat_pvp, solver_info,
                                       cm_ctx.F,
                                       Fs[TENSOR_Fn].m_pdata,
                                       Fs[TENSOR_pFnp1].m_pdata,
                                       Fs[TENSOR_pFn].m_pdata,
                                       vars[VAR_pc_np1],
                                       vars[VAR_pc_n],
                                       dt, nne, nsd);
  } else{
    double I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
    double pc_npa = 0.0;
    
    Tensor<2> Fnpa,pFnpa,hFnpa;

    mid_point_rule(Fnpa.data,  Fs[TENSOR_Fn].m_pdata,  cm_ctx.F,                   cm_ctx.alpha, DIM_3x3);
    mid_point_rule(pFnpa.data, Fs[TENSOR_pFn].m_pdata, Fs[TENSOR_pFnp1].m_pdata, cm_ctx.alpha, DIM_3x3);
    mid_point_rule(&pc_npa,    vars+VAR_pc_n,          vars+VAR_pc_np1,          cm_ctx.alpha, 1);
            
    if(cm_ctx.is_coulpled_with_thermal)
      mid_point_rule(hFnpa.data, cm_ctx.hFn, cm_ctx.hFnp1, cm_ctx.alpha, DIM_3x3);
    
    poro_visco_plasticity_compute_dMdu(dM_du, Grad_op, mat_pvp, solver_info,
                                       Fnpa.data,I,pFnpa.data,I,pc_npa,vars[VAR_pc_n],
                                       dt, nne, nsd);
  }
  
  return err;
}

/* THIS IS A FUNCTION STUB. */
int CM_PVP_PARAM::set_init_vals(Constitutive_model *m)
const
{
  // inital values are set in the more convoluted
  // read_constitutive_model_parameters->plasticity_model_read_parameters
  //   calling sequence
  
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *param      = (m->param)->model_param;
  double *vars       = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  double pJ   = param[PARAM_pJ];
  
  double pF11 = pow(pJ, 1.0/3.0);
  double F0[9] = {0.0,0.0,0.0,
                  0.0,0.0,0.0,
                  0.0,0.0,0.0};
  F0[0] = F0[4] = F0[8] =  pF11;

  vars[VAR_pc_nm1] = vars[VAR_pc_n] = vars[VAR_pc_np1] = param[PARAM_pc_0];
  if((m->param)->pF != NULL)
  {
    double pJ = det3x3((m->param)->pF);
    MaterialPoroViscoPlasticity *mat_pvp = (m->param)->cm_mat->mat_pvp;
    GcmSolverInfo *solver_info = (m->param)->gcm_solver_info;
    
    double pc = poro_visco_plasticity_compute_pc(pJ, param[PARAM_pc_0], mat_pvp, solver_info);
    vars[VAR_pc_nm1] = vars[VAR_pc_n] = vars[VAR_pc_np1] = pc;      
    memcpy(Fs[TENSOR_Fnm1 ].m_pdata, (m->param)->pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_Fn   ].m_pdata, (m->param)->pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_Fnp1 ].m_pdata, (m->param)->pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFnm1].m_pdata, (m->param)->pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFn  ].m_pdata, (m->param)->pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFnp1].m_pdata, (m->param)->pF, sizeof(double)*DIM_3x3);
  }
  else
  {
    memcpy(Fs[TENSOR_Fnm1 ].m_pdata, F0, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_Fn   ].m_pdata, F0, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_Fnp1 ].m_pdata, F0, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFnm1].m_pdata, F0, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFn  ].m_pdata, F0, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFnp1].m_pdata, F0, sizeof(double)*DIM_3x3);
  }    
  return 0;
}

int CM_PVP_PARAM::read_param(FILE *in)
const
{
  int err = 0;

  /* get pointer to parameter data */
  double *param     = this->model_param;
  int    *param_idx = this->model_param_index;
  assert(param     != NULL); // check the pointer
  assert(param_idx != NULL); // check the pointer

  /* scan to non-blank/comment line */
  err += scan_for_valid_line(in);

  /* READ PROPERTIES IN ALPHABETICAL ORDER */  
  int match =0;
  for(int ia=0; ia<PARAM_NO-4; ia++)
   match += fscanf(in, "%lf", param + ia);

  err += scan_for_valid_line(in);

  for(int ia=0; ia<2; ia++)
   match += fscanf(in, "%lf", param + PARAM_NO-4 + ia);
  
  err += scan_for_valid_line(in);
  
  match += fscanf(in, "%d %d %d %lf %lf", param_idx + PARAM_max_itr_stag,
                                          param_idx + PARAM_max_itr_M,
                                          param_idx + PARAM_max_subdivision,
                                          param     + PARAM_tol,
                                          param     + PARAM_computer_zero);


  set_gcm_solver_info(this->gcm_solver_info, 
                      param_idx[PARAM_max_itr_stag],
                      param_idx[PARAM_max_itr_M],
                      param_idx[PARAM_max_itr_M],
                      param[PARAM_tol],
                      param[PARAM_tol],
                      param[PARAM_computer_zero],
                      param_idx[PARAM_max_subdivision]);
                                                
  if(match != (PARAM_NO + PARAM_INX_NO))
  {   
    err++;
    assert(match == (PARAM_NO + PARAM_INX_NO) && "Did not read expected number of parameters"); 
  }
  err += scan_for_valid_line(in);
  
  /* not expecting EOF, check and return error if encountered */
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");

  this->cm_mat->mat_pvp = new MaterialPoroViscoPlasticity;
  set_properties_poro_visco_plasticity(this->cm_mat->mat_pvp,
                                       param[PARAM_yf_M],
                                       param[PARAM_yf_alpha],
                                       param[PARAM_flr_m],
                                       param[PARAM_flr_gamma0],
                                       param[PARAM_hr_a1],
                                       param[PARAM_hr_a2],
                                       param[PARAM_hr_Lambda1],
                                       param[PARAM_hr_Lambda2],
                                       param[PARAM_c_inf],
                                       param[PARAM_c_Gamma],
                                       param[PARAM_d_B],
                                       param[PARAM_d_pcb],
                                       param[PARAM_mu_0],
                                       param[PARAM_mu_1],
                                       param[PARAM_K_p0],
                                       param[PARAM_K_kappa],
                                       param[PARAM_pl_n],
                                       param[PARAM_cf_g0],
                                       param[PARAM_cf_pcinf]);
  
  return err;
}


int CM_PVP_PARAM::model_dependent_initialization(void)
{
  int err = 0;

  this->type              = POROVISCO_PLASTICITY;
  this->n_param           = PARAM_NO;
  this->model_param       = new double[PARAM_NO]();
  this->n_param_index     = PARAM_INX_NO;
  this->model_param_index = new int [PARAM_INX_NO]();
  this->gcm_solver_info   = new GcmSolverInfo;  
  
  return err;
}

int CM_PVP_PARAM::model_dependent_finalization(void)
{
  int err = 0;
  delete (this->cm_mat)->mat_pvp;
  delete this->gcm_solver_info;
  return err;
}
