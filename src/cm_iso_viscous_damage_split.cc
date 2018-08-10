/// This file defines the interface for the splited isotropic viscous
/// damage model.
/// 
/// REFERENCES:
/// 
/// Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
/// finite strains." Computational Mechanics 5.5 (1989): 375-400.
/// 
/// Authors:
///  Sangmin Lee, University of Notre Dame, <slee43@nd.edu>

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "allocation.h"
#include "cm_iso_viscous_damage_split.h"
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
  
  #define MIN(a,b) ((a)>(b)?(b):(a))
  #define MAX(a,b) ((a)>(b)?(a):(b))
  
  enum {TENSOR_Fnm1,
        TENSOR_Fn, 
        TENSOR_Fnp1, 
        TENSOR_end};
  
  enum {VAR_dw_np1,
        VAR_vw_np1,
        VAR_dH_np1,
        VAR_vH_np1,
        VAR_dX_np1,
        VAR_vX_np1,      
        VAR_dw_n,
        VAR_vw_n,
        VAR_dX_n,
        VAR_vX_n,
        VAR_dH_n,
        VAR_vH_n,
        VAR_dw_nm1,
        VAR_vw_nm1,
        VAR_dX_nm1,
        VAR_vX_nm1,            
        VAR_dH_nm1,
        VAR_vH_nm1,      
        VAR_end};      
  
  enum {FLAG_damaged_d_n,
        FLAG_damaged_v_n, 
        FLAG_damaged_d,
        FLAG_damaged_v, 
        FLAG_end};
  
  enum {PARAM_mu, 
        PARAM_w_max, 
        PARAM_P1, 
        PARAM_P2, 
        PARAM_Yin,
        PARAM_da,
        PARAM_db,
        PARAM_va,
        PARAM_vb, 
        PARAM_NO};
        
  // private context structure
  typedef struct {
    double *F;
    double dt;    // time increment
    double alpha; // mid point alpha
    double *eFnpa;
    int is_coulpled_with_thermal;
    double *hFn;
    double *hFnp1;
    int npa;
  } IvdsCtx;
}

/// Construct and initialize the poro-viscoplasticity model context 
/// for calling functions through the constitutive modeling interface
/// 
/// \param[in,out] ctx - handle to an opaque model context object.
/// \param[in] F The total deformation gradient.
/// \param[in] dt time increment
/// \param[in] alpha mid-point alpha
/// \param[in] eFnpa elastic deformation gradient at t = n + alpha
/// \param[in] hFn thermal part deformation gradient at t = n
/// \param[in] hFnp1 thermal part deformation gradient at t = n + 1
/// \param[in] is_coulpled_with_thermal flag for coupling with thermal
/// \return non-zero on internal error.
int iso_viscous_damage_model_split_ctx_build(void **ctx,
                                             double *F,
                                             const double dt,
                                             const double alpha,
                                             double *eFnpa,
                                             double *hFn,
                                             double *hFnp1,
                                             const int is_coulpled_with_thermal,
                                             const int npa)
{
  int err = 0;
  IvdsCtx *t_ctx = (IvdsCtx *) malloc(sizeof(IvdsCtx));

  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;
  t_ctx->hFn   = NULL;
  t_ctx->hFnp1 = NULL;  

  t_ctx->F = F;
  t_ctx->eFnpa = eFnpa;
  t_ctx->npa   = npa;  
  
  t_ctx->dt = dt;
  t_ctx->alpha = alpha;
  
  t_ctx->is_coulpled_with_thermal = is_coulpled_with_thermal;
  t_ctx->hFn  = hFn;
  t_ctx->hFnp1= hFnp1;  

  // assign handle
  *ctx = t_ctx;
  return err;
}

double ivds_compute_vw_npa(const Constitutive_model *m,
                           const int npa,
                           const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_vw_nm1, VAR_vw_n, VAR_vw_np1, npa, alpha);
}

double ivds_compute_dw_npa(const Constitutive_model *m,
                          const int npa,
                          const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_dw_nm1, VAR_dw_n, VAR_dw_np1, npa, alpha);
}

double ivds_compute_vH_npa(const Constitutive_model *m,
                           const int npa,
                           const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_vH_nm1, VAR_vH_n, VAR_vH_np1, npa, alpha);
}

double ivds_compute_dH_npa(const Constitutive_model *m,
                           const int npa,
                           const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_dH_nm1, VAR_dH_n, VAR_dH_np1, npa, alpha);
}

double ivds_compute_vX_npa(const Constitutive_model *m,
                           const int npa,
                           const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_vX_nm1, VAR_vX_n, VAR_vX_np1, npa, alpha);
}

double ivds_compute_dX_npa(const Constitutive_model *m,
                           const int npa,
                           const double alpha)
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->compute_state_vars_npa(VAR_dX_nm1, VAR_dX_n, VAR_dX_np1, npa, alpha);
}

int determine_damaged_npa(const int id_n,
                          const int id_np1,
                          const int *flags,
                          const int npa)
{
  int flag_npa = flags[id_np1];
  switch(npa)
  {
    case 0:
      flag_npa = flags[id_n];
      break;
    case 1:
      flag_npa = flags[id_np1];
      break;
  }
  return flag_npa;
} 

int determine_damaged_d_npa(const Constitutive_model *m,
                            const int npa)
{
  int *flags = m->vars_list[0][m->model_id].flags;
  return determine_damaged_npa(FLAG_damaged_d_n, FLAG_damaged_d, flags, npa);
}

int determine_damaged_v_npa(const Constitutive_model *m,
                            const int npa)
{
  int *flags = m->vars_list[0][m->model_id].flags;
  return determine_damaged_npa(FLAG_damaged_v_n, FLAG_damaged_v, flags, npa);
}

int CM_IVDS_PARAM::integration_algorithm(Constitutive_model *m,
                                         const void *ctx_in)
const
{
  
  int err = 0;
  auto ctx = (IvdsCtx *) ctx_in;

  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  int    *flags      = m->vars_list[0][m->model_id].flags;

  // integration algorithm
  ELASTICITY *elasticity = this->cm_elast;
  MATERIAL_CONTINUUM_DAMAGE *mat_d = this->cm_mat->mat_d;
      
  // store the deformation gradient
  memcpy(Fs[TENSOR_Fnp1].m_pdata, ctx->F, DIM_3x3 * sizeof(double));
  
  err += continuum_split_damage_integration_alg(mat_d,elasticity,
                                                vars+VAR_dw_np1,vars+VAR_vw_np1,
                                                vars+VAR_dX_np1,vars+VAR_vX_np1,
                                                vars+VAR_dH_np1,vars+VAR_vH_np1,
                                                flags+FLAG_damaged_d,
                                                flags+FLAG_damaged_v,
                                                vars[VAR_dw_n],vars[VAR_vw_n],
                                                vars[VAR_dX_n],vars[VAR_vX_n],
                                                ctx->dt, Fs[TENSOR_Fnp1].m_pdata); 

  return err;
}

int CM_IVDS_PARAM::compute_dev_stress(const Constitutive_model *m,
                                      const void *ctx_in,
                                      double *stress)
const
{
  return 0;
}

int CM_IVDS_PARAM::compute_dudj(const Constitutive_model *m,
                                const void *ctx_in,
                                double *dudj)
const
{
  int err = 0;
  auto ctx = (IvdsCtx *) ctx_in;
  double eJ = det3x3(ctx->F);
  
  double *vars  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;

  // scale by damage variable 
  ELASTICITY *elasticity = this->cm_elast;
  *dudj = damage_compute_dudj(elasticity, eJ, vars[VAR_vw_np1]);
  
  return err;
}

double CM_IVDS_PARAM::compute_dudj(const Constitutive_model *m,
                                   double theta_e,
                                   const int npa,
                                   const double alpha)
const 
{
  double w_npa = ivds_compute_vw_npa(m,npa,alpha);
  ELASTICITY *elasticity = this->cm_elast;
  return damage_compute_dudj(elasticity, theta_e, w_npa);          
}

int CM_IVDS_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                       const void *ctx_in,
                                       double *L)
const
{
  return 0;
}

int CM_IVDS_PARAM::compute_d2udj2(const Constitutive_model *m,
                                  const void *ctx_in,
                                  double *d2udj2)
const
{
  int err = 0;
  auto ctx = (IvdsCtx *) ctx_in;
  double eJ = det3x3(ctx->F);
  
  double *vars  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;

  // scale by damage variable 
  ELASTICITY *elasticity = this->cm_elast;
  *d2udj2 = damage_compute_d2udj2(elasticity, eJ, vars[VAR_vw_np1]);
  
  return err;
}

double CM_IVDS_PARAM::compute_d2udj2(const Constitutive_model *m,
                                     double theta_e,
                                     const int npa,
                                     const double alpha)
const 
{
  double w_npa = ivds_compute_vw_npa(m,npa,alpha);
  ELASTICITY *elasticity = this->cm_elast;
  return damage_compute_d2udj2(elasticity, theta_e, w_npa);           
}

int CM_IVDS_PARAM::update_elasticity(const Constitutive_model *m,
                                     const void *ctx_in,
                                     double *L,
                                     double *S,
                                     const int compute_stiffness)
const
{
  int err = 0;
  auto ctx = (IvdsCtx *) ctx_in;
  
  //double *vars  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  //int *flags = m->vars_list[0][m->model_id].flags;
  
  ELASTICITY *elasticity = this->cm_elast;
  MATERIAL_CONTINUUM_DAMAGE *mat_d = this->cm_mat->mat_d;
  
  double *S_tmp = elasticity->S;
  double *L_tmp = elasticity->L;
  
  elasticity->S = S;
  elasticity->L = L;
  
  double dw_npa = ivds_compute_dw_npa(m, ctx->npa, ctx->alpha);
  double vw_npa = ivds_compute_vw_npa(m, ctx->npa, ctx->alpha);
  double dH_npa = ivds_compute_dH_npa(m, ctx->npa, ctx->alpha);
  double vH_npa = ivds_compute_vH_npa(m, ctx->npa, ctx->alpha);
  int damage_d  = determine_damaged_d_npa(m, ctx->npa);
  int damage_v  = determine_damaged_v_npa(m, ctx->npa);
    
  if(ctx->eFnpa)
  {
    err += update_split_damage_elasticity(mat_d, elasticity, 
                                          dw_npa, vw_npa,
                                          dH_npa, vH_npa,
                                          damage_d,
                                          damage_v,
                                          ctx->dt, 
                                          ctx->eFnpa, 
                                          compute_stiffness);
  }    
  else
  {
    // shorthand of deformation gradients
    Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;  
    Tensor<2> eF = {};  
    
    if(ctx->is_coulpled_with_thermal){
      TensorA<2> Fnp1(Fs[TENSOR_Fnp1].m_pdata);
      TensorA<2> hFnp1(ctx->hFnp1);
      Tensor<2> hFnp1_I;

      err += inv(hFnp1, hFnp1_I);
      eF = Fnp1(i,k)*hFnp1_I(k,j);
  
      err += update_split_damage_elasticity(mat_d, elasticity, 
                                            dw_npa, vw_npa,
                                            dH_npa, vH_npa,
                                            damage_d, damage_v,
                                            ctx->dt, 
                                            eF.data, 
                                            compute_stiffness);
    }
    else{
      err += update_split_damage_elasticity(mat_d, elasticity, 
                                            dw_npa, vw_npa,
                                            dH_npa, vH_npa,
                                            damage_d, damage_v,
                                            ctx->dt, 
                                            Fs[TENSOR_Fnp1].m_pdata, 
                                            compute_stiffness);
    }      
  }
  
  elasticity->S = S_tmp;
  elasticity->L = L_tmp;   
  return err;
}

int CM_IVDS_PARAM::update_elasticity_dev(const Constitutive_model *m,
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
  
  double dw_npa = ivds_compute_dw_npa(m, npa, alpha);
  double dH_npa = ivds_compute_dH_npa(m, npa, alpha);
  int damage_d  = determine_damaged_d_npa(m, npa);
    
  err += update_damage_elasticity_dev(mat_d, elasticity, 
                                      dw_npa,
                                      dH_npa,
                                      damage_d,
                                      dt, 
                                      eFnpa, 
                                      compute_stiffness);
  elasticity->S = S_tmp;
  elasticity->L = L_tmp;                                      
  return err; 
}

int CM_IVDS_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  
  // deformation gradients
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;  
  Fs[TENSOR_Fnm1] = Fs[TENSOR_Fn];
  Fs[TENSOR_Fn]   = Fs[TENSOR_Fnp1];

  // state variables 
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;  
  vars[VAR_dw_nm1] = vars[VAR_dw_n];
  vars[VAR_dH_nm1] = vars[VAR_dH_n];
  vars[VAR_dX_nm1] = vars[VAR_dX_n];
  vars[VAR_vw_nm1] = vars[VAR_vw_n];
  vars[VAR_vH_nm1] = vars[VAR_vH_n];
  vars[VAR_vX_nm1] = vars[VAR_vX_n];  
    
  vars[VAR_dw_n] = vars[VAR_dw_np1];
  vars[VAR_dH_n] = vars[VAR_dH_np1];
  vars[VAR_dX_n] = vars[VAR_dX_np1];
  vars[VAR_vw_n] = vars[VAR_vw_np1];
  vars[VAR_vH_n] = vars[VAR_vH_np1];
  vars[VAR_vX_n] = vars[VAR_vX_np1];  

  // flags 
  int *flags = m->vars_list[0][m->model_id].flags;  
  flags[FLAG_damaged_d_n] = flags[FLAG_damaged_d];
  flags[FLAG_damaged_v_n] = flags[FLAG_damaged_v];
  return err;
}

int CM_IVDS_PARAM::reset_state_vars(Constitutive_model *m)
const 
{
  int err = 0;

  // deformation gradients 
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[TENSOR_Fnp1]  = Fs[TENSOR_Fn];

  // state variables
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;  
  vars[VAR_dw_np1] = vars[VAR_dw_n];
  vars[VAR_dH_np1] = vars[VAR_dH_n];
  vars[VAR_dX_np1] = vars[VAR_dX_n];
  vars[VAR_vw_np1] = vars[VAR_vw_n];
  vars[VAR_vH_np1] = vars[VAR_vH_n];
  vars[VAR_vX_np1] = vars[VAR_vX_n];  

  // flags
  int *flags = m->vars_list[0][m->model_id].flags;  
  flags[FLAG_damaged_d] = flags[FLAG_damaged_d_n];
  flags[FLAG_damaged_v] = flags[FLAG_damaged_v_n];
  
  return err;
}

int CM_IVDS_PARAM::get_subdiv_param(const Constitutive_model *m,
                                   double *subdiv_param,
                                   double dt)
const
{
  int err = 0;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  *subdiv_param = compute_subdiv_param4split_damage(vars[VAR_dw_n], vars[VAR_dw_np1],
                                                    vars[VAR_vw_n], vars[VAR_vw_np1]);
  return err;
}

int CM_IVDS_PARAM::get_var_info(Model_var_info &info)
const
{
  const CMVariableNames variable_names[] = {{VAR_dw_np1,"dw_np1"},
                                            {VAR_vw_np1,"vw_np1"},
                                            {VAR_dH_np1,"dH_np1"},
                                            {VAR_vH_np1,"vH_np1"},
                                            {VAR_dX_np1,"dX_np1"},
                                            {VAR_vX_np1,"vX_np1"},
                                            {VAR_dw_n,  "dw_n"},
                                            {VAR_vw_n,  "vw_n"},
                                            {VAR_dX_n,  "dX_n"},
                                            {VAR_vX_n,  "vX_n"},
                                            {VAR_dH_n,  "dH_n"},
                                            {VAR_vH_n,  "vH_n"},
                                            {VAR_dw_nm1,"dw_nm1"},
                                            {VAR_vw_nm1,"vw_nm1"},
                                            {VAR_dX_nm1,"dX_nm1"},
                                            {VAR_vX_nm1,"vX_nm1"},
                                            {VAR_dH_nm1,"dH_nm1"},
                                            {VAR_vH_nm1,"vH_nm1"},
                                           };
    
  const CMVariableNames tensor_names[] = {{TENSOR_Fnm1, "Fnm1"},
                                          {TENSOR_Fn,   "Fn"},
                                          {TENSOR_Fnp1, "Fnp1"} 
                                         };
  const CMVariableNames flag_names[] = {{FLAG_damaged_d_n, "damaged_d_n"},
                                        {FLAG_damaged_v_n, "damaged_v_n"},
                                        {FLAG_damaged_d,   "damaged_d"  },
                                        {FLAG_damaged_v,   "damaged_v"  },
                                       };
  
  int err = constitutive_model_info(info, VAR_end,    variable_names,
                                          TENSOR_end, tensor_names,
                                          FLAG_end,   flag_names);

  return err;
}

int ivds_get_F(const Constitutive_model *m,
               double *F,
               const int stepno)
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
      memcpy(F,Fs[TENSOR_Fnm1].m_pdata,DIM_3x3*sizeof(double));
      break;
    case 1: // n
      memcpy(F,Fs[TENSOR_Fn].m_pdata,  DIM_3x3*sizeof(double));
      break;
    case 2: // n+1
      memcpy(F,Fs[TENSOR_Fnp1].m_pdata,DIM_3x3*sizeof(double));
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;
}


int CM_IVDS_PARAM::get_F(const Constitutive_model *m,
                        double *F_out,
                        const int stepno)
const
{
  return ivds_get_F(m, F_out, stepno);
}


int CM_IVDS_PARAM::get_eF(const Constitutive_model *m,
                          double *eF_in,
                          const int stepno)
const
{
  return ivds_get_F(m, eF_in, stepno);
}

int CM_IVDS_PARAM::get_pF(const Constitutive_model *m,
                          double *F,
                          const int stepno)
const
{
  F[0] = F[4] = F[8] = 1.0;
  F[1] = F[2] = F[3] = F[5] = F[6] = F[7] = 0.0;
  
  return 0;
}

int CM_IVDS_PARAM::get_eF_of_hF(const Constitutive_model *m,
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

int CM_IVDS_PARAM::reset_state_vars_using_temporal(const Constitutive_model *m, 
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
  vars[VAR_dw_n] = vars_in[VAR_dw_n];
  vars[VAR_vw_n] = vars_in[VAR_vw_n];
  vars[VAR_dH_n] = vars_in[VAR_dH_n];  
  vars[VAR_vH_n] = vars_in[VAR_vH_n];
  vars[VAR_dX_n] = vars_in[VAR_dX_n];
  vars[VAR_vX_n] = vars_in[VAR_vX_n];  
  
  vars[VAR_dw_nm1] = vars_in[VAR_dw_nm1];
  vars[VAR_vw_nm1] = vars_in[VAR_vw_nm1];
  vars[VAR_dH_nm1] = vars_in[VAR_dH_nm1];  
  vars[VAR_vH_nm1] = vars_in[VAR_vH_nm1];
  vars[VAR_dX_nm1] = vars_in[VAR_dX_nm1];
  vars[VAR_vX_nm1] = vars_in[VAR_vX_nm1];  
  
  
  int *flags    = m->vars_list[0][m->model_id].flags;
  int *flags_in = var->flags;
  
  flags[FLAG_damaged_d_n] = flags_in[FLAG_damaged_d_n];
  flags[FLAG_damaged_v_n] = flags_in[FLAG_damaged_v_n];  
  
  return err;
}

int CM_IVDS_PARAM::update_np1_state_vars_to_temporal(const Constitutive_model *m, 
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
  vars[VAR_dw_np1] = vars_in[VAR_dw_np1];
  vars[VAR_vw_np1] = vars_in[VAR_vw_np1];
  vars[VAR_dH_np1] = vars_in[VAR_dH_np1];
  vars[VAR_vH_np1] = vars_in[VAR_vH_np1];
  vars[VAR_dX_np1] = vars_in[VAR_dX_np1];  
  vars[VAR_vX_np1] = vars_in[VAR_vX_np1];  
  
  int *flags_in = m->vars_list[0][m->model_id].flags;
  int *flags    = var->flags;
  
  flags[FLAG_damaged_d] = flags_in[FLAG_damaged_d];
  flags[FLAG_damaged_v] = flags_in[FLAG_damaged_v];
  
  return err;
}

int CM_IVDS_PARAM::save_state_vars_to_temporal(const Constitutive_model *m, 
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

int CM_IVDS_PARAM::get_hardening(const Constitutive_model *m,
                                double *var,
                                const int stepno)
const
{
  int err = 0;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  switch(stepno)
  {
    case 0: // n-1
      *var = 0.5*(vars[VAR_dw_nm1] + vars[VAR_vw_nm1]);
      break;
    case 1: // n
     *var = 0.5*(vars[VAR_dw_n] + vars[VAR_vw_n]);
      break;
    case 2: // n+1
     *var = 0.5*(vars[VAR_dw_np1] + vars[VAR_vw_np1]);
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;  
}

int CM_IVDS_PARAM::get_plast_strain_var(const Constitutive_model *m,
                                       double *chi)
const                                   
{
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  *chi = vars[VAR_dX_n] + vars[VAR_vX_n];
  return 0;
}

int CM_IVDS_PARAM::write_restart(FILE *out, const Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *vars  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  const int *flags = (m->vars_list[0][m->model_id]).flags;
  
  err += cm_write_tensor_restart(out, Fs[TENSOR_Fn].m_pdata);
  err += cm_write_tensor_restart(out, Fs[TENSOR_Fnm1].m_pdata);
    
  fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e\n", 
                             vars[VAR_dw_n], vars[VAR_vw_n], 
                             vars[VAR_dX_n], vars[VAR_vX_n], 
                             vars[VAR_dH_n], vars[VAR_vH_n]);
  fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e\n", 
                             vars[VAR_dw_nm1], vars[VAR_vw_nm1], 
                             vars[VAR_dX_nm1], vars[VAR_vX_nm1], 
                             vars[VAR_dH_nm1], vars[VAR_vH_nm1]);                             
                                                          
  fprintf(out, "%d %d\n", flags[FLAG_damaged_d_n], flags[FLAG_damaged_v_n]);

  return err;
}

int CM_IVDS_PARAM::read_restart(FILE *in, Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *vars  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  int *flags = (m->vars_list[0][m->model_id]).flags;
  
  cm_read_tensor_restart(in, Fs[TENSOR_Fn].m_pdata);
  cm_read_tensor_restart(in, Fs[TENSOR_Fnm1].m_pdata);
    
  fscanf(in, "%lf %lf %lf %lf %lf %lf", vars+VAR_dw_n, vars+VAR_vw_n, 
                                        vars+VAR_dX_n, vars+VAR_vX_n, 
                                        vars+VAR_dH_n, vars+VAR_vH_n);

  fscanf(in, "%lf %lf %lf %lf %lf %lf", vars+VAR_dw_nm1, vars+VAR_vw_nm1, 
                                        vars+VAR_dX_nm1, vars+VAR_vX_nm1, 
                                        vars+VAR_dH_nm1, vars+VAR_vH_nm1);
                                        
  fscanf(in, "%d %d", flags + FLAG_damaged_d_n, flags + FLAG_damaged_v_n);

  err += this->reset_state_vars(m);
  return err;
}

int CM_IVDS_PARAM::destroy_ctx(void **ctx)
const
{
  int err = 0;
  IvdsCtx *t_ctx = (IvdsCtx *) *ctx;
  // invalidate handle 
  *ctx = NULL;

  // no memory was created
  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;
  t_ctx->hFn   = NULL;
  t_ctx->hFnp1 = NULL; 

  free(t_ctx);
  return err;
}

int CM_IVDS_PARAM::set_init_vals(Constitutive_model *m)
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

int CM_IVDS_PARAM::read_param(FILE *in)
const
{
  int err = 0;
  // get pointer to parameter data
  double *param     = this->model_param;
  int    *param_idx = this->model_param_index;
  assert(param     != NULL); // check the pointer
  assert(param_idx != NULL); // check the pointer

  // scan to non-blank/comment line 
  err += scan_for_valid_line(in);

  // read properties
  int match = fscanf(in, "%lf %lf %lf %lf %lf",
                     param + PARAM_mu, param + PARAM_w_max, param + PARAM_P1,
                     param + PARAM_P2, param + PARAM_Yin);

  err += scan_for_valid_line(in);
  
  match += fscanf(in, "%lf %lf %lf %lf",
                  param + PARAM_da, param + PARAM_db, 
                  param + PARAM_va, param + PARAM_vb);                       
                     
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

int CM_IVDS_PARAM::model_dependent_initialization(void)
{
  int err = 0;
  this->type              = ISO_VISCOUS_SPLIT_DAMAGE;
  this->n_param           = PARAM_NO;
  this->model_param       = new double[PARAM_NO]();
  this->n_param_index     = 0;
  this->model_param_index = NULL;
  this->gcm_solver_info   = NULL;

  return err;
}

int CM_IVDS_PARAM::model_dependent_finalization(void)
{
  int err = 0;
  delete (this->cm_mat)->mat_d;
  return err;
}