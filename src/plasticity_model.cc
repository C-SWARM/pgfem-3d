/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 *  Sangmin Lee, University of Notre Dame, Notre Dame, IN, <slee43@nd.edu>
 *  Aaron Howell, University of Notre Dame, Notre Dame, IN, <ahowell3@nd.edu>
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "plasticity_model.h"
#include "allocation.h"
#include "cm_placeholder_functions.h"
#include "constitutive_model.h"
#include "crystal_plasticity_integration.h"
#include "elem3d.h"
#include "flowlaw.h"
#include "gen_path.h"
#include "hyperelasticity.h"
#include "material_properties.h"
#include "new_potentials.h"
#include "index_macros.h"
#include "utils.h"
#include <ttl/ttl.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <xmmintrin.h>

namespace {
const constexpr int          DIM_3 = 3;
const constexpr int        DIM_3x3 = 9;
const constexpr int      DIM_3x3x3 = 27;
const constexpr int    DIM_3x3x3x3 = 81;
const constexpr double MAX_D_GAMMA = 0.005;
const constexpr int       FLAG_end = 0;

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
  VAR_L_n,
  VAR_L_np1,
  VAR_g_n,
  VAR_g_np1,
  VAR_L_nm1,
  VAR_g_nm1,
  VAR_end
};

enum tensor_names {
  TENSOR_Fn,
  TENSOR_pFn,
  TENSOR_Fnp1,
  TENSOR_pFnp1,
  TENSOR_tau,
  TENSOR_R,
  TENSOR_gamma_dot,
  TENSOR_Fnm1,
  TENSOR_pFnm1,
  TENSOR_tau_n,
  TENSOR_gamma_dot_n,
  TENSOR_end
};

enum param_names {
  PARAM_gamma_dot_0,
  PARAM_m,
  PARAM_G0,
  PARAM_g0,
  PARAM_gs_0,
  PARAM_gamma_dot_s,
  PARAM_w,
  PARAM_tol_hardening,
  PARAM_tol_M,
  PARAM_computer_zero,
  PARAM_NO
};

enum param_index_names {
  PARAM_max_itr_stag,
  PARAM_max_itr_hardening,
  PARAM_max_itr_M,
  PARAM_max_subdivision,
  PARAM_INX_NO
};

typedef struct {
  int e_id;
  int n_ip;
  int mat_id;
  Matrix<int> ip_ids;
} IP_ID_LIST;

int plasticity_model_construct_elem_ip_map(IP_ID_LIST *elm_ip_map, EPS *eps, const Element *elem, int ne)
{
  int cnt = 0;
  for(int a=0; a<ne; a++)
  {
    long n_ip = 0;
    int_point(elem[a].toe,&n_ip);
    elm_ip_map[a].ip_ids.initialization(n_ip, 1);
    elm_ip_map[a].e_id = a;
    elm_ip_map[a].n_ip = n_ip;
    elm_ip_map[a].mat_id = elem[a].mat[0];

    for(int b=0; b<n_ip; b++)
    {
      elm_ip_map[a].ip_ids(b) = cnt;
      cnt++;
    }
  }
  return cnt;
}

/// Private structure for use exclusively with this model and
/// associated functions.
typedef struct plasticity_ctx {
  double *F;
  double dt;    // time increment
  double alpha; // mid point alpha
  double *eFnpa;
  int is_coulpled_with_thermal;
  double *hFn;
  double *hFnp1;
} plasticity_ctx;

int compute_M(double *M_in,
              double *pFn_in,
              double *N_in,
              double *pFnp1_I_in,
              const plasticity_ctx *ctx)
{
  int err = 0;

  TensorA<2> M(M_in), pFn(pFn_in), N(N_in), pFnp1_I(pFnp1_I_in);

  if(ctx->is_coulpled_with_thermal)
    M = pFn(i,k)*N(k,l)*pFnp1_I(l,j);
  else
    M = pFn(i,k)*pFnp1_I(k,j);
  return err;
}

int compute_eF(double *eF_in,
               double *F_in,
               double *hF_I_in,
               double *pF_I_in,
               const plasticity_ctx *ctx)
{
  int err = 0;
  TensorA<2> eF(eF_in), F(F_in), hF_I(hF_I_in), pF_I(pF_I_in);
  if(ctx->is_coulpled_with_thermal)
    eF = F(i,k)*hF_I(k,l)*pF_I(l,j);
  else
    eF = F(i,k)*pF_I(k,j);

  return err;
}

int CP_PARAM::update_state_vars(Constitutive_model *m)
  const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Fs[TENSOR_Fnm1] = Fs[TENSOR_Fn];
  Fs[TENSOR_Fn]   = Fs[TENSOR_Fnp1];
  Fs[TENSOR_pFnm1]= Fs[TENSOR_pFn];
  Fs[TENSOR_pFn]  = Fs[TENSOR_pFnp1];

  Fs[TENSOR_gamma_dot_n] = Fs[TENSOR_gamma_dot];
  Fs[TENSOR_tau_n]       = Fs[TENSOR_tau];
  state_var[VAR_g_nm1] = state_var[VAR_g_n];
  state_var[VAR_g_n]   = state_var[VAR_g_np1];
  state_var[VAR_L_nm1] = state_var[VAR_L_n];
  state_var[VAR_L_n]   = state_var[VAR_L_np1];
  return err;
}

int CP_PARAM::reset_state_vars(Constitutive_model *m)
  const 
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Fs[TENSOR_Fnp1]  = Fs[TENSOR_Fn];
  Fs[TENSOR_pFnp1] = Fs[TENSOR_pFn];
  state_var[VAR_g_np1] = state_var[VAR_g_n];
  state_var[VAR_L_np1] = state_var[VAR_L_n];
  return err;
}

int CP_PARAM::reset_state_vars_using_temporal(const Constitutive_model *m,
                                              State_variables *var)
  const 
{
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

  SLIP_SYSTEM *slip = (((m->param)->cm_mat)->mat_p)->slip;

  int N_SYS = slip->N_SYS;

  for(int ia=0; ia<N_SYS; ia++)
  {
    Fs[TENSOR_tau_n      ].m_pdata[ia] = Fs_in[TENSOR_tau_n      ].m_pdata[ia];
    Fs[TENSOR_gamma_dot_n].m_pdata[ia] = Fs_in[TENSOR_gamma_dot_n].m_pdata[ia];
  }

  state_var[VAR_g_n]   = state_var_in[VAR_g_n];
  state_var[VAR_L_n]   = state_var_in[VAR_L_n];
  state_var[VAR_g_nm1] = state_var_in[VAR_g_nm1];
  state_var[VAR_L_nm1] = state_var_in[VAR_L_nm1];

  return err;
}

int CP_PARAM::update_np1_state_vars_to_temporal(const Constitutive_model *m,
                                                State_variables *var)
  const 
{
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

  SLIP_SYSTEM *slip = (((m->param)->cm_mat)->mat_p)->slip;

  int N_SYS = slip->N_SYS;

  for(int ia=0; ia<N_SYS; ia++)
  {
    Fs[TENSOR_tau      ].m_pdata[ia] = Fs_in[TENSOR_tau      ].m_pdata[ia];
    Fs[TENSOR_gamma_dot].m_pdata[ia] = Fs_in[TENSOR_gamma_dot].m_pdata[ia];
  }

  state_var[VAR_g_np1]   = state_var_in[VAR_g_np1];
  state_var[VAR_L_np1]   = state_var_in[VAR_L_np1];


  return err;
}

int CP_PARAM::save_state_vars_to_temporal(const Constitutive_model *m,
                                          State_variables *var)
  const 
{
  int err = 0;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs    = var->Fs;
  double *state_var     = var->state_vars[0].m_pdata;
  double *state_var_in  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;

  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fn   ].m_pdata[ia] = Fs_in[TENSOR_Fn   ].m_pdata[ia];
    Fs[TENSOR_pFn  ].m_pdata[ia] = Fs_in[TENSOR_pFn  ].m_pdata[ia];
    Fs[TENSOR_Fnp1 ].m_pdata[ia] = Fs_in[TENSOR_Fnp1 ].m_pdata[ia];
    Fs[TENSOR_pFnp1].m_pdata[ia] = Fs_in[TENSOR_pFnp1].m_pdata[ia];
    Fs[TENSOR_R    ].m_pdata[ia] = Fs_in[TENSOR_R    ].m_pdata[ia];
    Fs[TENSOR_Fnm1 ].m_pdata[ia] = Fs_in[TENSOR_Fnm1 ].m_pdata[ia];
    Fs[TENSOR_pFnm1].m_pdata[ia] = Fs_in[TENSOR_pFnm1].m_pdata[ia];
  }

  // if size is not same redim size
  Fs[TENSOR_tau        ] = Fs_in[TENSOR_tau        ];
  Fs[TENSOR_gamma_dot  ] = Fs_in[TENSOR_gamma_dot  ];
  Fs[TENSOR_tau_n      ] = Fs_in[TENSOR_tau_n      ];
  Fs[TENSOR_gamma_dot_n] = Fs_in[TENSOR_gamma_dot_n];

  for(int ia=0; ia<VAR_end; ia++)
    state_var[ia] = state_var_in[ia];

  return err;
}

int CP_PARAM::get_var_info(Model_var_info &info)
  const 
{
  int Fno = TENSOR_end;

  info.n_Fs = Fno;
  info.F_names = (char **)malloc(sizeof(char*)*Fno);
  for(int a=0; a<Fno; a++)
    info.F_names[a] = (char *)malloc(sizeof(char)*1024);

  sprintf(info.F_names[TENSOR_Fn],        "Fn");
  sprintf(info.F_names[TENSOR_pFn],       "pFn");
  sprintf(info.F_names[TENSOR_Fnp1],      "Fnp1");
  sprintf(info.F_names[TENSOR_pFnp1],     "pFnp1");
  sprintf(info.F_names[TENSOR_tau],       "tau");
  sprintf(info.F_names[TENSOR_R],         "R");
  sprintf(info.F_names[TENSOR_gamma_dot], "gamma_dot");
  sprintf(info.F_names[TENSOR_Fnm1],      "Fnm1");
  sprintf(info.F_names[TENSOR_pFnm1],     "pFnm1");

  int varno = VAR_end;
  info.n_vars = varno;
  info.var_names = (char **)malloc(sizeof(char*)*varno);
  for(int a=0; a<varno; a++)
    info.var_names[a] = (char *)malloc(sizeof(char)*1024);

  sprintf(info.var_names[VAR_L_n],        "L_n");
  sprintf(info.var_names[VAR_L_np1],      "L_np1");
  sprintf(info.var_names[VAR_g_n],        "g_n");
  sprintf(info.var_names[VAR_g_np1],      "g_np1");

  sprintf(info.var_names[VAR_L_nm1],      "L_nm1");
  sprintf(info.var_names[VAR_g_nm1],      "g_nm1");

  info.n_flags = FLAG_end;
  info.flag_names = (char **) malloc(FLAG_end * sizeof(char*));  
  for(int a=0; a<FLAG_end; a++)
    info.flag_names[a] = (char *)malloc(sizeof(char)*1024);

  return 0;
}

int CP_PARAM::get_eF_of_hF(const Constitutive_model *m,
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

int CP_PARAM::get_pF(const Constitutive_model *m,
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

int CP_PARAM::get_F(const Constitutive_model *m,
                    double *F,
                    const int stepno)
  const 
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CP_PARAM::set_F(const Constitutive_model *m,
                    double *F,
                    const int stepno)
const
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CP_PARAM::get_eF(const Constitutive_model *m,
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

int CP_PARAM::get_hardening(const Constitutive_model *m,
                            double *var,
                            const int stepno)
  const 
{
  int err = 0;
  double *s_var = m->vars_list[0][m->model_id].state_vars->m_pdata;
  switch(stepno)
  {
   case 0: // n-1
    *var = s_var[VAR_g_nm1];
    break;
   case 1: // n
    *var = s_var[VAR_g_n];
    break;
   case 2: // n+1
    *var = s_var[VAR_g_np1];
    break;
   default:
    PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
    err++;
  }
  assert(err == 0);
  return err;
}

int CP_PARAM::compute_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 double *stress)
  const 
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;

  TensorA<2> F(CTX->F);
  Tensor<2> eC = F(k,i)*F(k,j);

  devStressFuncPtr Stress = getDevStressFunc(-1, m->param->p_hmat);
  Stress(eC.data, m->param->p_hmat, stress);
  return err;
}

int CP_PARAM::compute_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *dudj)
  const 
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
  TensorA<2> F(CTX->F);
  double J = ttl::det(F);
  Pressure(J,m->param->p_hmat,dudj);
  return err;
}

int CP_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  double *L)
  const 
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  TensorA<2> F(CTX->F);
  Tensor<2> eC = F(k,i)*F(k,j);

  matStiffFuncPtr Tangent = getMatStiffFunc(-1, m->param->p_hmat);
  Tangent(eC.data, m->param->p_hmat, L);
  return err;
}

int CP_PARAM::compute_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *d2udj2)
  const 
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  TensorA<2> F(CTX->F);
  double J = ttl::det(F);

  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
  D_Pressure(J,m->param->p_hmat,d2udj2);
  return err;
}

/// If drdtau is denormal it is set to zero. Any double precision number (64 bits)
/// is considered as denormal if it is smaller than ±2.23×10^308 
/// https://software.intel.com/en-us/node/523328 
inline double compute_drdtau(const double gamma_dot_0,
                             const double mm,
                             const double g_np1,
                             const double tau)
{
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  double drdtau = gamma_dot_0/mm/g_np1*pow(fabs(tau/g_np1), 1.0/mm - 1.0);
  return drdtau;
}

int plasticity_compute_dMdu(const Constitutive_model *con,
                            double *dMdu_in,
                            double *Grad_du_in,
                            double *eFn_in,
                            double *eFnp1_in,
                            double *M_in,
                            double *S_in,
                            double *L_in,
                            const double g_n,
                            const double g_np1,
                            const double *tau,
                            const double *gamma_dots,
                            double *Psys,
                            const double dt)
{
  // compute dMdu:U = -grad(du):B
  // Grad_du = Grad(du)

  MATERIAL_CRYSTAL_PLASTICITY *mat_p = ((con->param)->cm_mat)->mat_p;

  const int N_SYS          = (mat_p->slip)->N_SYS;
  const double gamma_dot_0 = mat_p->gamma_dot_0;
  const double gamma_dot_s = mat_p->gamma_dot_s;
  const double mm          = mat_p->m;
  const double g0          = mat_p->g0;
  const double G0          = mat_p->G0;
  const double gs_0        = mat_p->gs_0;
  const double w           = mat_p->w;

  double gamma_dot = 0.0;
  for(int a = 0; a<N_SYS; a++)
    gamma_dot += fabs(gamma_dots[a]);

  double gm_gms   = gamma_dot/gamma_dot_s;
  double sign_gm_gms = (gm_gms < 0) ? -1.0 : 1.0;

  double R3 = 0.0;
  double gs_np1 = 0.0;
  if(fabs(gm_gms)>1.0e-15)
  {
    R3 = gs_0*w/gamma_dot_s*sign_gm_gms*pow(fabs(gm_gms), w-1.0);
    gs_np1 = gs_0*pow(fabs(gm_gms),w);
  }

  double AA = R3*gamma_dot*(g_n - g0 + dt*G0*gamma_dot) +
              gs_np1 * (gs_np1 - g0 - g_n) + g0*g_n;
  double BB = gs_np1 - g0  - dt*G0*gamma_dot;
  double R4 = dt*G0*AA/BB/BB;

  double sum_1gm1gm = 0.0;

  // convert C struct matrices into ttl tensors
  const TensorA<2> Grad_du(Grad_du_in);
  const TensorA<2> eFn(eFn_in);
  const TensorA<2> eFnp1(eFnp1_in);
  const TensorA<2> M(M_in);
  const TensorA<2> S(S_in);
  const TensorA<4> L(L_in);

  const Tensor<2> C = eFnp1(k,i).to(i,k) * eFnp1(k,j);             // eFnp1' * eFnp1
  Tensor<2> MI = {};
  int err = inv(M, MI);
  const Tensor<2> MI_x_C = MI(j,i) * C(j,k);                       // M^{-T} * C
  const Tensor<2> M_x_eFn = M(j,i).to(i,j) * eFn(k,j).to(j,k);     // M' * eFn'

  // tensors are initialized to 0
  Tensor<2> sum_aC = {};
  Tensor<2> sum_Pa = {};
  Tensor<2> sum_aD = {};
  Tensor<4> U = {};
  Tensor<4> B = {};

  // set U to the 9x9 identity scaled by 1.0/dt
  U(i,j,k,l) =  ttl::identity(i,j,k,l);

  for(int a = 0; a<N_SYS; a++)
  {
    double drdtau = compute_drdtau(gamma_dot_0, mm, g_np1, tau[a]);
    double drdg   = -drdtau*tau[a]/g_np1;

    double R2_a = ((gamma_dots[a] < 0) ? -1.0 : 1.0)*drdtau;
    sum_1gm1gm += ((gamma_dots[a] < 0) ? -1.0 : 1.0)*drdg;

    // compute P alpha of Psys
    const Tensor<2, const double*> Pa(Psys + DIM_3x3 * a);

    // compute C alpha and D alpha using ttl operations
    auto t0 = Pa(i,k) * S(k,j);                                    // Pa * S
    auto t1 = S(i,k) * Pa(j,k).to(k,j);                            // S * Pa'
    auto t2 = L(i,k,n,l) * C(n,l) * Pa(k,j);                       // L:C * Pa
    const Tensor<2> AA = t0 + t1 + t2;
    const Tensor<2> aC = MI_x_C(i,j) * AA(j,k);
    const Tensor<2> aD = eFnp1(i,j) * AA(j,k) * M_x_eFn(k,l);

    sum_Pa(i,j) += drdg * Pa(i,j);
    sum_aC(i,j) += R2_a * aC(i,j);
    sum_aD(i,j) += R2_a * aD(i,j);

    // perform the Kronecker product using ttl and scales it by drdtau
    U(i,j,k,l) += dt*drdtau * aC(i,j) * Pa(k,l);
    B(i,j,k,l) += dt*drdtau * aD(i,j) * Pa(k,l);
  }

  double R1 = R4/(1.0-R4*sum_1gm1gm);

  // perform the Kronecker product using ttl and scales it by R1
  U(i,j,k,l) += dt*R1*sum_aC(i,j)*sum_Pa(k,l);
  B(i,j,k,l) += dt*R1*sum_aD(i,j)*sum_Pa(k,l);

  // cast _dmdu as a ttl tensor and compute its value
  // -1 * (inverse(U) * B:Grad_du)
  TensorA<2> dMdu(dMdu_in);
  try {
    dMdu(i,j) = -(ttl::inverse(U)(i,j,k,l)*B(k,l,m,n)*Grad_du(m,n));
  }
  catch (const int inverseException){
    dMdu(i,j) = 0.0*ttl::identity(i,j);
    err++;
  }

  return err;
}

static int plasticity_compute_dMdu_np1(const Constitutive_model *m,
                                       const void *ctx,
                                       const double *Grad_op,
                                       const int nne,
                                       const int ndofn,
                                       double *dM_du)
{
  int err = 0;
  // The existing function compute_dMdu takes Grad_op for a given
  // alpha,beta pair and computes the corresponding dMdu. We will
  // generate the required inputs for this function and call it
  // alpha*beta times. This should be improved

  auto CTX = (plasticity_ctx *) ctx;
  // shorthand of deformation gradients
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;

  // compute M at n+1, again, this information is in CM
  Tensor<2> M = {}, eFn = {}, eFnp1 = {};
  {
    Tensor<2> pFnp1_I;
    TensorA<2> pFnp1(Fs[TENSOR_pFnp1].m_pdata), pFn(Fs[TENSOR_pFn].m_pdata);
    err += inv(pFnp1,pFnp1_I);

    if(CTX->is_coulpled_with_thermal)
    {
      TensorA<2> hFn(CTX->hFn), hFnp1(CTX->hFnp1);
      Tensor<2> hFn_I, hFnp1_I, N = {}, pFn_I;

      err += inv(pFn,   pFn_I);
      err += inv(hFn,   hFn_I);
      err += inv(hFnp1, hFnp1_I);

      N = hFn(i,k)*hFnp1_I(k,j);
      err += compute_M(M.data, Fs[TENSOR_Fn].m_pdata, N.data, pFnp1_I.data, CTX);

      err += compute_eF(eFn.data,  Fs[TENSOR_Fn].m_pdata,  hFn_I.data,   pFn_I.data,   CTX);
      err += compute_eF(eFnp1.data,Fs[TENSOR_Fnp1].m_pdata,hFnp1_I.data, pFnp1_I.data, CTX);
    }
    else
    {
      M = pFn(i,k)*pFnp1_I(k,j);
      err += m->param->get_eF(m,eFn.data,1);
      err += m->param->get_eF(m,eFnp1.data,2);
    }
  }

  ELASTICITY *elast = (m->param)->cm_elast;
  elast->update_elasticity(elast,eFnp1.data, 1);

  // compute slip system
  SLIP_SYSTEM *slip_in = (((m->param)->cm_mat)->mat_p)->slip;
  SLIP_SYSTEM slip;
  construct_slip_system(&slip,slip_in->unit_cell);

  rotate_crystal_orientation(&slip, m->vars_list[0][m->model_id].Fs[TENSOR_R].m_pdata, slip_in);

  const double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  const double g_n         = state_var[VAR_g_n];
  const double g_np1       = state_var[VAR_g_np1];

  const double *tau        = (m->vars_list[0][m->model_id]).Fs[TENSOR_tau].m_pdata;
  const double *gamma_dots = (m->vars_list[0][m->model_id]).Fs[TENSOR_gamma_dot].m_pdata;

  double Grad_op_ab[DIM_3x3];
  // make successive calls to compute_dMdu for each node/dof
  for (int a = 0; a < nne; a++) {
    for(int b = 0; b < ndofn; b++) {
      int idx_ab = idx_4_gen(a,b,0,0,nne,ndofn,DIM_3,DIM_3);
      // call to compute_dMdu
      memcpy(Grad_op_ab,Grad_op+idx_ab,DIM_3x3*sizeof(double));
      err += plasticity_compute_dMdu(m, dM_du+idx_ab,Grad_op_ab, eFn.data, eFnp1.data, M.data,
                                     elast->S, elast->L, g_n, g_np1, tau, gamma_dots, slip.p_sys, CTX->dt);
    }
  }
  destruct_slip_system(&slip);
  return err;
}

static int plasticity_compute_dMdu_npa(const Constitutive_model *m,
                                       const void *ctx,
                                       const double *Grad_op,
                                       const int nne,
                                       const int ndofn,
                                       double *dM_du,
                                       double alpha)
{
  // Total Lagrangian based
  int err = 0;
  // The existing function compute_dMdu takes Grad_op for a given
  // alpha,beta pair and computes the corresponding dMdu. We will
  // generate the required inputs for this function and call it
  // alpha*beta times. This should be improved
  auto CTX = (plasticity_ctx *) ctx;

  Tensor<2> Fnpa,Fn,eFn = {},eFnpa = {},pFnpa,pFnpa_I,pFnp1,pFn,Mnpa = {},hFnpaI;
  TensorA<2> Fnp1(CTX->F);

  err += m->param->get_F(m,Fn.data,1);
  eFn = ttl::identity(i,j); // Total Lagrangian
  // if updated: err += m->param->get_eF(m,eFn.data,1);

  err += m->param->get_pF(m,pFnp1.data,2);
  err += m->param->get_pF(m,pFn.data,1);

  mid_point_rule(Fnpa.data, Fn.data, Fnp1.data, alpha, DIM_3x3);
  mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, DIM_3x3);

  err += inv(pFnpa, pFnpa_I); // Total Lagrangian

  hFnpaI = ttl::identity(i,j);
  if(CTX->is_coulpled_with_thermal)
  {
    // double hFnpa[9];
    // mid_point_rule(hFnpa, CTX->hFn, CTX->hFnp1, alpha, DIM_3x3);
    TensorA<2> hFnpa(CTX->hFnp1);
    err += inv(hFnpa, hFnpaI);
  }

  Mnpa = hFnpaI(i,j)*pFnpa_I(j,k);
  eFnpa = Fnpa(i,j)*Mnpa(j,k);

  ELASTICITY *elast = (m->param)->cm_elast;
  elast->update_elasticity(elast,eFnpa.data, 1);

  // compute slip system
  SLIP_SYSTEM *slip_in = (((m->param)->cm_mat)->mat_p)->slip;
  SLIP_SYSTEM slip;
  construct_slip_system(&slip,slip_in->unit_cell);

  rotate_crystal_orientation(&slip, m->vars_list[0][m->model_id].Fs[TENSOR_R].m_pdata, slip_in);

  // update plasticty variables
  int N_SYS = slip.N_SYS;

  const double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;
  double g_n         = state_var[VAR_g_n];
  double g_np1       = state_var[VAR_g_np1];
  double *tau_np1        = (m->vars_list[0][m->model_id]).Fs[TENSOR_tau].m_pdata;
  double *tau_n          = (m->vars_list[0][m->model_id]).Fs[TENSOR_tau_n].m_pdata;
  double *gamma_dots_np1 = (m->vars_list[0][m->model_id]).Fs[TENSOR_gamma_dot].m_pdata;
  double *gamma_dots_n   = (m->vars_list[0][m->model_id]).Fs[TENSOR_gamma_dot_n].m_pdata;

  double *tau        = (double *) malloc(sizeof(double)*N_SYS);
  double *gamma_dots = (double *) malloc(sizeof(double)*N_SYS);
  double g_npa = 0.0;

  mid_point_rule(tau, tau_n, tau_np1, alpha, N_SYS);
  mid_point_rule(gamma_dots, gamma_dots_n, gamma_dots_np1, alpha, N_SYS);
  mid_point_rule(&g_npa, &g_n, &g_np1, alpha, 1);

  // make successive calls to compute_dMdu for each node/dof
  double Grad_op_ab[DIM_3x3];
  for (int a = 0; a < nne; a++) {
    for(int b = 0; b < ndofn; b++) {
      int idx_ab = idx_4_gen(a,b,0,0,nne,ndofn,DIM_3,DIM_3);
      memcpy(Grad_op_ab,Grad_op+idx_ab,DIM_3x3*sizeof(double));
      err += plasticity_compute_dMdu(m, dM_du + idx_ab, Grad_op_ab, eFn.data, eFnpa.data, Mnpa.data,
                                     elast->S, elast->L, g_n, g_npa, tau, gamma_dots, slip.p_sys, CTX->dt);
    }
  }

  destruct_slip_system(&slip);
  free(tau);
  free(gamma_dots);
  return err;
}

int CP_PARAM::compute_dMdu(const Constitutive_model *m,
                           const void *ctx,
                           double *Grad_op,
                           const int nne,
                           const int ndofn,
                           double *dM_du)
  const
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  if(CTX->alpha<0)
    err += plasticity_compute_dMdu_np1(m,ctx,Grad_op,nne,ndofn,dM_du);
  else
    err += plasticity_compute_dMdu_npa(m,ctx,Grad_op,nne,ndofn,dM_du,CTX->alpha);
  return err;
}

int CP_PARAM::write_restart(FILE *fp, const Constitutive_model *m)
  const
{

  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);

  const int N_SYS = ((((m->param)->cm_mat)->mat_p)->slip)->N_SYS;
  for(int a=0; a<N_SYS; a++)
    fprintf(fp, "%.17e ", Fs[TENSOR_tau_n].m_pdata[a]);

  fprintf(fp, "\n");

  for(int a=0; a<N_SYS; a++)
    fprintf(fp, "%.17e ", Fs[TENSOR_gamma_dot_n].m_pdata[a]);

  fprintf(fp, "\n");

  fprintf(fp, "%.17e %.17e %.17e %.17e\n", state_var[VAR_g_n],
          state_var[VAR_g_nm1], state_var[VAR_L_n], state_var[VAR_L_nm1]);

  return err;
}

int CP_PARAM::read_restart(FILE *fp, Constitutive_model *m)
  const
{
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  cm_read_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);

  const int N_SYS = ((((m->param)->cm_mat)->mat_p)->slip)->N_SYS;
  for(int a=0; a<N_SYS; a++)
    CHECK_SCANF(fp, "%lf ", Fs[TENSOR_tau_n].m_pdata + a);

  for(int a=0; a<N_SYS; a++)
    CHECK_SCANF(fp, "%lf ", Fs[TENSOR_gamma_dot_n].m_pdata + a);

  CHECK_SCANF(fp, "%lf %lf %lf %lf\n", state_var+VAR_g_n, state_var+VAR_g_nm1, state_var+VAR_L_n, state_var+VAR_L_nm1);

  // set values at n+1

  Fs[TENSOR_Fnp1]  = Fs[TENSOR_Fn];
  Fs[TENSOR_pFnp1] = Fs[TENSOR_pFn];
  Fs[TENSOR_tau]   = Fs[TENSOR_tau_n];
  Fs[TENSOR_gamma_dot] = Fs[TENSOR_gamma_dot_n];
  state_var[VAR_g_np1] = state_var[VAR_g_n];
  state_var[VAR_L_np1] = state_var[VAR_L_n];
  return 0;
}

// THIS IS A FUNCTION STUB.
int CP_PARAM::read_param(FILE *in)
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

  // READ PROPERTIES IN ALPHABETICAL ORDER
  int param_in = PARAM_NO-3;
  int match = fscanf(in, "%lf %lf %lf %lf %lf %lf %lf",
                     param + PARAM_gamma_dot_0,
                     param + PARAM_m,
                     param + PARAM_G0,
                     param + PARAM_g0,
                     param + PARAM_gs_0,
                     param + PARAM_gamma_dot_s,
                     param + PARAM_w);

  err += scan_for_valid_line(in);

  SLIP_SYSTEM *slip = (SLIP_SYSTEM *) malloc(sizeof(SLIP_SYSTEM));

  int unit_cell = -1;
  match += fscanf(in, "%d", &unit_cell);
  param_in ++;

  construct_slip_system(slip,unit_cell);
  match += fscanf(in, "%d", (slip->ort_option)+0);
  param_in++;

  if(slip->ort_option[0] == 0)
  {
    match += fscanf(in, "%d", (slip->ort_option)+1);
    param_in++;
  }

  if(slip->ort_option[0] == 2)
  {
    match += fscanf(in, "%s", slip->ort_file_in);
    param_in++;
  }

  if(slip->ort_option[0] == 3)
  {
    match += fscanf(in, "%lf %lf %lf",
                    slip->ort_angles,
                    slip->ort_angles + 1,
                    slip->ort_angles + 2);
    param_in += 3;
  }

  if (match != param_in) err++;
  assert(match == param_in && "Did not read expected number of parameters");

  // scan past any other comment/blank lines in the block
  err += scan_for_valid_line(in);

  int set_cp_solver = 0;
  int param_read_no = fscanf(in, "%d", &set_cp_solver);

  int read_solver_info = 0;
  if(param_read_no==1)
  {
    if(set_cp_solver)
      read_solver_info = 1;
  }

  if(read_solver_info)
  {
    match = fscanf(in, "%d %d %d %d %lf %lf %lf",
                   param_idx + PARAM_max_itr_stag,
                   param_idx + PARAM_max_itr_hardening,
                   param_idx + PARAM_max_itr_M,
                   param_idx + PARAM_max_subdivision,
                   param     + PARAM_tol_hardening,
                   param     + PARAM_tol_M,
                   param     + PARAM_computer_zero);
    if (match != 7) err++;
    assert(match == 7 && "Did not read expected number of parameters");
  }
  else // use default
  {
    param_idx[PARAM_max_itr_stag]      = 15;
    param_idx[PARAM_max_itr_hardening] = 1;
    param_idx[PARAM_max_itr_M]         = 50;
    param_idx[PARAM_max_subdivision]   = -1;
    param[PARAM_tol_hardening]     = 1.0e-6;
    param[PARAM_tol_M]             = 1.0e-6;
    param[PARAM_computer_zero]     = 1.0e-15;
  }

  err += scan_for_valid_line(in);

  // not expecting EOF, check and return error if encountered
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");

  MATERIAL_CRYSTAL_PLASTICITY  *mat_p = (MATERIAL_CRYSTAL_PLASTICITY  *) malloc(sizeof(MATERIAL_CRYSTAL_PLASTICITY));
  set_properties_crystal_plasticity(mat_p,slip,param[PARAM_gamma_dot_0],param[PARAM_gamma_dot_s],
                                    param[PARAM_m],          param[PARAM_g0],
                                    param[PARAM_G0],         param[PARAM_gs_0],
                                    param[PARAM_w]);

  (this->cm_mat)->mat_p = mat_p;

  return err;
}

int CP_PARAM::update_elasticity(const Constitutive_model *m,
                                const void *ctx_in,
                                double *L,
                                double *S,
                                const int compute_stiffness)
  const
{
  int err = 0;
  auto ctx = (plasticity_ctx *) ctx_in;

  // if transient cases,
  // get_eF is not working because eF needs to be updated using mid-point alpha
  // below checks whether to use get_eF or give eFnpa in ctx

  if(ctx->eFnpa)
    err += constitutive_model_default_update_elasticity(m, ctx->eFnpa, L, S, compute_stiffness);
  else
  {
    // shorthand of deformation gradients
    Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
    Tensor<2> eF = {}, pFnp1_I;
    TensorA<2> pFnp1(Fs[TENSOR_pFnp1].m_pdata), Fnp1(Fs[TENSOR_Fnp1].m_pdata);

    err += inv(pFnp1, pFnp1_I);
    if(ctx->is_coulpled_with_thermal)
    {
      TensorA<2> hFnp1(ctx->hFnp1);
      Tensor<2> hFnp1_I;

      err += inv(hFnp1, hFnp1_I);
      err += compute_eF(eF.data, Fnp1.data,hFnp1_I.data, pFnp1_I.data, ctx);
    }
    else
      eF = Fnp1(i,k)*pFnp1_I(k,j);

    err += constitutive_model_default_update_elasticity(m, eF.data, L, S, compute_stiffness);
  }

  return err;
}

int CP_PARAM::get_subdiv_param(const Constitutive_model *m,
                               double *subdiv_param,
                               double dt)
  const
{
  int err = 0;
  SLIP_SYSTEM *slip = (((m->param)->cm_mat)->mat_p)->slip;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;

  double gamma_np1 = 0.0;
  for(int a=0; a<slip->N_SYS; a++)
    gamma_np1 += fabs(Fs[TENSOR_gamma_dot].m_pdata[a]);

  *subdiv_param = dt*gamma_np1/MAX_D_GAMMA;
  return err;
}

int CP_PARAM::model_dependent_initialization(void)
{
  this->type              = CRYSTAL_PLASTICITY;
  this->n_param           = PARAM_NO;
  this->model_param       = new double[PARAM_NO]();
  this->n_param_index     = PARAM_INX_NO;
  this->model_param_index = new int [PARAM_INX_NO]();
  return 0;
}

int CP_PARAM::model_dependent_finalization(void)
{
  int err = 0;
  destruct_slip_system(((this->cm_mat)->mat_p)->slip);
  free(((this->cm_mat)->mat_p)->slip);
  free((this->cm_mat)->mat_p);
  return err;
}

int plasticity_model_ctx_build(void **ctx,
                               double *F,
                               const double dt,
                               const double alpha,
                               double *eFnpa,
                               double *hFn,
                               double *hFnp1,
                               const int is_coulpled_with_thermal)
{
  int err = 0;
  plasticity_ctx *t_ctx = (plasticity_ctx *) malloc(sizeof(plasticity_ctx));

  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;
  t_ctx->hFn   = NULL;
  t_ctx->hFnp1 = NULL;

  t_ctx->F = F;
  t_ctx->eFnpa = eFnpa;


  t_ctx->dt = dt;
  t_ctx->alpha = alpha;

  t_ctx->is_coulpled_with_thermal = is_coulpled_with_thermal;
  t_ctx->hFn  = hFn;
  t_ctx->hFnp1= hFnp1;

  // assign handle
  *ctx = t_ctx;
  return err;
}

int CP_PARAM::destroy_ctx(void **ctx)
  const
{
  int err = 0;
  plasticity_ctx *t_ctx = (plasticity_ctx *) *ctx;
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

int CP_PARAM::integration_algorithm(Constitutive_model *m,
                                    const void *ctx)
  const
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  memcpy(m->vars_list[0][m->model_id].Fs[TENSOR_Fnp1].m_pdata, CTX->F, DIM_3x3 * sizeof(*(CTX->F)));

  const double dt = CTX->dt;
  double *param     = (m->param)->model_param;
  int    *param_idx = (m->param)->model_param_index;
    
  CRYSTAL_PLASTICITY_SOLVER_INFO solver_info;
  set_crystal_plasticity_solver_info(&solver_info,param_idx[PARAM_max_itr_stag],
                                     param_idx[PARAM_max_itr_hardening],
                                     param_idx[PARAM_max_itr_M],
                                     param[PARAM_tol_hardening],
                                     param[PARAM_tol_M],
                                     param[PARAM_computer_zero]);
  solver_info.max_subdivision = param_idx[PARAM_max_subdivision];

  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  double g_n   = state_var[VAR_g_n];
  double g_np1 = state_var[VAR_g_np1];
  double L_np1 = state_var[VAR_L_np1];

  Tensor<2> eFnp1,C,pFnp1_I;
  TensorA<2> pFnp1(Fs[TENSOR_pFnp1].m_pdata), pFn(Fs[TENSOR_pFn].m_pdata),
  Fn( Fs[TENSOR_Fn].m_pdata), Fnp1(Fs[TENSOR_Fnp1].m_pdata);

  MATERIAL_CONSTITUTIVE_MODEL *cm_mat = (m->param)->cm_mat;
  ELASTICITY *elasticity = (m->param)->cm_elast;

  // compute slip system and rotate orientation
  SLIP_SYSTEM *slip_in = (cm_mat->mat_p)->slip;
  SLIP_SYSTEM slip;
  construct_slip_system(&slip,slip_in->unit_cell);
  rotate_crystal_orientation(&slip, Fs[TENSOR_R].m_pdata, slip_in);
  (cm_mat->mat_p)->slip = &slip;
  
  
//  GcmCpIntegrator gcm_cp;
/*  gcm_cp.set_tensors(Fs[TENSOR_Fnp1].m_pdata,
                     Fs[TENSOR_Fn].m_pdata,
                     Fs[TENSOR_Fnm1].m_pdata,
                     Fs[TENSOR_pFnp1].m_pdata,
                     Fs[TENSOR_pFn].m_pdata,
                     CTX->hFnp1,
                     CTX->hFn);
                     
  */

  // perform integration algorithm for the crystal plasticity
  if(CTX->is_coulpled_with_thermal)
  {
    err += staggered_Newton_Rapson_generalized(pFnp1.data,
                                               &g_np1, &L_np1,
                                               pFn.data, Fn.data, Fnp1.data,
                                               CTX->hFn, CTX->hFnp1,
                                               g_n, dt, cm_mat, elasticity, &solver_info);
  }
  else
  {
    err += staggered_Newton_Rapson(pFnp1.data,
                                   &g_np1, &L_np1,
                                   pFn.data, Fn.data, Fnp1.data,
                                   g_n, dt, cm_mat, elasticity, &solver_info);
  }

  // update compute values from integration algorithm
  state_var[VAR_g_np1] =  g_np1;
  state_var[VAR_L_np1] =  L_np1;

  err += inv(pFnp1,pFnp1_I);

  if(CTX->is_coulpled_with_thermal)
  {
    TensorA<2> hFnp1(CTX->hFnp1);
    Tensor<2> hFnp1_I;
    inv(hFnp1,hFnp1_I);
    err += compute_eF(eFnp1.data, Fnp1.data, hFnp1_I.data, pFnp1_I.data, CTX);
  }
  else
    eFnp1 = Fnp1(i,k)*pFnp1_I(k,j);

  C = eFnp1(k,i)*eFnp1(k,j);

  elasticity->update_elasticity(elasticity,eFnp1.data, 0);
  err += compute_tau_alphas(Fs[TENSOR_tau].m_pdata,C.data, elasticity->S, &slip);
  err += compute_gamma_dots(Fs[TENSOR_gamma_dot].m_pdata, Fs[TENSOR_tau].m_pdata, g_np1, cm_mat->mat_p);

  (cm_mat->mat_p)->slip = slip_in;
  destruct_slip_system(&slip);

  return err;
}

int plasticity_model_construct_rotation(EPS *eps, Matrix<int> &e_ids, Matrix<double> &angles)
{
  int err = 0;
  double Ax[DIM_3x3], Az[DIM_3x3], Ay[DIM_3x3];

  for(int a=0; a<e_ids.m_row; a++)
  {
    int id = e_ids(a, 0);
    if(id<0)
      continue;

    int ip = e_ids(a, 1);
    Constitutive_model *m = &(eps[id].model[ip]);
    
    MATERIAL_CONSTITUTIVE_MODEL *cm_mat = (m->param)->cm_mat;
    SLIP_SYSTEM *slip = (cm_mat->mat_p)->slip;
    
    
    if((m->param)->type!=CRYSTAL_PLASTICITY)
      continue;
      
    Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
    double phi   = angles.m_pdata[a*DIM_3+0];     // NOTE: phi = Mat_v(*angles, a+1, 1) is not working, do not know why.
    double theta = angles.m_pdata[a*DIM_3+1];
    double psi   = angles.m_pdata[a*DIM_3+2];

    // compute rotation matrix of Euler angles
    // Fs[TENSOR_R] = Az(psi)*Ay(theta)*Ax(phi)
    err += rotation_matrix_of_Euler_angles(Fs[TENSOR_R].m_pdata,
                                           Ax, Ay, Az, phi, theta, psi, slip->ort_option[2]);
  }
  return err;
}

/// return Euler angle type
int plasticity_model_read_orientations(Matrix<int> &e_ids, Matrix<double> &angles, IP_ID_LIST *elm_ip_map, char *fn_in, int myrank, int ne)
{
  int EulerAngleType = 0;
  char fn[1024], line[1024];
  const char EulerAngleTypeName[] = " EulerAngleType:";
  int cno0 = strlen(EulerAngleTypeName);
  
  sprintf(fn, "%s_%d.in", fn_in, myrank);
  FILE *fp = fopen(fn, "r");
  if(fp==NULL)
  {
    PGFEM_printf("fail to read [%s]\n", fn);
    PGFEM_printf("set default onrientation [R=I]\n");
    return EulerAngleType;
  }

  
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
    {
      bool SetEulerAngleType = true;
      int cno = strlen(line);
      if(cno>cno0+2)
      {
        for(int ia=0; ia<cno0; ia++)
        {
          if(line[ia+1]!=EulerAngleTypeName[ia])
            SetEulerAngleType = false;
        }
      }
      if(SetEulerAngleType)
      {  
        sscanf(line+cno0+2, "%d", &EulerAngleType);
        PGFEM_printf("EulerAngleType: %d\n", EulerAngleType);
      }
      continue;
    }
    
    int e, ip;
    double x1, x2, x3;
    sscanf(line, "%d %d %lf %lf %lf", &e, &ip, &x1, &x2, &x3);

    if(e<0)
    {
      for(int a=0; a<ne; a++)
      {
        int n_ip = elm_ip_map[a].n_ip;
        int mat = elm_ip_map[a].mat_id;
        if(mat!=ip)
          continue;

        for(int b=0; b<n_ip; b++)
        {
          int ip_id = elm_ip_map[a].ip_ids.m_pdata[b];
          e_ids(ip_id, 0) = a;  // +1 is needed because ip_id starts from 0
          e_ids(ip_id, 1) = b;
          angles(ip_id, 0) = x1;
          angles(ip_id, 1) = x2;
          angles(ip_id, 2) = x3;
        }
      }
    }
    else
    {
      int ip_id = elm_ip_map[e].ip_ids.m_pdata[ip];
      e_ids(ip_id, 0) = e;  // +1 is needed because ip_id starts from 0
      e_ids(ip_id, 1) = ip;
      angles(ip_id, 0) = x1;
      angles(ip_id, 1) = x2;
      angles(ip_id, 2) = x3;
    }
  }

  fclose(fp);
  return EulerAngleType;
}

int plasticity_model_generate_random_orientation_element(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix<int> &e_ids, Matrix<double> &angles, int diff_ort_at_ip)
{
  int err = 0;

  double *angle_temp = (double *) malloc(sizeof(double)*ne*DIM_3);
  err += generate_random_crystal_orientation(angle_temp, ne);

  for(int a=0; a<ne; a++)
  {
    int n_ip = elm_ip_map[a].n_ip;
    int mat = elm_ip_map[a].mat_id;
    if(mat!=mat_id)
      continue;

    double   phi = angle_temp[a*DIM_3 + 0];
    double theta = angle_temp[a*DIM_3 + 1];
    double   psi = angle_temp[a*DIM_3 + 2];

    for(int ip=0; ip<n_ip; ip++)
    {
      int ip_id = elm_ip_map[a].ip_ids.m_pdata[ip];
      if(diff_ort_at_ip && ip>0)
      {
        double temp[3];
        err += generate_random_crystal_orientation(temp, 1);
        phi = temp[0];
        theta = temp[1];
        psi = temp[2];
      }
      angles(ip_id,0) = phi;
      angles(ip_id,1) = theta;
      angles(ip_id,2) = psi;
      e_ids(ip_id,0) = a;
      e_ids(ip_id,1) = ip;
    }
  }
  free(angle_temp);
  return err;
}

int plasticity_model_generate_random_orientation_crystal(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix<int> &e_ids, Matrix<double> &angles)
{
  int err = 0;

  double temp[3];
  err += generate_random_crystal_orientation(temp, 1);

  for(int a=0; a<ne; a++)
  {
    int n_ip = elm_ip_map[a].n_ip;
    int mat = elm_ip_map[a].mat_id;
    if(mat!=mat_id)
      continue;

    for(int ip=0; ip<n_ip; ip++)
    {
      int ip_id = elm_ip_map[a].ip_ids.m_pdata[ip];
      angles(ip_id,0) = temp[0];
      angles(ip_id,1) = temp[1];
      angles(ip_id,2) = temp[2];
      e_ids(ip_id,0) = a;
      e_ids(ip_id,1) = ip;
    }
  }
  return err;
}

int plasticity_model_set_given_orientation_crystal(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix<int> &e_ids, Matrix<double> &angles, double *angle_in)
{
  int err = 0;

  for(int a=0; a<ne; a++)
  {
    int n_ip = elm_ip_map[a].n_ip;
    int mat = elm_ip_map[a].mat_id;
    if(mat!=mat_id)
      continue;

    for(int ip=0; ip<n_ip; ip++)
    {
      int ip_id = elm_ip_map[a].ip_ids.m_pdata[ip];
      angles(ip_id,0) = angle_in[0];
      angles(ip_id,1) = angle_in[1];
      angles(ip_id,2) = angle_in[2];
      e_ids(ip_id,0) = a;
      e_ids(ip_id,1) = ip;
    }
  }
  return err;
}


int plasticity_model_set_zero_angles(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix<int> &e_ids, Matrix<double> &angles)
{
  int err = 0;
  for(int a=0; a<ne; a++)
  {
    int n_ip = elm_ip_map[a].n_ip;
    int mat = elm_ip_map[a].mat_id;
    if(mat!=mat_id)
      continue;

    for(int ip=0; ip<n_ip; ip++)
    {
      int ip_id = elm_ip_map[a].ip_ids.m_pdata[ip];

      angles(ip_id,0) = 0.0;
      angles(ip_id,1) = 0.0;
      angles(ip_id,2) = 0.0;
      e_ids(ip_id,0) = a;
      e_ids(ip_id,1) = ip;
    }
  }
  return err;
}

int plasticity_model_set_orientations(EPS *eps,
                                      const int ne,
                                      const Element *elem,
                                      const int n_mat,
                                      const HOMMAT *hmat_list,
				      int myrank)
{
  int err = 0;
  int crystal_plasticity_included = 0;
  for(int i = 0; i < n_mat; i++)
  {
    if(hmat_list[i].param->type==CRYSTAL_PLASTICITY)
      crystal_plasticity_included++;
  }

  if(crystal_plasticity_included==0)
    return err;

  // build element ip ids that will be used to assign element orientation
  IP_ID_LIST *elm_ip_map = new IP_ID_LIST[ne];
  int cnt_of_ips = plasticity_model_construct_elem_ip_map(elm_ip_map, eps, elem, ne);

  Matrix<int> e_ids(cnt_of_ips, 2, -1);
  Matrix<double> angles(cnt_of_ips,DIM_3,0.0);

  for(int i = 0; i < n_mat; i++)
  {
    if(hmat_list[i].param->type==CRYSTAL_PLASTICITY)
    {
      SLIP_SYSTEM *slip = ((hmat_list[i].param->cm_mat)->mat_p)->slip;
      if(slip->ort_option[0]==2)
      {
        slip->ort_option[2] = plasticity_model_read_orientations(e_ids, angles, elm_ip_map, slip->ort_file_in, myrank, ne);
        break;
      }
    }
  }

  int save_orientations = 0;
  for(int i = 0; i < n_mat; i++)
  {
    if(hmat_list[i].param->type==CRYSTAL_PLASTICITY)
    {
      SLIP_SYSTEM *slip = ((hmat_list[i].param->cm_mat)->mat_p)->slip;
      switch(slip->ort_option[0])
      {
       case -1:
        plasticity_model_set_zero_angles(ne, elm_ip_map, hmat_list[i].param->mat_id, e_ids, angles);
        break;
       case 0:
         {
           err += plasticity_model_generate_random_orientation_element(ne, elm_ip_map, hmat_list[i].param->mat_id, e_ids, angles, slip->ort_option[1]);
           save_orientations++;
           break;
         }
       case 1:
         {
           err += plasticity_model_generate_random_orientation_crystal(ne, elm_ip_map, hmat_list[i].param->mat_id, e_ids, angles);
           save_orientations++;
           break;
         }
       case 3:
         {
           err += plasticity_model_set_given_orientation_crystal(ne, elm_ip_map, hmat_list[i].param->mat_id, e_ids, angles, slip->ort_angles);
           break;
         }
       default:
        break;
      }
    }
  }
  plasticity_model_construct_rotation(eps, e_ids, angles);

  char default_ort_dir[1024];
  sprintf(default_ort_dir, "CRYSTAL_ORIENTATION");

  if(save_orientations)
  {
    // check default CRYSTAL_ORIENTATION directory exists
    if(make_path(default_ort_dir,DIR_MODE) != 0)
    {
      PGFEM_printf("Directory [%s] not created!\n",default_ort_dir);
      abort();
    }

    // write crystal orientations
    char fn_orientation[2048];
    sprintf(fn_orientation, "%s/orientation_%d.in", default_ort_dir, myrank);
    FILE *fp_ort = fopen(fn_orientation, "w");
    fprintf(fp_ort, "# Element (crystal) orientations are generated randomly\n");
    fprintf(fp_ort, "# element_ID, Integration_point_ID, phi [radian], theta [radian], psi [radian]\n");
    for(int a=0; a<e_ids.m_row; a++)
    {
      if(e_ids(a, 0)<0)
        continue;

      fprintf(fp_ort, "%d %d %e %e %e\n", e_ids(a,0), e_ids(a,1), angles(a,0), angles(a,1), angles(a,2));
    }
    fclose(fp_ort);
  }

  for(int a=0; a<ne; a++) {
    long n_ip = 0;
    int_point(elem[a].toe,&n_ip);
    for(int ip=0; ip<n_ip; ip++)
    {
      Constitutive_model *m = &(eps[a].model[ip]);

      if(m->param->type != CRYSTAL_PLASTICITY)
        continue;

      double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
      MATERIAL_CRYSTAL_PLASTICITY *mat_p = ((m->param)->cm_mat)->mat_p;

      // set state variables to initial values
      state_var[VAR_g_nm1] = mat_p->g0;
      state_var[VAR_g_n]   = mat_p->g0;
      state_var[VAR_g_np1] = mat_p->g0;
      state_var[VAR_L_nm1] = 0.0;
      state_var[VAR_L_n]   = 0.0;
      state_var[VAR_L_np1] = 0.0;

      // set the dimensions for the tau and gamma_dot varables
      Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
      const int N_SYS = (mat_p->slip)->N_SYS;
      Fs[TENSOR_tau].initialization(N_SYS, 1);
      Fs[TENSOR_tau_n].initialization(N_SYS, 1);
      Fs[TENSOR_gamma_dot].initialization(N_SYS, 1);
      Fs[TENSOR_gamma_dot_n].initialization(N_SYS, 1);

      // intitialize to zeros
      for(int ia=0; ia<N_SYS; ia++)
      {
        Fs[TENSOR_tau].m_pdata[ia] = 0.0;
        Fs[TENSOR_tau_n].m_pdata[ia] = 0.0;
        Fs[TENSOR_gamma_dot].m_pdata[ia] = 0.0;
        Fs[TENSOR_gamma_dot_n].m_pdata[ia] = 0.0;
      }
    }
  }

  delete[] elm_ip_map;
  return err;
}

int cm3f_plasticity_compute_dM(const Constitutive_model *con,
                               double *dMdu_in,
                               double *dMdt_in,
                               double *Grad_du_in,
                               double d_theta_r,
                               double *eFn_in,
                               double *eFnp1_in,
                               double *M_in,
                               double *S_in,
                               double *L_in,
                               const double g_n,
                               const double g_np1,
                               const double *tau,
                               const double *gamma_dots,
                               double *Psys,
                               double *tFr_in,
                               double theta_r_in,
                               const double dt)
{
  int err = 0;
  // compute dMdu:U = -grad(du):B
  // Grad_du = Grad(du)

  MATERIAL_CRYSTAL_PLASTICITY *mat_p = ((con->param)->cm_mat)->mat_p;

  const int N_SYS          = (mat_p->slip)->N_SYS;
  const double gamma_dot_0 = mat_p->gamma_dot_0;
  const double gamma_dot_s = mat_p->gamma_dot_s;
  const double mm          = mat_p->m;
  const double g0          = mat_p->g0;
  const double G0          = mat_p->G0;
  const double gs_0        = mat_p->gs_0;
  const double w           = mat_p->w;

  double Jtheta  = 1.0;
  double theta_r = 1.0;
  Tensor<2> tFrI = {};

  if(con->param->cm3f)
  {
    TensorA<2> tFr(tFr_in);

    double tJr     = ttl::det(tFr);
    theta_r = theta_r_in;
    Jtheta = pow(theta_r/tJr, 1.0/3.0);
    err += inv(tFr, tFrI);
  }

  double gamma_dot = 0.0;
  for(int a = 0; a<N_SYS; a++)
    gamma_dot += fabs(gamma_dots[a]);

  double gm_gms   = gamma_dot/gamma_dot_s;
  double sign_gm_gms = (gm_gms < 0) ? -1.0 : 1.0;

  double R3 = 0.0;
  double gs_np1 = 0.0;
  if(fabs(gm_gms)>1.0e-15)
  {
    R3 = gs_0*w/gamma_dot_s*sign_gm_gms*pow(fabs(gm_gms), w-1.0);
    gs_np1 = gs_0*pow(fabs(gm_gms),w);
  }

  double AA = R3*gamma_dot*(g_n - g0 + dt*G0*gamma_dot) +
              gs_np1 * (gs_np1 - g0 - g_n) + g0*g_n;
  double BB = gs_np1 - g0  - dt*G0*gamma_dot;
  double R4 = dt*G0*AA/BB/BB;

  double sum_1gm1gm = 0.0;

  // convert C struct matrices into ttl tensors
  const TensorA<2> Grad_du(Grad_du_in);
  const TensorA<2> eFn(eFn_in);
  const TensorA<2> eFnp1(eFnp1_in);
  const TensorA<2> M(M_in);
  const TensorA<2> S(S_in);
  const TensorA<4> L(L_in);

  const Tensor<2> C = eFnp1(k,i).to(i,k) * eFnp1(k,j);             // eFnp1' * eFnp1
  Tensor<2> MI = {};
  err += inv(M, MI);
  const Tensor<2> MI_x_C = MI(j,i) * C(j,k);                       // M^{-T} * C
  const Tensor<2> M_x_eFn = M(j,i).to(i,j) * eFn(k,j).to(j,k);     // M' * eFn'

  // tensors are initialized to 0
  Tensor<2> sum_aC = {};
  Tensor<2> sum_Pa = {};
  Tensor<2> sum_aD = {};
  Tensor<2> sum_CP = {};
  Tensor<4> U = {};
  Tensor<4> B = {};
  double sum_CY = 0.0;

  // set U to the 9x9 identity scaled by 1.0/dt
  U(i,j,k,l) =  ttl::identity(i,j,k,l);

  for(int a = 0; a<N_SYS; a++)
  {
    double drdtau = gamma_dot_0/mm/g_np1*pow(fabs(tau[a]/g_np1), 1.0/mm - 1.0);
    double drdg   = -drdtau*tau[a]/g_np1;

    double R2_a = ((gamma_dots[a] < 0) ? -1.0 : 1.0)*drdtau;
    sum_1gm1gm += ((gamma_dots[a] < 0) ? -1.0 : 1.0)*drdg;

    // compute P alpha of Psys
    const Tensor<2, const double*> Pa(Psys + DIM_3x3 * a);

    // compute C alpha and D alpha using ttl operations
    auto t0 = Pa(i,k) * S(k,j);                                    // Pa * S
    auto t1 = S(i,k) * Pa(j,k).to(k,j);                            // S * Pa'
    auto t2 = L(i,k,n,l) * C(n,l) * Pa(k,j);                       // L:C * Pa
    const Tensor<2> AA = t0 + t1 + t2;
    const Tensor<2> aC = MI_x_C(i,j) * AA(j,k);
    const Tensor<2> aD = Jtheta*eFnp1(i,j)*AA(j,k)*M_x_eFn(k,l) - 1.0/3.0*C(m,n)*AA(m,n)*tFrI(l,i).to(i,l);

    sum_Pa(i,j) += drdg * Pa(i,j);
    sum_aC(i,j) += R2_a * aC(i,j);
    sum_aD(i,j) += R2_a * aD(i,j);
    if(con->param->cm3f)
    {
      double CY = C(i,j)*AA(i,j);
      sum_CP(i,j) += CY*Pa(i,j);
      sum_CY      += R2_a * C(i,j)*AA(i,j);
    }

    // perform the Kronecker product using ttl and scales it by drdtau
    U(i,j,k,l) += dt*drdtau * aC(i,j) * Pa(k,l);
    B(i,j,k,l) += dt*drdtau * aD(i,j) * Pa(k,l);
  }

  double R1 = R4/(1.0-R4*sum_1gm1gm);

  // perform the Kronecker product using ttl and scales it by R1
  U(i,j,k,l) += dt*R1*sum_aC(i,j)*sum_Pa(k,l);
  B(i,j,k,l) += dt*R1*sum_aD(i,j)*sum_Pa(k,l);

  Tensor<4> UI = {};

  try {
    UI = ttl::inverse(U);
  }
  catch (const int inverseException){
    err++;
  }

  // cast _dmdu as a ttl tensor and compute its value
  // -1 * (inverse(U) * B:Grad_du)
  TensorA<2> dMdu(dMdu_in);

  if(err==0)
    dMdu(i,j) = -(UI(i,j,k,l)*B(k,l,m,n)*Grad_du(m,n));
  else
    dMdu(i,j) = 0.0*ttl::identity(i,j);

  if(con->param->cm3f)
  {
    TensorA<2> dMdt(dMdt_in);
    if(err==0)
    {
      Tensor<2> Z = 1.0/3.0/theta_r*(sum_CP(i,j) + R1*sum_CY*sum_Pa(i,j));
      dMdt(i,j) = -(UI(i,j,k,l)*Z(k,l)*d_theta_r);
    }
    else
      dMdt(i,j) = 0.0*ttl::identity(i,j);
  }

  return err;
}

int CP_PARAM::set_init_vals(Constitutive_model *m)
const
{
  // inital values are set in the more convoluted
  // read_constitutive_model_parameters->plasticity_model_read_parameters
  //   calling sequence

  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;

  if((m->param)->pF != NULL)
  {
    double pF[DIM_3x3];
    double pJ = det3x3((m->param)->pF);
    pJ = pow(pJ, 1.0/3.0);

    for(int ia=0; ia<DIM_3x3; ia++)
      pF[ia] = (m->param)->pF[ia]/pJ;

    memcpy(Fs[TENSOR_Fnm1 ].m_pdata, pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_Fn   ].m_pdata, pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_Fnp1 ].m_pdata, pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFnm1].m_pdata, pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFn  ].m_pdata, pF, sizeof(double)*DIM_3x3);
    memcpy(Fs[TENSOR_pFnp1].m_pdata, pF, sizeof(double)*DIM_3x3);
  }  
  return 0;
}
