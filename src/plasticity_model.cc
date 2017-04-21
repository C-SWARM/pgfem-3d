/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 *  Sangmin Lee, University of Notre Dame, Notre Dame, IN, <slee43@nd.edu>
 *  Aaron Howell, University of Notre Dame, Notre Dame, IN, <ahowell3@nd.edu>
 */

#include "plasticity_model.h"
#include "constitutive_model.h"
#include "cm_placeholder_functions.h"
#include "new_potentials.h"
#include "data_structure_c.h"
#include "elem3d.h"
#include "gen_path.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"
#include <ttl/ttl.h>

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX_D_GAMMA 0.005


static const int FLAG_end = 0;

// ttl declarations
namespace {
/// Only dealing with 3 dimensional double data, but it is sometimes const.
template <int R, class S = double>
using Tensor = ttl::Tensor<R, 3, S>;

template <int R, class S = double>
using Tensor9x9 = ttl::Tensor<R, 9, S>;

template <int R>
using Delta = ttl::Tensor<R, 3, double>;

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

Define_Matrix(double);
Define_Matrix(int);

typedef struct {
  int e_id;
  int n_ip;
  int mat_id;
  Matrix(int) ip_ids;
} IP_ID_LIST;

static size_t cp_get_size(const Constitutive_model *m)
{
  return state_variables_get_packed_size(m->vars_list[0]+m->model_id);
}

static int cp_pack(const Constitutive_model *m,
                    char *buffer,
                    size_t *pos)
{
  return state_variables_pack(m->vars_list[0]+m->model_id, buffer, pos);
}

static int cp_unpack(Constitutive_model *m,
                     const char *buffer,
                     size_t *pos)
{
  return state_variables_unpack(m->vars_list[0]+m->model_id, buffer, pos);
}

int compute_dMdu(const Constitutive_model *m,
                 Matrix(double) *dMdu,
                 const Matrix(double) *Grad_du,
                 const Matrix(double) *eFn,
                 const Matrix(double) *eFnp1,
                 const Matrix(double) *M,
                 const Matrix(double) *S,
                 const Matrix(double) *L,
                 const double g_n,
                 const double g_np1,
                 const double *tau,
                 const double *gamma_dots,
                 double *Psys,
                 const double dt);
int plasticity_model_construct_elem_ip_map(IP_ID_LIST *elm_ip_map, EPS *eps, const ELEMENT *elem, int ne)
{
  int cnt = 0;
  for(int a=0; a<ne; a++)
  {
    long n_ip = 0;
    int_point(elem[a].toe,&n_ip);         
    Matrix_construct_redim(int,elm_ip_map[a].ip_ids,n_ip, 1);
    elm_ip_map[a].e_id = a;
    elm_ip_map[a].n_ip = n_ip;
    elm_ip_map[a].mat_id = elem[a].mat[0];
      
    for(int b=1; b<=n_ip; b++)
    {
      Mat_v(elm_ip_map[a].ip_ids, b, 1) = cnt;
      cnt++;
    }    
  }
  return cnt;
}

int plasticity_model_cleanup_elem_ip_map(IP_ID_LIST *elm_ip_map, int ne)
{
  int err = 0;
  for(int a=0; a<ne; a++)
    Matrix_cleanup(elm_ip_map[a].ip_ids);

  return err;
} 

/**
 * Private structure for use exclusively with this model and
 * associated functions.
 */
typedef struct plasticity_ctx {
  double *F;
  double dt;    // time increment
  double alpha; // mid point alpha
  double *eFnpa;
  int is_coulpled_with_thermal;
  double *hFn;
  double *hFnp1;  
} plasticity_ctx;

 int compute_M(Matrix(double) *M, 
                      const Matrix(double) *pFn, 
                      const Matrix(double) *N,
                      const Matrix(double) *pFnp1_I,                      
                      const plasticity_ctx *ctx)
{
  int err = 0;
  if(ctx->is_coulpled_with_thermal)
    Matrix_Tns2_AxBxC(*M,1.0,0.0,*pFn,*N,*pFnp1_I);
  else
    Matrix_AxB(*M,1.0,0.0,*pFn,0,*pFnp1_I,0);  
  return err;
}

 int compute_eF(Matrix(double) *eF, 
                      const Matrix(double) *F,
                      const Matrix(double) *hF_I,
                      const Matrix(double) *pF_I,
                      const plasticity_ctx *ctx)
{
  int err = 0;
  if(ctx->is_coulpled_with_thermal)
    Matrix_Tns2_AxBxC(*eF,1.0,0.0,*F,*hF_I,*pF_I);
  else
    Matrix_AxB(*eF,1.0,0.0,*F,0,*pF_I,0);
    
  return err;  
}               
               
static double compute_bulk_mod(const HOMMAT *mat)
{
  return hommat_get_kappa(mat);
}

static int plasticity_int_alg(Constitutive_model *m,
                              const void *ctx);

static int cp_compute_eC(const double * restrict eF,
                         double * restrict eC)
{
  int err = 0;
  memset(eC, 0, DIM_3x3 * sizeof(*eC));
  for (int i = 0; i < DIM_3; i++) {
    for (int j = 0; j < DIM_3; j++) {
      for (int k = 0; k < DIM_3; k++) {
        eC[idx_2(i,j)] += eF[idx_2(k,i)] * eF[idx_2(k,j)];
      }
    }
  }
  return err;
}

static int cp_dev_stress(const double *eC,
                         const HOMMAT *p_hmat,
                         double *Sdev)
{
  int err = 0;
  devStressFuncPtr Stress = getDevStressFunc(-1, p_hmat);
  Stress(eC, p_hmat, Sdev);
  return err;
}

static int cp_dev_tangent(const double *eC,
                          const HOMMAT *p_hmat,
                          double *Ldev)
{
  int err = 0;
  matStiffFuncPtr Tangent = getMatStiffFunc(-1, p_hmat);
  Tangent(eC, p_hmat, Ldev);
  return err;
}

static int plasticity_update(Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Matrix_AeqB(Fs[TENSOR_Fnm1], 1.0,Fs[TENSOR_Fn]);    
  Matrix_AeqB(Fs[TENSOR_Fn], 1.0,Fs[TENSOR_Fnp1]);
  Matrix_AeqB(Fs[TENSOR_pFnm1],1.0,Fs[TENSOR_pFn]);
  Matrix_AeqB(Fs[TENSOR_pFn],1.0,Fs[TENSOR_pFnp1]);  
  Matrix_AeqB(Fs[TENSOR_gamma_dot_n],1.0,Fs[TENSOR_gamma_dot]);  
  Matrix_AeqB(Fs[TENSOR_tau_n],1.0,Fs[TENSOR_tau]);      
  state_var[VAR_g_nm1] = state_var[VAR_g_n];
  state_var[VAR_g_n] = state_var[VAR_g_np1];  
  state_var[VAR_L_nm1] = state_var[VAR_L_n];
  state_var[VAR_L_n] = state_var[VAR_L_np1];  
  return err;
}

static int plasticity_reset(Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Matrix_AeqB(Fs[TENSOR_Fnp1], 1.0,Fs[TENSOR_Fn]);
  Matrix_AeqB(Fs[TENSOR_pFnp1],1.0,Fs[TENSOR_pFn]);
  state_var[VAR_g_np1] = state_var[VAR_g_n];
  state_var[VAR_L_np1] = state_var[VAR_L_n];
  return err;
}

static int plasticity_reset_using_temporal(const Constitutive_model *m, State_variables *var)
{
  int err = 0;
  Matrix(double) *Fs    = m->vars_list[0][m->model_id].Fs;
  Matrix(double) *Fs_in = var->Fs;
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

static int plasticity_update_np1_to_temporal(const Constitutive_model *m, State_variables *var)
{
  int err = 0;
  Matrix(double) *Fs    = var->Fs;
  Matrix(double) *Fs_in = m->vars_list[0][m->model_id].Fs;
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

static int plasticity_save_to_temporal(const Constitutive_model *m, State_variables *var)
{
  int err = 0;
  Matrix(double) *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix(double) *Fs    = var->Fs;
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
  Matrix_AeqB(Fs[TENSOR_tau        ], 1.0, Fs_in[TENSOR_tau     ]);
  Matrix_AeqB(Fs[TENSOR_gamma_dot  ], 1.0, Fs_in[TENSOR_gamma_dot  ]);
  Matrix_AeqB(Fs[TENSOR_tau_n      ], 1.0, Fs_in[TENSOR_tau_n      ]);
  Matrix_AeqB(Fs[TENSOR_gamma_dot_n], 1.0, Fs_in[TENSOR_gamma_dot_n]);
  
  for(int ia=0; ia<VAR_end; ia++)
    state_var[ia] = state_var_in[ia];
    
  return err;
}

static int plasticity_info(Model_var_info **info)
{
  *info = malloc(sizeof(**info));
  int Fno = TENSOR_end;
  
  (*info)->n_Fs = Fno;
  (*info)->F_names = (char **)malloc(sizeof(char*)*Fno);
  for(int a=0; a<Fno; a++)
    (*info)->F_names[a] = (char *)malloc(sizeof(char)*1024);

  sprintf((*info)->F_names[TENSOR_Fn],        "Fn");
  sprintf((*info)->F_names[TENSOR_pFn],       "pFn");
  sprintf((*info)->F_names[TENSOR_Fnp1],      "Fnp1");
  sprintf((*info)->F_names[TENSOR_pFnp1],     "pFnp1");
  sprintf((*info)->F_names[TENSOR_tau],       "tau");
  sprintf((*info)->F_names[TENSOR_R],         "R");
  sprintf((*info)->F_names[TENSOR_gamma_dot], "gamma_dot");
  sprintf((*info)->F_names[TENSOR_Fnm1],      "Fnm1");
  sprintf((*info)->F_names[TENSOR_pFnm1],     "pFnm1");
  
  int varno = VAR_end;
  (*info)->n_vars = varno;
  (*info)->var_names = (char **)malloc(sizeof(char*)*varno);
  for(int a=0; a<varno; a++)
    (*info)->var_names[a] = (char *)malloc(sizeof(char)*1024);
  
  sprintf((*info)->var_names[VAR_L_n],        "L_n");
  sprintf((*info)->var_names[VAR_L_np1],      "L_np1");
  sprintf((*info)->var_names[VAR_g_n],        "g_n");
  sprintf((*info)->var_names[VAR_g_np1],      "g_np1");

  sprintf((*info)->var_names[VAR_L_nm1],      "L_nm1");
  sprintf((*info)->var_names[VAR_g_nm1],      "g_nm1");

  (*info)->n_flags = FLAG_end;
  (*info)->flag_names = malloc(FLAG_end * sizeof( ((*info)->flag_names) ));

  return 0;
}

static int plasticity_get_eF_with_thermal(const Constitutive_model *m,
                                          Matrix(double) *eF,
                                          const Matrix(double) *hFI,
                                          const int stepno)
{
  int err = 0;
  double temp[DIM_3x3];
  Matrix(double) pFI;
  pFI.m_pdata = temp; pFI.m_row = pFI.m_col = DIM_3;
  
  switch(stepno)
  {
    case 0: // n-1
      Matrix_inv(m->vars_list[0][m->model_id].Fs[TENSOR_pFnm1], pFI);
      Matrix_Tns2_AxBxC(*eF,1.0,0.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fnm1],*hFI,pFI);
      break;
    case 1: // n
      Matrix_inv(m->vars_list[0][m->model_id].Fs[TENSOR_pFn], pFI);
      Matrix_Tns2_AxBxC(*eF,1.0,0.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fn],*hFI,pFI);      
      break;
    case 2: // n+1
      Matrix_inv(m->vars_list[0][m->model_id].Fs[TENSOR_pFnp1], pFI);
      Matrix_Tns2_AxBxC(*eF,1.0,0.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fnp1],*hFI,pFI);      
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);

  return err;      
}                          
                                  
static int plasticity_get_pF(const Constitutive_model *m,
                             Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_pFnp1]);
  return err;
}

static int plasticity_get_F(const Constitutive_model *m,
                            Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fnp1]);
  return err;
}

static int plasticity_get_Fn(const Constitutive_model *m,
                             Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fn]);
  return err;
}

static int plasticity_get_pFn(const Constitutive_model *m,
                              Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_pFn]);
  return err;
}

static int plasticity_get_eF(const Constitutive_model *m,
                             Matrix(double) *F)
{
  int err = 0;
  Matrix(double) invFp;
  Matrix_construct_redim(double,invFp,DIM_3,DIM_3);
  err += inv3x3(m->vars_list[0][m->model_id].Fs[TENSOR_pFnp1].m_pdata,invFp.m_pdata);
  Matrix_AxB(*F, 1.0, 0.0, m->vars_list[0][m->model_id].Fs[TENSOR_Fnp1], 0, invFp, 0);
  Matrix_cleanup(invFp);
  return err;
}

static int plasticity_get_eFn(const Constitutive_model *m,
                              Matrix(double) *F)
{
  int err = 0;
  Matrix(double) invFp;
  Matrix_construct_redim(double,invFp,DIM_3,DIM_3);
  err += inv3x3(m->vars_list[0][m->model_id].Fs[TENSOR_pFn].m_pdata, invFp.m_pdata);
  Matrix_AxB(*F, 1.0, 0.0, m->vars_list[0][m->model_id].Fs[TENSOR_Fn], 0, invFp, 0);
  Matrix_cleanup(invFp);
  return err;
}

static int plasticity_get_pFnm1(const Constitutive_model *m,
                              Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_pFnm1]);
  return err;
}

static int plasticity_get_Fnm1(const Constitutive_model *m,
                             Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fnm1]);
  return err;
}

static int plasticity_get_eFnm1(const Constitutive_model *m,
                              Matrix(double) *F)
{
  int err = 0;
  Matrix(double) invFp;
  Matrix_construct_redim(double,invFp,DIM_3,DIM_3);  
  err += inv3x3(m->vars_list[0][m->model_id].Fs[TENSOR_pFnm1].m_pdata,invFp.m_pdata);
  Matrix_AxB(*F, 1.0, 0.0, m->vars_list[0][m->model_id].Fs[TENSOR_Fnm1], 0, invFp, 0);
  Matrix_cleanup(invFp);
  return err;
}

static int plasticity_get_hardening(const Constitutive_model *m,
                                    double *var)
{
  int err = 0;
  *var = m->vars_list[0][m->model_id].state_vars->m_pdata[VAR_g_n];
  return err;
}
static int plasticity_get_hardening_nm1(const Constitutive_model *m,
                                    double *var)
{
  int err = 0;
  *var = m->vars_list[0][m->model_id].state_vars->m_pdata[VAR_g_nm1];
  return err;
}

static int plasticity_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 Matrix(double) *stress)
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  double eC[DIM_3x3] = {};
      
  err += cp_compute_eC(CTX->F,eC);
  err += cp_dev_stress(eC,m->param->p_hmat,stress->m_pdata);

  return err;
}

static int plasticity_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *dudj)
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
  double J = det3x3(CTX->F);
  Pressure(J,m->param->p_hmat,dudj);
  return err;
}

static int plasticity_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  Matrix(double) *tangent)
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  double eC[DIM_3x3] = {};
  
  err += cp_compute_eC(CTX->F,eC);
  err += cp_dev_tangent(eC,m->param->p_hmat,tangent->m_pdata);
  return err;
}

static int plasticity_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *d2udj2)
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  double J = det3x3(CTX->F);

  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
  D_Pressure(J,m->param->p_hmat,d2udj2);
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
  /* The existing function compute_dMdu takes Grad_op for a given
     alpha,beta pair and computes the corresponding dMdu. We will
     generate the required inputs for this function and call it
     alpha*beta times. This should be improved */
 
  auto CTX = (plasticity_ctx *) ctx;
  // shorthand of deformation gradients
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  
  /* compute M at n+1, again, this information is in CM */
  Matrix(double) M, eFn, eFnp1;
  {
    Matrix(double) pFnp1_I;
    Matrix_construct_redim(double, pFnp1_I, DIM_3, DIM_3);   
    Matrix_construct_redim(double, M,       DIM_3, DIM_3);
    Matrix_construct_redim(double, eFn,     DIM_3, DIM_3);
    Matrix_construct_redim(double, eFnp1,   DIM_3, DIM_3);
    
    for(int ia=0; ia<DIM_3x3; ia++)
    {
          M.m_pdata[ia] = 0.0;
        eFn.m_pdata[ia] = 0.0;
      eFnp1.m_pdata[ia] = 0.0;
    }

    err += inv3x3(Fs[TENSOR_pFnp1].m_pdata, pFnp1_I.m_pdata);

    if(CTX->is_coulpled_with_thermal)
    {
      Matrix(double) hFn, hFnp1, hFn_I, hFnp1_I, N, pFn_I;
        hFn.m_row =   hFn.m_col = DIM_3;   hFn.m_pdata = CTX->hFn;
      hFnp1.m_row = hFnp1.m_col = DIM_3; hFnp1.m_pdata = CTX->hFnp1;

      Matrix_construct_redim(double, hFnp1_I, DIM_3, DIM_3);
      Matrix_construct_redim(double, hFn_I,   DIM_3, DIM_3);
      Matrix_construct_redim(double, N,       DIM_3, DIM_3);
      Matrix_construct_redim(double, pFn_I,   DIM_3, DIM_3);     
        
      err += inv3x3(Fs[TENSOR_pFn].m_pdata,   pFn_I.m_pdata);
      err += inv3x3(CTX->hFn,     hFn_I.m_pdata);
      err += inv3x3(CTX->hFnp1, hFnp1_I.m_pdata);    
    
      Matrix_init(N,0.0);
      Matrix_AxB(N,1.0,0.0,hFn,0,hFnp1_I,0);
        
      err += compute_M(&M, Fs+TENSOR_Fn, &N, &pFnp1_I, CTX);

      err += compute_eF(&eFn,   Fs+TENSOR_Fn,  &hFn_I,   &pFn_I,   CTX);
      err += compute_eF(&eFnp1, Fs+TENSOR_Fnp1,&hFnp1_I, &pFnp1_I, CTX);
    
      Matrix_cleanup(pFn_I);
      Matrix_cleanup(hFn_I);
      Matrix_cleanup(hFnp1_I);
      Matrix_cleanup(N);
    }
    else
    {
      Matrix_AxB(M, 1.0, 0.0, m->vars_list[0][m->model_id].Fs[TENSOR_pFn], 0, pFnp1_I, 0);
      err += m->param->get_eFn(m,&eFn);
      err += m->param->get_eF(m,&eFnp1);      
    }
    Matrix_cleanup(pFnp1_I);
  }

  ELASTICITY *elast = (m->param)->cm_elast;
  elast->update_elasticity(elast,eFnp1.m_pdata, 1);

  Matrix(double) L,S;
  S.m_row = S.m_col = DIM_3; S.m_pdata = elast->S;
  L.m_row = DIM_3x3x3x3;
  L.m_col = 1; L.m_pdata = elast->L;  
      
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

  /* make successive calls to compute_dMdu for each node/dof. To avoid
     copying, I am abusing access to the internal data structure of
     the Matrix structure. */
  Matrix(double) dMdu_ab, Grad_op_ab; // no memory created for this, no need Matrix_cleanup
  Matrix_construct(double, Grad_op_ab);
  
  for (int a = 0; a < nne; a++) {
    for(int b = 0; b < ndofn; b++) {
      int idx_ab = idx_4_gen(a,b,0,0,nne,ndofn,DIM_3,DIM_3);
      /* reset dimensions of the matrix objects and set pointer */
      dMdu_ab.m_row = DIM_3;
      dMdu_ab.m_col = DIM_3;
      dMdu_ab.m_pdata = dM_du + idx_ab;      

      /* need to copy Grad_op due to const qualifier */
      Matrix_init_w_array(Grad_op_ab, DIM_3, DIM_3, Grad_op + idx_ab);

      /* call to compute_dMdu */
      
      err += compute_dMdu(m, &dMdu_ab, &Grad_op_ab, &eFn, &eFnp1, &M,
                          &S, &L, g_n, g_np1, tau, gamma_dots, slip.p_sys, CTX->dt);
    }
  }

  /* clean up */
  Matrix_cleanup(M);
  Matrix_cleanup(eFn);
  Matrix_cleanup(eFnp1);
  Matrix_cleanup(Grad_op_ab);
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
  /* The existing function compute_dMdu takes Grad_op for a given
     alpha,beta pair and computes the corresponding dMdu. We will
     generate the required inputs for this function and call it
     alpha*beta times. This should be improved */
  auto CTX = (plasticity_ctx *) ctx;
  
  enum {Fnpa,Fn,Fnp1,eFn,eFnpa,pFnpa,pFnpa_I,pFnp1,pFn,Mnpa,hFnpaI,Fend}; 
  
  // second-order tensors
  Matrix(double) *F2 = malloc(Fend*sizeof(Matrix(double)));
  for (int a = 0; a < Fend; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
  }  
  
  Matrix_init_w_array(F2[Fnp1], DIM_3, DIM_3, CTX->F);  
  err += m->param->get_Fn(m,&F2[Fn]); 
  Matrix_init(F2[eFn],0.0); 
  err += m->param->get_eFn(m,&F2[eFn]);
  Matrix_eye(F2[eFn],DIM_3); // <== recompute F2[eFn] to make Total Lagrangian
    
  err += m->param->get_pF(m,&F2[pFnp1]);
  err += m->param->get_pFn(m,&F2[pFn]);
 
  mid_point_rule(F2[Fnpa].m_pdata, F2[Fn].m_pdata, F2[Fnp1].m_pdata, alpha, DIM_3x3);  
  mid_point_rule(F2[pFnpa].m_pdata, F2[pFn].m_pdata, F2[pFnp1].m_pdata, alpha, DIM_3x3);  

  err += inv3x3(F2[pFnpa].m_pdata, F2[pFnpa_I].m_pdata); // Total Lagrangian

  Matrix_eye(F2[hFnpaI],DIM_3);
  if(CTX->is_coulpled_with_thermal)
  {
    double hFnpa[9];
    //mid_point_rule(hFnpa, CTX->hFn, CTX->hFnp1, alpha, DIM_3x3);  
    err += inv3x3(CTX->hFnp1, F2[hFnpaI].m_pdata);
  }
  
  Matrix_init(F2[Mnpa], 0.0);
  Matrix_init(F2[eFnpa], 0.0);
  Matrix_AxB(F2[Mnpa],1.0,0.0,F2[hFnpaI],0,F2[pFnpa_I],0);  
  Matrix_AxB(F2[eFnpa], 1.0,0.0,F2[Fnpa],0,F2[Mnpa],0);
    
  ELASTICITY *elast = (m->param)->cm_elast;
  elast->update_elasticity(elast,F2[eFnpa].m_pdata, 1);

  Matrix(double) L,S;
  S.m_row = S.m_col = DIM_3; S.m_pdata = elast->S;
  L.m_row = DIM_3x3x3x3;
  L.m_col = 1; L.m_pdata = elast->L;  
      
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
  
  /* make successive calls to compute_dMdu for each node/dof. To avoid
     copying, I am abusing access to the internal data structure of
     the Matrix structure. */
  Matrix(double) dMdu_ab, Grad_op_ab;
  Matrix_construct(double, dMdu_ab); // no memory created for this, no need Matrix_cleanup
  Matrix_construct(double, Grad_op_ab);
  for (int a = 0; a < nne; a++) {
    for(int b = 0; b < ndofn; b++) {
      int idx_ab = idx_4_gen(a,b,0,0,nne,ndofn,DIM_3,DIM_3);
      /* reset dimensions of the matrix objects and set pointer */
      dMdu_ab.m_row = DIM_3;
      dMdu_ab.m_col = DIM_3;
      dMdu_ab.m_pdata = dM_du + idx_ab;

      /* need to copy Grad_op due to const qualifier */
      Matrix_init_w_array(Grad_op_ab, DIM_3, DIM_3, Grad_op + idx_ab);

      /* call to compute_dMdu */
      
      err += compute_dMdu(m, &dMdu_ab, &Grad_op_ab, &F2[eFn], &F2[eFnpa], &F2[Mnpa],
                          &S, &L, g_n, g_npa, tau, gamma_dots, slip.p_sys, CTX->dt);
    }
  }

  // destroyu list of second-order tenosors
  for(int a = 0; a < Fend; a++) {
    Matrix_cleanup(F2[a]);
  }
  free(F2);

  /* clean up */
  Matrix_cleanup(Grad_op_ab);
  destruct_slip_system(&slip);
  free(tau);
  free(gamma_dots);
  return err;
}


static int plasticity_compute_dMdu(const Constitutive_model *m,
                                   const void *ctx,
                                   const double *Grad_op,
                                   const int nne,
                                   const int ndofn,
                                   double *dM_du)
{
  int err = 0;
  auto CTX = (plasticity_ctx *) ctx;
  if(CTX->alpha<0)
    err += plasticity_compute_dMdu_np1(m,ctx,Grad_op,nne,ndofn,dM_du);
  else
    err += plasticity_compute_dMdu_npa(m,ctx,Grad_op,nne,ndofn,dM_du,CTX->alpha);
  return err;
}            

static int cp_write_tensor_restart(FILE *fp, const double *tensor)
{
  int err = 0;
  fprintf(fp, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          tensor[0], tensor[1], tensor[2],
          tensor[3], tensor[4], tensor[5],
          tensor[6], tensor[7], tensor[8]);
  return err;
}

static int cp_read_tensor_restart(FILE *fp, double *tensor)
{
  int err = 0;
  fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &tensor[0], &tensor[1], &tensor[2],
         &tensor[3], &tensor[4], &tensor[5],
         &tensor[6], &tensor[7], &tensor[8]);
  return err;
}

static int plasticity_write_restart(FILE *fp, const Constitutive_model *m)
{

  int err = 0;
  Matrix(double) *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  err += cp_write_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  err += cp_write_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  err += cp_write_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  err += cp_write_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);

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

static int plasticity_read_restart(FILE *fp, Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  err += cp_read_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  err += cp_read_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  err += cp_read_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  err += cp_read_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);

  const int N_SYS = ((((m->param)->cm_mat)->mat_p)->slip)->N_SYS;                                                 
  for(int a=0; a<N_SYS; a++)
    fscanf(fp, "%lf ", Fs[TENSOR_tau_n].m_pdata + a);
  
  for(int a=0; a<N_SYS; a++)
    fscanf(fp, "%lf ", Fs[TENSOR_gamma_dot_n].m_pdata + a);
                                                        
  fscanf(fp, "%lf %lf %lf %lf\n", state_var+VAR_g_n, state_var+VAR_g_nm1, state_var+VAR_L_n, state_var+VAR_L_nm1);

  /* set values at n+1 */
  Matrix_AeqB(Fs[TENSOR_Fnp1],     1.0,Fs[TENSOR_Fn]);                                                   
  Matrix_AeqB(Fs[TENSOR_pFnp1],    1.0,Fs[TENSOR_pFn]);  
  Matrix_AeqB(Fs[TENSOR_tau],      1.0,Fs[TENSOR_tau_n]);
  Matrix_AeqB(Fs[TENSOR_gamma_dot],1.0,Fs[TENSOR_gamma_dot_n]);                                                   
  state_var[VAR_g_np1] = state_var[VAR_g_n];
  state_var[VAR_L_np1] = state_var[VAR_L_n];  
  return 0;  
}

/* THIS IS A FUNCTION STUB. */
static int cp_set_init_vals(Constitutive_model *m)
{
  /* inital values are set in the more convoluted
     read_constitutive_model_parameters->plasticity_model_read_parameters
     calling sequence
  */
  return 0;
}

/* THIS IS A FUNCTION STUB. */
static int cp_read(Model_parameters *p,
                   FILE *in)
{
  int err = 0;

  /* get pointer to parameter data */
  double *param     = p->model_param;
  int    *param_idx = p->model_param_index;
  assert(param     != NULL); // check the pointer
  assert(param_idx != NULL); // check the pointer

  /* scan to non-blank/comment line */
  err += scan_for_valid_line(in);

  /* READ PROPERTIES IN ALPHABETICAL ORDER */  
  int param_in = PARAM_NO-3;
  int match = fscanf(in, "%lf %lf %lf %lf %lf %lf %lf",
                     param + PARAM_gamma_dot_0, param + PARAM_m,    param + PARAM_G0,
                     param + PARAM_g0,          param + PARAM_gs_0, param + PARAM_gamma_dot_s,
                     param + PARAM_w);            

  err += scan_for_valid_line(in);
  
  SLIP_SYSTEM *slip = malloc(sizeof(SLIP_SYSTEM));

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
    match += fscanf(in, "%lf %lf %lf", (slip->ort_angles)+0, (slip->ort_angles)+1, (slip->ort_angles)+2);
    param_in += 3;
  }    
  
  if (match != param_in) err++;
  assert(match == param_in && "Did not read expected number of parameters");

  /* scan past any other comment/blank lines in the block */
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
    match = fscanf(in, "%d %d %d %d %lf %lf %lf", param_idx + PARAM_max_itr_stag,
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
  
  /* not expecting EOF, check and return error if encountered */
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");
    
  MATERIAL_CRYSTAL_PLASTICITY  *mat_p = malloc(sizeof(MATERIAL_CRYSTAL_PLASTICITY));  
  set_properties_crystal_plasticity(mat_p,slip,param[PARAM_gamma_dot_0],param[PARAM_gamma_dot_s], 
                                               param[PARAM_m],          param[PARAM_g0],
                                               param[PARAM_G0],         param[PARAM_gs_0],
                                               param[PARAM_w]);
  
  (p->cm_mat)->mat_p = mat_p;

  return err;
}

int plasticity_model_update_elasticity(const Constitutive_model *m,
                                       const void *ctx_in,
                                       Matrix_double *L,
                                       Matrix_double *S,
                                       const int compute_stiffness)
{
  int err = 0;
  auto ctx = (plasticity_ctx *) ctx_in;
  
  // if transient cases, 
  // get_eF is not working because eF needs to be updated using mid-point alpha
  // below checks whether to use get_eF or give eFnpa in ctx

  if(ctx->eFnpa)
  {
    Matrix(double) eF;
    eF.m_row = eF.m_col = DIM_3; eF.m_pdata = ctx->eFnpa;
    err += constitutive_model_default_update_elasticity(m, &eF, L, S, compute_stiffness);  
  }    
  else
  {
    // shorthand of deformation gradients
    Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;  
    Matrix(double) eF, pFnp1_I;  
    Matrix_construct_init( double,eF,      DIM_3,DIM_3, 0.0); 
    Matrix_construct_redim(double,pFnp1_I, DIM_3, DIM_3);  
  
    err += inv3x3(Fs[TENSOR_pFnp1].m_pdata, pFnp1_I.m_pdata);    
    if(ctx->is_coulpled_with_thermal)
    {
      Matrix(double) hFnp1, hFnp1_I;
      hFnp1.m_row = hFnp1.m_col = DIM_3; hFnp1.m_pdata = ctx->hFnp1;
      Matrix_construct_redim(double, hFnp1_I, DIM_3, DIM_3);

      err += inv3x3(ctx->hFnp1,   hFnp1_I.m_pdata);    
      err += compute_eF(&eF, Fs+TENSOR_Fnp1,&hFnp1_I, &pFnp1_I, ctx);
      Matrix_cleanup(hFnp1_I);
    }
    else
      Matrix_AxB(eF,1.0,0.0,Fs[TENSOR_Fnp1],0,pFnp1_I,0);
      
    err += constitutive_model_default_update_elasticity(m, &eF, L, S, compute_stiffness);  

    Matrix_cleanup(eF);
    Matrix_cleanup(pFnp1_I);
  }
      
  return err;
}

int plasticity_model_get_subdiv_param(const Constitutive_model *m,
                                double *subdiv_param, double dt)
{
  int err = 0;
  SLIP_SYSTEM *slip = (((m->param)->cm_mat)->mat_p)->slip;  
  Matrix(double) *Fs = (m->vars_list[0][m->model_id]).Fs;
  
  double gamma_np1 = 0.0;
  for(int a=0; a<slip->N_SYS; a++)
    gamma_np1 += fabs(Fs[TENSOR_gamma_dot].m_pdata[a]);
  
  *subdiv_param = dt*gamma_np1/MAX_D_GAMMA;
  return err;
}

int plasticity_model_initialize(Model_parameters *p)
{
  int err = 0;

  /* set functions */
  p->integration_algorithm = plasticity_int_alg;
  p->compute_dev_stress = plasticity_dev_stress;
  p->compute_dudj = plasticity_dudj;
  p->compute_dev_tangent = plasticity_dev_tangent;
  p->update_elasticity = plasticity_model_update_elasticity;
  p->compute_d2udj2 = plasticity_d2udj2;
  p->update_state_vars = plasticity_update;
  p->reset_state_vars = plasticity_reset;
  p->reset_state_vars_using_temporal = plasticity_reset_using_temporal;
  p->update_np1_state_vars_to_temporal = plasticity_update_np1_to_temporal;
  p->save_state_vars_to_temporal = plasticity_save_to_temporal;
  p->get_subdiv_param = plasticity_model_get_subdiv_param;
  p->get_var_info = plasticity_info;
  p->get_F     = plasticity_get_F;
  p->get_Fn    = plasticity_get_Fn;
  p->get_Fnm1  = plasticity_get_Fnm1;  
  p->get_pF    = plasticity_get_pF;
  p->get_pFn   = plasticity_get_pFn;
  p->get_pFnm1 = plasticity_get_pFnm1;    
  p->get_eF    = plasticity_get_eF;
  p->get_eFn   = plasticity_get_eFn;
  p->get_eFnm1 = plasticity_get_eFnm1;
  p->get_eF_of_hF = plasticity_get_eF_with_thermal;
    
  p->get_hardening     = plasticity_get_hardening;
  p->get_hardening_nm1 = plasticity_get_hardening_nm1;  
  p->get_plast_strain_var = cm_get_lam_p;
  p->write_restart = plasticity_write_restart;
  p->read_restart  = plasticity_read_restart;
  
  p->destroy_ctx   = plasticity_model_ctx_destroy;
  p->compute_dMdu  = plasticity_compute_dMdu;

  p->set_init_vals = cp_set_init_vals;
  p->read_param = cp_read;

  p->get_size = cp_get_size;
  p->pack = cp_pack;
  p->unpack = cp_unpack;

  p->type = CRYSTAL_PLASTICITY;

  p->n_param = PARAM_NO;
  p->model_param = calloc(PARAM_NO, sizeof(*(p->model_param)));
  p->n_param_index = PARAM_INX_NO;
  p->model_param_index = calloc(PARAM_INX_NO, sizeof(*(p->model_param_index)));
  return err;
}

int plasticity_model_destory(Model_parameters *p)
{
  int err = 0;
  destruct_slip_system(((p->cm_mat)->mat_p)->slip);
  free(((p->cm_mat)->mat_p)->slip);  
  free((p->cm_mat)->mat_p); 
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
  plasticity_ctx *t_ctx = malloc(sizeof(plasticity_ctx));

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

  /* assign handle */
  *ctx = t_ctx;
  return err;
}

int plasticity_model_ctx_destroy(void **ctx)
{
  int err = 0;
  plasticity_ctx *t_ctx = *ctx;
  /* invalidate handle */
  *ctx = NULL;

  // no memory was created
  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;
  t_ctx->hFn   = NULL;
  t_ctx->hFnp1 = NULL; 

  free(t_ctx);
  return err;
}


int compute_dMdu(const Constitutive_model *con,
                 Matrix(double) *_dMdu,
                 const Matrix(double) *_Grad_du,
                 const Matrix(double) *_eFn,
                 const Matrix(double) *_eFnp1,
                 const Matrix(double) *_M,
                 const Matrix(double) *_S,
                 const Matrix(double) *_L,
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
  const double *state_var = con->vars_list[0][con->model_id].state_vars[0].m_pdata;
  
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
  const Tensor<2, const double*> Grad_du(_Grad_du->m_pdata);
  const Tensor<2, const double*> eFn(_eFn->m_pdata);
  const Tensor<2, const double*> eFnp1(_eFnp1->m_pdata);
  const Tensor<2, const double*> M(_M->m_pdata);
  const Tensor<2, const double*> S(_S->m_pdata);
  const Tensor<4, const double*> L(_L->m_pdata);

  const Tensor<2> C = eFnp1(k,i).to(i,k) * eFnp1(k,j);             // eFnp1' * eFnp1
  Tensor<2> MI = {};
  int err = inv3x3(_M->m_pdata, MI.data);
  const Tensor<2> MI_x_C = MI(j,i) * C(j,k);                       // M^{-T} * C
  const Tensor<2> M_x_eFn = M(j,i).to(i,j) * eFn(k,j).to(j,k);     // M' * eFn'

  // tensors are initialized to 0
  Tensor<2> sum_aC = {};
  Tensor<2> sum_Pa = {};
  Tensor<2> sum_aD = {};
  Tensor<4> U = {};
  Tensor<4> B = {};

  // set U to the 9x9 identity scaled by 1.0/dt
  U(i,j,k,l) = (1.0/dt) * ttl::identity(i,j,k,l);                    

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
    const Tensor<2> aD = eFnp1(i,j) * AA(j,k) * M_x_eFn(k,l);
    
    sum_Pa(i,j) += drdg * Pa(i,j);
    sum_aC(i,j) += R2_a * aC(i,j);
    sum_aD(i,j) += R2_a * aD(i,j);

    // perform the Kronecker product using ttl and scales it by drdtau
    U(i,j,k,l) += drdtau * aC(i,j) * Pa(k,l);
    B(i,j,k,l) += drdtau * aD(i,j) * Pa(k,l);
  }

  double R1 = R4/(1.0-R4*sum_1gm1gm);

  // perform the Kronecker product using ttl and scales it by R1
  U(i,j,k,l) += R1 * sum_aC(i,j) * sum_Pa(k,l);
  B(i,j,k,l) += R1 * sum_aD(i,j) * sum_Pa(k,l);

  // cast _dmdu as a ttl tensor and compute its value
  // -1 * (inverse(U) * B:Grad_du)
  Tensor<2,double*>(_dMdu->m_pdata)(i,j) = -(ttl::inverse(U)(i,j,k,l) * 
					     B(k,l,m,n) * Grad_du(m,n));  

  return err;
}


static int plasticity_int_alg(Constitutive_model *m,
                              const void *ctx)
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
    
  enum {M,eFnp1,C,pFnp1_I,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3 ,0.0);
  }
    
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double g_n   = state_var[VAR_g_n];
  double g_np1 = state_var[VAR_g_np1];
  double L_np1 = state_var[VAR_L_np1];
  
  MATERIAL_CONSTITUTIVE_MODEL *cm_mat = (m->param)->cm_mat;
  ELASTICITY *elasticity = (m->param)->cm_elast;

  // compute slip system and rotate orientation
  SLIP_SYSTEM *slip_in = (cm_mat->mat_p)->slip;
  SLIP_SYSTEM slip;
  construct_slip_system(&slip,slip_in->unit_cell);
  rotate_crystal_orientation(&slip, Fs[TENSOR_R].m_pdata, slip_in);
  (cm_mat->mat_p)->slip = &slip;
 
  // perform integration algorithm for the crystal plasticity 
  if(CTX->is_coulpled_with_thermal)
  {    
    err += staggered_Newton_Rapson_generalized(Fs[TENSOR_pFnp1].m_pdata,
                                               F2[M].m_pdata, 
                                               &g_np1, &L_np1, 
                                               Fs[TENSOR_pFn].m_pdata, 
                                               Fs[TENSOR_Fn].m_pdata, 
                                               Fs[TENSOR_Fnp1].m_pdata,
                                               CTX->hFn, CTX->hFnp1,
                                               g_n, dt, cm_mat, elasticity, &solver_info);  
  }
  else
  {
    err += staggered_Newton_Rapson(Fs[TENSOR_pFnp1].m_pdata,
                                   F2[M].m_pdata, &g_np1, &L_np1, 
                                   Fs[TENSOR_pFn].m_pdata, 
                                   Fs[TENSOR_Fn].m_pdata, 
                                   Fs[TENSOR_Fnp1].m_pdata,
                                   g_n, dt, cm_mat, elasticity, &solver_info);     
  }    
 
  // update compute values from integration algorithm  
  state_var[VAR_g_np1] =  g_np1;
  state_var[VAR_L_np1] =  L_np1;
  
  err += inv3x3(Fs[TENSOR_pFnp1].m_pdata,F2[pFnp1_I].m_pdata);

  if(CTX->is_coulpled_with_thermal)
  {
    Matrix(double) hFnp1,hFnp1_I;
    hFnp1.m_row = hFnp1.m_col = DIM_3; hFnp1.m_pdata = CTX->hFnp1;        
    Matrix_construct_redim(double,hFnp1_I,DIM_3,DIM_3);
    
    Matrix_inv(hFnp1,hFnp1_I);
    err += compute_eF(F2+eFnp1, Fs+TENSOR_Fnp1, &hFnp1_I, F2+pFnp1_I, CTX); 
    Matrix_cleanup(hFnp1_I);
  }
  else
    Matrix_AxB(F2[eFnp1],1.0,0.0,Fs[TENSOR_Fnp1],0,F2[pFnp1_I],0);  
  
  Matrix_AxB(F2[C], 1.0, 0.0, F2[eFnp1],1,F2[eFnp1],0);
  
  elasticity->update_elasticity(elasticity,F2[eFnp1].m_pdata, 0);
  err += compute_tau_alphas(Fs[TENSOR_tau].m_pdata,F2[C].m_pdata, elasticity->S, &slip);
  err += compute_gamma_dots(Fs[TENSOR_gamma_dot].m_pdata, Fs[TENSOR_tau].m_pdata, g_np1, cm_mat->mat_p);

  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);   
  free(F2);
  
  (cm_mat->mat_p)->slip = slip_in;
  destruct_slip_system(&slip);

  return err;
}

int plasticity_model_construct_rotation(EPS *eps, Matrix(int) *e_ids, Matrix(double) *angles)
{
  int err = 0;
  double Ax[DIM_3x3], Az[DIM_3x3], Ay[DIM_3x3];   
  
  for(int a=0; a<e_ids->m_row; a++)
  {
    int id = Mat_v(*e_ids, a+1, 1);
    if(id<0)
      continue;        
    
    int ip = Mat_v(*e_ids, a+1, 2);
    Constitutive_model *m = &(eps[id].model[ip]);
    Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;    
    double phi   = angles->m_pdata[a*DIM_3+0]; // NOTE: phi = Mat_v(*angles, a+1, 1) is not working, do not know why.
    double theta = angles->m_pdata[a*DIM_3+1];
    double psi   = angles->m_pdata[a*DIM_3+2];

    // compute rotation matrix of Euler angles
    // Fs[TENSOR_R] = Az(psi)*Ay(theta)*Ax(phi)
    err += rotation_matrix_of_Euler_angles(Fs[TENSOR_R].m_pdata, 
                                           Ax, Ay, Az, phi, theta, psi);     
  }
  return err;
}

int plasticity_model_read_orientations(Matrix(int) *e_ids, Matrix(double) *angles, IP_ID_LIST *elm_ip_map, char *fn_in, int myrank, int ne)
{
  int err = 0;
  char fn[1024], line[1024];
  sprintf(fn, "%s_%d.in", fn_in, myrank);
  FILE *fp = fopen(fn, "r");
  if(fp==NULL)
  { 
    printf("fail to read [%s]\n", fn);
    printf("set default onrientation [R=I]\n");     
    return err;
  }
  
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
	    continue;
        
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
           Mat_v(*e_ids,  ip_id+1, 1) = a;  // +1 is needed because ip_id starts from 0
           Mat_v(*e_ids,  ip_id+1, 2) = b;
           Mat_v(*angles, ip_id+1, 1) = x1;    
           Mat_v(*angles, ip_id+1, 2) = x2;        
           Mat_v(*angles, ip_id+1, 3) = x3;          
        }
      }
    }  
    else
    {  
      int ip_id = elm_ip_map[e].ip_ids.m_pdata[ip];
      Mat_v(*e_ids,  ip_id+1, 1) = e;  // +1 is needed because ip_id starts from 0
      Mat_v(*e_ids,  ip_id+1, 2) = ip;
      Mat_v(*angles, ip_id+1, 1) = x1;    
      Mat_v(*angles, ip_id+1, 2) = x2;        
      Mat_v(*angles, ip_id+1, 3) = x3;
    }
  } 
  
  fclose(fp);
  return err;
}

int plasticity_model_generate_random_orientation_element(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix(int) *e_ids, Matrix(double) *angles, int diff_ort_at_ip)
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
      Mat_v(*angles,ip_id+1,1) = phi;
      Mat_v(*angles,ip_id+1,2) = theta;
      Mat_v(*angles,ip_id+1,3) = psi;            
      Mat_v(*e_ids, ip_id+1,1) = a;
      Mat_v(*e_ids, ip_id+1,2) = ip;                     
    }  
  }
  free(angle_temp);
  return err;
}

int plasticity_model_generate_random_orientation_crystal(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix(int) *e_ids, Matrix(double) *angles)
{
  int err = 0;
  
  double temp[3];
  err += generate_random_crystal_orientation(temp, 1);
      
  int cnt = 0;          
  for(int a=0; a<ne; a++)
  {
    int n_ip = elm_ip_map[a].n_ip;
    int mat = elm_ip_map[a].mat_id;
    if(mat!=mat_id)
      continue;    
    
    for(int ip=0; ip<n_ip; ip++) 
    {
      int ip_id = elm_ip_map[a].ip_ids.m_pdata[ip];
      Mat_v(*angles,ip_id+1,1) = temp[0];
      Mat_v(*angles,ip_id+1,2) = temp[1];
      Mat_v(*angles,ip_id+1,3) = temp[2];            
      Mat_v(*e_ids, ip_id+1,1) = a;
      Mat_v(*e_ids, ip_id+1,2) = ip;                      
    } 
  } 
  return err;
}

int plasticity_model_set_given_orientation_crystal(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix(int) *e_ids, Matrix(double) *angles, double *angle_in)
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
      Mat_v(*angles,ip_id+1,1) = angle_in[0];
      Mat_v(*angles,ip_id+1,2) = angle_in[1];
      Mat_v(*angles,ip_id+1,3) = angle_in[2];            
      Mat_v(*e_ids, ip_id+1,1) = a;
      Mat_v(*e_ids, ip_id+1,2) = ip;                      
    } 
  } 
  return err;
}  


int plasticity_model_set_zero_angles(const int ne, const IP_ID_LIST *elm_ip_map, int mat_id, Matrix(int) *e_ids, Matrix(double) *angles)
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

      Mat_v(*angles,ip_id+1,1) = 0.0;
      Mat_v(*angles,ip_id+1,2) = 0.0;
      Mat_v(*angles,ip_id+1,3) = 0.0;            
      Mat_v(*e_ids, ip_id+1,1) = a;
      Mat_v(*e_ids, ip_id+1,2) = ip;                     
    }  
  }
  return err;
}

int plasticity_model_set_orientations(EPS *eps,
                                const int ne,
                                const ELEMENT *elem,
                                const int n_mat,
                                const Model_parameters *param_list)
{
  int err = 0;
  int crystal_plasticity_included = 0;
  for(int i = 0; i < n_mat; i++)
  {
    if(param_list[i].type==CRYSTAL_PLASTICITY)
      crystal_plasticity_included++;
  }
  
  if(crystal_plasticity_included==0)
    return err;  
  
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int myrank = 0;  
  MPI_Comm_rank (mpi_comm,&myrank);

  // build element ip ids that will be used to assign element orientation
  IP_ID_LIST *elm_ip_map = malloc(sizeof(IP_ID_LIST)*ne);
  int cnt_of_ips = plasticity_model_construct_elem_ip_map(elm_ip_map, eps, elem, ne);
    
  Matrix(int) e_ids;
  Matrix(double) angles;
  Matrix_construct_init(int,e_ids,cnt_of_ips, 2, -1);
  Matrix_construct_init(double,angles,cnt_of_ips,DIM_3,0.0);  
  
  for(int i = 0; i < n_mat; i++)
  {
    if(param_list[i].type==CRYSTAL_PLASTICITY)
    {
      SLIP_SYSTEM *slip = ((param_list[i].cm_mat)->mat_p)->slip;
      if(slip->ort_option[0]==2)
      {
        err += plasticity_model_read_orientations(&e_ids, &angles, elm_ip_map, slip->ort_file_in, myrank, ne);
        break;
      }
    }
  }
  
  int save_orientations = 0;
  for(int i = 0; i < n_mat; i++)
  {
    if(param_list[i].type==CRYSTAL_PLASTICITY)
    {
      SLIP_SYSTEM *slip = ((param_list[i].cm_mat)->mat_p)->slip;     
      switch(slip->ort_option[0])
      {
        case -1:
          plasticity_model_set_zero_angles(ne, elm_ip_map, param_list[i].mat_id, &e_ids, &angles);
          break;
        case 0:
        {  
          err += plasticity_model_generate_random_orientation_element(ne, elm_ip_map, param_list[i].mat_id, &e_ids, &angles, slip->ort_option[1]);
          save_orientations++;
          break;
        }
        case 1:
        {
          err += plasticity_model_generate_random_orientation_crystal(ne, elm_ip_map, param_list[i].mat_id, &e_ids, &angles);
          save_orientations++;
          break;
        }  
        case 3:
        {
          err += plasticity_model_set_given_orientation_crystal(ne, elm_ip_map, param_list[i].mat_id, &e_ids, &angles, slip->ort_angles);
          break;
        }
        default:
          break;
      }
    }  
  }
  plasticity_model_construct_rotation(eps, &e_ids, &angles); 
  
  char file_in_ort[1024], default_ort_dir[1024];    
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
    char fn_orientation[1024];
    sprintf(fn_orientation, "%s/orientation_%d.in", default_ort_dir, myrank);  
    FILE *fp_ort = fopen(fn_orientation, "w");
    fprintf(fp_ort, "# Element (crystal) orientations are generated randomly\n");
    fprintf(fp_ort, "# element_ID, Integration_point_ID, phi [radian], theta [radian], psi [radian]\n");  
    for(int a=1; a<=e_ids.m_row; a++)
    {
      if(Mat_v(e_ids, a, 1)<0)
        continue;
    
      fprintf(fp_ort, "%d %d %e %e %e\n", Mat_v(e_ids, a, 1), Mat_v(e_ids, a, 2),
                                        Mat_v(angles, a, 1), Mat_v(angles, a, 2), Mat_v(angles, a, 3));             
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
            
      /* set state variables to initial values */
      state_var[VAR_g_nm1] = mat_p->g0;      
      state_var[VAR_g_n]   = mat_p->g0;      
      state_var[VAR_g_np1] = mat_p->g0;    
      state_var[VAR_L_nm1] = 0.0;
      state_var[VAR_L_n]   = 0.0;  
      state_var[VAR_L_np1] = 0.0;
      
      /* set the dimensions for the tau and gamma_dot varables */
      Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
      const int N_SYS = (mat_p->slip)->N_SYS;
      Matrix_redim(Fs[TENSOR_tau],N_SYS, 1);
      Matrix_redim(Fs[TENSOR_tau_n],N_SYS, 1);
      Matrix_redim(Fs[TENSOR_gamma_dot],N_SYS, 1);      
      Matrix_redim(Fs[TENSOR_gamma_dot_n],N_SYS, 1);
      
      /* intitialize to zeros */
      for(int ia=0; ia<N_SYS; ia++)
      {
                Fs[TENSOR_tau].m_pdata[ia] = 0.0;
              Fs[TENSOR_tau_n].m_pdata[ia] = 0.0;
          Fs[TENSOR_gamma_dot].m_pdata[ia] = 0.0;
        Fs[TENSOR_gamma_dot_n].m_pdata[ia] = 0.0;
      }  
    }
  }

  Matrix_cleanup(e_ids);
  Matrix_cleanup(angles);
  err += plasticity_model_cleanup_elem_ip_map(elm_ip_map, ne);
  free(elm_ip_map);  
  return err;      
}

void test_crystal_plasticity_single_crystal(void)
{
  // test for defined F
  // F = [1 - t,       0,       0
  //          0, 1 + t/2,       0
  //          0,       0, 1 + t/2];   
  double lame1 = 75600.0;
  double lame2     = 26100.0;
  double E = 70.0e+3;
  double nu = 0.25;
  
  double gamma_dot_0 = 1.0;
  double gamma_dot_s = 50.0e+9;
  double m           = 0.1;  
  double g0          = 210.0;
  double G0          = 200.0;
  double gs_0        = 330.0;
  double w           = 0.005;
  
  int max_itr_stag      = 100;
  int max_itr_hardening = 5;
  int max_itr_M         = 100;
  double tol_hardening  = 1.0e-6;
  double tol_M          = 1.0e-6;
  double computer_zero  = 1.0e-15;

  // create material properties: Elasticity
  MATERIAL_ELASTICITY mat_e;  
  set_properties_using_E_and_nu(&mat_e,E,nu); 
  // or you can use : set_properties_using_Lame_constants(&mat_e,lame1,lame2);
  //print_material_property_elasticity(&mat_e); // <= this is optional
  
  // create slip system : 0 for FCC
  // it is needed to setup plasticity
  SLIP_SYSTEM slip;
  construct_slip_system(&slip,0); 
  
  // create material properties: Plasticity
  MATERIAL_CRYSTAL_PLASTICITY mat_p;
  set_properties_crystal_plasticity(&mat_p,&slip,gamma_dot_0,gamma_dot_s, 
                                     m,g0,G0,gs_0,w);
  //print_material_property_crystal_plasticity(&mat_p);  // <= this is optional 

  // create material plasticity: it needs material properties for elasticity and plasticity
  MATERIAL_CONSTITUTIVE_MODEL mat;
  set_properties_constitutive_model(&mat,&mat_e,&mat_p);
  
  // create solver info: criteria for numerical iterations
  CRYSTAL_PLASTICITY_SOLVER_INFO solver_info;
  set_crystal_plasticity_solver_info(&solver_info,max_itr_stag,
                                                  max_itr_hardening,
                                                  max_itr_M,
                                                  tol_hardening,
                                                  tol_M,
                                                  computer_zero);  
  //print_crystal_plasticity_solver_info(&solver_info); // <= this is optional
  
  // create elasticity object for integration
  // this creates memory for stress and elasticity tensor s.t. requres destructor
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);  

  // set variables for integration
  enum {M,MI,pFn,pFnp1,Fn,Fnp1,eFnp1,eFPK2,pFnp1_I,sigma,PK2dev,sigma_dev,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3,0.0);
    Matrix_eye(F2[a],DIM_3);
  } 
  
  double g_n,g_np1;
  g_n = g_np1 = mat_p.g0;
  
  double dt = 0.1;
    
  // start integration  
  Matrix(double) PK2;
  PK2.m_row = PK2.m_col = DIM_3; PK2.m_pdata = elast.S;
  
  FILE *fp = fopen("single_crystal_results.txt", "w");
  
  for(int a = 1; a<=1000; a++)
  {
    double lambda = 0.0;
    double t = a*dt;
    
    // compute total deformation gradient using velocity gradient
    // Fnp1_Implicit(F2[Fnp1].m_pdata, F2[Fn].m_pdata, F2[L].m_pdata, dt); 
    //define Fnp1    
    Mat_v(F2[Fnp1],1,1) = 1.0 - t*1.0e-3;
    Mat_v(F2[Fnp1],2,2) = Mat_v(F2[Fnp1],3,3) = 1.0 + t*0.5*1.0e-3;
    
    staggered_Newton_Rapson(F2[pFnp1].m_pdata,F2[M].m_pdata, &g_np1, &lambda, 
                            F2[pFn].m_pdata, F2[Fn].m_pdata,F2[Fnp1].m_pdata, 
                            g_n, dt, &mat, &elast, &solver_info);
    Matrix_AeqB(F2[pFn],1.0,F2[pFnp1]);
    Matrix_AeqB(F2[Fn],1.0,F2[Fnp1]);  
    
    g_n = g_np1;
    
    
    // print result at time t
    double det_eF;
    
    // compute Caush stress
    Matrix_inv(F2[pFnp1], F2[pFnp1_I]);
    Matrix_AxB(F2[eFnp1],1.0,0.0,F2[Fnp1],0,F2[pFnp1_I],0);
    Matrix_det(F2[eFnp1], det_eF);
    Matrix_AxB(F2[eFPK2],1.0,0.0,F2[eFnp1],0,PK2,0);
    Matrix_AxB(F2[sigma],1.0/det_eF,0.0,F2[eFPK2],0,F2[eFnp1],1);  
    
    double trPK2, tr_sigma;
    Matrix_trace(PK2,trPK2);
    Matrix_trace(F2[sigma],tr_sigma);
    Matrix_eye(F2[PK2dev], DIM_3);
    Matrix_eye(F2[sigma_dev], DIM_3);
        
    Matrix_AplusB(F2[PK2dev],    1.0, PK2,      -trPK2/3.0, F2[PK2dev]);
    Matrix_AplusB(F2[sigma_dev], 1.0, F2[sigma], -tr_sigma/3.0, F2[sigma_dev]);    
    
    double norm_sigma, norm_PK2;
    Matrix_ddot(F2[PK2dev],F2[PK2dev],norm_PK2);    
    Matrix_ddot(F2[sigma_dev],F2[sigma_dev],norm_sigma);
    
    double sigma_eff=sqrt(3.0/2.0*norm_sigma);
    double PK2_eff = sqrt(3.0/2.0*norm_PK2);    

    fprintf(fp, "%e %e %e %e %e %e\n",t,sigma_eff,PK2_eff, g_np1, 0.0, Mat_v(PK2,1,1));
  }    
  
  fclose(fp);  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  free(F2);    
  destruct_elasticity(&elast);
  destruct_slip_system(&slip);
}


