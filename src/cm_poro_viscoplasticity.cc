
/// Define poro-visco plasticity model functions for the constitutive model interface
/// 
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// Alberto Salvadori, [1], <asalvad2@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

//#include <ttl/ttl.h>
#include "cm_poro_viscoplasticity.h"

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

#include <iostream>
#include <iomanip>
#include <time.h>
#include <ttl/ttl.h>

#include "KMS-IJSS2017.h"
#include "pvp_interface.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX_D_ALPHA 0.005


static const int FLAG_end = 0;

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

Define_Matrix(double);
Define_Matrix(int);

void matrix_print_name(Matrix(double) *A, char name[])
{
  printf("double %s[9] = {\t", name);
  for(int ia=0; ia<9; ia++)
  {
    if(ia==3 || ia==6)
      printf("\n\t\t\t");
    printf("%.17e", A->m_pdata[ia]);
    if(ia<8)
      printf(",");
  }
  printf("};\n");
}

int test_cm_poro_viscoplasticity_model(Constitutive_model *m);

/// Private structure for use exclusively with this model and
// associated functions.
typedef struct poro_viscoplasticity_ctx {
  double *F;
  double dt;    // time increment
  double alpha; // mid point alpha
  double *eFnpa;
  int is_coulpled_with_thermal;
  double *hFn;
  double *hFnp1;
  int npa;  
} poro_viscoplasticity_ctx;


int poro_viscoplasticity_model_ctx_build(void **ctx,
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
  poro_viscoplasticity_ctx *t_ctx = (poro_viscoplasticity_ctx *) malloc(sizeof(poro_viscoplasticity_ctx));

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

  /* assign handle */
  *ctx = t_ctx;
  return err;
}

static int cm_pvp_int_alg(Constitutive_model *m,
                              const void *ctx)
{
  int err = 0;
  
//  err += test_cm_poro_viscoplasticity_model(m);


  auto CTX = (poro_viscoplasticity_ctx *) ctx;

  const double dt = CTX->dt;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *param     = (m->param)->model_param;
  int    *param_idx = (m->param)->model_param_index;  
  
  memcpy(Fs[TENSOR_Fnp1].m_pdata, CTX->F, DIM_3x3*sizeof(double));
  KMS_IJSS2017_Parameters *mat_pvp = (m->param)->cm_mat_2->mat_pvp;
  
  int myrank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  err += pvp_intf_perform_integration_algorithm(Fs[TENSOR_Fnp1].m_pdata,
                                                Fs[TENSOR_Fn].m_pdata,
                                                Fs[TENSOR_pFnp1].m_pdata,
                                                Fs[TENSOR_pFn].m_pdata,
                                                vars+VAR_pc_np1,
                                                vars[VAR_pc_n],
                                                mat_pvp,
                                                dt);                                                                                                  
  return err;
}

static int cm_pvp_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 Matrix(double) *stress)
{
  int err = 0;
  auto CTX = (poro_viscoplasticity_ctx *) ctx;
  double eC[DIM_3x3] = {};
      
//  err += cp_compute_eC(CTX->F,eC);
//  err += cp_dev_stress(eC,m->param->p_hmat,stress->m_pdata);

  return err;
}

static int cm_pvp_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *dudj)
{
  int err = 0;
  auto CTX = (poro_viscoplasticity_ctx *) ctx;
//  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
//  double J = det3x3(CTX->F);
//  Pressure(J,m->param->p_hmat,dudj);
  return err;
}

static int cm_pvp_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  Matrix(double) *tangent)
{
  int err = 0;
  auto CTX = (poro_viscoplasticity_ctx *) ctx;
  double eC[DIM_3x3] = {};
  
//  err += cp_compute_eC(CTX->F,eC);
//  err += cp_dev_tangent(eC,m->param->p_hmat,tangent->m_pdata);
  return err;
}


int cm_pvp_update_elasticity(const Constitutive_model *m,
                                       const void *ctx_in,
                                       Matrix_double *L_in,
                                       Matrix_double *S_in,
                                       const int compute_stiffness)
{
  int err = 0;
  auto ctx = (poro_viscoplasticity_ctx *) ctx_in;
  
  double *vars = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  
  double pc = vars[VAR_pc_np1];
  switch(ctx->npa)
  {
    case 0:
      mid_point_rule(&pc, vars + VAR_pc_nm1, vars + VAR_pc_n, ctx->alpha, 1);
      break;
    case 1:  
      mid_point_rule(&pc, vars + VAR_pc_n, vars + VAR_pc_np1, ctx->alpha, 1);
      break;
  }
      
  KMS_IJSS2017_Parameters *mat_pvp = (m->param)->cm_mat_2->mat_pvp;    
  
  double *L = NULL;
  if(compute_stiffness)
    L = L_in->m_pdata;
  
  if(ctx->eFnpa)
  {
    err += pvp_intf_update_elasticity(ctx->eFnpa, 
                                      pc, 
                                      S_in->m_pdata, 
                                      L, 
                                      mat_pvp,
                                      compute_stiffness);
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
      Matrix_Tns2_AxBxC(eF,1.0,0.0,Fs[TENSOR_Fnp1],hFnp1_I,pFnp1_I);
      Matrix_cleanup(hFnp1_I);
    }
    else      
      Matrix_AxB(eF,1.0,0.0,Fs[TENSOR_Fnp1],0,pFnp1_I,0);
    
    err += pvp_intf_update_elasticity(eF.m_pdata, 
                                      pc, 
                                      S_in->m_pdata, 
                                      L, 
                                      mat_pvp,
                                      compute_stiffness);
    Matrix_cleanup(eF);
    Matrix_cleanup(pFnp1_I);
  }
  if(err!=0)
    printf("model id = %e\n", m->model_id);
     
  return err;
}

static int cm_pvp_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *d2udj2)
{
  int err = 0;
  auto CTX = (poro_viscoplasticity_ctx *) ctx;
  double J = det3x3(CTX->F);

//  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
//  D_Pressure(J,m->param->p_hmat,d2udj2);
  return err;
}

static int cm_pvp_update(Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var  = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fnm1].m_pdata[ia]  = Fs[TENSOR_Fn].m_pdata[ia];
    Fs[TENSOR_Fn].m_pdata[ia]    = Fs[TENSOR_Fnp1].m_pdata[ia];
    Fs[TENSOR_pFnm1].m_pdata[ia] = Fs[TENSOR_pFn].m_pdata[ia];
    Fs[TENSOR_pFn].m_pdata[ia]   = Fs[TENSOR_pFnp1].m_pdata[ia];    
  }
  state_var[VAR_pc_nm1] = state_var[VAR_pc_n];
  state_var[VAR_pc_n]   = state_var[VAR_pc_np1];
  return err;
}

static int cm_pvp_reset(Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double *state_var = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fnp1].m_pdata[ia]  = Fs[TENSOR_Fn].m_pdata[ia];
    Fs[TENSOR_pFnp1].m_pdata[ia] = Fs[TENSOR_pFn].m_pdata[ia];
  }  
  state_var[VAR_pc_np1] = state_var[VAR_pc_n];
  return err;
}

int cm_pvp_get_subdiv_param(const Constitutive_model *m,
                                double *subdiv_param, double dt)
{
  int err = 0;
  Matrix(double) *Fs = (m->vars_list[0][m->model_id]).Fs;
  
  double alpha = 0.0;
  *subdiv_param = alpha/MAX_D_ALPHA;
  return err;
}

static int cm_pvp_info(Model_var_info **info)
{
  *info = malloc(sizeof(**info));
  int Fno = TENSOR_end;
  
  (*info)->n_Fs = Fno;
  (*info)->F_names = (char **)malloc(sizeof(char*)*Fno);
  for(int a=0; a<Fno; a++)
    (*info)->F_names[a] = (char *)malloc(sizeof(char)*1024);

  sprintf((*info)->F_names[TENSOR_Fn],    "Fn");
  sprintf((*info)->F_names[TENSOR_pFn],   "pFn");
  sprintf((*info)->F_names[TENSOR_Fnp1],  "Fnp1");
  sprintf((*info)->F_names[TENSOR_pFnp1], "pFnp1");
  sprintf((*info)->F_names[TENSOR_Fnm1],  "Fnm1");
  sprintf((*info)->F_names[TENSOR_pFnm1], "pFnm1");
  
  int varno = VAR_end;
  (*info)->n_vars = varno;
  (*info)->var_names = (char **)malloc(sizeof(char*)*varno);
  for(int a=0; a<varno; a++)
    (*info)->var_names[a] = (char *)malloc(sizeof(char)*1024);
  
  sprintf((*info)->var_names[VAR_pc_n],   "pc_n");
  sprintf((*info)->var_names[VAR_pc_np1], "pc_np1");
  sprintf((*info)->var_names[VAR_pc_nm1], "pc_nm1");

  (*info)->n_flags = FLAG_end;
  (*info)->flag_names = malloc(FLAG_end * sizeof( ((*info)->flag_names) ));

  return 0;
}


static int cm_pvp_get_F(const Constitutive_model *m,
                            Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fnp1]);
  return err;
}

static int cm_pvp_get_Fn(const Constitutive_model *m,
                             Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fn]);
  return err;
}

static int cm_pvp_get_Fnm1(const Constitutive_model *m,
                             Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_Fnm1]);
  return err;
}

static int cm_pvp_get_pF(const Constitutive_model *m,
                             Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_pFnp1]);
  return err;
}

static int cm_pvp_get_pFn(const Constitutive_model *m,
                              Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_pFn]);
  return err;
}

static int cm_pvp_get_pFnm1(const Constitutive_model *m,
                              Matrix(double) *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars_list[0][m->model_id].Fs[TENSOR_pFnm1]);
  return err;
}

static int cm_pvp_get_eF(const Constitutive_model *m,
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

static int cm_pvp_get_eFn(const Constitutive_model *m,
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

static int cm_pvp_get_eFnm1(const Constitutive_model *m,
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

static int cm_pvp_get_eF_with_thermal(const Constitutive_model *m,
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

static int cm_pvp_reset_using_temporal(const Constitutive_model *m, State_variables *var)
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
  
  state_var[VAR_pc_n]   = state_var_in[VAR_pc_n];
  state_var[VAR_pc_nm1] = state_var_in[VAR_pc_nm1];  
    
  return err;
}

static int cm_pvp_update_np1_to_temporal(const Constitutive_model *m, State_variables *var)
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
  
  state_var[VAR_pc_np1]   = state_var_in[VAR_pc_np1];
  

  return err;
}

static int cm_pvp_save_to_temporal(const Constitutive_model *m, State_variables *var)
{
  int err = 0;
  Matrix(double) *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix(double) *Fs    = var->Fs;
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

static int cm_pvp_get_hardening(const Constitutive_model *m,
                                    double *var)
{
  int err = 0;
  *var = m->vars_list[0][m->model_id].state_vars->m_pdata[VAR_pc_np1];
        
  return err;
}
static int cm_pvp_get_hardening_nm1(const Constitutive_model *m,
                                    double *var)
{
  int err = 0;
  *var = m->vars_list[0][m->model_id].state_vars->m_pdata[VAR_pc_nm1];
  return err;
}

static int cm_pvp_write_restart(FILE *fp, const Constitutive_model *m)
{

  int err = 0;
  Matrix(double) *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var  = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);
                                                   
  fprintf(fp, "%.17e %.17e\n", state_var[VAR_pc_n], state_var[VAR_pc_nm1]);
  
  return err;
}

static int cm_pvp_read_restart(FILE *fp, Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *state_var = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  err += cm_read_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  err += cm_read_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  err += cm_read_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  err += cm_read_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);
                                                        
  fscanf(fp, "%lf %lf\n", state_var+VAR_pc_n, state_var+VAR_pc_nm1);
  
  Matrix_AeqB(Fs[TENSOR_Fnp1],     1.0,Fs[TENSOR_Fn]);
  Matrix_AeqB(Fs[TENSOR_pFnp1],    1.0,Fs[TENSOR_pFn]);
  state_var[VAR_pc_np1] = state_var[VAR_pc_n];

  err += cm_pvp_reset(m);
  return 0;  
}

int cm_pvp_ctx_destroy(void **ctx)
{
  int err = 0;
  poro_viscoplasticity_ctx *t_ctx = *ctx;
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

static int cm_pvp_compute_dMdu(const Constitutive_model *m,
                                   const void *ctx,
                                   const double *Grad_op,
                                   const int nne,
                                   const int nsd,
                                   double *dM_du)
{
  int err = 0;
  auto CTX = (poro_viscoplasticity_ctx *) ctx;
  const double dt = CTX->dt;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *param     = (m->param)->model_param;
  int    *param_idx = (m->param)->model_param_index;  

  KMS_IJSS2017_Parameters *mat_pvp = (m->param)->cm_mat_2->mat_pvp;
  
  double dMdF_in[DIM_3x3x3x3];
  Matrix(double) dMdF;
  
  err += pvp_intf_compute_dMdF(dMdF_in,
                               Fs[TENSOR_Fnp1].m_pdata,
                               Fs[TENSOR_Fn].m_pdata,
                               Fs[TENSOR_pFnp1].m_pdata,
                               Fs[TENSOR_pFn].m_pdata,
                               vars[VAR_pc_np1],
                               vars[VAR_pc_n],
                               mat_pvp,
                               dt);

  dMdF.m_pdata = dMdF_in;
  dMdF.m_row = DIM_3x3x3x3;
  dMdF.m_col = 1;

  Matrix(double) dMdu_ab, Grad_op_ab;
     dMdu_ab.m_row =    dMdu_ab.m_col = DIM_3;
  Grad_op_ab.m_row = Grad_op_ab.m_col = DIM_3;
  

    
  for (int a = 0; a<nne; a++)
  {
    for(int b = 0; b<nsd; b++)
    {
      int idx_ab = idx_4_gen(a,b,0,0,nne,nsd,DIM_3,DIM_3);
         dMdu_ab.m_pdata = dM_du + idx_ab;
      Grad_op_ab.m_pdata = Grad_op + idx_ab;
      
      for(int w = 1; w<=DIM_3; w++)
      {
        for(int x = 1; x<=DIM_3; x++)
        {
          Mat_v(dMdu_ab,w,x) = 0.0;
/*          for(int y = 1; y<=DIM_3; y++)
          {
            for(int z = 1; z<=DIM_3; z++)
              Mat_v(dMdu_ab,w,x) += Tns4_v(dMdF,w,x,y,z)*Mat_v(Grad_op_ab,y,z);
          }
*/          
        }
      }
    }
  }        
  return err;
}

/* THIS IS A FUNCTION STUB. */
static int cm_pvp_set_init_vals(Constitutive_model *m)
{
  // inital values are set in the more convoluted
  // read_constitutive_model_parameters->plasticity_model_read_parameters
  //   calling sequence
  
  Matrix(double) *Fs = (m->vars_list[0][m->model_id]).Fs;
  double *param      = (m->param)->model_param;
  double *vars       = (m->vars_list[0][m->model_id]).state_vars[0].m_pdata;

  KMS_IJSS2017_Parameters *mat_pvp = (m->param)->cm_mat_2->mat_pvp;

  double pJ   = param[PARAM_pJ];
  
  double pF11 = pow(pJ, 1.0/3.0);
  double F0[9] = {0.0,0.0,0.0,
                  0.0,0.0,0.0,
                  0.0,0.0,0.0};
  F0[0] = F0[4] = F0[8] =  pF11;

  vars[VAR_pc_nm1] = vars[VAR_pc_n] = vars[VAR_pc_np1] = param[PARAM_pc_0];
  if((m->param)->pF != NULL)
  {
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

/* THIS IS A FUNCTION STUB. */
static int cm_pvp_read(Model_parameters *p,
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
  if(match != (PARAM_NO + PARAM_INX_NO))
  {   
    err++;
    assert(match == (PARAM_NO + PARAM_INX_NO) && "Did not read expected number of parameters"); 
  }
  err += scan_for_valid_line(in);
  
  /* not expecting EOF, check and return error if encountered */
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");

  bool usingSmoothMacauleyBrackets = true;  
  p->cm_mat_2->mat_pvp = new KMS_IJSS2017_Parameters;
  (p->cm_mat_2->mat_pvp)->set_parameters(param[PARAM_yf_M],
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
                                         param[PARAM_cf_pcinf],
                                         usingSmoothMacauleyBrackets);
  
  return err;
}

static size_t cm_pvp_get_size(const Constitutive_model *m)
{
  return state_variables_get_packed_size(m->vars_list[0]+m->model_id);
}

static int cm_pvp_pack(const Constitutive_model *m,
                    char *buffer,
                    size_t *pos)
{
  return state_variables_pack(m->vars_list[0]+m->model_id, buffer, pos);
}

static int cm_pvp_unpack(Constitutive_model *m,
                     const char *buffer,
                     size_t *pos)
{
  return state_variables_unpack(m->vars_list[0]+m->model_id, buffer, pos);
}
                                        
int poro_viscoplasticity_model_initialize(Model_parameters *p)
{
  int err = 0;

  /* set functions */
  p->integration_algorithm = cm_pvp_int_alg;
  p->compute_dev_stress    = cm_pvp_dev_stress;
  p->compute_dudj          = cm_pvp_dudj;
  p->compute_dev_tangent   = cm_pvp_dev_tangent;
  p->update_elasticity     = cm_pvp_update_elasticity;
  p->compute_d2udj2        = cm_pvp_d2udj2;
  p->update_state_vars     = cm_pvp_update;
  p->reset_state_vars      = cm_pvp_reset;
  p->get_subdiv_param      = cm_pvp_get_subdiv_param;
  p->get_var_info          = cm_pvp_info;
  p->get_F                 = cm_pvp_get_F;
  p->get_Fn                = cm_pvp_get_Fn;
  p->get_Fnm1              = cm_pvp_get_Fnm1;  
  p->get_pF                = cm_pvp_get_pF;
  p->get_pFn               = cm_pvp_get_pFn;
  p->get_pFnm1             = cm_pvp_get_pFnm1;    
  p->get_eF                = cm_pvp_get_eF;
  p->get_eFn               = cm_pvp_get_eFn;
  p->get_eFnm1             = cm_pvp_get_eFnm1;
  p->get_eF_of_hF          = cm_pvp_get_eF_with_thermal;

  p->reset_state_vars_using_temporal   = cm_pvp_reset_using_temporal;
  p->update_np1_state_vars_to_temporal = cm_pvp_update_np1_to_temporal;
  p->save_state_vars_to_temporal       = cm_pvp_save_to_temporal;
    
  p->get_hardening         = cm_pvp_get_hardening;
  p->get_hardening_nm1     = cm_pvp_get_hardening_nm1;  
  p->get_plast_strain_var  = cm_get_lam_p;
  p->write_restart         = cm_pvp_write_restart;
  p->read_restart          = cm_pvp_read_restart;
  p->destroy_ctx           = cm_pvp_ctx_destroy;
  p->compute_dMdu          = cm_pvp_compute_dMdu;
  p->set_init_vals         = cm_pvp_set_init_vals;
  p->read_param            = cm_pvp_read;
  p->get_size              = cm_pvp_get_size;
  p->pack                  = cm_pvp_pack;
  p->unpack                = cm_pvp_unpack;
  p->type                  = POROVISCO_PLASTICITY;

  p->n_param           = PARAM_NO;
  p->model_param       = (double *) calloc(PARAM_NO, sizeof(*(p->model_param)));
  p->n_param_index     = PARAM_INX_NO;
  p->model_param_index = (int *) calloc(PARAM_INX_NO, sizeof(*(p->model_param_index)));
  return err;
}

int poro_viscoplasticity_model_destroy(Model_parameters *p)
{
  int err = 0;
  delete (p->cm_mat_2)->mat_pvp;  
  return err;
}

int test_cm_poro_viscoplasticity_model(Constitutive_model *m)
{
  int err = 0;
  // simulation parameters
  double dt = 0.01;
  double strainrate = -0.005;
  int stepno = 13200;

  double *param = (m->param)->model_param;
  
  KMS_IJSS2017_Parameters *mat_pvp = (m->param)->cm_mat_2->mat_pvp;
  
  double p0 = param[PARAM_K_p0];
  double h  = pvp_intf_hardening_law(p0, mat_pvp);
  double HardLawJp0Coeff = pow(exp(h), 1.0/3.0);

  double Fnp1[9], Fn[9], pFnp1[9], pFn[9];
  double  F0[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double   I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};                  
  double Sh2[9] = {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};                

  F0[0] = F0[4] = F0[8] =  HardLawJp0Coeff;
                 
  double pc_n, pc_np1;
  pc_n = pc_np1 = p0;
  
  memcpy(pFn,   F0, sizeof(double)*DIM_3x3);
  memcpy(pFnp1, F0, sizeof(double)*DIM_3x3);
  memcpy( Fn,    I, sizeof(double)*DIM_3x3);
  memcpy( Fnp1,  I, sizeof(double)*DIM_3x3);   
  
  
  double t = 0;
  for(int iA=1; iA<=stepno; iA++)
  {
    t += dt;
    // update total deformation gradient
    double lambdac = (t > 46.0) ?strainrate * 46.0*sin(46.0*M_PI/66.0):strainrate*t*sin(t*M_PI/66.0);
    double lambdas = (t > 42.6) ? ( 42.6 - t ) * strainrate : 0.0;
                   
    for(int ia=0; ia<9; ia++)
      Fnp1[ia] = F0[ia] + lambdac*I[ia] + lambdas*Sh2[ia];    

    err += pvp_intf_perform_integration_algorithm(Fnp1,Fn,pFnp1,pFn,&pc_np1,pc_n,mat_pvp,dt);

    memcpy(pFn,pFnp1,sizeof(double)*DIM_3x3);
    memcpy( Fn, Fnp1,sizeof(double)*DIM_3x3);
    pc_n = pc_np1;
  }
    
  return err;
}
