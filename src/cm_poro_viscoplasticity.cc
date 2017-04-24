
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
#include "ttl-tools.h"

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
} poro_viscoplasticity_ctx;


int poro_viscoplasticity_model_ctx_build(void **ctx,
                                         double *F,
                                         const double dt,
                                         const double alpha,
                                         double *eFnpa,
                                         double *hFn,
                                         double *hFnp1,
                                         const int is_coulpled_with_thermal)
{
  int err = 0;
  poro_viscoplasticity_ctx *t_ctx = malloc(sizeof(poro_viscoplasticity_ctx));

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

static int cm_pvp_int_alg(Constitutive_model *m,
                              const void *ctx)
{
  int err = 0;

  auto CTX = (poro_viscoplasticity_ctx *) ctx;
  memcpy(m->vars_list[0][m->model_id].Fs[TENSOR_Fnp1].m_pdata, CTX->F, DIM_3x3 * sizeof(*(CTX->F)));

  const double dt = CTX->dt;
  Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs;
  double *vars       = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  double *param     = (m->param)->model_param;
  int    *param_idx = (m->param)->model_param_index;  

  KMS_IJSS2017_Parameters *mat_pvp = (m->param)->cm_mat_2->mat_pvp;
  
  // data for KMS_IJSS2017 definition

  unsigned PrintEveryNSteps= 100;

  std::string pathstr = "./out/KMS_IJSS2017";
  std::string logstr = pathstr + ".implicit.cst.smb.log";
  std::string Fstr = pathstr + ".implicit.cst.smb.F.txt";
  std::string pFstr = pathstr + ".implicit.cst.smb.Fp.txt";
  std::string Sstr = pathstr + ".implicit.cst.smb.S.txt";
  std::string KSstr = pathstr + ".implicit.cst.smb.KS.txt";
  std::string sigmastr = pathstr + ".implicit.cst.smb.sigma.txt";
  std::string pcstr = pathstr + ".implicit.cst.smb.pc.txt";
  std::string matstr = pathstr + ".implicit.cst.smb.mat.txt";
    
  KMS_IJSS2017_IO IO(logstr, Fstr, pFstr, Sstr, KSstr, sigmastr, pcstr, matstr, PrintEveryNSteps);
    
  
  // Time integrator
  double InitialTime = 0.0; 
  double FinalTime   = dt; 
  double CurrentTime = 0.0;
  size_t TimeStep = 0;
  
  TimeIntegrationDataManager Time_intg( InitialTime, FinalTime, dt, CurrentTime, TimeStep);
  
  

  // Parameters and Model definition
  bool usingSmoothMacauleyBrackets = true;
  bool Verbose = false;
     
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, &IO, &Time_intg );

  ttl::Tensor<2, DIM_3, double*>  Fn(  Fs[TENSOR_Fn   ].m_pdata);
  ttl::Tensor<2, DIM_3, double*>  Fnp1(Fs[TENSOR_Fnp1 ].m_pdata);
  ttl::Tensor<2, DIM_3, double*> pFn(  Fs[TENSOR_pFn  ].m_pdata);
  ttl::Tensor<2, DIM_3, double*> pFnp1(Fs[TENSOR_pFnp1].m_pdata);
    
  ImplicitModelInstance.Set_Fs_and_state_variables(&Fnp1,&Fn,&pFnp1,&pFn,vars[VAR_pc_np1],vars[VAR_pc_n]);
  ImplicitModelInstance.StepUpdate(Fnp1,dt, Verbose );            
  
  Matrix_print(Fs[TENSOR_pFnp1]);
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
                                       Matrix_double *L,
                                       Matrix_double *S,
                                       const int compute_stiffness)
{
  int err = 0;
  auto ctx = (poro_viscoplasticity_ctx *) ctx_in;
  
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
      Matrix_Tns2_AxBxC(eF,1.0,0.0,Fs[TENSOR_Fnp1],hFnp1_I,pFnp1_I);
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
  state_var[VAR_pc_n] = state_var[VAR_pc_np1];  
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
  
  sprintf((*info)->var_names[VAR_pc_n],   "L_n");
  sprintf((*info)->var_names[VAR_pc_np1], "L_np1");
  sprintf((*info)->var_names[VAR_pc_nm1], "L_nm1");

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
  *var = 0.0;
//  *var = m->vars_list[0][m->model_id].state_vars->m_pdata[VAR_g_n];
  return err;
}
static int cm_pvp_get_hardening_nm1(const Constitutive_model *m,
                                    double *var)
{
  int err = 0;
  *var = 0.0;
//*var = m->vars_list[0][m->model_id].state_vars->m_pdata[VAR_g_nm1];
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
                                   const int ndofn,
                                   double *dM_du)
{
  int err = 0;
  auto CTX = (poro_viscoplasticity_ctx *) ctx;
  for(int ia=0; ia<nne*ndofn*DIM_3x3; ia++)
    dM_du[ia] = 0.0;
/*  if(CTX->alpha<0)
    err += plasticity_compute_dMdu_np1(m,ctx,Grad_op,nne,ndofn,dM_du);
  else
    err += plasticity_compute_dMdu_npa(m,ctx,Grad_op,nne,ndofn,dM_du,CTX->alpha);*/
  return err;
}

/* THIS IS A FUNCTION STUB. */
static int cm_pvp_set_init_vals(Constitutive_model *m)
{
  /* inital values are set in the more convoluted
     read_constitutive_model_parameters->plasticity_model_read_parameters
     calling sequence
  */
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
  for(int ia=0; ia<PARAM_NO; ia++)
   match += fscanf(in, "%lf", param + ia);

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