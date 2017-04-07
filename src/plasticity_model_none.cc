/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */

#include "plasticity_model_none.h"
#include "constitutive_model.h"
#include "cm_placeholder_functions.h"
#include "new_potentials.h"
#include "data_structure_c.h"

Define_Matrix(double);

#define tensor 9
#define dim 3

static const int g_n_Fs = 3;
static const int g_n_vars = 0;
static const int g_n_flags = 0;
enum {Fnm1, Fn, Fnp1, F};

/**
 * Private structure for use exclusively with this model and
 * associated functions.
 */
typedef struct none_ctx {
  double *F;
  double *eFnpa;
  int is_coulpled_with_thermal;
  double *hFn;
  double *hFnp1;   
} none_ctx;

static size_t he_get_size(const Constitutive_model *m)
{
  return ((g_n_Fs * tensor + g_n_vars) * sizeof(double)
          + g_n_flags * sizeof(int));
}

static int he_pack(const Constitutive_model *m,
                   char *buffer,
                   size_t *pos)
{
  /* pack/unpack Fs */
  const Matrix_double *Fs = m->vars_list[0][m->model_id].Fs;
  pack_data(Fs[Fn].m_pdata, buffer, pos, tensor, sizeof(double));
  pack_data(Fs[Fnp1].m_pdata, buffer, pos, tensor, sizeof(double));
  return 0;
}

static int he_unpack(Constitutive_model *m,
                     const char *buffer,
                     size_t *pos)
{
  Matrix_double *Fs = m->vars_list[0][m->model_id].Fs;
  unpack_data(buffer, Fs[Fn].m_pdata, pos, tensor, sizeof(double));
  unpack_data(buffer, Fs[Fnp1].m_pdata, pos, tensor, sizeof(double));
  return 0;
}

static void he_compute_C(double * restrict C,
                         const double * restrict F)
{
  memset(C, 0, tensor * sizeof(*C));
  for (int i = 0; i < dim; i++){
    for (int j = 0; j < dim; j++){
      for (int k = 0; k < dim; k++){
        C[idx_2(i,j)] += F[idx_2(k,i)] * F[idx_2(k,j)];
      }
    }
  }
}

static double compute_bulk_mod(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static int plasticity_none_int_alg(Constitutive_model *m,
                                   const void *ctx)
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  memcpy(m->vars_list[0][m->model_id].Fs[Fnp1].m_pdata, CTX->F, tensor * sizeof(*CTX->F));
  return err;
}

static int plasticity_none_dev_stress(const Constitutive_model *m,
                                      const void *ctx,
                                      Matrix_double *stress)
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  devStressFuncPtr Stress = getDevStressFunc(-1,m->param->p_hmat);
  double C[tensor] = {};  
  he_compute_C(C,CTX->F);
  Stress(C,m->param->p_hmat,stress->m_pdata);
  return err;
}

static int plasticity_none_dudj(const Constitutive_model *m,
                                const void *ctx,
                                double *dudj)
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
  const double J = det3x3(CTX->F);
  Pressure(J,m->param->p_hmat,dudj);
  return err;
}

static int plasticity_none_dev_tangent(const Constitutive_model *m,
                                       const void *ctx,
                                       Matrix_double *tangent)
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,m->param->p_hmat);
  double C[tensor] = {};
  he_compute_C(C,CTX->F);
  Tangent(C,m->param->p_hmat,tangent->m_pdata);
  return err;
}

static int plasticity_none_d2udj2(const Constitutive_model *m,
                                  const void *ctx,
                                  double *d2udj2)
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
  const double J = det3x3(CTX->F);
  D_Pressure(J,m->param->p_hmat,d2udj2);
  return err;
}

static int plasticity_none_update(Constitutive_model *m)
{
  int err = 0;
  Matrix_AeqB(m->vars_list[0][m->model_id].Fs[Fnm1], 1.0, m->vars_list[0][m->model_id].Fs[Fn]);
  Matrix_AeqB(m->vars_list[0][m->model_id].Fs[Fn],   1.0, m->vars_list[0][m->model_id].Fs[Fnp1]);
  return err;
}

static int plasticity_none_reset(Constitutive_model *m)
{
  int err = 0;
  Matrix_AeqB(m->vars_list[0][m->model_id].Fs[Fnp1], 1.0, m->vars_list[0][m->model_id].Fs[Fn]);
  return err;
}


static int plasticity_none_reset_using_temporal(const Constitutive_model *m, State_variables *var)
{
  int err = 0;
  Matrix(double) *Fs    = (m->vars_list[0][m->model_id]).Fs;
  Matrix(double) *Fs_in = var->Fs;
  Matrix_AeqB(Fs[Fn],    1.0,Fs_in[Fn]);
  Matrix_AeqB(Fs[Fnm1],  1.0,Fs_in[Fnm1]);

  return err;
}

static int plasticity_none_update_np1_to_temporal(const Constitutive_model *m, State_variables *var)
{
  int err = 0;
  Matrix(double) *Fs    = var->Fs;
  Matrix(double) *Fs_in = (m->vars_list[0][m->model_id]).Fs;
  Matrix_AeqB(Fs[Fnp1], 1.0,Fs_in[Fnp1]);

  return err;
}

static int plasticity_none_save_to_temporal(const Constitutive_model *m, State_variables *var)
{
  int err = 0;
  Matrix(double) *Fs_in = (m->vars_list[0][m->model_id]).Fs;
  Matrix(double) *Fs    = var->Fs;
  Matrix_AeqB(Fs[Fn],    1.0,Fs_in[Fn]);
  Matrix_AeqB(Fs[Fnm1],  1.0,Fs_in[Fnm1]);
  
  return err;
}

static int plasticity_none_info(Model_var_info **info)
{
  int err = 0;

  /* make sure I don't leak memory */
  if( *info != NULL) err += model_var_info_destroy(info);

  /* allocate pointers */
  (*info) = malloc(sizeof(**info));
  (*info)->n_Fs = g_n_Fs;
  (*info)->n_vars = g_n_vars;
  (*info)->n_flags = g_n_flags;
  (*info)->F_names = malloc(g_n_Fs * sizeof( ((*info)->F_names) ));
  (*info)->var_names = malloc( g_n_vars * sizeof( ((*info)->var_names) ));
  (*info)->flag_names = malloc( g_n_flags * sizeof( ((*info)->flag_names) ));

  /* allocate/copy strings */
  (*info)->F_names[Fnm1] = strdup("Fnm1");  
  (*info)->F_names[Fn]   = strdup("Fn");
  (*info)->F_names[Fnp1] = strdup("F");
  
  return err;
}

static int he_get_eF_with_thermal(const Constitutive_model *m,
                                  Matrix(double) *eF,
                                  const Matrix(double) *hFI,
                                  const int stepno)
{
  int err = 0;
  
  switch(stepno)
  {
    case 0: // n-1
      Matrix_AxB(*eF,1.0,0.0,m->vars_list[0][m->model_id].Fs[Fnm1],0,*hFI,0);
      break;
    case 1: // n
      Matrix_AxB(*eF,1.0,0.0,m->vars_list[0][m->model_id].Fs[Fn],0,*hFI,0);
      break;
    case 2: // n+1
      Matrix_AxB(*eF,1.0,0.0,m->vars_list[0][m->model_id].Fs[Fnp1],0,*hFI,0);
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);

  return err;      
}

static int he_get_F(const Constitutive_model *m,
                    Matrix_double *F)
{
  int err = 0;
  Matrix_AeqB(*F, 1.0, m->vars_list[0][m->model_id].Fs[Fnp1]);
  return err;
}

static int he_get_Fn(const Constitutive_model *m,
                     Matrix_double *F)
{
  int err = 0;
  Matrix_AeqB(*F, 1.0, m->vars_list[0][m->model_id].Fs[Fn]);
  return err;
}

static int he_get_Fnm1(const Constitutive_model *m,
                     Matrix_double *F)
{
  int err = 0;
  Matrix_AeqB(*F, 1.0, m->vars_list[0][m->model_id].Fs[Fnm1]);
  return err;
}

static int he_get_eye(const Constitutive_model *m,
                     Matrix_double *F)
{
  int err = 0;
  Matrix_eye(*F,3);
  return err;
}

static int he_get_eF(const Constitutive_model *m,
                     Matrix_double *F)
{
  int err = 0;
  Matrix_AeqB(*F, 1.0, m->vars_list[0][m->model_id].Fs[Fnp1]);
  return err;
}

static int he_get_eFn(const Constitutive_model *m,
                      Matrix_double *F)
{
  int err = 0;
  err += he_get_Fn(m,F);
  return err;
}

static int he_get_eFnm1(const Constitutive_model *m,
                      Matrix_double *F)
{
  int err = 0;
  err += he_get_Fnm1(m,F);
  return err;
}

static int he_compute_dMdu(const Constitutive_model *m,
                           const void *ctx,
                           const double *Grad_op,
                           const int nne,
                           const int ndofn,
                           double *dM_du)
{
  int err = 0;
  /* there is no plastic deformation in this formulation, return zeros
     in dM_du */
  memset(dM_du, 0, nne * ndofn * tensor * sizeof(*dM_du));
  return err;
}

static int he_read(Model_parameters *p,
                   FILE *in)
{
  /* there are no parameters to read */
  return scan_for_valid_line(in);
}

static int he_set_initial_vals(Constitutive_model *m)
{
  /* do nothing */
  return 0;
}

static int he_write_restart(FILE *out,
                            const Constitutive_model *m)
{
  /* write Fn to file */
  int err = 0;
  const double *F = m->vars_list[0][m->model_id].Fs[Fn].m_pdata;
  if (fprintf(out,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
              F[0], F[1], F[2],
              F[3], F[4], F[5],
              F[6], F[7], F[8]) < 0) err ++;
  return err;
}

static int he_read_restart(FILE *in,
                           Constitutive_model *m)
{
  /* read Fn from file and set Fnp1 = Fn */
  int err = 0;
  double *FN = m->vars_list[0][m->model_id].Fs[Fn].m_pdata;
  double *FNP1 = m->vars_list[0][m->model_id].Fs[Fnp1].m_pdata;
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            FN, FN + 1, FN + 2,
            FN + 3, FN + 4, FN + 5,
            FN + 6, FN + 7, FN + 8) != tensor) err++;
  err += plasticity_none_reset(m);
  return err;
}

int plasticity_model_none_elasticity(const Constitutive_model *m,
                                     const void *ctx_in,
                                     Matrix_double *L,
                                     Matrix_double *S,
                                     const int compute_stiffness)
{
  int err = 0;
  auto ctx = (none_ctx *) ctx_in;
  
  // if transient cases, 
  // get_eF is not working because eF needs to be updated using mid-point alpha
  // below checks whether to use get_eF or give eFnpa in ctx

  if(ctx->eFnpa)
  {
    Matrix(double) eF;
    eF.m_row = eF.m_col = dim; eF.m_pdata = ctx->eFnpa;
    err += constitutive_model_default_update_elasticity(m, &eF, L, S, compute_stiffness);
  }
  else
  {
  	Matrix(double) *Fs = m->vars_list[0][m->model_id].Fs; 

    if(ctx->is_coulpled_with_thermal)
    {
      Matrix(double) hFnp1, hFnp1_I, eF;    
      Matrix_construct_redim(double,eF,dim,dim);
          
      hFnp1.m_row = hFnp1.m_col = dim; hFnp1.m_pdata = ctx->hFnp1;
      Matrix_construct_redim(double, hFnp1_I, dim, dim);

      err += inv3x3(ctx->hFnp1,hFnp1_I.m_pdata);
      Matrix_AxB(eF, 1.0,0.0,Fs[Fnp1],0,hFnp1_I,0);
      err += constitutive_model_default_update_elasticity(m, &eF, L, S, compute_stiffness);
      Matrix_cleanup(hFnp1_I);
      Matrix_cleanup(eF);      
    }
    else
    	err += constitutive_model_default_update_elasticity(m, Fs+Fnp1, L, S, compute_stiffness);
  }
      
  return err;
}

int plasticity_model_none_initialize(Model_parameters *p)
{
  int err = 0;

  /* set functions */
  p->integration_algorithm = plasticity_none_int_alg;
  p->compute_dev_stress = plasticity_none_dev_stress;
  p->compute_dudj = plasticity_none_dudj;
  p->compute_dev_tangent = plasticity_none_dev_tangent;
  p->update_elasticity = plasticity_model_none_elasticity;
  p->compute_d2udj2 = plasticity_none_d2udj2;
  p->update_state_vars = plasticity_none_update;
  p->reset_state_vars = plasticity_none_reset;
  p->reset_state_vars_using_temporal = plasticity_none_reset_using_temporal;
  p->update_np1_state_vars_to_temporal = plasticity_none_update_np1_to_temporal;
  p->save_state_vars_to_temporal = plasticity_none_save_to_temporal;  
  p->get_var_info = plasticity_none_info;
  p->get_F     = he_get_F;  
  p->get_Fn    = he_get_Fn;
  p->get_Fnm1  = he_get_Fnm1;  
  p->get_pF    = he_get_eye;
  p->get_pFn   = he_get_eye;
  p->get_pFnm1 = he_get_eye;      
  p->get_eF    = he_get_eF;
  p->get_eFn   = he_get_eFn;
  p->get_eFnm1 = he_get_eFnm1;
  p->get_eF_of_hF = he_get_eF_with_thermal;

  p->get_hardening = cm_get_var_zero;
  p->get_plast_strain_var = cm_get_var_zero;

  p->write_restart = he_write_restart;
  p->read_restart = he_read_restart;

  p->destroy_ctx = plasticity_model_none_ctx_destroy;
  p->compute_dMdu = he_compute_dMdu;

  p->set_init_vals = he_set_initial_vals;
  p->read_param = he_read;

  p->get_size = he_get_size;
  p->pack = he_pack;
  p->unpack = he_unpack;

  p->type = HYPER_ELASTICITY;

  p->n_param = 0;
  p->model_param = NULL;

  return err;
}

int plasticity_model_none_ctx_build(void **ctx,
                                    double *F,
                                    double *eFnpa,
                                    double *hFn,
                                    double *hFnp1,
                                    const int is_coulpled_with_thermal)
{
  int err = 0;
  none_ctx *t_ctx = malloc(sizeof(none_ctx));

  /* assign internal pointers. NOTE: We are copying the pointer NOT
     the value. No additional memory is allocated. */

  t_ctx->F     = NULL;  
  t_ctx->eFnpa = NULL;
  t_ctx->hFn   = NULL;
  t_ctx->hFnp1 = NULL;  

  t_ctx->F = F;    
  t_ctx->eFnpa = eFnpa;
  
  t_ctx->is_coulpled_with_thermal = is_coulpled_with_thermal;
  t_ctx->hFn  = hFn;
  t_ctx->hFnp1= hFnp1;  

  /* assign handle */
  *ctx = t_ctx;
  return err;
}

int plasticity_model_none_ctx_destroy(void **ctx)
{
  int err = 0;
  none_ctx *t_ctx = *ctx;
  /* invalidate handle */
  *ctx = NULL;

  /* we do not control memory for internal pointers */

  // no memory was created
  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;
  t_ctx->hFn   = NULL;
  t_ctx->hFnp1 = NULL;    
    
  free(t_ctx);
  return err;
}
