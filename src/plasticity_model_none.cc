/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 *  Sangmin Lee, University of Notre Dame, Notre Dame, IN, <slee43@nd.edu> 
 */

#include "plasticity_model_none.h"
#include "allocation.h"
#include "constitutive_model.h"
#include "cm_placeholder_functions.h"
#include "new_potentials.h"
#include <ttl/ttl.h>
#include "utils.h"

#include "crystal_plasticity_integration.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

enum tensor_names {
  TENSOR_Fnm1, 
  TENSOR_Fn,
  TENSOR_Fnp1,
  TENSOR_end
};

static const int VAR_no = 0;
static const int FLAG_no = 0;

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

/// Private structure for use exclusively with this model and
/// associated functions.
typedef struct none_ctx {
  double *F;
  double *eFnpa;
  int is_coulpled_with_thermal;
  double *hFn;
  double *hFnp1;
} none_ctx;

static void he_compute_C(double *C_in,
                         double *F_in)
{
  TensorA<2> C(C_in), F(F_in);
  C = F(k,i)*F(k,j);
}

int HE_PARAM::integration_algorithm(Constitutive_model *m,
                                    const void *ctx)
const
{
  int err = 0;
  ++perIter_ODE_EXA_metric;      //accumulate the EXA metric for hyperelastic model (1 per integration point)
  
  auto CTX = (none_ctx *) ctx;
  memcpy(m->vars_list[0][m->model_id].Fs[TENSOR_Fnp1].m_pdata, CTX->F, DIM_3x3 * sizeof(*CTX->F));
  return err;
}

int HE_PARAM::compute_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 double *stress)
const
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  devStressFuncPtr Stress = getDevStressFunc(-1,m->param->p_hmat);
  double C[DIM_3x3] = {};  
  he_compute_C(C,CTX->F);
  Stress(C,m->param->p_hmat,stress);
  return err;
}

int HE_PARAM::compute_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *dudj)
const
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
  TensorA<2> F(CTX->F);
  const double J = ttl::det(F);
  Pressure(J,m->param->p_hmat,dudj);
  return err;
}

int HE_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  double *L)
const
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,m->param->p_hmat);
  double C[DIM_3x3] = {};
  he_compute_C(C,CTX->F);
  Tangent(C,m->param->p_hmat,L);
  return err;
}

int HE_PARAM::compute_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *d2udj2)
const
{
  int err = 0;
  auto CTX = (none_ctx *) ctx;
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
  TensorA<2> F(CTX->F);
  const double J = ttl::det(F);
  D_Pressure(J,m->param->p_hmat,d2udj2);
  return err;
}

int HE_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[TENSOR_Fnm1] = Fs[TENSOR_Fn];
  Fs[TENSOR_Fn]   = Fs[TENSOR_Fnp1];
  return err;
}

int HE_PARAM::reset_state_vars(Constitutive_model *m)
const 
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[TENSOR_Fnp1] = Fs[TENSOR_Fn];
  return err;
}


int HE_PARAM::reset_state_vars_using_temporal(const Constitutive_model *m, 
                                              State_variables *var)
const
{
  int err = 0;
  Matrix<double> *Fs    = (m->vars_list[0][m->model_id]).Fs;
  Matrix<double> *Fs_in = var->Fs;
  Fs[TENSOR_Fn]   = Fs_in[TENSOR_Fn];
  Fs[TENSOR_Fnm1] = Fs_in[TENSOR_Fnm1];

  return err;
}

int HE_PARAM::update_np1_state_vars_to_temporal(const Constitutive_model *m, 
                                                State_variables *var)
const{
  int err = 0;
  Matrix<double> *Fs    = var->Fs;
  Matrix<double> *Fs_in = (m->vars_list[0][m->model_id]).Fs;
  Fs[TENSOR_Fnp1] = Fs_in[TENSOR_Fnp1];

  return err;
}

int HE_PARAM::save_state_vars_to_temporal(const Constitutive_model *m, 
                                          State_variables *var)
const {
  int err = 0;
  Matrix<double> *Fs_in = (m->vars_list[0][m->model_id]).Fs;
  Matrix<double> *Fs    = var->Fs;
  Fs[TENSOR_Fn]   = Fs_in[TENSOR_Fn];
  Fs[TENSOR_Fnm1] = Fs_in[TENSOR_Fnm1];
  
  return err;
}

int HE_PARAM::get_var_info(Model_var_info &info)
const
{
  int err = 0;
  int Fno = TENSOR_end;
  
  info.n_Fs = Fno;
  info.F_names = (char **)malloc(sizeof(char*)*Fno);
  for(int a=0; a<Fno; a++)
    info.F_names[a] = (char *)malloc(sizeof(char)*1024);

  info.var_names =  NULL;
  info.flag_names = NULL;

  // allocate/copy strings 
  sprintf(info.F_names[TENSOR_Fnm1], "Fnm1");
  sprintf(info.F_names[TENSOR_Fn],   "Fn");
  sprintf(info.F_names[TENSOR_Fnp1], "Fnp1");
  
  return err;
}

int HE_PARAM::get_eF_of_hF(const Constitutive_model *m,
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

int HE_get_F(const Constitutive_model *m,
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

int HE_PARAM::get_F(const Constitutive_model *m,
                    double *F,
                    const int stepno)
const
{
  return HE_get_F(m,F,stepno);
}

int HE_PARAM::get_eF(const Constitutive_model *m,
                     double *F,
                     const int stepno)
const
{
  return HE_get_F(m,F,stepno);
}

int HE_PARAM::get_pF(const Constitutive_model *m,
                     double *F,
                     const int stepno)
const 
{
  F[0] = F[4] = F[8] = 1.0;
  F[1] = F[2] = F[3] = F[5] = F[6] = F[7] = 0.0;
  
  return 0;
}

int HE_PARAM::read_param(FILE *in)
const
{
  // there are no parameters to read 
  return scan_for_valid_line(in);
}

int HE_PARAM::write_restart(FILE *out,
                            const Constitutive_model *m)
const
{
  // write Fn to file 
  int err = 0;
  const double *F = m->vars_list[0][m->model_id].Fs[TENSOR_Fn].m_pdata;
  if (fprintf(out,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
              F[0], F[1], F[2],
              F[3], F[4], F[5],
              F[6], F[7], F[8]) < 0) err ++;
  return err;
}

int HE_PARAM::read_restart(FILE *in,
                           Constitutive_model *m)
const
{
  // read Fn from file and set Fnp1 = Fn 
  int err = 0;
  double *FN = m->vars_list[0][m->model_id].Fs[TENSOR_Fn].m_pdata;
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            FN, FN + 1, FN + 2,
            FN + 3, FN + 4, FN + 5,
            FN + 6, FN + 7, FN + 8) != DIM_3x3) err++;
  err += this->reset_state_vars(m);
  return err;
}

int HE_PARAM::update_elasticity(const Constitutive_model *m,
                                const void *ctx_in,
                                double *L,
                                double *S,
                                const int compute_stiffness)
const
{
  int err = 0;
  auto ctx = (none_ctx *) ctx_in;

  // if transient cases,
  // get_eF is not working because eF needs to be updated using mid-point alpha
  // below checks whether to use get_eF or give eFnpa in ctx

  if(ctx->eFnpa)
    err += constitutive_model_default_update_elasticity(m, ctx->eFnpa, L, S, compute_stiffness);  
  else
  {
  	Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs; 
    
    if(ctx->is_coulpled_with_thermal)
    {
      Tensor<2> hFnp1_I, eF;    
      TensorA<2> Fnp1(Fs[TENSOR_Fnp1].m_pdata), hFnp1(ctx->hFnp1);  
      err += inv(hFnp1, hFnp1_I);
      eF = Fnp1(i,k)*hFnp1_I(k,j);

      err += constitutive_model_default_update_elasticity(m, eF.data, L, S, compute_stiffness);     
    }
    else
      err += constitutive_model_default_update_elasticity(m, Fs[TENSOR_Fnp1].m_pdata, L, S, compute_stiffness);
  }

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
  none_ctx *t_ctx = (none_ctx *) malloc(sizeof(none_ctx));

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

int HE_PARAM::destroy_ctx(void **ctx)
const
{
  int err = 0;
  none_ctx *t_ctx = (none_ctx *)*ctx;
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
