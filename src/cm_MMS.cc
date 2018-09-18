/// Define poro-visco plasticity model functions for the constitutive model interface
/// 
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#include "cm_MMS.h"

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
#include"poro_visco_plasticity.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX_D_ALPHA 0.0002

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static constexpr double MMS4cm_c4 = 0.7;

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


static const int VAR_end = 0;
static const int FLAG_no = 0;
static const int PARAM_NO = 0;
static const int PARAM_INX_NO = 0;

enum tensor_names {
  TENSOR_Fn,
  TENSOR_pFn,
  TENSOR_Fnp1,
  TENSOR_pFnp1,
  TENSOR_Fnm1,
  TENSOR_pFnm1,
  TENSOR_end
};

/// Private structure for use exclusively with this model and
// associated functions.
typedef struct cm_mms_ctx {
  double *F;
  double t;
  double x;
  double y;
  double z;
  double alpha; // mid point alpha
  double *eFnpa;
  int is_coulpled_with_thermal;
  int npa;  
} cm_mms_ctx;

/// Construct and initialize the poro-viscoplasticity model context 
/// for calling functions through the constitutive modeling interface
/// 
/// \param[in,out] ctx - handle to an opaque model context object.
/// \param[in] F The total deformation gradient.
/// \param[in] alpha mid-point alpha
/// \param[in] eFnpa elastic deformation gradient at t = n + alpha
/// \param[in] npa   mid-point rule
/// \param[in] x     coordinate. x = x[0], y = x[1], z = x[2]
/// \param[in] t     time
/// \return non-zero on internal error.
int cm_mms_ctx_build(void **ctx,
                     double *F,
                     const double alpha,
                     double *eFnpa,
                     const int npa,
                     const double *x,
                     const double t)
{
  int err = 0;
  cm_mms_ctx *t_ctx = (cm_mms_ctx *) malloc(sizeof(cm_mms_ctx));

  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;

  t_ctx->F = F;
  t_ctx->eFnpa = eFnpa;
  t_ctx->npa   = npa;  
  
  t_ctx->t  = t;
  t_ctx->x  = x[0];
  t_ctx->y  = x[1];
  t_ctx->z  = x[2];
  t_ctx->alpha = alpha;

  /* assign handle */
  *ctx = t_ctx;
  return err;
}            

void MMS4cm_pF(double *pF, double t, double X, double Y, double Z);

int CM_MMS_PARAM::integration_algorithm(Constitutive_model *m,
                                        const void *ctx)
const
{
  int err = 0;
  auto CTX = (cm_mms_ctx *) ctx;

  const double t = CTX->t;
  const double x = CTX->x;
  const double y = CTX->y;
  const double z = CTX->z;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  
  MMS4cm_pF(Fs[TENSOR_pFnp1].m_pdata, t, x, y, z);
  
  memcpy(Fs[TENSOR_Fnp1].m_pdata, CTX->F, DIM_3x3*sizeof(double));
  return err;
}

int CM_MMS_PARAM::update_elasticity(const Constitutive_model *m,
                                    const void *ctx_in,
                                    double *L_in,
                                    double *S,
                                    const int compute_stiffness)
const
{
  int err = 0;
  auto ctx = (cm_mms_ctx *) ctx_in;
  
  double *L = NULL;
  if(compute_stiffness)
    L = L_in;
  
  if(ctx->eFnpa)
    err += constitutive_model_default_update_elasticity(m, ctx->eFnpa, L, S, compute_stiffness);
  else
  {
    // shorthand of deformation gradients
    Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;  

    TensorA<2> pFnp1(Fs[TENSOR_pFnp1].m_pdata), Fnp1(Fs[TENSOR_Fnp1].m_pdata);

    Tensor<2> pFnp1_I, eF;    
    err += inv(pFnp1, pFnp1_I);    
    eF(i,j) = Fnp1(i,k)*pFnp1_I(k,j);

    err += constitutive_model_default_update_elasticity(m, eF.data, L, S, compute_stiffness);          
  }     
  return err;
}


int CM_MMS_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[TENSOR_Fnm1] = Fs[TENSOR_Fn];    
  Fs[TENSOR_Fn]   = Fs[TENSOR_Fnp1];
  Fs[TENSOR_pFnm1]= Fs[TENSOR_pFn];
  Fs[TENSOR_pFn]  = Fs[TENSOR_pFnp1];

  return err;
}

int CM_MMS_PARAM::reset_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[TENSOR_Fnp1]  = Fs[TENSOR_Fn];
  Fs[TENSOR_pFnp1] = Fs[TENSOR_pFn];  

  return err;
}

int CM_MMS_PARAM::get_var_info(Model_var_info &info)
const 
{
  const CMVariableNames variable_names[] = {};

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

int CM_MMS_PARAM::get_pF(const Constitutive_model *m,
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

int CM_MMS_PARAM::get_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const 
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_MMS_PARAM::set_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const
{
  State_variables *sv = m->vars_list[0] + m->model_id;
  return sv->get_F(F, TENSOR_Fnm1, TENSOR_Fn, TENSOR_Fnp1, stepno);
}

int CM_MMS_PARAM::get_eF(const Constitutive_model *m,
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

int CM_MMS_PARAM::get_eF_of_hF(const Constitutive_model *m,
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

int CM_MMS_PARAM::reset_state_vars_using_temporal(const Constitutive_model *m, 
                                                  State_variables *var)
const {
  int err = 0;
  Matrix<double> *Fs    = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs_in = var->Fs;

  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fn   ].m_pdata[ia] = Fs_in[TENSOR_Fn   ].m_pdata[ia];
    Fs[TENSOR_pFn  ].m_pdata[ia] = Fs_in[TENSOR_pFn  ].m_pdata[ia];
    Fs[TENSOR_Fnm1 ].m_pdata[ia] = Fs_in[TENSOR_Fnm1 ].m_pdata[ia];
    Fs[TENSOR_pFnm1].m_pdata[ia] = Fs_in[TENSOR_pFnm1].m_pdata[ia];
  }
    
  return err;
}

int CM_MMS_PARAM::update_np1_state_vars_to_temporal(const Constitutive_model *m, 
                                                    State_variables *var)
const{
  int err = 0;
  Matrix<double> *Fs    = var->Fs;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
    
  for(int ia=0; ia<DIM_3x3; ia++)
  {
    Fs[TENSOR_Fnp1 ].m_pdata[ia] = Fs_in[TENSOR_Fnp1 ].m_pdata[ia];
    Fs[TENSOR_pFnp1].m_pdata[ia] = Fs_in[TENSOR_pFnp1].m_pdata[ia];
  }

  return err;
}

int CM_MMS_PARAM::save_state_vars_to_temporal(const Constitutive_model *m, 
                                              State_variables *var)
const{
  int err = 0;
  Matrix<double> *Fs_in = m->vars_list[0][m->model_id].Fs;
  Matrix<double> *Fs    = var->Fs;
  
  for(int ia=0; ia<TENSOR_end; ia++){
    for(int ib=0; ib<DIM_3x3; ib++)
      Fs[ia].m_pdata[ib] = Fs_in[ia].m_pdata[ib];
  }
    
  return err;
}


int CM_MMS_PARAM::write_restart(FILE *fp, const Constitutive_model *m)
const
{

  int err = 0;
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;

  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  err += cm_write_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);
  
  return err;
}

int CM_MMS_PARAM::read_restart(FILE *fp, Constitutive_model *m)
const
{
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  
  cm_read_tensor_restart(fp, Fs[TENSOR_Fn].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_Fnm1].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_pFn].m_pdata);
  cm_read_tensor_restart(fp, Fs[TENSOR_pFnm1].m_pdata);

  this->reset_state_vars(m);
  return 0;  
}

int CM_MMS_PARAM::destroy_ctx(void **ctx)
const
{
  int err = 0;
  cm_mms_ctx *t_ctx = (cm_mms_ctx *) *ctx;
  /* invalidate handle */
  *ctx = NULL;

  // no memory was created
  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;

  free(t_ctx);
  return err;
}


int CM_MMS_PARAM::compute_dMdu(const Constitutive_model *m,
                               const void *ctx,
                               double *Grad_op,
                               const int nne,
                               const int nsd,
                               double *dM_du)
const
{
  memset(dM_du, 0, nne*nsd*DIM_3x3*sizeof(double));
  return 0;
}

/* THIS IS A FUNCTION STUB. */
int CM_MMS_PARAM::set_init_vals(Constitutive_model *m)
const
{
  // inital values are set in the more convoluted
  // read_constitutive_model_parameters->plasticity_model_read_parameters
  //   calling sequence
  
  Matrix<double> *Fs = (m->vars_list[0][m->model_id]).Fs;
  double F0[DIM_3x3] = {1.0,0.0,0.0,
                        0.0,1.0,0.0,
                        0.0,0.0,1.0};
  memcpy(Fs[TENSOR_Fnm1 ].m_pdata, F0, sizeof(double)*DIM_3x3);
  memcpy(Fs[TENSOR_Fn   ].m_pdata, F0, sizeof(double)*DIM_3x3);
  memcpy(Fs[TENSOR_Fnp1 ].m_pdata, F0, sizeof(double)*DIM_3x3);
  memcpy(Fs[TENSOR_pFnm1].m_pdata, F0, sizeof(double)*DIM_3x3);
  memcpy(Fs[TENSOR_pFn  ].m_pdata, F0, sizeof(double)*DIM_3x3);
  memcpy(Fs[TENSOR_pFnp1].m_pdata, F0, sizeof(double)*DIM_3x3);

  return 0;
}

int CM_MMS_PARAM::read_param(FILE *in)
const
{
  // there are no parameters to read 
  return scan_for_valid_line(in);
}

void MMS4cm_pF(double *pF, double t, double X, double Y, double Z)
{
  double c4 = MMS4cm_c4;
  double pU = -X*Y*Z*c4*(exp(t*-5.0E1)-1.0)*(X-1.0)*(Y-1.0)*(Z-1.0)+1.0;
  pF[0] = 1.0;
  pF[1] = 0.0;
  pF[2] = 0.0;
  pF[3] = 0.0;
  pF[4] = 1.0;
  pF[5] = 0.0;
  pF[6] = 0.0;
  pF[7] = 0.0;
  pF[8] = pU;
}

void MMS4cm_displacement(double *u, double t, double X, double Y, double Z)
{
  double xU1 = (X*X)*Y*Z*(Y*2.0-1.0)*pow(X-1.0,2.0)*(Y-1.0)*(Z-1.0)*-1.0E1;
  double xU2 = X*(Y*Y)*Z*(X*2.0-1.0)*(X-1.0)*pow(Y-1.0,2.0)*(Z-1.0)*1.0E1;
  double xU3 = X*Y*Z*(X-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/2.0);
  double tU = cos(M_PI*(1.0/4.0)-t*1.0E3)/(t*1.0E3+1.0)-sqrt(2.0)*(1.0/2.0);

  u[0] = tU*xU1;
  u[1] = tU*xU2;
  u[2] = tU*xU3;
}

void MMS4cm_velocity(double *v, double t, double X, double Y, double Z)
{
  double xU1 = (X*X)*Y*Z*(Y*2.0-1.0)*pow(X-1.0,2.0)*(Y-1.0)*(Z-1.0)*-1.0E1;
  double xU2 = X*(Y*Y)*Z*(X*2.0-1.0)*(X-1.0)*pow(Y-1.0,2.0)*(Z-1.0)*1.0E1;
  double xU3 = X*Y*Z*(X-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/2.0);
  double tUt = cos(M_PI*(1.0/4.0)-t*1.0E3)*1.0/pow(t*1.0E3+1.0,2.0)*-1.0E3+(sin(M_PI*(1.0/4.0)-t*1.0E3)*1.0E3)/(t*1.0E3+1.0);

  v[0] = tUt*xU1;
  v[1] = tUt*xU2;
  v[2] = tUt*xU3;
}

void MMS4cm_initial_velocity(double *v, double X, double Y, double Z)
{
  MMS4cm_velocity(v, 0.0, X, Y, Z);
}

void MMS4cm_pressure_volume(double *P, double *V, ELASTICITY *elast, double t, double X, double Y, double Z)
{
  double c4 = MMS4cm_c4;
  double xUx11 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double xUx21 = (Y*Y)*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double xUx31 = Y*Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/2.0);
  double xUx12 = (X*X)*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1;
  double xUx22 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double xUx32 = X*Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/2.0);
  double xUx13 = (X*X)*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1;
  double xUx23 = X*(Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double xUx33 = X*Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/2.0);
  double tU = cos(M_PI*(1.0/4.0)-t*1.0E3)/(t*1.0E3+1.0)-sqrt(2.0)*(1.0/2.0);
  double pU = -X*Y*Z*c4*(exp(t*-5.0E1)-1.0)*(X-1.0)*(Y-1.0)*(Z-1.0)+1.0;
  double eJ = (tU*xUx11+tU*xUx22+tU*xUx33+(tU*tU)*xUx11*xUx22-(tU*tU)*xUx12*xUx21+(tU*tU)*xUx11*xUx33-(tU*tU)*xUx13*xUx31+(tU*tU)*xUx22*xUx33-(tU*tU)*xUx23*xUx32+(tU*tU*tU)*xUx11*xUx22*xUx33-(tU*tU*tU)*xUx11*xUx23*xUx32-(tU*tU*tU)*xUx12*xUx21*xUx33+(tU*tU*tU)*xUx12*xUx23*xUx31+(tU*tU*tU)*xUx13*xUx21*xUx32-(tU*tU*tU)*xUx13*xUx22*xUx31+1.0)/pU;
  double pJ = pU;

  double J = pJ*eJ;
  double dudJ = 0.0;
  elast->compute_dudj(&dudJ, eJ);
  *P = dudJ*elast->mat->kappa;
  *V = J;
}

void MMS4cm_body_force(double *b, const HOMMAT *hommat, ELASTICITY *elast, double t, double X, double Y, double Z)
{
  double rho = hommat->density;
  double c4 = MMS4cm_c4;
  double xU1 = (X*X)*Y*Z*(Y*2.0-1.0)*pow(X-1.0,2.0)*(Y-1.0)*(Z-1.0)*-1.0E1;
  double xU2 = X*(Y*Y)*Z*(X*2.0-1.0)*(X-1.0)*pow(Y-1.0,2.0)*(Z-1.0)*1.0E1;
  double xU3 = X*Y*Z*(X-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/2.0);
  double xUx11 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double xUx21 = (Y*Y)*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double xUx31 = Y*Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/2.0);
  double xUx12 = (X*X)*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1;
  double xUx22 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double xUx32 = X*Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/2.0);
  double xUx13 = (X*X)*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1;
  double xUx23 = X*(Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double xUx33 = X*Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/2.0);
  double tU = cos(M_PI*(1.0/4.0)-t*1.0E3)/(t*1.0E3+1.0)-sqrt(2.0)*(1.0/2.0);
  double tUtt = (cos(M_PI*(1.0/4.0)-t*1.0E3)*-1.0E6)/(t*1.0E3+1.0)+cos(M_PI*(1.0/4.0)-t*1.0E3)*1.0/pow(t*1.0E3+1.0,3.0)*2.0E6-sin(M_PI*(1.0/4.0)-t*1.0E3)*1.0/pow(t*1.0E3+1.0,2.0)*2.0E6;

  double pU = -X*Y*Z*c4*(exp(t*-5.0E1)-1.0)*(X-1.0)*(Y-1.0)*(Z-1.0)+1.0;
  double pUx1 = Y*Z*c4*exp(t*-5.0E1)*(X*2.0-1.0)*(exp(t*5.0E1)-1.0)*(Y-1.0)*(Z-1.0);
  double pUx2 = X*Z*c4*exp(t*-5.0E1)*(Y*2.0-1.0)*(exp(t*5.0E1)-1.0)*(X-1.0)*(Z-1.0);
  double pUx3 = X*Y*c4*exp(t*-5.0E1)*(Z*2.0-1.0)*(exp(t*5.0E1)-1.0)*(X-1.0)*(Y-1.0);
  double pJx1 = Y*Z*c4*exp(t*-5.0E1)*(X*2.0-1.0)*(exp(t*5.0E1)-1.0)*(Y-1.0)*(Z-1.0);
  double pJx2 = X*Z*c4*exp(t*-5.0E1)*(Y*2.0-1.0)*(exp(t*5.0E1)-1.0)*(X-1.0)*(Z-1.0);
  double pJx3 = X*Y*c4*exp(t*-5.0E1)*(Z*2.0-1.0)*(exp(t*5.0E1)-1.0)*(X-1.0)*(Y-1.0);
  double Utt1 = tUtt*xU1;
  double Utt2 = tUtt*xU2;
  double Utt3 = tUtt*xU3;

  double eF[9] = {};
  eF[0] = tU*xUx11+1.0;
  eF[1] = tU*xUx12;
  eF[2] = (tU*xUx13)/pU;
  eF[3] = tU*xUx21;
  eF[4] = tU*xUx22+1.0;
  eF[5] = (tU*xUx23)/pU;
  eF[6] = tU*xUx31;
  eF[7] = tU*xUx32;
  eF[8] = (tU*xUx33+1.0)/pU;

  double *tempS = elast->S;
  double *tempL = elast->L;
  double eS[9] = {}, L[81]={};
  elast->S = eS;
  elast->L = L;
  elast->update_elasticity(elast,eF,1);
  elast->S = tempS;
  elast->L = tempL;
  double pJ = pU;
  double xUxx11 = Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-X*Y*Z*(X*4.0-3.0)*(Z-1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double xUxx21 = (Y*Y)*Z*(X*1.2E1-6.0)*pow(Y-1.0,2.0)*(Z-1.0)*1.0E1;
  double xUxx31 = -Y*Z*(Y-1.0)*(Z-1.0);
  double xUxx12 = X*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-2.0E1-(X*X)*Z*(X*2.0-2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*1.0E1;
  double xUxx22 = Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1+X*Y*Z*(X*4.0-3.0)*(Z-1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double xUxx32 = Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/2.0)-X*Z*(Y*2.0-1.0)*(Z-1.0)*(1.0/2.0);
  double xUxx13 = X*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-(X*X)*Y*(X*2.0-2.0)*(Z*2.0-1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*1.0E1;
  double xUxx23 = (Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1+X*(Y*Y)*(X*4.0-3.0)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*1.0E1;
  double xUxx33 = Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/2.0)-X*Y*(Z*2.0-1.0)*(Y-1.0)*(1.0/2.0);
  double xUxy11 = X*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-X*Y*Z*(Y*4.0-3.0)*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1;
  double xUxy21 = Y*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*2.0E1+(Y*Y)*Z*(Y*2.0-2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double xUxy31 = Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/2.0)-Y*Z*(X*2.0-1.0)*(Z-1.0)*(1.0/2.0);
  double xUxy12 = (X*X)*Z*(Y*1.2E1-6.0)*pow(X-1.0,2.0)*(Z-1.0)*-1.0E1;
  double xUxy22 = X*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1+X*Y*Z*(Y*4.0-3.0)*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1;
  double xUxy32 = -X*Z*(X-1.0)*(Z-1.0);
  double xUxy13 = (X*X)*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1-(X*X)*Y*(Y*4.0-3.0)*(Z*2.0-1.0)*pow(X-1.0,2.0)*1.0E1;
  double xUxy23 = X*Y*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1+X*(Y*Y)*(Y*2.0-2.0)*(Z*2.0-1.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double xUxy33 = X*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/2.0)-X*Y*(Z*2.0-1.0)*(X-1.0)*(1.0/2.0);
  double xUxz11 = X*Y*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-X*Y*Z*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double xUxz21 = (Y*Y)*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1+(Y*Y)*Z*pow(Y-1.0,2.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double xUxz31 = Y*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/2.0)-Y*Z*(X*2.0-1.0)*(Y-1.0)*(1.0/2.0);
  double xUxz12 = (X*X)*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1-(X*X)*Z*pow(X-1.0,2.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*1.0E1;
  double xUxz22 = X*Y*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1+X*Y*Z*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double xUxz32 = X*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/2.0)-X*Z*(Y*2.0-1.0)*(X-1.0)*(1.0/2.0);
  double xUxz13 = (X*X)*Y*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double xUxz23 = X*(Y*Y)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1;
  double xUxz33 = -X*Y*(X-1.0)*(Y-1.0);

  b[0] = Utt1*rho-(eF[0]*eS[2]*pJx3+eF[1]*eS[5]*pJx3+eF[2]*eS[8]*pJx3)/pU-(-eS[8]*pJ*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eS[2]*pJ*tU*xUxz11+eS[5]*pJ*tU*xUxz12)/pU+(eF[0]*pJ*(-L[18]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[22]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[19]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[21]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[26]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[20]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[23]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[24]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[25]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0)))+eF[1]*pJ*(-L[45]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[49]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[46]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[48]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[53]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[47]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[50]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[51]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[52]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0)))+eF[2]*pJ*(-L[72]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[76]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[73]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[75]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[80]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[74]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[77]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[78]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[79]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))))/pU-eF[0]*eS[0]*pJx1-eF[0]*eS[1]*pJx2-eF[1]*eS[3]*pJx1-eF[1]*eS[4]*pJx2-eF[2]*eS[6]*pJx1-eF[2]*eS[7]*pJx2+eS[6]*pJ*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eS[7]*pJ*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[0]*pJ*(-L[0]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[4]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[1]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[3]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[8]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[2]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[5]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[6]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[7]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[0]*pJ*(-L[9]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[13]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[10]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[12]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[17]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[11]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[14]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[15]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[16]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+eF[1]*pJ*(-L[27]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[31]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[28]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[30]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[35]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[29]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[32]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[33]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[34]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[1]*pJ*(-L[36]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[40]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[37]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[39]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[44]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[38]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[41]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[42]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[43]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+eF[2]*pJ*(-L[54]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[58]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[55]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[57]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[62]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[56]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[59]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[60]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[61]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[2]*pJ*(-L[63]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[67]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[64]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[66]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[71]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[65]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[68]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[69]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[70]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+1.0/(pU*pU)*pUx3*(eF[0]*eS[2]*pJ+eF[1]*eS[5]*pJ+eF[2]*eS[8]*pJ)-eS[0]*pJ*tU*xUxx11-eS[1]*pJ*tU*xUxy11-eS[3]*pJ*tU*xUxx12-eS[4]*pJ*tU*xUxy12;
  b[1] = Utt2*rho-(eF[3]*eS[2]*pJx3+eF[4]*eS[5]*pJx3+eF[5]*eS[8]*pJx3)/pU-(-eS[8]*pJ*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eS[2]*pJ*tU*xUxz21+eS[5]*pJ*tU*xUxz22)/pU+(eF[3]*pJ*(-L[18]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[22]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[19]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[21]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[26]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[20]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[23]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[24]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[25]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0)))+eF[4]*pJ*(-L[45]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[49]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[46]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[48]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[53]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[47]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[50]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[51]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[52]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0)))+eF[5]*pJ*(-L[72]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[76]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[73]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[75]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[80]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[74]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[77]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[78]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[79]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))))/pU-eF[3]*eS[0]*pJx1-eF[3]*eS[1]*pJx2-eF[4]*eS[3]*pJx1-eF[4]*eS[4]*pJx2-eF[5]*eS[6]*pJx1-eF[5]*eS[7]*pJx2+eS[6]*pJ*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eS[7]*pJ*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[3]*pJ*(-L[0]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[4]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[1]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[3]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[8]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[2]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[5]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[6]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[7]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[3]*pJ*(-L[9]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[13]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[10]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[12]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[17]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[11]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[14]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[15]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[16]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+eF[4]*pJ*(-L[27]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[31]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[28]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[30]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[35]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[29]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[32]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[33]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[34]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[4]*pJ*(-L[36]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[40]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[37]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[39]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[44]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[38]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[41]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[42]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[43]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+eF[5]*pJ*(-L[54]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[58]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[55]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[57]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[62]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[56]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[59]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[60]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[61]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[5]*pJ*(-L[63]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[67]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[64]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[66]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[71]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[65]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[68]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[69]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[70]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+1.0/(pU*pU)*pUx3*(eF[3]*eS[2]*pJ+eF[4]*eS[5]*pJ+eF[5]*eS[8]*pJ)-eS[0]*pJ*tU*xUxx21-eS[1]*pJ*tU*xUxy21-eS[3]*pJ*tU*xUxx22-eS[4]*pJ*tU*xUxy22;
  b[2] = Utt3*rho-(eF[6]*eS[2]*pJx3+eF[7]*eS[5]*pJx3+eF[8]*eS[8]*pJx3)/pU-(-eS[8]*pJ*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)+eS[2]*pJ*tU*xUxz31+eS[5]*pJ*tU*xUxz32)/pU+(eF[6]*pJ*(-L[18]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[22]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[19]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[21]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[26]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[20]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[23]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[24]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[25]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0)))+eF[7]*pJ*(-L[45]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[49]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[46]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[48]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[53]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[47]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[50]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[51]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[52]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0)))+eF[8]*pJ*(-L[72]*(eF[0]*tU*xUxz11+eF[3]*tU*xUxz21+eF[6]*tU*xUxz31)-L[76]*(eF[1]*tU*xUxz12+eF[4]*tU*xUxz22+eF[7]*tU*xUxz32)-L[73]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))-L[75]*(eF[0]*tU*xUxz12*(1.0/2.0)+eF[1]*tU*xUxz11*(1.0/2.0)+eF[3]*tU*xUxz22*(1.0/2.0)+eF[4]*tU*xUxz21*(1.0/2.0)+eF[6]*tU*xUxz32*(1.0/2.0)+eF[7]*tU*xUxz31*(1.0/2.0))+L[80]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU))+L[74]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[77]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))+L[78]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz11*(1.0/2.0)-eF[5]*tU*xUxz21*(1.0/2.0)-eF[8]*tU*xUxz31*(1.0/2.0))+L[79]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx3-(tU*xUxz13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx3-(tU*xUxz23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx3-(tU*xUxz33)/pU)*(1.0/2.0)-eF[2]*tU*xUxz12*(1.0/2.0)-eF[5]*tU*xUxz22*(1.0/2.0)-eF[8]*tU*xUxz32*(1.0/2.0))))/pU-eF[6]*eS[0]*pJx1-eF[6]*eS[1]*pJx2-eF[7]*eS[3]*pJx1-eF[7]*eS[4]*pJx2-eF[8]*eS[6]*pJx1-eF[8]*eS[7]*pJx2+eS[6]*pJ*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)+eS[7]*pJ*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)+eF[6]*pJ*(-L[0]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[4]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[1]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[3]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[8]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[2]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[5]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[6]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[7]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[6]*pJ*(-L[9]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[13]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[10]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[12]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[17]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[11]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[14]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[15]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[16]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+eF[7]*pJ*(-L[27]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[31]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[28]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[30]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[35]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[29]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[32]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[33]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[34]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[7]*pJ*(-L[36]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[40]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[37]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[39]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[44]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[38]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[41]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[42]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[43]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+eF[8]*pJ*(-L[54]*(eF[0]*tU*xUxx11+eF[3]*tU*xUxx21+eF[6]*tU*xUxx31)-L[58]*(eF[1]*tU*xUxx12+eF[4]*tU*xUxx22+eF[7]*tU*xUxx32)-L[55]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))-L[57]*(eF[0]*tU*xUxx12*(1.0/2.0)+eF[1]*tU*xUxx11*(1.0/2.0)+eF[3]*tU*xUxx22*(1.0/2.0)+eF[4]*tU*xUxx21*(1.0/2.0)+eF[6]*tU*xUxx32*(1.0/2.0)+eF[7]*tU*xUxx31*(1.0/2.0))+L[62]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU))+L[56]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[59]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0))+L[60]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx11*(1.0/2.0)-eF[5]*tU*xUxx21*(1.0/2.0)-eF[8]*tU*xUxx31*(1.0/2.0))+L[61]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx1-(tU*xUxx13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx1-(tU*xUxx23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx1-(tU*xUxx33)/pU)*(1.0/2.0)-eF[2]*tU*xUxx12*(1.0/2.0)-eF[5]*tU*xUxx22*(1.0/2.0)-eF[8]*tU*xUxx32*(1.0/2.0)))+eF[8]*pJ*(-L[63]*(eF[0]*tU*xUxy11+eF[3]*tU*xUxy21+eF[6]*tU*xUxy31)-L[67]*(eF[1]*tU*xUxy12+eF[4]*tU*xUxy22+eF[7]*tU*xUxy32)-L[64]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))-L[66]*(eF[0]*tU*xUxy12*(1.0/2.0)+eF[1]*tU*xUxy11*(1.0/2.0)+eF[3]*tU*xUxy22*(1.0/2.0)+eF[4]*tU*xUxy21*(1.0/2.0)+eF[6]*tU*xUxy32*(1.0/2.0)+eF[7]*tU*xUxy31*(1.0/2.0))+L[71]*(eF[2]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)+eF[5]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)+eF[8]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU))+L[65]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[68]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0))+L[69]*(eF[0]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[3]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[6]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy11*(1.0/2.0)-eF[5]*tU*xUxy21*(1.0/2.0)-eF[8]*tU*xUxy31*(1.0/2.0))+L[70]*(eF[1]*(eF[2]*1.0/(pU*pU)*pUx2-(tU*xUxy13)/pU)*(1.0/2.0)+eF[4]*(eF[5]*1.0/(pU*pU)*pUx2-(tU*xUxy23)/pU)*(1.0/2.0)+eF[7]*(eF[8]*1.0/(pU*pU)*pUx2-(tU*xUxy33)/pU)*(1.0/2.0)-eF[2]*tU*xUxy12*(1.0/2.0)-eF[5]*tU*xUxy22*(1.0/2.0)-eF[8]*tU*xUxy32*(1.0/2.0)))+1.0/(pU*pU)*pUx3*(eF[6]*eS[2]*pJ+eF[7]*eS[5]*pJ+eF[8]*eS[8]*pJ)-eS[0]*pJ*tU*xUxx31-eS[1]*pJ*tU*xUxy31-eS[3]*pJ*tU*xUxx32-eS[4]*pJ*tU*xUxy32;
}

