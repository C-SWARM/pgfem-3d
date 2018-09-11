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


int cm_mms_ctx_build(void **ctx,
                     double *F,
                     const double t,
                     const double x,
                     const double y,
                     const double z,
                     const double alpha,
                     double *eFnpa,
                     const int npa)
{
  int err = 0;
  cm_mms_ctx *t_ctx = (cm_mms_ctx *) malloc(sizeof(cm_mms_ctx));

  t_ctx->F     = NULL;
  t_ctx->eFnpa = NULL;

  t_ctx->F = F;
  t_ctx->eFnpa = eFnpa;
  t_ctx->npa   = npa;  
  
  t_ctx->t  = t;
  t_ctx->x  = x;
  t_ctx->y  = y;
  t_ctx->z  = z;
  t_ctx->alpha = alpha;

  /* assign handle */
  *ctx = t_ctx;
  return err;
}            

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

  double c4 = 0.7;
  double c5 = -50.0;
  Fs[TENSOR_pFnp1].m_pdata[0] = Fs[TENSOR_pFnp1].m_pdata[4] = 1.0;
  Fs[TENSOR_pFnp1].m_pdata[1] = Fs[TENSOR_pFnp1].m_pdata[2] = 0.0;
  Fs[TENSOR_pFnp1].m_pdata[3] = Fs[TENSOR_pFnp1].m_pdata[5] = 0.0;
  Fs[TENSOR_pFnp1].m_pdata[6] = Fs[TENSOR_pFnp1].m_pdata[7] = 0.0;  
  Fs[TENSOR_pFnp1].m_pdata[8]  = 1.0 - c4*(1.0-exp(c5*t))*x*(1.0-x)*y*(1.0-y)*z*(1.0-z);  
  
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

void MMS4cm_displacement(double *u, double t, double X, double Y, double Z)
{
  double exUx11 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double exUx21 = (Y*Y)*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double exUx31 = Y*Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx12 = (X*X)*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1;
  double exUx22 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUx32 = X*Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx13 = (X*X)*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1;
  double exUx23 = X*(Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double exUx33 = X*Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/1.0E2);
  double etU = cos(M_PI*(1.0/4.0)-t*1.0E3)/(t*1.0E3+1.0)-sqrt(2.0)*(1.0/2.0);
  double pU = -X*Y*Z*(X-1.0)*(Y-1.0)*(Z-1.0)*(exp(t*-5.0E1)*(7.0/1.0E1)-7.0/1.0E1)+1.0;

  u[0] = X*etU*exUx11+Y*etU*exUx12+Z*etU*exUx13*pU;
  u[1] = X*etU*exUx21+Y*etU*exUx22+Z*etU*exUx23*pU;
  u[2] = Z*(pU*(etU*exUx33+1.0)-1.0)+X*etU*exUx31+Y*etU*exUx32;
}

void MMS4cm_velocity(double *v, double t, double X, double Y, double Z)
{
  double exUx11 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double exUx21 = (Y*Y)*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double exUx31 = Y*Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx12 = (X*X)*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1;
  double exUx22 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUx32 = X*Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx13 = (X*X)*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1;
  double exUx23 = X*(Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double exUx33 = X*Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/1.0E2);
  double etU = cos(M_PI*(1.0/4.0)-t*1.0E3)/(t*1.0E3+1.0)-sqrt(2.0)*(1.0/2.0);
  double etUt = cos(M_PI*(1.0/4.0)-t*1.0E3)*1.0/pow(t*1.0E3+1.0,2.0)*-1.0E3+(sin(M_PI*(1.0/4.0)-t*1.0E3)*1.0E3)/(t*1.0E3+1.0);
  double pU = -X*Y*Z*(X-1.0)*(Y-1.0)*(Z-1.0)*(exp(t*-5.0E1)*(7.0/1.0E1)-7.0/1.0E1)+1.0;
  double pUt = X*Y*Z*exp(t*-5.0E1)*(X-1.0)*(Y-1.0)*(Z-1.0)*3.5E1;

  v[0] = Z*(etUt*exUx13*pU+etU*exUx13*pUt)+X*etUt*exUx11+Y*etUt*exUx12;
  v[1] = Z*(etUt*exUx23*pU+etU*exUx23*pUt)+X*etUt*exUx21+Y*etUt*exUx22;
  v[2] = Z*(pUt*(etU*exUx33+1.0)+etUt*exUx33*pU)+X*etUt*exUx31+Y*etUt*exUx32;
}

void MMS4cm_initial_velocity(double *v, double X, double Y, double Z)
{
  MMS4cm_velocity(v, 0.0, X, Y, Z);
}

void MMS4cm_pressure_volume(double *P, double *V, ELASTICITY *elast, double t, double X, double Y, double Z)
{
  double exUx11 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double exUx21 = (Y*Y)*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double exUx31 = Y*Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx12 = (X*X)*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1;
  double exUx22 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUx32 = X*Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx13 = (X*X)*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1;
  double exUx23 = X*(Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double exUx33 = X*Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/1.0E2);
  double etU = cos(M_PI*(1.0/4.0)-t*1.0E3)/(t*1.0E3+1.0)-sqrt(2.0)*(1.0/2.0);
  double pU = -X*Y*Z*(X-1.0)*(Y-1.0)*(Z-1.0)*(exp(t*-5.0E1)*(7.0/1.0E1)-7.0/1.0E1)+1.0;
  double eJ = etU*exUx11+etU*exUx22+etU*exUx33+(etU*etU)*exUx11*exUx22-(etU*etU)*exUx12*exUx21+(etU*etU)*exUx11*exUx33-(etU*etU)*exUx13*exUx31+(etU*etU)*exUx22*exUx33-(etU*etU)*exUx23*exUx32+(etU*etU*etU)*exUx11*exUx22*exUx33-(etU*etU*etU)*exUx11*exUx23*exUx32-(etU*etU*etU)*exUx12*exUx21*exUx33+(etU*etU*etU)*exUx12*exUx23*exUx31+(etU*etU*etU)*exUx13*exUx21*exUx32-(etU*etU*etU)*exUx13*exUx22*exUx31+1.0;
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
  double exUx11 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double exUx21 = (Y*Y)*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double exUx31 = Y*Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx12 = (X*X)*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1;
  double exUx22 = X*Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUx32 = X*Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/1.0E2);
  double exUx13 = (X*X)*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1;
  double exUx23 = X*(Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double exUx33 = X*Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/1.0E2);
  double etU = cos(M_PI*(1.0/4.0)-t*1.0E3)/(t*1.0E3+1.0)-sqrt(2.0)*(1.0/2.0);
  double etUt = cos(M_PI*(1.0/4.0)-t*1.0E3)*1.0/pow(t*1.0E3+1.0,2.0)*-1.0E3+(sin(M_PI*(1.0/4.0)-t*1.0E3)*1.0E3)/(t*1.0E3+1.0);
  double etUtt = (cos(M_PI*(1.0/4.0)-t*1.0E3)*-1.0E6)/(t*1.0E3+1.0)+cos(M_PI*(1.0/4.0)-t*1.0E3)*1.0/pow(t*1.0E3+1.0,3.0)*2.0E6-sin(M_PI*(1.0/4.0)-t*1.0E3)*1.0/pow(t*1.0E3+1.0,2.0)*2.0E6;
  double pU = -X*Y*Z*(X-1.0)*(Y-1.0)*(Z-1.0)*(exp(t*-5.0E1)*(7.0/1.0E1)-7.0/1.0E1)+1.0;
  double pUt = X*Y*Z*exp(t*-5.0E1)*(X-1.0)*(Y-1.0)*(Z-1.0)*3.5E1;
  double pUx3 = X*Y*exp(t*-5.0E1)*(Z*2.0-1.0)*(exp(t*5.0E1)-1.0)*(X-1.0)*(Y-1.0)*(7.0/1.0E1);

  double pJx1 = Y*Z*exp(t*-5.0E1)*(X*2.0-1.0)*(exp(t*5.0E1)-1.0)*(Y-1.0)*(Z-1.0)*(7.0/1.0E1);
  double pJx2 = X*Z*exp(t*-5.0E1)*(Y*2.0-1.0)*(exp(t*5.0E1)-1.0)*(X-1.0)*(Z-1.0)*(7.0/1.0E1);
  double pJx3 = X*Y*exp(t*-5.0E1)*(Z*2.0-1.0)*(exp(t*5.0E1)-1.0)*(X-1.0)*(Y-1.0)*(7.0/1.0E1);
  double Utt1 = Z*(etUtt*exUx13*pU+etU*exUx13*pUt+etUt*exUx13*pUt*2.0)+X*etUtt*exUx11+Y*etUtt*exUx12;
  double Utt2 = Z*(etUtt*exUx23*pU+etU*exUx23*pUt+etUt*exUx23*pUt*2.0)+X*etUtt*exUx21+Y*etUtt*exUx22;
  double Utt3 = Z*(pUt*(etU*exUx33+1.0)+etUtt*exUx33*pU+etUt*exUx33*pUt*2.0)+X*etUtt*exUx31+Y*etUtt*exUx32;

  double eF[9] = {};
  eF[0] = etU*exUx11+1.0;
  eF[1] = etU*exUx12;
  eF[2] = etU*exUx13;
  eF[3] = etU*exUx21;
  eF[4] = etU*exUx22+1.0;
  eF[5] = etU*exUx23;
  eF[6] = etU*exUx31;
  eF[7] = etU*exUx32;
  eF[8] = etU*exUx33+1.0;

  double *tempS = elast->S;
  double *tempL = elast->L;
  double eS[9] = {}, L[81]={};
  elast->S = eS;
  elast->L = L;
  elast->update_elasticity(elast,eF,1);
  elast->S = tempS;
  elast->L = tempL;
  double pJ = pU;
  double exUxx11 = Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-X*Y*Z*(X*4.0-3.0)*(Z-1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUxx21 = (Y*Y)*Z*(X*1.2E1-6.0)*pow(Y-1.0,2.0)*(Z-1.0)*1.0E1;
  double exUxx31 = Y*Z*(Y-1.0)*(Z-1.0)*(-1.0/5.0E1);
  double exUxx12 = X*Z*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-2.0E1-(X*X)*Z*(X*2.0-2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*1.0E1;
  double exUxx22 = Y*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1+X*Y*Z*(X*4.0-3.0)*(Z-1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUxx32 = Z*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/1.0E2)-X*Z*(Y*2.0-1.0)*(Z-1.0)*(1.0/1.0E2);
  double exUxx13 = X*Y*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-(X*X)*Y*(X*2.0-2.0)*(Z*2.0-1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*1.0E1;
  double exUxx23 = (Y*Y)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1+X*(Y*Y)*(X*4.0-3.0)*(Z*2.0-1.0)*pow(Y-1.0,2.0)*1.0E1;
  double exUxx33 = Y*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/1.0E2)-X*Y*(Z*2.0-1.0)*(Y-1.0)*(1.0/1.0E2);
  double exUxy11 = X*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-X*Y*Z*(Y*4.0-3.0)*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1;
  double exUxy21 = Y*Z*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*2.0E1+(Y*Y)*Z*(Y*2.0-2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double exUxy31 = Z*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/1.0E2)-Y*Z*(X*2.0-1.0)*(Z-1.0)*(1.0/1.0E2);
  double exUxy12 = (X*X)*Z*(Y*1.2E1-6.0)*pow(X-1.0,2.0)*(Z-1.0)*-1.0E1;
  double exUxy22 = X*Z*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1+X*Y*Z*(Y*4.0-3.0)*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1;
  double exUxy32 = X*Z*(X-1.0)*(Z-1.0)*(-1.0/5.0E1);
  double exUxy13 = (X*X)*(Z*2.0-1.0)*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-1.0E1-(X*X)*Y*(Y*4.0-3.0)*(Z*2.0-1.0)*pow(X-1.0,2.0)*1.0E1;
  double exUxy23 = X*Y*(Z*2.0-1.0)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1+X*(Y*Y)*(Y*2.0-2.0)*(Z*2.0-1.0)*(X*-3.0+(X*X)*2.0+1.0)*1.0E1;
  double exUxy33 = X*(Z*2.0-1.0)*(X-1.0)*(Y-1.0)*(-1.0/1.0E2)-X*Y*(Z*2.0-1.0)*(X-1.0)*(1.0/1.0E2);
  double exUxz11 = X*Y*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1-X*Y*Z*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUxz21 = (Y*Y)*pow(Y-1.0,2.0)*(Z-1.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1+(Y*Y)*Z*pow(Y-1.0,2.0)*(X*-6.0+(X*X)*6.0+1.0)*1.0E1;
  double exUxz31 = Y*(X*2.0-1.0)*(Y-1.0)*(Z-1.0)*(-1.0/1.0E2)-Y*Z*(X*2.0-1.0)*(Y-1.0)*(1.0/1.0E2);
  double exUxz12 = (X*X)*pow(X-1.0,2.0)*(Z-1.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*-1.0E1-(X*X)*Z*pow(X-1.0,2.0)*(Y*-6.0+(Y*Y)*6.0+1.0)*1.0E1;
  double exUxz22 = X*Y*(Z-1.0)*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1+X*Y*Z*(X*-3.0+(X*X)*2.0+1.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*2.0E1;
  double exUxz32 = X*(Y*2.0-1.0)*(X-1.0)*(Z-1.0)*(-1.0/1.0E2)-X*Z*(Y*2.0-1.0)*(X-1.0)*(1.0/1.0E2);
  double exUxz13 = (X*X)*Y*pow(X-1.0,2.0)*(Y*-3.0+(Y*Y)*2.0+1.0)*-2.0E1;
  double exUxz23 = X*(Y*Y)*pow(Y-1.0,2.0)*(X*-3.0+(X*X)*2.0+1.0)*2.0E1;
  double exUxz33 = X*Y*(X-1.0)*(Y-1.0)*(-1.0/5.0E1);
  double eCx11 = etUt*exUxx11*(etU*exUx11+1.0)*2.0+etU*etUt*exUx21*exUxx21*2.0+etU*etUt*exUx31*exUxx31*2.0;
  double eCx21 = etUt*exUxx12*(etU*exUx11+1.0)+etUt*exUxx21*(etU*exUx22+1.0)+etU*etUt*exUx12*exUxx11+etU*etUt*exUx21*exUxx22+etU*etUt*exUx31*exUxx32+etU*etUt*exUx32*exUxx31;
  double eCx31 = etUt*exUxx13*(etU*exUx11+1.0)+etUt*exUxx31*(etU*exUx33+1.0)+etU*etUt*exUx13*exUxx11+etU*etUt*exUx21*exUxx23+etU*etUt*exUx23*exUxx21+etU*etUt*exUx31*exUxx33;
  double eCx12 = etUt*exUxx12*(etU*exUx11+1.0)+etUt*exUxx21*(etU*exUx22+1.0)+etU*etUt*exUx12*exUxx11+etU*etUt*exUx21*exUxx22+etU*etUt*exUx31*exUxx32+etU*etUt*exUx32*exUxx31;
  double eCx22 = etUt*exUxx22*(etU*exUx22+1.0)*2.0+etU*etUt*exUx12*exUxx12*2.0+etU*etUt*exUx32*exUxx32*2.0;
  double eCx32 = etUt*exUxx23*(etU*exUx22+1.0)+etUt*exUxx32*(etU*exUx33+1.0)+etU*etUt*exUx12*exUxx13+etU*etUt*exUx13*exUxx12+etU*etUt*exUx23*exUxx22+etU*etUt*exUx32*exUxx33;
  double eCx13 = etUt*exUxx13*(etU*exUx11+1.0)+etUt*exUxx31*(etU*exUx33+1.0)+etU*etUt*exUx13*exUxx11+etU*etUt*exUx21*exUxx23+etU*etUt*exUx23*exUxx21+etU*etUt*exUx31*exUxx33;
  double eCx23 = etUt*exUxx23*(etU*exUx22+1.0)+etUt*exUxx32*(etU*exUx33+1.0)+etU*etUt*exUx12*exUxx13+etU*etUt*exUx13*exUxx12+etU*etUt*exUx23*exUxx22+etU*etUt*exUx32*exUxx33;
  double eCx33 = etUt*exUxx33*(etU*exUx33+1.0)*2.0+etU*etUt*exUx13*exUxx13*2.0+etU*etUt*exUx23*exUxx23*2.0;
  double eCy11 = etUt*exUxy11*(etU*exUx11+1.0)*2.0+etU*etUt*exUx21*exUxy21*2.0+etU*etUt*exUx31*exUxy31*2.0;
  double eCy21 = etUt*exUxy12*(etU*exUx11+1.0)+etUt*exUxy21*(etU*exUx22+1.0)+etU*etUt*exUx12*exUxy11+etU*etUt*exUx21*exUxy22+etU*etUt*exUx31*exUxy32+etU*etUt*exUx32*exUxy31;
  double eCy31 = etUt*exUxy13*(etU*exUx11+1.0)+etUt*exUxy31*(etU*exUx33+1.0)+etU*etUt*exUx13*exUxy11+etU*etUt*exUx21*exUxy23+etU*etUt*exUx23*exUxy21+etU*etUt*exUx31*exUxy33;
  double eCy12 = etUt*exUxy12*(etU*exUx11+1.0)+etUt*exUxy21*(etU*exUx22+1.0)+etU*etUt*exUx12*exUxy11+etU*etUt*exUx21*exUxy22+etU*etUt*exUx31*exUxy32+etU*etUt*exUx32*exUxy31;
  double eCy22 = etUt*exUxy22*(etU*exUx22+1.0)*2.0+etU*etUt*exUx12*exUxy12*2.0+etU*etUt*exUx32*exUxy32*2.0;
  double eCy32 = etUt*exUxy23*(etU*exUx22+1.0)+etUt*exUxy32*(etU*exUx33+1.0)+etU*etUt*exUx12*exUxy13+etU*etUt*exUx13*exUxy12+etU*etUt*exUx23*exUxy22+etU*etUt*exUx32*exUxy33;
  double eCy13 = etUt*exUxy13*(etU*exUx11+1.0)+etUt*exUxy31*(etU*exUx33+1.0)+etU*etUt*exUx13*exUxy11+etU*etUt*exUx21*exUxy23+etU*etUt*exUx23*exUxy21+etU*etUt*exUx31*exUxy33;
  double eCy23 = etUt*exUxy23*(etU*exUx22+1.0)+etUt*exUxy32*(etU*exUx33+1.0)+etU*etUt*exUx12*exUxy13+etU*etUt*exUx13*exUxy12+etU*etUt*exUx23*exUxy22+etU*etUt*exUx32*exUxy33;
  double eCy33 = etUt*exUxy33*(etU*exUx33+1.0)*2.0+etU*etUt*exUx13*exUxy13*2.0+etU*etUt*exUx23*exUxy23*2.0;
  double eCz11 = etUt*exUxz11*(etU*exUx11+1.0)*2.0+etU*etUt*exUx21*exUxz21*2.0+etU*etUt*exUx31*exUxz31*2.0;
  double eCz21 = etUt*exUxz12*(etU*exUx11+1.0)+etUt*exUxz21*(etU*exUx22+1.0)+etU*etUt*exUx12*exUxz11+etU*etUt*exUx21*exUxz22+etU*etUt*exUx31*exUxz32+etU*etUt*exUx32*exUxz31;
  double eCz31 = etUt*exUxz13*(etU*exUx11+1.0)+etUt*exUxz31*(etU*exUx33+1.0)+etU*etUt*exUx13*exUxz11+etU*etUt*exUx21*exUxz23+etU*etUt*exUx23*exUxz21+etU*etUt*exUx31*exUxz33;
  double eCz12 = etUt*exUxz12*(etU*exUx11+1.0)+etUt*exUxz21*(etU*exUx22+1.0)+etU*etUt*exUx12*exUxz11+etU*etUt*exUx21*exUxz22+etU*etUt*exUx31*exUxz32+etU*etUt*exUx32*exUxz31;
  double eCz22 = etUt*exUxz22*(etU*exUx22+1.0)*2.0+etU*etUt*exUx12*exUxz12*2.0+etU*etUt*exUx32*exUxz32*2.0;
  double eCz32 = etUt*exUxz23*(etU*exUx22+1.0)+etUt*exUxz32*(etU*exUx33+1.0)+etU*etUt*exUx12*exUxz13+etU*etUt*exUx13*exUxz12+etU*etUt*exUx23*exUxz22+etU*etUt*exUx32*exUxz33;
  double eCz13 = etUt*exUxz13*(etU*exUx11+1.0)+etUt*exUxz31*(etU*exUx33+1.0)+etU*etUt*exUx13*exUxz11+etU*etUt*exUx21*exUxz23+etU*etUt*exUx23*exUxz21+etU*etUt*exUx31*exUxz33;
  double eCz23 = etUt*exUxz23*(etU*exUx22+1.0)+etUt*exUxz32*(etU*exUx33+1.0)+etU*etUt*exUx12*exUxz13+etU*etUt*exUx13*exUxz12+etU*etUt*exUx23*exUxz22+etU*etUt*exUx32*exUxz33;
  double eCz33 = etUt*exUxz33*(etU*exUx33+1.0)*2.0+etU*etUt*exUx13*exUxz13*2.0+etU*etUt*exUx23*exUxz23*2.0;
  double L1111 = L[0];
  double L1112 = L[1];
  double L1113 = L[2];
  double L1121 = L[3];
  double L1122 = L[4];
  double L1123 = L[5];
  double L1131 = L[6];
  double L1132 = L[7];
  double L1133 = L[8];
  double L1211 = L[9];
  double L1212 = L[10];
  double L1213 = L[11];
  double L1221 = L[12];
  double L1222 = L[13];
  double L1223 = L[14];
  double L1231 = L[15];
  double L1232 = L[16];
  double L1233 = L[17];
  double L1311 = L[18];
  double L1312 = L[19];
  double L1313 = L[20];
  double L1321 = L[21];
  double L1322 = L[22];
  double L1323 = L[23];
  double L1331 = L[24];
  double L1332 = L[25];
  double L1333 = L[26];
  double L2111 = L[27];
  double L2112 = L[28];
  double L2113 = L[29];
  double L2121 = L[30];
  double L2122 = L[31];
  double L2123 = L[32];
  double L2131 = L[33];
  double L2132 = L[34];
  double L2133 = L[35];
  double L2211 = L[36];
  double L2212 = L[37];
  double L2213 = L[38];
  double L2221 = L[39];
  double L2222 = L[40];
  double L2223 = L[41];
  double L2231 = L[42];
  double L2232 = L[43];
  double L2233 = L[44];
  double L2311 = L[45];
  double L2312 = L[46];
  double L2313 = L[47];
  double L2321 = L[48];
  double L2322 = L[49];
  double L2323 = L[50];
  double L2331 = L[51];
  double L2332 = L[52];
  double L2333 = L[53];
  double L3111 = L[54];
  double L3112 = L[55];
  double L3113 = L[56];
  double L3121 = L[57];
  double L3122 = L[58];
  double L3123 = L[59];
  double L3131 = L[60];
  double L3132 = L[61];
  double L3133 = L[62];
  double L3211 = L[63];
  double L3212 = L[64];
  double L3213 = L[65];
  double L3221 = L[66];
  double L3222 = L[67];
  double L3223 = L[68];
  double L3231 = L[69];
  double L3232 = L[70];
  double L3233 = L[71];
  double L3311 = L[72];
  double L3312 = L[73];
  double L3313 = L[74];
  double L3321 = L[75];
  double L3322 = L[76];
  double L3323 = L[77];
  double L3331 = L[78];
  double L3332 = L[79];
  double L3333 = L[80];
  double eF11 = eF[0];
  double eF12 = eF[1];
  double eF13 = eF[2];
  double eF21 = eF[3];
  double eF22 = eF[4];
  double eF23 = eF[5];
  double eF31 = eF[6];
  double eF32 = eF[7];
  double eF33 = eF[8];
  double eS11 = eS[0];
  double eS12 = eS[1];
  double eS13 = eS[2];
  double eS21 = eS[3];
  double eS22 = eS[4];
  double eS23 = eS[5];
  double eS31 = eS[6];
  double eS32 = eS[7];
  double eS33 = eS[8];

  b[0] = -(eF11*pJ*(L1311*eCz11*(1.0/2.0)+L1312*eCz12*(1.0/2.0)+L1313*eCz13*(1.0/2.0)+L1321*eCz21*(1.0/2.0)+L1322*eCz22*(1.0/2.0)+L1323*eCz23*(1.0/2.0)+L1331*eCz31*(1.0/2.0)+L1332*eCz32*(1.0/2.0)+L1333*eCz33*(1.0/2.0))+eF12*pJ*(L2311*eCz11*(1.0/2.0)+L2312*eCz12*(1.0/2.0)+L2313*eCz13*(1.0/2.0)+L2321*eCz21*(1.0/2.0)+L2322*eCz22*(1.0/2.0)+L2323*eCz23*(1.0/2.0)+L2331*eCz31*(1.0/2.0)+L2332*eCz32*(1.0/2.0)+L2333*eCz33*(1.0/2.0))+eF13*pJ*(L3311*eCz11*(1.0/2.0)+L3312*eCz12*(1.0/2.0)+L3313*eCz13*(1.0/2.0)+L3321*eCz21*(1.0/2.0)+L3322*eCz22*(1.0/2.0)+L3323*eCz23*(1.0/2.0)+L3331*eCz31*(1.0/2.0)+L3332*eCz32*(1.0/2.0)+L3333*eCz33*(1.0/2.0)))/pU+Utt1*rho-(eF11*eS13*pJx3+eF12*eS23*pJx3+eF13*eS33*pJx3)/pU-(eS13*etU*exUxz11*pJ+eS23*etU*exUxz12*pJ+eS33*etU*exUxz13*pJ)/pU-eF11*pJ*(L1111*eCx11*(1.0/2.0)+L1112*eCx12*(1.0/2.0)+L1113*eCx13*(1.0/2.0)+L1121*eCx21*(1.0/2.0)+L1122*eCx22*(1.0/2.0)+L1123*eCx23*(1.0/2.0)+L1131*eCx31*(1.0/2.0)+L1132*eCx32*(1.0/2.0)+L1133*eCx33*(1.0/2.0))-eF11*pJ*(L1211*eCy11*(1.0/2.0)+L1212*eCy12*(1.0/2.0)+L1213*eCy13*(1.0/2.0)+L1221*eCy21*(1.0/2.0)+L1222*eCy22*(1.0/2.0)+L1223*eCy23*(1.0/2.0)+L1231*eCy31*(1.0/2.0)+L1232*eCy32*(1.0/2.0)+L1233*eCy33*(1.0/2.0))-eF12*pJ*(L2111*eCx11*(1.0/2.0)+L2112*eCx12*(1.0/2.0)+L2113*eCx13*(1.0/2.0)+L2121*eCx21*(1.0/2.0)+L2122*eCx22*(1.0/2.0)+L2123*eCx23*(1.0/2.0)+L2131*eCx31*(1.0/2.0)+L2132*eCx32*(1.0/2.0)+L2133*eCx33*(1.0/2.0))-eF12*pJ*(L2211*eCy11*(1.0/2.0)+L2212*eCy12*(1.0/2.0)+L2213*eCy13*(1.0/2.0)+L2221*eCy21*(1.0/2.0)+L2222*eCy22*(1.0/2.0)+L2223*eCy23*(1.0/2.0)+L2231*eCy31*(1.0/2.0)+L2232*eCy32*(1.0/2.0)+L2233*eCy33*(1.0/2.0))-eF13*pJ*(L3111*eCx11*(1.0/2.0)+L3112*eCx12*(1.0/2.0)+L3113*eCx13*(1.0/2.0)+L3121*eCx21*(1.0/2.0)+L3122*eCx22*(1.0/2.0)+L3123*eCx23*(1.0/2.0)+L3131*eCx31*(1.0/2.0)+L3132*eCx32*(1.0/2.0)+L3133*eCx33*(1.0/2.0))-eF13*pJ*(L3211*eCy11*(1.0/2.0)+L3212*eCy12*(1.0/2.0)+L3213*eCy13*(1.0/2.0)+L3221*eCy21*(1.0/2.0)+L3222*eCy22*(1.0/2.0)+L3223*eCy23*(1.0/2.0)+L3231*eCy31*(1.0/2.0)+L3232*eCy32*(1.0/2.0)+L3233*eCy33*(1.0/2.0))-eF11*eS11*pJx1-eF11*eS12*pJx2-eF12*eS21*pJx1-eF12*eS22*pJx2-eF13*eS31*pJx1-eF13*eS32*pJx2+1.0/(pU*pU)*pUx3*(eF11*eS13*pJ+eF12*eS23*pJ+eF13*eS33*pJ)-eS11*etU*exUxx11*pJ-eS12*etU*exUxy11*pJ-eS21*etU*exUxx12*pJ-eS22*etU*exUxy12*pJ-eS31*etU*exUxx13*pJ-eS32*etU*exUxy13*pJ;
  b[1] = -(eF21*pJ*(L1311*eCz11*(1.0/2.0)+L1312*eCz12*(1.0/2.0)+L1313*eCz13*(1.0/2.0)+L1321*eCz21*(1.0/2.0)+L1322*eCz22*(1.0/2.0)+L1323*eCz23*(1.0/2.0)+L1331*eCz31*(1.0/2.0)+L1332*eCz32*(1.0/2.0)+L1333*eCz33*(1.0/2.0))+eF22*pJ*(L2311*eCz11*(1.0/2.0)+L2312*eCz12*(1.0/2.0)+L2313*eCz13*(1.0/2.0)+L2321*eCz21*(1.0/2.0)+L2322*eCz22*(1.0/2.0)+L2323*eCz23*(1.0/2.0)+L2331*eCz31*(1.0/2.0)+L2332*eCz32*(1.0/2.0)+L2333*eCz33*(1.0/2.0))+eF23*pJ*(L3311*eCz11*(1.0/2.0)+L3312*eCz12*(1.0/2.0)+L3313*eCz13*(1.0/2.0)+L3321*eCz21*(1.0/2.0)+L3322*eCz22*(1.0/2.0)+L3323*eCz23*(1.0/2.0)+L3331*eCz31*(1.0/2.0)+L3332*eCz32*(1.0/2.0)+L3333*eCz33*(1.0/2.0)))/pU+Utt2*rho-(eF21*eS13*pJx3+eF22*eS23*pJx3+eF23*eS33*pJx3)/pU-(eS13*etU*exUxz21*pJ+eS23*etU*exUxz22*pJ+eS33*etU*exUxz23*pJ)/pU-eF21*pJ*(L1111*eCx11*(1.0/2.0)+L1112*eCx12*(1.0/2.0)+L1113*eCx13*(1.0/2.0)+L1121*eCx21*(1.0/2.0)+L1122*eCx22*(1.0/2.0)+L1123*eCx23*(1.0/2.0)+L1131*eCx31*(1.0/2.0)+L1132*eCx32*(1.0/2.0)+L1133*eCx33*(1.0/2.0))-eF21*pJ*(L1211*eCy11*(1.0/2.0)+L1212*eCy12*(1.0/2.0)+L1213*eCy13*(1.0/2.0)+L1221*eCy21*(1.0/2.0)+L1222*eCy22*(1.0/2.0)+L1223*eCy23*(1.0/2.0)+L1231*eCy31*(1.0/2.0)+L1232*eCy32*(1.0/2.0)+L1233*eCy33*(1.0/2.0))-eF22*pJ*(L2111*eCx11*(1.0/2.0)+L2112*eCx12*(1.0/2.0)+L2113*eCx13*(1.0/2.0)+L2121*eCx21*(1.0/2.0)+L2122*eCx22*(1.0/2.0)+L2123*eCx23*(1.0/2.0)+L2131*eCx31*(1.0/2.0)+L2132*eCx32*(1.0/2.0)+L2133*eCx33*(1.0/2.0))-eF22*pJ*(L2211*eCy11*(1.0/2.0)+L2212*eCy12*(1.0/2.0)+L2213*eCy13*(1.0/2.0)+L2221*eCy21*(1.0/2.0)+L2222*eCy22*(1.0/2.0)+L2223*eCy23*(1.0/2.0)+L2231*eCy31*(1.0/2.0)+L2232*eCy32*(1.0/2.0)+L2233*eCy33*(1.0/2.0))-eF23*pJ*(L3111*eCx11*(1.0/2.0)+L3112*eCx12*(1.0/2.0)+L3113*eCx13*(1.0/2.0)+L3121*eCx21*(1.0/2.0)+L3122*eCx22*(1.0/2.0)+L3123*eCx23*(1.0/2.0)+L3131*eCx31*(1.0/2.0)+L3132*eCx32*(1.0/2.0)+L3133*eCx33*(1.0/2.0))-eF23*pJ*(L3211*eCy11*(1.0/2.0)+L3212*eCy12*(1.0/2.0)+L3213*eCy13*(1.0/2.0)+L3221*eCy21*(1.0/2.0)+L3222*eCy22*(1.0/2.0)+L3223*eCy23*(1.0/2.0)+L3231*eCy31*(1.0/2.0)+L3232*eCy32*(1.0/2.0)+L3233*eCy33*(1.0/2.0))-eF21*eS11*pJx1-eF21*eS12*pJx2-eF22*eS21*pJx1-eF22*eS22*pJx2-eF23*eS31*pJx1-eF23*eS32*pJx2+1.0/(pU*pU)*pUx3*(eF21*eS13*pJ+eF22*eS23*pJ+eF23*eS33*pJ)-eS11*etU*exUxx21*pJ-eS12*etU*exUxy21*pJ-eS21*etU*exUxx22*pJ-eS22*etU*exUxy22*pJ-eS31*etU*exUxx23*pJ-eS32*etU*exUxy23*pJ;
  b[2] = -(eF31*pJ*(L1311*eCz11*(1.0/2.0)+L1312*eCz12*(1.0/2.0)+L1313*eCz13*(1.0/2.0)+L1321*eCz21*(1.0/2.0)+L1322*eCz22*(1.0/2.0)+L1323*eCz23*(1.0/2.0)+L1331*eCz31*(1.0/2.0)+L1332*eCz32*(1.0/2.0)+L1333*eCz33*(1.0/2.0))+eF32*pJ*(L2311*eCz11*(1.0/2.0)+L2312*eCz12*(1.0/2.0)+L2313*eCz13*(1.0/2.0)+L2321*eCz21*(1.0/2.0)+L2322*eCz22*(1.0/2.0)+L2323*eCz23*(1.0/2.0)+L2331*eCz31*(1.0/2.0)+L2332*eCz32*(1.0/2.0)+L2333*eCz33*(1.0/2.0))+eF33*pJ*(L3311*eCz11*(1.0/2.0)+L3312*eCz12*(1.0/2.0)+L3313*eCz13*(1.0/2.0)+L3321*eCz21*(1.0/2.0)+L3322*eCz22*(1.0/2.0)+L3323*eCz23*(1.0/2.0)+L3331*eCz31*(1.0/2.0)+L3332*eCz32*(1.0/2.0)+L3333*eCz33*(1.0/2.0)))/pU+Utt3*rho-(eF31*eS13*pJx3+eF32*eS23*pJx3+eF33*eS33*pJx3)/pU-(eS13*etU*exUxz31*pJ+eS23*etU*exUxz32*pJ+eS33*etU*exUxz33*pJ)/pU-eF31*pJ*(L1111*eCx11*(1.0/2.0)+L1112*eCx12*(1.0/2.0)+L1113*eCx13*(1.0/2.0)+L1121*eCx21*(1.0/2.0)+L1122*eCx22*(1.0/2.0)+L1123*eCx23*(1.0/2.0)+L1131*eCx31*(1.0/2.0)+L1132*eCx32*(1.0/2.0)+L1133*eCx33*(1.0/2.0))-eF31*pJ*(L1211*eCy11*(1.0/2.0)+L1212*eCy12*(1.0/2.0)+L1213*eCy13*(1.0/2.0)+L1221*eCy21*(1.0/2.0)+L1222*eCy22*(1.0/2.0)+L1223*eCy23*(1.0/2.0)+L1231*eCy31*(1.0/2.0)+L1232*eCy32*(1.0/2.0)+L1233*eCy33*(1.0/2.0))-eF32*pJ*(L2111*eCx11*(1.0/2.0)+L2112*eCx12*(1.0/2.0)+L2113*eCx13*(1.0/2.0)+L2121*eCx21*(1.0/2.0)+L2122*eCx22*(1.0/2.0)+L2123*eCx23*(1.0/2.0)+L2131*eCx31*(1.0/2.0)+L2132*eCx32*(1.0/2.0)+L2133*eCx33*(1.0/2.0))-eF32*pJ*(L2211*eCy11*(1.0/2.0)+L2212*eCy12*(1.0/2.0)+L2213*eCy13*(1.0/2.0)+L2221*eCy21*(1.0/2.0)+L2222*eCy22*(1.0/2.0)+L2223*eCy23*(1.0/2.0)+L2231*eCy31*(1.0/2.0)+L2232*eCy32*(1.0/2.0)+L2233*eCy33*(1.0/2.0))-eF33*pJ*(L3111*eCx11*(1.0/2.0)+L3112*eCx12*(1.0/2.0)+L3113*eCx13*(1.0/2.0)+L3121*eCx21*(1.0/2.0)+L3122*eCx22*(1.0/2.0)+L3123*eCx23*(1.0/2.0)+L3131*eCx31*(1.0/2.0)+L3132*eCx32*(1.0/2.0)+L3133*eCx33*(1.0/2.0))-eF33*pJ*(L3211*eCy11*(1.0/2.0)+L3212*eCy12*(1.0/2.0)+L3213*eCy13*(1.0/2.0)+L3221*eCy21*(1.0/2.0)+L3222*eCy22*(1.0/2.0)+L3223*eCy23*(1.0/2.0)+L3231*eCy31*(1.0/2.0)+L3232*eCy32*(1.0/2.0)+L3233*eCy33*(1.0/2.0))-eF31*eS11*pJx1-eF31*eS12*pJx2-eF32*eS21*pJx1-eF32*eS22*pJx2-eF33*eS31*pJx1-eF33*eS32*pJx2+1.0/(pU*pU)*pUx3*(eF31*eS13*pJ+eF32*eS23*pJ+eF33*eS33*pJ)-eS11*etU*exUxx31*pJ-eS12*etU*exUxy31*pJ-eS21*etU*exUxx32*pJ-eS22*etU*exUxy32*pJ-eS31*etU*exUxx33*pJ-eS32*etU*exUxy33*pJ;
}

