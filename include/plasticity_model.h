/**
 * Declare/define functions that integrate with the plasticty model
 * interface (see constitutive_model.h). 
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 *  Sangmin Lee, University of Notre Dame, Notre Dame, IN, <slee43@nd.edu>
 */
#pragma once
#ifndef PLASTICITY_MODEL_H
#define PLASTICITY_MODEL_H

struct Model_parameters;
#ifndef TYPE_MODEL_PARAMETERS
#define TYPE_MODEL_PARAMETERS
typedef struct Model_parameters Model_parameters;
#endif

#ifndef TYPE_CONSTITUTIVE_MODEL
#define TYPE_CONSTITUTIVE_MODEL
typedef struct Constitutive_model Constitutive_model;
#endif

struct Matrix_double;
#ifndef TYPE_MATRIX_DOUBLE
#define TYPE_MATRIX_DOUBLE
typedef struct Matrix_double Matrix_double;
#endif

enum variable_names {
  VAR_L_n,
  VAR_L_np1,  
  VAR_g_n,
  VAR_g_np1,  
  VAR_gamma_dot_0,
  VAR_gamma_dot_s, 
  VAR_m,
  VAR_g0,  
  VAR_G0,    
  VAR_gs_0,
  VAR_w
};
  
enum tensor_names {
  TENSOR_Fn,
  TENSOR_pFn,
  TENSOR_Fnp1,
  TENSOR_pFnp1,  
  TENSOR_tau,
  TENSOR_gamma_dot
};

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Constitutive_model object
 * \return non-zero on internal error
 */
int plasticity_model_initialize(Model_parameters *p);

/**
 * Construct and initialize the model context for calling functions
 * through the plasticity interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] C - The total left Cauchy-Green deformation tensor.
 * \param[in] J_or_Theta - The Jacobian of the deformation -OR- the
 *   volume field (depending on FE formulation)
 * \return non-zero on internal error.
 *
 * CAVEATES: The addresses of C and J_or_Theta must remain valid
 * throughout the existence of the ctx. Destroying these memory
 * locations before calling plasticity_model_none_ctx_destroy may
 * invalidate ctx.
 */
int plasticity_model_ctx_build(void **ctx,
                               const double *C,
                               const double *J_or_Theta,
                               const double dt);

/**
 * Destroy the model context and invalidate the handle.
 *
 * \param[in,out] ctx - handle to context object. ctx = NULL on exit.
 * \return non-zero on internal error.
 */
int plasticity_model_ctx_destroy(void **ctx);

int compute_dMdu(const Constitutive_model *m,
                 Matrix_double *dMdu,
                 const Matrix_double *Grad_du,
                 const Matrix_double *eFn,
                 const Matrix_double *eFnp1,
                 const Matrix_double *M,
                 const Matrix_double *S,
                 const Matrix_double *L,
                 const double dt);
/**
 * compute tangent of the plasticity part of deformation gradient with respect to the deformation
 *
 * \param[in,out] dMdu - computed tangent of plasticity part of deformation gradient
 * \param[in] du - Grad(du)
 * \param[in,out] m - pointer to a Constitutive_model object.
 * \param[in] dt - time step size
 * \return non-zero on internal error.
 */
int plasticity_model_slip_system(Matrix_double *P);

int plasticity_model_integration_ip(Matrix_double *pFnp1,
                                    Constitutive_model *m,
                                    const Matrix_double *Fnp1,
                                    const Matrix_double *Fe_n,
                                    const double dt);

int plasticity_model_read_parameters(Constitutive_model *m);

typedef struct HOMMAT HOMMAT;

int plasticity_model_test(const HOMMAT *hmat, Matrix_double *L_in, int Load_Type);
 
#endif
