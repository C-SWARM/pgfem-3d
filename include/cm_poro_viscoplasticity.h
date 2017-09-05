/// Declare/define functions for constitutive model interface performing
/// poro-viscoplasticity integration algorithm
/// 
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// Alberto Salvadori, [1], <asalvad2@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#pragma once
#ifndef CM_PORO_VISCO_PLASTICITY_H
#define CM_PORO_VISCO_PLASTICITY_H

#ifndef TYPE_CONSTITUTIVE_MODEL
#define TYPE_CONSTITUTIVE_MODEL
typedef struct Constitutive_model Constitutive_model;
#endif

#ifndef TYPE_MODEL_PARAMETERS
#define TYPE_MODEL_PARAMETERS
typedef struct Model_parameters Model_parameters;
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

/// Initialize the Model_parameters object for this particular model.
/// 
/// \param[in,out] p - pointer to a Model_parameters object
/// \return non-zero on internal error
int poro_viscoplasticity_model_initialize(Model_parameters *p);
 
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
int poro_viscoplasticity_model_ctx_build(void **ctx,
                                         double *F,
                                         const double dt,
                                         const double alpha,
                                         double *eFnpa,
                                         double *hFn,
                                         double *hFnp1,
                                         const int is_coulpled_with_thermal,
                                         const int npa);
                                         
int poro_viscoplasticity_model_destroy(Model_parameters *p);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */


#endif