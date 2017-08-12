/**
 * Declare/define functions that integrate with the plasticty model
 * interface (see constitutive_model.h). This file defines the
 * interface for the purely hyperelastic model.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */
#pragma once
#ifndef PLASTICITY_MODEL_NONE_H
#define PLASTICITY_MODEL_NONE_H

struct Model_parameters;

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Constitutive_model object
 * \return non-zero on internal error
 */
int plasticity_model_none_initialize(Model_parameters *p);

/// Construct and initialize the model context for calling functions
/// through the constitutive model interface.
///
/// \param[in,out] ctx - handle to an opaque model context object.
/// \param[in] F The total deformation gradient.
/// \param[in] eFnpa elastic deformation gradient at t = n + alpha
/// \param[in] hFn thermal part deformation gradient at t = n
/// \param[in] hFnp1 thermal part deformation gradient at t = n + 1
/// \param[in] is_coulpled_with_thermal flag for coupling with thermal
/// \return non-zero on internal error.
int plasticity_model_none_ctx_build(void **ctx,
                                    double *F,
                                    double *eFnpa,
                                    double *hFn,
                                    double *hFnp1,
                                    const int is_coulpled_with_thermal);

/**
 * Destroy the model context and invalidate the handle.
 *
 * \param[in,out] ctx - handle to context object. ctx = NULL on exit.
 * \return non-zero on internal error.
 */
int plasticity_model_none_ctx_destroy(void **ctx);

#endif
