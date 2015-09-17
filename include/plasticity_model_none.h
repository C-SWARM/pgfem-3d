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
#ifndef TYPE_MODEL_PARAMETERS
#define TYPE_MODEL_PARAMETERS
typedef struct Model_parameters Model_parameters;
#endif

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Constitutive_model object
 * \return non-zero on internal error
 */
int plasticity_model_none_initialize(Model_parameters *p);

/**
 * Construct and initialize the model context for calling functions
 * through the plasticity interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F - The total deformation gradient.
 * \return non-zero on internal error.
 */
int plasticity_model_none_ctx_build(void **ctx,
                                    const double *F);

/**
 * Destroy the model context and invalidate the handle.
 *
 * \param[in,out] ctx - handle to context object. ctx = NULL on exit.
 * \return non-zero on internal error.
 */
int plasticity_model_none_ctx_destroy(void **ctx);

#endif
