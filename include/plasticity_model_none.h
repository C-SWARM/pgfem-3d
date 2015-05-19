/**
 * Declare/define functions that integrate with the plasticty model
 * interface (see plasticity.h). This file defines the interface for
 * the purely hyperelastic model.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */
#pragma once
#ifndef PLASTICITY_MODEL_NONE_H
#define PLASTICITY_MODEL_NONE_H

/** Pre-declare and typedef Plasticity structure */
struct Plasticity;
#ifndef TYPE_PLASTICITY
#define TYPE_PLASTICITY
typedef struct Plasticity Plasticity;
#endif

/**
 * Initialize the plasticity object for this particular model.
 *
 * \param[in,out] p - pointer to a Plasticity object
 * \return non-zero on internal error
 */
int plasticity_model_none_initialize(Plasticity *p);

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
int plasticity_model_none_ctx_build(void **ctx,
                                    const double *C,
                                    const double *J_or_Theta);

/**
 * Destroy the model context and invalidate the handle.
 *
 * \param[in,out] ctx - handle to context object. ctx = NULL on exit.
 * \return non-zero on internal error.
 */
int plasticity_model_none_ctx_destroy(void **ctx);

#endif
