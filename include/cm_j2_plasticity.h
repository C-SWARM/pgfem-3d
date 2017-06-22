/* HEADER */
/**
 * This file defines the interface for the J2 plasticity + damage
 * model.
 *
 * REFERENCES:
 *
 * Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
 * finite strains." Computational Mechanics 5.5 (1989): 375-400.
 *
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#pragma once
#ifndef CM_J2_PLASTICITY_H
#define CM_J2_PLASTICITY_H

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

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Model_parameters object
 * \return non-zero on internal error
 */
int j2d_plasticity_model_initialize(Model_parameters *p){return 0;};

/**
 * Construct and initialize the model context for calling functions
 * through the constitutive modeling interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F, _total_ deformation gradient
 * \param[in] dt, time increment
 * \return non-zero on internal error.
 */
int j2d_plasticity_model_ctx_build(void **ctx,
                                   const double *F,
                                   const double dt) {return 0;};


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
