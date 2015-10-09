/**
 * This file defines the interface to the BPA plasticity model.
 * REFERENCES:
 *
 * Mary C. BOYCE, David M. PARKS, Ali S. ARGON (1988). LARGE INELASTIC
 * DEFORMATION OF GLASSY POLYMERS.PART I: RATE DEPENDENT CONSTITUTIVE
 * MODEL. Mechanics of Materials, 7:15-33.
 *
 * Holopainen, S. (2013). Modeling of the mechanical behavior of
 * amorphous glassy polymers under variable loadings and comparison
 * with state-of-the-art model predictions. Mechanics of Materials,
 * 66:35–58.
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#pragma once
#ifndef PLASTICITY_MODEL_BPA_H
#define PLASTICITY_MODEL_BPA_H

#include <stdio.h>

#ifndef TYPE_CONSTITUTIVE_MODEL
#define TYPE_CONSTITUTIVE_MODEL
typedef struct Constitutive_model Constitutive_model;
#endif

#ifndef TYPE_MODEL_PARAMETERS
#define TYPE_MODEL_PARAMETERS
typedef struct Model_parameters Model_parameters;
#endif

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Model_parameters object
 * \return non-zero on internal error
 */
int plasticity_model_BPA_initialize(Model_parameters *p);

/**
 * Read the model parameters from a file.
 *
 * \param[in,out] p - pointer to an inialized Model_parameters object
 * \return non-zero on internal error
 */
int plasticity_model_BPA_read(Model_parameters *p,
                              FILE *in);

/**
 * Set the initial values for the state variables.
 *
 * \return non-zero on error.
 */
int plasticity_model_BPA_set_initial_values(Constitutive_model *m);

/**
 * Construct and initialize the model context for calling functions
 * through the plasticity interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F, _total_ deformation gradient
 * \param[in] dt, time increment
 * \return non-zero on internal error.
 */
int plasticity_model_BPA_ctx_build(void **ctx,
                                   const double *F,
                                   const double dt);

/**
 * Destroy the model context and invalidate the handle.
 *
 * \param[in,out] ctx - handle to context object. ctx = NULL on exit.
 * \return non-zero on internal error.
 */
int plasticity_model_BPA_ctx_destroy(void **ctx);


#endif
