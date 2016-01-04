/* HEADER */
/**
 * This file defines the interface for the isotropic viscous
 * damage model.
 *
 * REFERENCES:
 *
 * Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
 * finite strains." Computational Mechanics 5.5 (1989): 375-400.
 *
 * Mosby, Matthew and K. Matous. "On mechanics and material length scales of failure
 * in heterogeneous interfaces using a finite strain high performance
 * solver." Modelling and Simulation in Materials Science and
 * Engineering 23.8 (2015): 085014.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#pragma once
#ifndef CM_ISO_VISCOUS_DAMAGE_H
#define CM_ISO_VISCOUS_DAMAGE_H

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
int iso_viscous_damage_model_initialize(Model_parameters *p);

/**
 * Construct and initialize the model context for calling functions
 * through the constitutive modeling interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F, _total_ deformation gradient
 * \param[in] dt, time increment
 * \return non-zero on internal error.
 */
int iso_viscous_damage_model_ctx_build(void **ctx,
                                       const double *F,
                                       const double dt);

/**
 * Destroy the model context and invalidate the handle.
 *
 * \param[in,out] ctx - handle to context object. ctx = NULL on exit.
 * \return non-zero on internal error.
 */
int iso_viscous_damage_model_ctx_destroy(void **ctx);

/**
 * Provide a public interface to the integration algorithm to allow
 * other constitutive models to include damage as well.
 */
int ivd_public_int_alg(double *var_w,
                       double *var_X,
                       double *var_H,
                       int *flag_damaged,
                       const double var_wn,
                       const double var_Xn,
                       const double dt,
                       const double Ybar,
                       const double param_mu,
                       const double param_ome_max,
                       const double param_p1,
                       const double param_p2,
                       const double param_Yin);

#endif
