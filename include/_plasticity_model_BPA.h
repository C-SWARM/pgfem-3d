/**
 * This file defines the _private_ interface to the BPA plasticity
 * model. It is intended to provide greater access to driver programs,
 * particularly during the testing/development stages. In general,
 * this file should only be included in the BPA plasticity
 * implementation code.
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */

#pragma once
#ifndef _PLASTICITY_MODEL_BPA_H
#define _PLASTICITY_MODEL_BPA_H

#ifndef TYPE_MATRIX_DOUBLE
#define TYPE_MATRIX_DOUBLE
typedef struct Matrix_double Matrix_double;
#endif

#ifndef TYPE_CONSTITUTIVE_MODEL
#define TYPE_CONSTITUTIVE_MODEL
typedef struct Constitutive_model Constitutive_model;
#endif

#ifndef TYPE_MODEL_VAR_INFO
#define TYPE_MODEL_VAR_INFO
typedef struct Model_var_info Model_var_info;
#endif

#ifndef TYPE_HOMMAT
#define TYPE_HOMMAT
typedef struct HOMMAT HOMMAT;
#endif

/**
 * Private structure used exclsively with this model and associated
 * functions.
 */
typedef struct {
  double F[9];
  double dt;
} BPA_ctx;


/* The following functions are accessed throught the constitutive
   modeling interface. Direct access is provided here for driver
   programs/debugging purposes. */
/**
 * The integration algorithm for the BPA plasticity model.
 */
int BPA_int_alg(Constitutive_model *m,
                const void *ctx);

/**
 * Compute the deviatoric PK2 stress for the BPA model.
 */
int BPA_dev_stress(const Constitutive_model *m,
                   const void *ctx,
                   Matrix_double *dev_stress);

/**
 * Compute the pressure term for the BPA model.
 */
int BPA_dudj(const Constitutive_model *m,
             const void *ctx,
             double *dudj);

/**
 * Compute the deviatoric material tangent for the BPA model.
 */
int BPA_dev_tangent(const Constitutive_model *m,
                    const void *ctx,
                    Matrix_double *dev_tangent);

/**
 * Compute the pressure tangent for the BPA model.
 */
int BPA_d2udj2(const Constitutive_model *m,
               const void *ctx,
               double *d2udj2);

/**
 * Update the BPA state variables (n <- n+1)
 */
int BPA_update_vars(Constitutive_model *m);

/**
 * reset the BPA state variables (n+1 <- n)
 */
int BPA_reset_vars(Constitutive_model *m);

/**
 * Construct and populate a model info object for the BPA model.
 */
int BPA_model_info(Model_var_info **info);

#endif
