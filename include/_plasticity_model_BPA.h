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

/**
 * Private structure used exclsively with this model and associated
 * functions.
 */
typedef struct {
  double F[9];
} BPA_ctx;

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

/**
 * Evaluate the inverse Langevin function using a Pade approximation
 * for y in [0,1).
 *
 * \param[in] y in [0,1)
 * \param[out] inv_lang, function value in [0, inf)
 * \return non-zero if y outside bounds (resulting in undefined
 * behavior).
 */
int BPA_inverse_langevin(const double y,
                         double *inv_lang);

/**
 * Evaluate the derivative of the Pade approximation for the inverse
 * Langevin function for y in [0,1).
 *
 * \param[in] y in [0,1)
 * \param[out] der_inv_lang, derivative of the function at y
 * \return non-zero if y outside bounds (resulting in undefined
 * behavior).
 */
int BPA_der_inverse_langevin(const double y,
                             double  *der_inv_lang);

/**
 * Compute the plastic slip rate. NOTE: makes use of global values at
 * file scope within plasticity_mode_BPA.c.
 *
 * \param[out] gdot, the plastic slip rate
 * \param[in] tau, the equivelent plastic stress
 * \param[in] s_s, the pressure-dependent athermal shear stress
 * \return non-zero on internal error.
 */
int BPA_compute_gdot(double *gdot,
                     const double tau,
                     const double s_s);

/**
 * Compute the athermal shear hardenening rate. NOTE: makes use of
 * global values at file scope within plasticity_mode_BPA.c.
 *
 * \param[out] sdot, the hardening rate
 * \param[in] s, the current athermal shear stress
 * \param[in] gdot, the current slip rate
 * \return non-zero on internal error
 */
int BPA_compute_sdot(double *sdot,
                     const double s,
                     const double gdot);

/**
 * Compute the deviatoric backstress in the intermediate
 * configuration. NOTE: makes use of global values at file scope
 * within plasticity_mode_BPA.c.
 *
 * \param[out] Bdev, the plastic slip rate
 * \param[in] Fp, the (total) plastic deformation gradient
 * \return non-zero on internal error.
 */
int BPA_compute_Bdev(double *Bdev,
                     const double *Fp);

/**
 * Compute the plastic loading direction. Note that all pointers are
 * _restrict_ qualified for performace. Therefore, the memory pointed
 * to by each argument must not overlap. Undefined behavior will occur
 * if they do.
 *
 * \param[out] normal, the direction of the plastic loading
 * \param[out] eq_sig_dev, the equivalent plastic stress tensor
 * \param[out] tau, the magnitude of eq_sig_dev
 * \param[in] Sdev, the (elastic) deviatoric PK2 stress in the
 *                  intermediate configuration
 * \param[in] Bdev, the deviatoric back stress in the intermediate
 *                  configuration
 * \return non-zero on internal error
 */
int BPA_compute_loading_dir(double *normal,
                            double *eq_sig_dev,
                            double *tau,
                            const double *Sdev,
                            const double *Bdev,
                            const double *Fe);

#endif
