/**
 * plasticity.h
 *
 * Define the interface to generalized integration algorithms with
 * their associated data structure(s). The user is responsible for
 * defining and linking the functions associated with a particular
 * model. NEED MORE USAGE DETAILS.
 *
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */
#pragma once
#ifndef PLASTICITY_H
#define PLASTICITY_H

#include <stdlib.h>

/** Pre-declare MATERIAL structure */
struct MATERIAL;
#ifndef TYPE_MATERIAL
#define TYPE_MATERIAL
typedef struct MATERIAL MATERIAL;
#endif

/** Pre-declare HOMMAT structure */
struct HOMMAT;
#ifndef TYPE_HOMMAT
#define TYPE_HOMMAT
typedef struct HOMMAT HOMMAT;
#endif

/** Handle to a matrix object */
struct Matrix_double;
#ifndef TYPE_MATRIX_HANDLE
#define TYPE_MATRIX_HANDLE
typedef struct Matrix_double Matrix_handle;
#endif

/** Handle to a vector object */
#ifndef TYPE_VECTOR_HANDLE
#define TYPE_VECTOR_HANDLE
typedef struct Matrix_double Vector_handle;
#endif

/** Pre-declare and typedef Plasticity structure */
struct Plasticity;
#ifndef TYPE_PLASTICITY
#define TYPE_PLASTICITY
typedef struct Plasticity Plasticity;
#endif

/**
 * User defined function to modify the Plasticity object.
 *
 * \param[in,out] p - Pointer to a plasticity object.
 * \param[in,out] usr_ctx - handle to a user defined structure that
 *   controls execution of the user-defined function. This is the
 *   mechanism for passing model-specific information into the general
 *   function interface.
 * 
 * \return non-zero on internal error that should be handled by the
 * calling function.
 */
typedef int (*usr_func)(Plasticity *p,
                        const void *usr_ctx);

/**
 * User defined function to compute constitutive tensors.
 *
 * \param[in] p - pointer to Plasticity object.
 * \param[in,out] usr_ctx - handle to a user defined structure that
 *   controls execution of the user-defined function. This is the
 *   mechanism for passing model-specific information into the general
 *   function interface.
 * \param[out] tensor - pointer to a Matrix object (treated as a tensor)
 *   that is populated with values from the constitutive model.
 *
 * \return non-zero value on internal error that should be handled by
 * the calling function.
 */
typedef int (*usr_tensor)(const Plasticity *p,
                          const void *usr_ctx,
                          Matrix_handle *tensor);

/**
 * User defined function to compute constitutive scalars.
 *
 * \param[in] p - pointer to Plasticity object.
 * \param[in,out] usr_ctx - handle to a user defined structure that
 *   controls execution of the user-defined function. This is the
 *   mechanism for passing model-specific information into the general
 *   function interface.
 * \param[out] scalar - scalar passed by reference
 *
 * \return non-zero value on internal error that should be handled by
 * the calling function.
 */
typedef int (*usr_scalar)(const Plasticity *p,
                          const void *usr_ctx,
                          double *scalar);

/**
 * User defined function that to increment/decrement the internally
 * stored state variables, i.e., advance: (n) <- (n + 1),
 * or reset: (n + 1) <- (n) the internal state variables.
 *
 * \param[in,out] p - pointer to Plasticity object.
 * \return non-zero on internal error that should be handled by the
 * calling function.
 */
typedef int (*usr_increment)(Plasticity *p);

/**
 * Generalized interface for integrating a constitutive model.
 */
struct Plasticity {
  Matrix_handle *Fs;        /*< Array of handles to deformation
                              gradients, e.g., Fp, Ft, Fe,... */
  Vector_handle *state_vars; /*< Handle to vector of state variables */
  const MATERIAL *p_mat;    /*< Pointer to a material (props,
                              orientation, etc.)*/
  const HOMMAT *p_hmat;     /*< Pointer to add'l material props */

  usr_func integration_algorithm;  /*< integrate the constitutive
                                     model (compute n + 1 state
                                     variables) */
  usr_tensor compute_dev_stress;   /*< compute the (modified)
                                     deviatoric stress */
  usr_scalar compute_pressure;   /*< compute the (modified)
                                     pressure */ 
  usr_tensor compute_dev_tangent;  /*< compute the (modified)
                                     deviatoric tangent */
  usr_scalar compute_pressure_tangent;  /*< compute the (modified)
                                          pressure tangent */
  usr_increment update_state_vars; /*< state variables
                                     (n) <- (n + 1) */
  usr_increment reset_state_vars; /*< state variables
                                    (n + 1) <- (n) */

  size_t n_Fs; /*< Number of handles to deformation gradients. */
  size_t type; /*< Model type, see enumeration @model_type */
};

/**
 * Enumeration for the model type
 */
enum model_type {
  NONE,
  CRYSTAL_PLASTICITY,
  BPA_PLASTICITY
};

/**
 * Construct a plasticity object. The object is left in an invalid
 * state until plasticity initialize is called.
 *
 * \return non-zero on error.
 */
int plasticity_construct(Plasticity *p);

/**
 * Initialize the plasticity object given the material type and
 * properties. The object may be used after calling this function.
 *
 * \return non-zero on error.
 */
int plasticity_initialize(Plasticity *p,
                          const MATERIAL *p_mat,
                          const HOMMAT *p_hmat,
                          const size_t type);

/**
 * Destroy a plasticity object.
 *
 * \return non-zero on error.
 */
int plasticity_destroy(Plasticity *p);

#endif
