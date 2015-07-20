/**
 * Define the interface to generalized constitutive models with
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
#ifndef CONSTITUTIVE_MODEL_H
#define CONSTITUTIVE_MODEL_H

#include <stdlib.h>
#include "state_variables.h" /* provides declaration of Matrix_double */

/**
 * Enumeration for the model type
 */
enum model_type {
  HYPER_ELASTICITY,
  CRYSTAL_PLASTICITY,
  BPA_PLASTICITY
};

/**
 * Pre-declare the Model_parameters structure
 */
struct Model_parameters;
#ifndef TYPE_MODEL_PARAMETERS
#define TYPE_MODEL_PARAMETERS
typedef struct Model_parameters Model_parameters;
#endif

/**
 * General interface to a constitutive model.
 *
 * Holds a Model_parameters object.
 * Has a State_variables object.
 */
struct Constitutive_model {
  const Model_parameters *param;
  State_variables vars;
};

#ifndef TYPE_CONSTITUTIVE_MODEL
#define TYPE_CONSTITUTIVE_MODEL
typedef struct Constitutive_model Constitutive_model;
#endif

/**
 * Construct a Constitutive_model object. The object is left in an
 * invalid state until constitutive_model initialize is called.
 *
 * \return non-zero on error.
 */
int constitutive_model_construct(Constitutive_model *m);

/**
 * Initialize the Constitutive_model object given the material type
 * and properties. The object may be used after calling this function.
 *
 * \return non-zero on error.
 */
int constitutive_model_initialize(Constitutive_model *m,
                                  const Model_parameters *param);

/**
 * Destroy a Constitutive_model object.
 *
 * \return non-zero on error.
 */
int constitutive_model_destroy(Constitutive_model *m);

/**
 * User defined function for the Constitutive_model integration algorithm.
 *
 * \param[in,out] m - pointer to a Constitutive_model object.
 * \param[in,out] usr_ctx - handle to a user defined structure that
 *   controls execution of the user-defined integration algorithm. This
 *   is the mechanism for passing model-specific information into the
 *   general function interface.
 * 
 * \return non-zero on internal error that should be handled by the
 * calling function.
 */
typedef int (*usr_int_alg)(Constitutive_model *m,
                             const void *usr_ctx);

/**
 * User defined function to compute constitutive tensors.
 *
 * \param[in] m - pointer to Constitutive_model object.
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
typedef int (*usr_tensor)(const Constitutive_model *m,
                          const void *usr_ctx,
                          Matrix_double *tensor);

/**
 * User defined function to compute constitutive scalars.
 *
 * \param[in] m - pointer to Constitutive_model object.
 * \param[in,out] usr_ctx - handle to a user defined structure that
 *   controls execution of the user-defined function. This is the
 *   mechanism for passing model-specific information into the general
 *   function interface.
 * \param[out] scalar - scalar passed by reference
 *
 * \return non-zero value on internal error that should be handled by
 * the calling function.
 */
typedef int (*usr_scalar)(const Constitutive_model *m,
                          const void *usr_ctx,
                          double *scalar);

/**
 * User defined function that to increment/decrement the internally
 * stored state variables, i.e., advance: (n) <- (n + 1),
 * or reset: (n + 1) <- (n) the internal state variables.
 *
 * \param[in,out] m - pointer to Constitutive_model object.
 * \return non-zero on internal error that should be handled by the
 * calling function.
 */
typedef int (*usr_increment)(Constitutive_model *m);

/** Pre-declare MATERIAL structure */
struct MATERIAL;
#ifndef TYPE_MATERIAL
#define TYPE_MATERIAL
typedef struct MATERIAL MATERIAL;
#endif

/** Pre-declare MATERIAL structure */
struct MATGEOM_1;
#ifndef TYPE_MATGEOM_1
#define TYPE_MATGEOM_1
typedef struct MATGEOM_1 MATGEOM_1;
#endif

/** Pre-declare HOMMAT structure */
struct HOMMAT;
#ifndef TYPE_HOMMAT
#define TYPE_HOMMAT
typedef struct HOMMAT HOMMAT;
#endif

/**
 * Object for querying/describing the state variables.
 */
struct Model_var_info {
  char **F_names;
  char **var_names;
  size_t n_Fs;
  size_t n_vars;
};

#ifndef TYPE_MODEL_VAR_INFO
#define TYPE_MODEL_VAR_INFO
typedef struct Model_var_info Model_var_info;
#endif

/**
 * destroy a Model_var_info object. Assumes full control of all
 * internal pointers.
 *
 * \return non-zero on error
 */
int model_var_info_destroy(Model_var_info **info);

/**
 * A user described function that allocates a Model_var_info object
 * and populates the internal structure. This function should fully
 * allocate the internals of Model_var_info and *copy* data rather
 * than hold references.
 */
typedef int (*usr_info)(Model_var_info **info);

/**
 * Interface for accessing model parameters and modifying/updating the
 * associated state variable(s) at integration points.
 */
struct Model_parameters {
  /** Pointer to anisotropic material properties (props,orientation, etc.) */
  const MATERIAL *p_mat;
  /** Pointer to anisotropic material properties (props,orientation, etc.) */
  const MATGEOM_1  *p_mgeom;  
  /** Pointer to isotropic material props */
  const HOMMAT *p_hmat;

  /** access to user-defined functions */
  usr_int_alg integration_algorithm;
  usr_tensor compute_dev_stress;
  usr_scalar compute_dudj;
  usr_tensor compute_dev_tangent;
  usr_scalar compute_d2udj2;
  usr_increment update_state_vars;
  usr_increment reset_state_vars;
  usr_info get_var_info;
  /** Model type, see enumeration @model_type */
  size_t type;
};

/**
 * Construct a Model_parameters object. The object is left in an
 * invalid state until model_parameters_initialize is called.
 *
 * \return non-zero on error.
 */
int model_parameters_construct(Model_parameters *p);

/**
 * Initialize the Model_parameters object. The object may be used
 * after calling this function. Calling this function on an already
 * initialized object is undefined.
 *
 * \return non-zero on error.
 */
int model_parameters_initialize(Model_parameters *p,
                                const MATERIAL *p_mat,
                                const MATGEOM_1 *p_mgeom,
                                const HOMMAT *p_hmat,
                                const size_t type);

/**
 * Destroy a Model_parameters object.
 *
 * \return non-zero on error.
 */
int model_parameters_destroy(Model_parameters *p);

int constitutive_model_update_elasticity(Constitutive_model *m, Matrix_double *Fe, double dt,
                               Matrix_double *L, Matrix_double *S, int compute_stiffness);
// if compute_stiffness == 1: compute stiffness (L)
// if compute_stiffness == 0: just stress (S) is updated and stiffness (L) is not computed

int constitutive_model_update_plasticity(Matrix_double *pFnp1,Matrix_double *eFn, Matrix_double *pFn, Constitutive_model *m, double dt);

int constitutive_model_update_dMdu(Constitutive_model *m, Matrix_double *dMdu, Matrix_double *Fe, Matrix_double *S, Matrix_double *L, Matrix_double *Grad_du, double dt);

#endif
