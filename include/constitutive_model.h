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
#include <stdio.h>
#include "state_variables.h" /* provides declaration of Matrix_double */
#include "sig.h"
#include "supp.h"
#include "PGFem3D_options.h"

typedef struct EPS EPS;
typedef struct ELEMENT ELEMENT;
typedef struct NODE NODE;


#ifndef PGFEM3D_DEV_TEST
#define PGFEM3D_DEV_TEST 1
#endif

/**
 * Enumeration for the model type
 */
enum model_type {
  HYPER_ELASTICITY,
  CRYSTAL_PLASTICITY,
  BPA_PLASTICITY,
  ISO_VISCOUS_DAMAGE,
  J2_PLASTICITY_DAMAGE,
  NUM_MODELS,
  TESTING=99,
  CM_UQCM=100
};

/**
 * Enumeration for the frame
 */
enum integration_frame {
  UPDATED_LAGRANGIAN,
  TOTAL_LAGRANGIAN,
  MIXED_ANALYSIS_MODE
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
 * Initialize the constitutive model object at each integration point.
 *
 * \return non-zero on error.
 */
int init_all_constitutive_model(EPS *eps,
                                const int ne,
                                const ELEMENT *elem,
                                const int n_mat,
                                const Model_parameters *param_list);

/**
 * Reset the cinstitutive model at each integration point in the
 * domain by calling its respective reset function.
 *
 * \return non-zero on error.
 */
int constitutive_model_reset_state(EPS *eps,
                                   const int ne,
                                   const ELEMENT *elem);

/**
 * Compute the elastic stress and tangent tensors. Note that the total
 * stress and tangent (not only the deviatoric parts) are returned.
 *
 * \param[in] F, *total* deformation gradient
 * \param[in] compute_stiffness, flag to comptue the stiffness tensor
 * (L) if non-zero. Otherwise L is unchanged.
 *
 * \return non-zero on internal error.
 */ 

int constitutive_model_defaut_update_elasticity(const Constitutive_model *m,
                                                Matrix_double *eF,
                                                Matrix_double *L,
                                                Matrix_double *S,
                                                const int compute_stiffness);

typedef int (*usr_update_elasticity) (const Constitutive_model *m,
                                      const void *ctx,
                                      Matrix_double *L,
                                      Matrix_double *S,
                                      const int compute_stiffness);


/**
 * User defined function for the Constitutive_model integration
 * algorithm. This function shall be implemented such that it modifies
 * the internal state to contain the updated values upon exit, i.e.,
 * subsequent calls to usr_get_F functions will return the correct
 * values without re-integrating the constitutive model.
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

/**
 * User defined function to return a state variable. Note that this
 * function *DOES NOT* modify the internal state variables. This
 * function shall be implemented such that destroying the returned
 * variable does not modify the internal state of the
 * Constitutive_model object.
 *
 * \param[in] m - constant reference to a Constitiutive_model object.
 * \param[out] var - contains the value of the state variable upon
 *                   exit.
 * \return non-zero on internal error
 */
typedef int (*usr_get_var)(const Constitutive_model *m,
                           double *var);

/**
 * User defined function to return a time dependent state variable. 
 * Note that this
 * function *DOES NOT* modify the internal state variables. This
 * function shall be implemented such that destroying the returned
 * variable does not modify the internal state of the
 * Constitutive_model object.
 *
 * \param[in] m - constant reference to a Constitiutive_model object.
 * \param[in] t - constant t (time). 
 * \param[out] var - contains the value of the state variable upon
 *                   exit.
 * \return non-zero on internal error
 */
typedef int (*usr_get_var_of_t)(const Constitutive_model *m,
                           double *var,
                           const double t);
                           
/**
 * User defined function to return the deformation gradient. Note that
 * this function *DOES NOT* modify the internal state variables. This
 * function shall be implemented such that destroying the returned
 * deformation gradient does not modify the internal state of the
 * Constitutive_model object.
 *
 * \param[in] m - constant reference to a Constitiutive_model object.
 * \param[out] F - reference to Matrix object that contains the
 *                 deformation gradient upon exit.
 * \return non-zero on internal error
 */
typedef int (*usr_get_F)(const Constitutive_model *m,
                         Matrix_double *F);

/**
 * User defined function to destroy a context for the model. This
 * function shall destroy any internally allocated data and
 * invalidate the handle, i.e., *ctx = NULL.
 */
typedef int (*usr_destroy_ctx)(void **ctx);

/**
 * User defined function to compute the linearization of the plastic
 * deformation w.r.t. the displacement variable. In the current
 * formulations, the formulation for the PDE is in terms of M =
 * inv(pFr). Thus the linearization is dM_du. ***This may be changed
 * in the future, deprecating this function***
 *
 * \param[in] m, Constitutive model object
 * \param[in] ctx, model specific user context
 * \param[in] Grad_op, 4-index FE gradient operator where the 1st two
 *                     indices are the node followed by the DOF.
 * \param[in] nne, number of nodes on the element (length of idx 1 in
 *                 Grad_op)
 * \param[in] ndofn, number of dofs on each node (length of idx 2 in
 *                  Grad_op)
 * \param[out] dM_du, 4-index FE linearization of M, same dimensions
 *                   as Grad_op
 *
 * \return non-zero on internal error
 */
typedef int (*usr_compute_dM_du)(const Constitutive_model *m,
                                 const void *ctx,
                                 const double *Grad_op,
                                 const int nne,
                                 const int ndofn,
                                 double *dM_du);

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
  char **flag_names;
  size_t n_Fs;
  size_t n_vars;
  size_t n_flags;
};

#ifndef TYPE_MODEL_VAR_INFO
#define TYPE_MODEL_VAR_INFO
typedef struct Model_var_info Model_var_info;
#endif

/**
 * Print the object to the specified file.
 */
int model_var_info_print(FILE *f,
                         const Model_var_info * info);

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
 * A user described function that writes and read restart file at
 * integration point
 */
typedef int (*usr_w_restart)(FILE *fp,
                             const Constitutive_model *m);
typedef int (*usr_r_restart)(FILE *fp,
                             Constitutive_model *m);

/**
 * User defined function to set the initial values of the state
 * variables for the particular model.
 */
typedef int (*usr_set_init_vals)(Constitutive_model *m);

/**
 * User defined function for reading in material parameters from a
 * file.
 */
typedef int (*usr_read_params)(Model_parameters *p,
                               FILE *in);

/**
 * User defined function that returns the size of the data to be
 * packed/unpacked.
 *
 * Does not modify the CM object or any of the data it holds.
 * \return size in bytes of the pack/unpack data
 */
typedef size_t (*usr_get_size)(const Constitutive_model *m);

/**
 * User defined function to pack the CM data into a buffer (see
 * pack_data).
 *
 * Does not modify the CM object or any of the data it holds.
 * \param[in,out] buffer, a buffer to insert data to
 *
 * \param[in,out] pos, insert position in the buffer. Upon exit - next
 *                     insertion position.
 * \return non-zero on error.
 */
typedef int (*usr_pack)(const Constitutive_model *m,
                        char *buffer,
                        size_t *pos);

/**
 * User defined function to unpack CM data from a buffer (see also
 * usr_pack, unpack_data).
 *
 * \param[out] m, CM object with internal data set from the buffer
 * \param[in] buffer, the buffer to read data from

 * \param[in,out] pos, the position in buffer to begin reading from.
 *                     Upon exit - position for next read.
 * \return non-zero on error.
 */
typedef int (*usr_unpack)(Constitutive_model *m,
                          const char *buffer,
                          size_t *pos);

/**
 * Interface for accessing model parameters and modifying/updating the
 * associated state variable(s) at integration points.
 */
 
struct MATERIAL_CONSTITUTIVE_MODEL;
#ifndef TYPE_MATERIAL_CONSTITUTIVE_MODEL
#define TYPE_MATERIAL_CONSTITUTIVE_MODEL
typedef struct MATERIAL_CONSTITUTIVE_MODEL MATERIAL_CONSTITUTIVE_MODEL;
#endif

struct ELASTICITY;
#ifndef TYPE_ELASTICITY
#define TYPE_ELASTICITY
typedef struct ELASTICITY ELASTICITY;
#endif
 
struct Model_parameters {
  /** Pointer to isotropic material props */
  const HOMMAT *p_hmat;
  int mat_id; // Global material id, mat_id may not be the same as the hommat id
  int uqcm;   // UQ study through constitutive model 0: no, or yes
  
  MATERIAL_CONSTITUTIVE_MODEL *cm_mat; 
  ELASTICITY *cm_elast;

  /** access to user-defined functions */
  usr_int_alg integration_algorithm;
  usr_tensor compute_dev_stress;
  usr_scalar compute_dudj;
  usr_tensor compute_dev_tangent;
  usr_scalar compute_d2udj2;

  /* compute the elastic algorithmic stiffness tangent */
  usr_tensor compute_AST;
  usr_update_elasticity update_elasticity;

  usr_increment update_state_vars;
  usr_increment reset_state_vars;
  usr_info get_var_info;
  usr_get_F get_Fn;
  usr_get_F get_Fnm1;  
  usr_get_F get_pF;
  usr_get_F get_pFn;
  usr_get_F get_pFnm1;
  usr_get_F get_eF;
  usr_get_F get_eFn;
  usr_get_F get_eFnm1;
    
  usr_get_var get_hardening;
  usr_get_var get_hardening_nm1;
  usr_get_var get_plast_strain_var;
  usr_get_var_of_t get_subdiv_param;
  
  usr_w_restart write_restart;
  usr_r_restart read_restart;

  usr_destroy_ctx destroy_ctx;
  usr_compute_dM_du compute_dMdu;

  usr_set_init_vals set_init_vals;
  usr_read_params read_param;

  usr_get_size get_size;
  usr_pack pack;
  usr_unpack unpack;

  /** Model type, see enumeration @model_type */
  size_t type;

  /* array for storing the model constants. */
  size_t n_param;
  double *model_param;
  /* array for storing the model integer type constant */
  size_t n_param_index;
  int *model_param_index;
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
                                const HOMMAT *p_hmat,
                                const size_t type);


/**
 * Destroy a Model_parameters object.
 *
 * \return non-zero on error.
 */
int model_parameters_destroy(Model_parameters *p);


/**
 * Allocate and populate a list of Model_parameters given the number
 * of materials.
 */
int read_model_parameters_list(Model_parameters **param_list,
                               const int n_mat,
                               const HOMMAT *hmat_list,
                               FILE *in);
/**
 * Free all of the memory assiciated with the list of model
 * parameters.
 *
 * \return non-zero on error.
 */
int destroy_model_parameters_list(const int n_mat,
                                  Model_parameters *param_list);

/**
 * update values for next time step: variables[tn] = variables[tn+1]
 * \return non-zero on error.
 */
int constitutive_model_update_time_steps_test(const ELEMENT *elem,
                                              NODE *node,
                                              EPS *eps,
                                              const int ne,
                                              const int nn,
                                              const int ndofn,
                                              const double* r,
                                              const double dt,
                                              const int total_Lagrangian);

int constitutive_model_test(const HOMMAT *hmat,
                            Matrix_double *L_in,
                            int Print_results);

int stiffness_el_hyper_elasticity(double *lk,
                                  const int ii,
                                  const int ndofn,
                                  const int nne,
                                  const int nsd,
                                  const ELEMENT *elem,
                                  const long *nod,
                                  const NODE *node,
                                  const double dt,
                                  EPS *eps,
                                  const SUPP sup,
                                  const double *r_e);
        
int residuals_el_hyper_elasticity(double *f,
                                  const int ii,
                                  const int ndofn,
                                  const int nne,
                                  const int nsd,
                                  const ELEMENT *elem,
                                  const long *nod,
                                  const NODE *node,
                                  const double dt,
                                  EPS *eps,
                                  const SUPP sup,
                                  const double *r_e);

int stiffness_el_crystal_plasticity(double *lk,
                                    const int ii,
                                    const int ndofn,
                                    const int nne,
                                    const int nsd,
                                    const ELEMENT *elem,
                                    const long *nod,
                                    const NODE *node,
                                    const double dt,
                                    EPS *eps,
                                    const SUPP sup,
                                    const double *r_e,
                                    const int total_Lagrangian);
        
int residuals_el_crystal_plasticity(double *f,
                                    const int ii,
                                    const int ndofn,
                                    const int nne,
                                    const int nsd,
                                    const ELEMENT *elem,
                                    const long *nod,
                                    const NODE *node,
                                    const double dt,
                                    EPS *eps,
                                    const SUPP sup,
                                    const double *r_e,
                                    const int total_Lagrangian);

int constitutive_model_update_output_variables(SIG *sig,
                                               EPS *eps,
                                               NODE *node,
                                               ELEMENT *elem,
                                               const int ne,
                                               const double dt,
                                               PGFem3D_opt *opts,
                                               double alpha);
                                               
int stiffness_el_crystal_plasticity_w_inertia(double *lk,
                                              const int ii,
                                              const int ndofn,
                                              const int nne,
                                              const int nsd,
                                              const ELEMENT *elem,
                                              const long *nod,
                                              const NODE *node,
                                              const double dt,
                                              EPS *eps,
                                              const SUPP sup,
                                              const double *r_e,
                                              double alpha);
                                    
int residuals_el_crystal_plasticity_w_inertia(double *f,
                                              const int ii,
                                              const int ndofn,
                                              const int nne,
                                              const int nsd,
                                              const ELEMENT *elem,
                                              const long *nod,
                                              const NODE *node,
                                              const double *dts,
                                              EPS *eps,
                                              const SUPP sup,
                                              const double *r_e,
                                              const double alpha);

/**
 * Compute the physics-based subdivision paramter for all integration
 * points on the domain.
 */
int cm_get_subdivision_parameter(double *subdiv_param,
                                 const int ne,
                                 const ELEMENT *elem,
                                 const EPS *eps,
                                 const double dt);

/**
 * Construct the model context for any model.
 */
int construct_model_context(void **ctx,
                                   const int type,
                                   const double *F,
                                   const double dt,
                                   const double alpha,
                                   const double *eFnpa);
#endif
