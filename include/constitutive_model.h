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
#include "matgeom.h"
#include "supp.h"

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
 * Initialize the constitutive model object at each integration point.
 *
 * \return non-zero on error.
 */
int init_all_constitutive_model(EPS *eps,
                                const int ne,
                                const ELEMENT *elem,
                                const Model_parameters *param_list);


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
int constitutive_model_update_elasticity(const Constitutive_model *m,
                                         const Matrix_double *F,
                                         const double dt,
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
  
  Matrix_double *Psys;
  int N_SYS;

  /** access to user-defined functions */
  usr_int_alg integration_algorithm;
  usr_tensor compute_dev_stress;
  usr_scalar compute_dudj;
  usr_tensor compute_dev_tangent;
  usr_scalar compute_d2udj2;
  usr_increment update_state_vars;
  usr_increment reset_state_vars;
  usr_info get_var_info;
  usr_get_F get_Fn;
  usr_get_F get_pF;
  usr_get_F get_pFn;
  usr_get_F get_eF;
  usr_get_F get_eFn;
  usr_destroy_ctx destroy_ctx;
  usr_compute_dM_du compute_dMdu;

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
 * function for reading extra prameters for plasticity, this is a temporal.
 * decision needs to be made how to pass material properties
 * \param[in/out] param_list, list of Model_parameters, unallocated on
 *                            entry -- allocated on exit.
 * \param[in] eps, pointer to strains EPS
 * \param[in] ne, number of elements
 * \param[in] elem, pointer to elements
 * \param[in] param_list, list of constitutive model parameters
 *
 * \return non-zero on error.
 */
int read_constitutive_model_parameters(EPS *eps, 
                                       const int ne, 
                                       const ELEMENT *elem, 
                                       const int n_mat, 
                                       Model_parameters *param_list,
                                       const int type); 

/**
 * Destroy a Model_parameters object.
 *
 * \return non-zero on error.
 */
int model_parameters_destroy(Model_parameters *p);

/**
 * Allocate and populate a list of Model_parameters given the number
 * of materials. It seems that only one MATGEOM_1 exists, no matter
 * how many orientations are present, and that MATERIAL is not used by
 * any of the current implementation. Therefore the length of the list
 * is n_hmat.
 *
 * \param[in/out] param_list, list of Model_parameters, unallocated on
 *                            entry -- allocated on exit.
 * \param[in] n_mat, length oc hmat_list
 * \param[in] p_mgeom, pointer to material geometry object
 * \param[in] hmat_list, list of homogenized material properies
 *
 * \return non-zero on error.
 */
int build_model_parameters_list(Model_parameters **param_list,
                                const int n_mat,
                                const MATGEOM_1 *p_mgeom,
                                const HOMMAT *hmat_list,
                                const int type);

/**
 * Free all of the memory assiciated with the list of model
 * parameters.
 *
 * \return non-zero on error.
 */
int destroy_model_parameters_list(const int n_mat,
                                  Model_parameters *param_list);

int constitutive_model_update_time_steps(EPS *eps, const int ne, const ELEMENT *elem);
/**
 * update values for next time step: variables[tn] = variables[tn+1]
 * \return non-zero on error.
 */
int constitutive_model_update_time_steps_test(ELEMENT *elem, NODE *node, HOMMAT *hommat, EPS *eps, 
                                        const int ne, const int nn, const int ndofn,
                                        double* r, double dt); 

int constitutive_model_test(const HOMMAT *hmat, Matrix_double *L_in, int Print_results);

int stiffness_el_hyper_elasticity(double *lk,const int ii,const int ndofn,const int nne,const int nsd,
        const ELEMENT *elem,const HOMMAT *hommat,MATGEOM matgeom,const long *nod,const NODE *node,
        double dt,SIG *sig,EPS *eps,const SUPP sup,double *r_e);
        
int residuals_el_hyper_elasticity(double *f,const int ii,const int ndofn,const int nne,const int nsd,
        const ELEMENT *elem,const HOMMAT *hommat,MATGEOM matgeom,const long *nod,const NODE *node,
        double dt,SIG *sig,EPS *eps,const SUPP sup,double *r_e);

int stiffness_el_crystal_plasticity(double *lk,const int ii,const int ndofn,const int nne,const int nsd,
        const ELEMENT *elem,const HOMMAT *hommat,MATGEOM matgeom,const long *nod,const NODE *node,
        double dt,SIG *sig,EPS *eps,const SUPP sup,double *r_e);
        
int residuals_el_crystal_plasticity(double *f,const int ii,const int ndofn,const int nne,const int nsd,
        const ELEMENT *elem,const HOMMAT *hommat,MATGEOM matgeom,const long *nod,const NODE *node,
        double dt,SIG *sig,EPS *eps,const SUPP sup,double *r_e);
        
#endif
