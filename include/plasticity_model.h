/**
 * Declare/define functions that integrate with the plasticty model
 * interface (see constitutive_model.h). 
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 *  Sangmin Lee, University of Notre Dame, Notre Dame, IN, <slee43@nd.edu>
 */
#pragma once
#ifndef PLASTICITY_MODEL_H
#define PLASTICITY_MODEL_H

struct Model_parameters;
#ifndef TYPE_MODEL_PARAMETERS
#define TYPE_MODEL_PARAMETERS
typedef struct Model_parameters Model_parameters;
#endif

#ifndef TYPE_CONSTITUTIVE_MODEL
#define TYPE_CONSTITUTIVE_MODEL
typedef struct Constitutive_model Constitutive_model;
#endif

struct Matrix_double;
#ifndef TYPE_MATRIX_DOUBLE
#define TYPE_MATRIX_DOUBLE
typedef struct Matrix_double Matrix_double;
#endif

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Constitutive_model object
 * \return non-zero on internal error
 */
int plasticity_model_initialize(Model_parameters *p);

int plasticity_model_destory(Model_parameters *p);

/**
 * Construct and initialize the model context for calling functions
 * through the plasticity interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F - The *total* deformation gradient.
 * \return non-zero on internal error.
 */
int plasticity_model_ctx_build(void **ctx,
                               const double *F,
                               const double dt,
                               const double alpha,
                               const double *eFnpa);

/**
 * Destroy the model context and invalidate the handle.
 *
 * \param[in,out] ctx - handle to context object. ctx = NULL on exit.
 * \return non-zero on internal error.
 */
int plasticity_model_ctx_destroy(void **ctx);

/**
 * compute tangent of the plasticity part of deformation gradient with respect to the deformation
 *
 * \param[in,out] dMdu - computed tangent of plasticity part of deformation gradient
 * \param[in] du - Grad(du)
 * \param[in,out] m - pointer to a Constitutive_model object.
 * \param[in] dt - time step size
 * \return non-zero on internal error.
 */
int plasticity_model_slip_system(Matrix_double *P);

#ifndef ELEMENT_H
typedef struct EPS EPS;
#endif

#ifndef EPS_H
typedef struct ELEMENT ELEMENT;
#endif

int plasticity_model_set_orientations(EPS *eps,
                                const int ne,
                                const ELEMENT *elem,
                                const int n_mat,
                                const Model_parameters *param_list);
/** read material properties for plasticity
 * need to provide model_params.in with format as below
 *
 *------------------------------------------------------------------------------------------------
 * model_params.in
 *------------------------------------------------------------------------------------------------ 
 * # <= this denotes this line is a comment. It will be ignored.
 * # Number of material
 * 2
 * ########################################################################
 * # material 0
 * ########################################################################
 * # material_ID meter_scaling_factor
 * 0 1.0
 * # properties
 * # num_of_properties gamma_dot_0    m    G0    g0  gs_0 gamma_dot_s     w
 *                   7         1.0 0.05 200.0 210.0 330.0     50.0e+9 0.005
 * #
 * ######################################################################## 
 * # unit_cell orientation
 * ########################################################################
 * # 0: FCC
 * # 1: BCC not implemented
 * # 2: HCP not implemented
 * #
 * # -1: no orientation is used
 * # 0: random - each element will have random orientation using built in function
 * #             if 0 is followed, integration points in a element will have same orientation
 * #             if 1 is followed, integration points in a element will have different orientations
 * # 1: crystals - each crystal will have random orientation using built in function
 * # 2: file - orientation is givne by a file, need to provide file path with part of file name
 * #           path/orientation where path has files with name as orientation_*.in
 * # 3: provide material orientation directly
 * #    e.g) 0 3 0.1 0.1 0.1
 * 0 2 CRYSTAL_ORIENTATION/orientation
 */

typedef struct HOMMAT HOMMAT;

int plasticity_model_test(const HOMMAT *hmat, Matrix_double *L_in, int Load_Type);
void test_crystal_plasticity_single_crystal(void); 
#endif
