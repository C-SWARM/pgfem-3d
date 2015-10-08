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

/**
 * Construct and initialize the model context for calling functions
 * through the plasticity interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F - The *total* deformation gradient.
 * \return non-zero on internal error.
 *
 * CAVEATES: The addresses of C and J_or_Theta must remain valid
 * throughout the existence of the ctx. Destroying these memory
 * locations before calling plasticity_model_none_ctx_destroy may
 * invalidate ctx.
 */
int plasticity_model_ctx_build(void **ctx,
                               const double *F,
                               const double dt,
                               const double alpha);

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

int plasticity_model_integration_ip(Matrix_double *pFnp1,
                                    Constitutive_model *m,
                                    const Matrix_double *Fnp1,
                                    const Matrix_double *Fe_n,
                                    const double dt);

#ifndef ELEMENT_H
typedef struct EPS EPS;
#endif

#ifndef EPS_H
typedef struct ELEMENT ELEMENT;
#endif

int plasticity_model_read_parameters(EPS *eps,
                                       const int ne,
                                       const ELEMENT *elem,
                                       const int n_mat,
                                       Model_parameters *param_list);
/** read material properties for plasticity
 * need to provide MATERIAL_PROPERTY.in with format as below
 *
 *------------------------------------------------------------------------------------------------
 * MATERIAL.in
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
 * # orientation
 * # -1: no orientation is used
 * # 0: random - each element will have random orientation using built in function
 * #             if 0 is followed, integration points in a element will have same orientation
 * #             if 1 is followed, integration points in a element will have different orientations
 * # 1: crystals - each crystal will have random orientation using built in function
 * # 2: file - orientation is givne by a file, need to provide file path with part of file name
 * #           path/orientation where path has files with name as orientation_*.in
 * # 3: provide material orientation directly
 * #    e.g) 3 0.1 0.1 0.1
 * 2 CRYSTAL_ORIENTATION/orientation
 * ########################################################################
 * # material 1
 * ########################################################################
 * # Material properties can be read from other material card as:
 * # # to read from other material card, use -1 followed by material file name
 * # 0.001 meter_scaling_factor is used in case if your domain unit is [mm], 
 * # because all material properties in MATERIAL_DATA are in [m]
 * 1 0.001
 * -1 /scratch365/cswarm/MATERIAL_DATA/MATERIAL_ALUMINUM-1100
 * -1
 */

typedef struct HOMMAT HOMMAT;

int plasticity_model_test(const HOMMAT *hmat, Matrix_double *L_in, int Load_Type);
 
#endif
