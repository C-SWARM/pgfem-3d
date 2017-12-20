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

#include "constitutive_model.h"
#include "cm_placeholder_functions.h"

class CP_PARAM: public Model_parameters
{
 public:

  virtual int model_dependent_initialization(void);
  virtual int model_dependent_finalization(void);

  virtual int integration_algorithm(Constitutive_model *m,
                                    const void *usr_ctx) const;
  virtual int compute_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 double *S) const;
  virtual int compute_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *value) const;
  virtual int compute_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  double *L) const;
  virtual int compute_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *value) const;
  virtual int update_elasticity(const Constitutive_model *m,
                                const void *ctx,
                                double *L,
                                double *S,
                                const int compute_stiffness) const;
  virtual int update_state_vars(Constitutive_model *m) const;
  virtual int reset_state_vars(Constitutive_model *m) const;
  virtual int reset_state_vars_using_temporal(const Constitutive_model *m,
                                              State_variables *var) const;
  virtual int update_np1_state_vars_to_temporal(const Constitutive_model *m,
                                                State_variables *var) const;
  virtual int save_state_vars_to_temporal(const Constitutive_model *m,
                                          State_variables *var) const;
  virtual int get_var_info(Model_var_info &info) const;
  virtual int get_F(const Constitutive_model *m,
                    double *F,
                    const int stepno) const;
  virtual int get_pF(const Constitutive_model *m,
                     double *F,
                     const int stepno) const;
  virtual int get_eF(const Constitutive_model *m,
                     double *F,
                     const int stepno) const;
  virtual int get_eF_of_hF(const Constitutive_model *m,
                           double *F,
                           double *hFI,
                           const int stepno) const;
  virtual int get_hardening(const Constitutive_model *m,
                            double *var,
                            const int stepno) const;
  virtual int get_plast_strain_var(const Constitutive_model *m,
                                   double *lam_p)
  const { return cm_get_lam_p(m, lam_p);};

  virtual int get_subdiv_param(const Constitutive_model *m,
                               double *var,
                               const double t) const;
  virtual int write_restart(FILE *fp,
                            const Constitutive_model *m) const;
  virtual int read_restart(FILE *fp,
                           Constitutive_model *m) const;
  virtual int destroy_ctx(void **ctx) const;
  virtual int compute_dMdu(const Constitutive_model *m,
                           const void *ctx,
                           const double *Grad_op,
                           const int nne,
                           const int ndofn,
                           double *dM_du) const;
  virtual int read_param(FILE *in) const;
  virtual int set_init_vals(Constitutive_model *m) const;
};

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Constitutive_model object
 * \return non-zero on internal error
 */
//int plasticity_model_initialize(Model_parameters *p);

//int plasticity_model_destory(Model_parameters *p);

/// Construct and initialize the model context for calling functions
/// through the constitutive model interface.
///
/// \param[in,out] ctx - handle to an opaque model context object.
/// \param[in] F The total deformation gradient.
/// \param[in] dt time increment
/// \param[in] alpha mid-point alpha
/// \param[in] eFnpa elastic deformation gradient at t = n + alpha
/// \param[in] hFn thermal part deformation gradient at t = n
/// \param[in] hFnp1 thermal part deformation gradient at t = n + 1
/// \param[in] is_coulpled_with_thermal flag for coupling with thermal
/// \return non-zero on internal error.
int plasticity_model_ctx_build(void **ctx,
                               double *F,
                               const double dt,
                               const double alpha,
                               double *eFnpa,
                               double *hFn,
                               double *hFnp1,
                               const int is_coulpled_with_thermal);

/**
 * compute tangent of the plasticity part of deformation gradient with respect to the deformation
 *
 * \param[in,out] dMdu - computed tangent of plasticity part of deformation gradient
 * \param[in] du - Grad(du)
 * \param[in,out] m - pointer to a Constitutive_model object.
 * \param[in] dt - time step size
 * \return non-zero on internal error.
 */
int plasticity_model_slip_system(double *P);

struct EPS;
struct Element;

int plasticity_model_set_orientations(EPS *eps,
                                      const int ne,
                                      const Element *elem,
                                      const int n_mat,
                                      const HOMMAT *hmat_list,
				      int myrank);
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

#endif
