/**
 * Declare/define functions that integrate with the plasticty model
 * interface (see constitutive_model.h). This file defines the
 * interface for the purely hyperelastic model.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */
#pragma once
#ifndef PLASTICITY_MODEL_NONE_H
#define PLASTICITY_MODEL_NONE_H

#include "constitutive_model.h"

class HE_PARAM: public Model_parameters
{
  public:

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
  virtual int set_F(const Constitutive_model *m,
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
  virtual int write_restart(FILE *fp,
                            const Constitutive_model *m) const;
  virtual int read_restart(FILE *fp,
                           Constitutive_model *m) const;
  virtual int destroy_ctx(void **ctx) const;
  virtual int read_param(FILE *in) const;
};


/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Constitutive_model object
 * \return non-zero on internal error
 */
//int plasticity_model_none_initialize(Model_parameters *p);

/// Construct and initialize the model context for calling functions
/// through the constitutive model interface.
///
/// \param[in,out] ctx - handle to an opaque model context object.
/// \param[in] F The total deformation gradient.
/// \param[in] eFnpa elastic deformation gradient at t = n + alpha
/// \param[in] hFn thermal part deformation gradient at t = n
/// \param[in] hFnp1 thermal part deformation gradient at t = n + 1
/// \param[in] is_coulpled_with_thermal flag for coupling with thermal
/// \return non-zero on internal error.
int plasticity_model_none_ctx_build(void **ctx,
                                    double *F,
                                    double *eFnpa,
                                    double *hFn,
                                    double *hFnp1,
                                    const int is_coulpled_with_thermal);


#endif
