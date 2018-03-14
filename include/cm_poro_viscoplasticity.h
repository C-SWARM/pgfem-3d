/// Declare/define functions for constitutive model interface performing
/// poro-viscoplasticity integration algorithm
/// 
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// Alberto Salvadori, [1], <asalvad2@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#pragma once
#ifndef CM_PORO_VISCO_PLASTICITY_H
#define CM_PORO_VISCO_PLASTICITY_H

#include "constitutive_model.h"
#include "cm_placeholder_functions.h"

class CM_PVP_PARAM: public Model_parameters
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
  virtual double compute_dudj(const Constitutive_model *m,
                              double theta_e,
                              const int npa,
                              const double alpha) const;                           
  virtual int compute_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  double *L) const;
  virtual int compute_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *value) const;
  virtual double compute_d2udj2(const Constitutive_model *m,
                                double theta_e,
                                const int npa,
                                const double alpha) const;                             
  virtual int update_elasticity(const Constitutive_model *m,
                                const void *ctx,
                                double *L,
                                double *S,
                                const int compute_stiffness) const;
  virtual int update_elasticity_dev(const Constitutive_model *m,
                                    double *eF,
                                    double *L,
                                    double *S,
                                    const int npa,
                                    const double alpha,
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

 
/// Construct and initialize the poro-viscoplasticity model context 
/// for calling functions through the constitutive modeling interface
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
int poro_viscoplasticity_model_ctx_build(void **ctx,
                                         double *F,
                                         const double dt,
                                         const double alpha,
                                         double *eFnpa,
                                         double *hFn,
                                         double *hFnp1,
                                         const int is_coulpled_with_thermal,
                                         const int npa);

#endif