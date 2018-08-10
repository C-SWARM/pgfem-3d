/// This file defines the implementation for finite strain J2
/// plasticity with damage employing the stress split. This
/// formulations was chosen as it allows for simple computation of the
/// exact algorithmic tangent for the coupled plasticity and
/// damage. This is because the evaluation of the damage parameter is
/// not a function of any of the plastic variables.
/// 
/// REFERENCES:
/// 
/// Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
/// finite strains." Computational Mechanics 5.5 (1989): 375-400.
/// 
/// Mosby, Matthew and K. Matous. "On mechanics and material length scales of failure
/// in heterogeneous interfaces using a finite strain high performance
/// solver." Modelling and Simulation in Materials Science and
/// Engineering 23.8 (2015): 085014.
/// 
/// Authors:
///  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
///  Sangmin Lee, University of Notre Dame, <slee43@nd.edu>

#pragma once
#ifndef CM_J2_PLASTICITY_H
#define CM_J2_PLASTICITY_H


#include "constitutive_model.h"
#include "cm_placeholder_functions.h"

class CM_J2P_PARAM: public Model_parameters
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
                                    double *eFnpa,
                                    double *L,
                                    double *S,
                                    const int npa,
                                    const double alpha,
                                    const double dt,
                                    const int compute_stiffness = 0) const;
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
                                   double *chi) const;

  virtual int get_subdiv_param(const Constitutive_model *m,
                               double *var,
                               const double t) const;
  virtual int write_restart(FILE *fp,
                            const Constitutive_model *m) const;
  virtual int read_restart(FILE *fp,
                           Constitutive_model *m) const;
  virtual int destroy_ctx(void **ctx) const;
  virtual int read_param(FILE *in) const;
  virtual int set_init_vals(Constitutive_model *m) const;
};

// Construct and initialize the model context for calling functions
// through the constitutive modeling interface.
int j2d_plasticity_model_ctx_build(void **ctx,
                                   double *F,
                                   const double dt,
                                   const double alpha,
                                   double *eFnpa,
                                   double *hFn,
                                   double *hFnp1,
                                   const int is_coulpled_with_thermal,
                                   const int npa);

#endif
