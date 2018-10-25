/// This file defines the interface for the splited isotropic viscous
/// damage model.
/// 
/// REFERENCES:
/// 
/// Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
/// finite strains." Computational Mechanics 5.5 (1989): 375-400.
/// 
/// Authors:
///  Sangmin Lee, University of Notre Dame, <slee43@nd.edu>

#pragma once
#ifndef CM_ISO_VISCOUS_DAMAGE_H
#define CM_ISO_VISCOUS_DAMAGE_H

#include "constitutive_model.h"
#include "cm_placeholder_functions.h"

class CM_IVD_PARAM: public Model_parameters
{
  public:

  virtual int model_dependent_initialization(void);
  virtual int model_dependent_finalization(void);

  virtual int integration_algorithm(Constitutive_model *m,
                                    CM_Ctx &cm_ctx) const;
  virtual int compute_dev_stress(const Constitutive_model *m,
                                 CM_Ctx &cm_ctx,
                                 double *S) const;
  virtual int compute_dudj(const Constitutive_model *m,
                           CM_Ctx &cm_ctx,
                           double *value) const;
  virtual double compute_dudj(const Constitutive_model *m,
                              double theta_e,
                              const int npa,
                              const double alpha) const;                           
  virtual int compute_dev_tangent(const Constitutive_model *m,
                                  CM_Ctx &cm_ctx,
                                  double *L) const;
  virtual int compute_d2udj2(const Constitutive_model *m,
                             CM_Ctx &cm_ctx,
                             double *value) const;
  virtual double compute_d2udj2(const Constitutive_model *m,
                                double theta_e,
                                const int npa,
                                const double alpha) const;                             
  virtual int update_elasticity(const Constitutive_model *m,
                                CM_Ctx &cm_ctx,
                                double *L,
                                double *S,
                                const bool compute_stiffness) const;
  virtual int update_elasticity_dev(const Constitutive_model *m,
                                    double *eFnpa,
                                    double *L,
                                    double *S,
                                    const int npa,
                                    const double alpha,
                                    const double dt,
                                    const bool compute_stiffness = 0) const;
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
  virtual int read_param(FILE *in) const;
  virtual int set_init_vals(Constitutive_model *m) const;
};

#endif
