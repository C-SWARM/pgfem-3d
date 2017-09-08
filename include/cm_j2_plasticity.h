/* HEADER */
/**
 * This file defines the interface for the J2 plasticity + damage
 * model.
 *
 * REFERENCES:
 *
 * Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
 * finite strains." Computational Mechanics 5.5 (1989): 375-400.
 *
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#pragma once
#ifndef CM_J2_PLASTICITY_H
#define CM_J2_PLASTICITY_H


#include "constitutive_model.h"
#include "cm_placeholder_functions.h"

class CM_J2P_PARAM: public Model_parameters
{
  public:

  virtual int model_dependent_initialization(void);

  virtual int integration_algorithm(Constitutive_model *m,
                                    const void *usr_ctx) const;
  virtual int compute_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 double *S) const;
  virtual int compute_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *value) const;
  virtual int update_elasticity(const Constitutive_model *m,
                                const void *ctx,
                                double *L,
                                double *S,
                                const int compute_stiffness) const;
  virtual int update_state_vars(Constitutive_model *m) const;  
  virtual int reset_state_vars(Constitutive_model *m) const;
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
  virtual int get_hardening(const Constitutive_model *m,
                            double *var,
                            const int stepno) const;
  virtual int get_plast_strain_var(const Constitutive_model *m,
                                   double *lam_p) const;                              
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

/**
 * Construct and initialize the model context for calling functions
 * through the constitutive modeling interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F, _total_ deformation gradient
 * \param[in] dt, time increment
 * \return non-zero on internal error.
 */
int j2d_plasticity_model_ctx_build(void **ctx,
                                   const double *F,
                                   const double dt);

#endif
