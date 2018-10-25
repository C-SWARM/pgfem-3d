/**
 * This file defines the interface to the BPA plasticity model.
 * REFERENCES:
 *
 * Mary C. BOYCE, David M. PARKS, Ali S. ARGON (1988). LARGE INELASTIC
 * DEFORMATION OF GLASSY POLYMERS.PART I: RATE DEPENDENT CONSTITUTIVE
 * MODEL. Mechanics of Materials, 7:15-33.
 *
 * Holopainen, S. (2013). Modeling of the mechanical behavior of
 * amorphous glassy polymers under variable loadings and comparison
 * with state-of-the-art model predictions. Mechanics of Materials,
 * 66:35–58.
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#pragma once
#ifndef PLASTICITY_MODEL_BPA_H
#define PLASTICITY_MODEL_BPA_H

#include "constitutive_model.h"
#include "cm_placeholder_functions.h"

class BPA_PARAM: public Model_parameters
{
  public:

  virtual int model_dependent_initialization(void);

  virtual int integration_algorithm(Constitutive_model *m,
                                    CM_Ctx &cm_ctx) const;
  virtual int compute_dev_stress(const Constitutive_model *m,
                                 CM_Ctx &cm_ctx,
                                 double *S) const;
  virtual int compute_dudj(const Constitutive_model *m,
                           CM_Ctx &cm_ctx,
                           double *value) const;
  virtual int compute_dev_tangent(const Constitutive_model *m,
                                  CM_Ctx &cm_ctx,
                                  double *L) const;
  virtual int compute_d2udj2(const Constitutive_model *m,
                             CM_Ctx &cm_ctx,
                             double *value) const;
  virtual int update_elasticity(const Constitutive_model *m,
                                CM_Ctx &cm_ctx,
                                double *L,
                                double *S,
                                const bool compute_stiffness) const;
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
                                   double *lam_p)
  const { return cm_get_lam_p(m, lam_p);};
                              
  virtual int write_restart(FILE *fp,
                            const Constitutive_model *m) const;
  virtual int read_restart(FILE *fp,
                           Constitutive_model *m) const;
  virtual int compute_dMdu(const Constitutive_model *m,
                           CM_Ctx &cm_ctx,
                           const double *Grad_op,
                           const int nne,
                           const int ndofn,
                           double *dM_du) const;
  virtual int read_param(FILE *in) const;
  virtual int set_init_vals(Constitutive_model *m) const;
};

#endif
