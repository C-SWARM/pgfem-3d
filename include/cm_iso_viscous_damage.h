/* HEADER */
/**
 * This file defines the interface for the isotropic viscous
 * damage model.
 *
 * REFERENCES:
 *
 * Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
 * finite strains." Computational Mechanics 5.5 (1989): 375-400.
 *
 * Mosby, Matthew and K. Matous. "On mechanics and material length scales of failure
 * in heterogeneous interfaces using a finite strain high performance
 * solver." Modelling and Simulation in Materials Science and
 * Engineering 23.8 (2015): 085014.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#pragma once
#ifndef CM_ISO_VISCOUS_DAMAGE_H
#define CM_ISO_VISCOUS_DAMAGE_H

#include "constitutive_model.h"
#include "cm_placeholder_functions.h"

class CM_IVD_PARAM: public Model_parameters
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

/**
 * Initialize the Model_parameters object for this particular model.
 *
 * \param[in,out] p - pointer to a Model_parameters object
 * \return non-zero on internal error
 */
int iso_viscous_damage_model_initialize(Model_parameters *p) {return 0;};

/**
 * Construct and initialize the model context for calling functions
 * through the constitutive modeling interface.
 *
 * \param[in,out] ctx - handle to an opaque model context object.
 * \param[in] F, _total_ deformation gradient
 * \param[in] dt, time increment
 * \return non-zero on internal error.
 */
int iso_viscous_damage_model_ctx_build(void **ctx,
                                       const double *F,
                                       const double dt);

/**
 * Provide a public interface to the integration algorithm to allow
 * other constitutive models to include damage as well.
 */
int ivd_public_int_alg(double *var_w,
                       double *var_X,
                       double *var_H,
                       int *flag_damaged,
                       const double var_wn,
                       const double var_Xn,
                       const double dt,
                       const double Ybar,
                       const double param_mu,
                       const double param_ome_max,
                       const double param_p1,
                       const double param_p2,
                       const double param_Yin);

/**
 * Provide a public interface to the subdivision parameter to allow
 * other constitutive models to subdivide in a consistent fashion when
 * using the damage integration algorithm.
 */
int ivd_public_subdiv_param(const double var_wn,
                            const double var_w,
                            double *subdiv_param);

#endif
