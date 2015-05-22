/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */

#include "plasticity_model_none.h"
#include "constitutive_model.h"
#include "new_potentials.h"

/**
 * Private structure for use exclusively with this model and
 * associated functions.
 */
typedef struct none_ctx {
  const double *C; /*< pointer to the left Cauchy-Green deformation tensor */
  const double *J; /*< pointer to the Jacobian -OR- Volume term */
} none_ctx;

static double compute_bulk_mod(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static int plasticity_none_int_alg(Constitutive_model *p,
                                   const void *ctx)
{
  int err = 0;
  /* hyperelastic, no integration algorithm */
  return err;
}

static int plasticity_none_dev_stress(const Constitutive_model *p,
                                      const void *ctx,
                                      Matrix_handle *stress)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  devStressFuncPtr Stress = getDevStressFunc(-1,p->p_hmat);
  /* Stress(CTX->C,p->p_hmat,stress->data); */
  return err;
}

static int plasticity_none_pressure(const Constitutive_model *p,
                                    const void *ctx,
                                    double *pressure)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  const double kappa = compute_bulk_mod(p->p_hmat);
  dUdJFuncPtr Pressure = getDUdJFunc(-1,p->p_hmat);
  Pressure(*(CTX->J),p->p_hmat,pressure);
  (*pressure) *= kappa / 2.0;
  return err;
}

static int plasticity_none_dev_tangent(const Constitutive_model *p,
                                       const void *ctx,
                                       Matrix_handle *tangent)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,p->p_hmat);
  /* Tangent(CTX->C,p->p_hmat,tangent->data); */
  return err;
}

static int plasticity_none_pressure_tangent(const Constitutive_model *p,
                                            const void *ctx,
                                            double *pres_tan)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  const double kappa = compute_bulk_mod(p->p_hmat);
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,p->p_hmat);
  D_Pressure(*(CTX->J),p->p_hmat,pres_tan);
  (*pres_tan) *= kappa / 2.0;
  return err;
}

static int plasticity_none_update(Constitutive_model *p)
{
  int err = 0;
  return err;
}

static int plasticity_none_reset(Constitutive_model *p)
{
  int err = 0;
  return err;
}

static int plasticity_none_info(Model_var_info **info)
{
  *info = malloc(sizeof(**info));
  (*info)->F_names = NULL;
  (*info)->var_names = NULL;
  (*info)->n_Fs = 0;
  (*info)->n_vars = 0;
  return 0;
}

int plasticity_model_none_initialize(Model_parameters *p)
{
  int err = 0;

  /* set functions */
  p->integration_algorithm = plasticity_none_int_alg;
  p->compute_dev_stress = plasticity_none_dev_stress;
  p->compute_pressure = plasticity_none_pressure;
  p->compute_dev_tangent = plasticity_none_dev_tangent;
  p->compute_pressure_tangent = plasticity_none_pressure_tangent;
  p->update_state_vars = plasticity_none_update;
  p->reset_state_vars = plasticity_none_reset;
  p->get_var_info = plasticity_none_info;

  return err;
}

int plasticity_model_none_ctx_build(void **ctx,
                                    const double *C,
                                    const double *J_or_Theta)
{
  int err = 0;
  none_ctx *t_ctx = malloc(sizeof(none_ctx));

  /* assign internal pointers. NOTE: We are copying the pointer NOT
     the value. No additional memory is allocated. */
  t_ctx->C = C;
  t_ctx->J = J_or_Theta;

  /* assign handle */
  *ctx = t_ctx;
  return err;
}

int plasticity_model_none_ctx_destroy(void **ctx)
{
  int err = 0;
  none_ctx *t_ctx = *ctx;
  /* invalidate handle */
  *ctx = NULL;

  /* we do not control memory for internal pointers */

  /* free object memory */
  free(t_ctx);
  return err;
}
