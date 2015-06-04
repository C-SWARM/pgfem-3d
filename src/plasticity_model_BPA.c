/**
 * This file defines the implementation for the BPA plasticity model.
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

#include "plasticity_model_BPA.h"
#include "constitutive_model.h"
#include "new_potentials.h"
#include "data_structure_c.h"

Define_Matrix(double);

/**
 * Private structure used exclsively with this model and associated
 * functions.
 */

typedef struct {
  int dummy;
} BPA_ctx;

static int BPA_int_alg(Constitutive_model *m,
                       const void *ctx)
{
  int err = 0;

  return err;
}

static int BPA_dev_stress(const Constitutive_model *m,
                          const void *ctx,
                          Matrix_double *dev_stress)
{
  int err = 0;
  return err;
}

static int BPA_dudj(const Constitutive_model *m,
                    const void *ctx,
                    double *dudj)
{
  int err = 0;
  return err;
}

static int BPA_dev_tangent(const Constitutive_model *m,
                           const void *ctx,
                           Matrix_double *dev_tangent)
{
  int err = 0;

  return err;
}

static int BPA_d2udj2(const Constitutive_model *m,
                      const void *ctx,
                      double *d2udj2)
{
  int err = 0;
  return err;
}

static int BPA_update_vars(Constitutive_model *m)
{
  int err = 0;
  return err;
}

static int BPA_reset_vars(Constitutive_model *m)
{
  int err = 0;
  return err;
}

static int BPA_model_info(Model_var_info **info)
{
  int err = 0;
  if (*info != NULL) err += model_var_info_destroy(info);
  (*info) = malloc(sizeof(**info));
  (*info)->F_names = NULL;
  (*info)->var_names = NULL;
  (*info)->n_Fs = 0;
  (*info)->n_vars = 0;
  return err;
}

int plasticity_model_BPA_initialize(Model_parameters *p)
{
  int err = 0;
  return err;
}

int plasticity_model_BPA_ctx_build(void **ctx
                                   //...
                                   )
{
  int err = 0;
  BPA_ctx *t_ctx = malloc(sizeof(*t_ctx));
  t_ctx->dummy = -1;
  *ctx = t_ctx;
  return err;
}

int plasticity_model_BPA_ctx_destroy(void **ctx)
{
  int err = 0;
  free(*ctx);
  *ctx = NULL;
  return err;
}
