/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */

#include "plasticity_model_none.h"
#include "constitutive_model.h"
#include "new_potentials.h"
#include "data_structure_c.h"

Define_Matrix(double);

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

static int plasticity_none_int_alg(Constitutive_model *m,
                                   const void *ctx)
{
  int err = 0;
  /* hyperelastic, no integration algorithm */
  return err;
}

static int plasticity_none_dev_stress(const Constitutive_model *m,
                                      const void *ctx,
                                      Matrix_double *stress)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  devStressFuncPtr Stress = getDevStressFunc(-1,m->param->p_hmat);
  Stress(CTX->C,m->param->p_hmat,stress->m_pdata);
  return err;
}

static int plasticity_none_dudj(const Constitutive_model *m,
                                const void *ctx,
                                double *dudj)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
  Pressure(*(CTX->J),m->param->p_hmat,dudj);
  return err;
}

static int plasticity_none_dev_tangent(const Constitutive_model *m,
                                       const void *ctx,
                                       Matrix_double *tangent)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,m->param->p_hmat);
  Tangent(CTX->C,m->param->p_hmat,tangent->m_pdata);
  return err;
}

static int plasticity_none_d2udj2(const Constitutive_model *m,
                                  const void *ctx,
                                  double *d2udj2)
{
  int err = 0;
  const none_ctx *CTX = ctx;
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
  D_Pressure(*(CTX->J),m->param->p_hmat,d2udj2);
  return err;
}

static int plasticity_none_update(Constitutive_model *m)
{
  int err = 0;
  /* there are no state variables to update */
  return err;
}

static int plasticity_none_reset(Constitutive_model *m)
{
  int err = 0;
  /* there are no state varables to update */
  return err;
}

static int plasticity_none_info(Model_var_info **info)
{
  int err = 0;

  /* make sure I don't leak memory */
  if( *info != NULL) err += model_var_info_destroy(info);

  /* construct valid info object */
  *info = malloc(sizeof(**info));
  (*info)->F_names = malloc(1*sizeof(char *));
  (*info)->F_names[0] = (char *) malloc(1024*sizeof(char));
  
  sprintf((*info)->F_names[0], "F");
  (*info)->var_names = NULL;
  (*info)->n_Fs = 1;
  (*info)->n_vars = 0;
  return err;
}

static int he_get_Fn(const Constitutive_model *m,
                     Matrix_double *F)
{
  assert(0 && "this function is not implemented");
  return 1;
}

static int he_get_pF(const Constitutive_model *m,
                     Matrix_double *F)
{
  int err = 0;
  Matrix_eye(*F,3);
  return err;
}

static int he_get_pFn(const Constitutive_model *m,
                      Matrix_double *F)
{
  return he_get_pF(m,F);
}

static int he_get_eF(const Constitutive_model *m,
                     Matrix_double *F)
{
  assert(0 && "this function is not implemented");
  return 1;
}

static int he_get_eFn(const Constitutive_model *m,
                      Matrix_double *F)
{
  assert(0 && "this function is not implemented");
  return 1;
}

static int he_get_hardening(const Constitutive_model *m,
                            double *var)
{
  *var = 0.0;
  return 0;
}

static int he_compute_dMdu(const Constitutive_model *m,
                           const void *ctx,
                           const double *Grad_op,
                           const int nne,
                           const int ndofn,
                           double *dM_du)
{
  int err = 0;
  /* there is no plastic deformation in this formulation, return zeros
     in dM_du */
  memset(dM_du, 0, nne * ndofn * 9 * nne * ndofn * 9 * sizeof(*dM_du));
  return err;
}

int plasticity_model_none_initialize(Model_parameters *p)
{
  int err = 0;

  /* set functions */
  p->integration_algorithm = plasticity_none_int_alg;
  p->compute_dev_stress = plasticity_none_dev_stress;
  p->compute_dudj = plasticity_none_dudj;
  p->compute_dev_tangent = plasticity_none_dev_tangent;
  p->compute_d2udj2 = plasticity_none_d2udj2;
  p->update_state_vars = plasticity_none_update;
  p->reset_state_vars = plasticity_none_reset;
  p->get_var_info = plasticity_none_info;
  p->get_Fn = he_get_Fn;
  p->get_pF = he_get_pF;
  p->get_pFn = he_get_pFn;
  p->get_eF = he_get_eF;
  p->get_eFn = he_get_eFn;
  p->get_hardening = he_get_hardening;
  p->destroy_ctx = plasticity_model_none_ctx_destroy;
  p->compute_dMdu = he_compute_dMdu;

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
