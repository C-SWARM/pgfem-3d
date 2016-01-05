/**
 * Function definitions of some common placeholder functions to
 * improve code reuse.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#include "cm_placeholder_functions.h"
#include "data_structure_c.h"

Define_Matrix(double);
static const int dim = 3;
static const int tensor = 9;

int cm_get_var_zero(const Constitutive_model *m,
                    double *var) { *var = 0.0; return 0; }

int cm_get_F_zero(const Constitutive_model *m,
                  Matrix_double *F) { Matrix_init(*F, 0.0); return 0; }

int cm_get_F_eye(const Constitutive_model *m,
                 Matrix_double *F) { Matrix_eye(*F, 3); return 0; }

int cm_get_lam_p(const Constitutive_model *m,
                 double *lam_p)
{
  int err = 0;
  Matrix_double Fp, Cp;
  Matrix_construct_redim(double, Fp, dim, dim);
  Matrix_construct_redim(double, Cp, dim, dim);

  /* get the plastic deformation gradient at (n + 1) */
  m->param->get_pFn(m, &Fp);

  /* compute lam_p = sqrt( tr(Cp) / 3 ) */
  /* tr(Cp) = Fp : Fp */
  *lam_p = 0.0;
  for (int i = 0; i < tensor; i++) {
    *lam_p += Cp.m_pdata[i] * Cp.m_pdata[i];
  }
  *lam_p = sqrt( *lam_p / 3.0 );

  Matrix_cleanup(Fp);
  Matrix_cleanup(Cp);
  return err;
}

int cm_compute_null_dMdu(const Constitutive_model *m,
                         const void *ctx,
                         const double *Grad_op,
                         const int nne,
                         const int ndofn,
                         double *dM_du)
{
  memset(dM_du, 0, tensor * nne *dim * sizeof(*dM_du));
  return 0;
}

