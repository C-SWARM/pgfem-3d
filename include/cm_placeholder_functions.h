/**
 * Function declarations of some common placeholder functions to
 * improve code reuse.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#pragma once
#ifndef CM_PLACEHOLDER_FUNCTIONS_H
#define CM_PLACEHOLDER_FUNCTIONS_H

#include "constitutive_model.h"

/**
 * Function to return 0 for a particular variable.
 */
int cm_get_var_zero(const Constitutive_model *m,
                    double *var);

/**
 * Function to return 2nd order tensor filled with zeros.
 */
int cm_get_F_zero(const Constitutive_model *m,
                  double *F);

/**
 * Funtion to return 2nd order identity tensor.
 */
int cm_get_F_eye(const Constitutive_model *m,
                 double *F);

/**
 * Function to compute the plastic stretch \lambda^p used by several
 * models for the plastic strain measure.
 *
 * lam_p = sqrt( tr(Fp Fp') / 3 )
 */
int cm_get_lam_p(const Constitutive_model *m,
                 double *lam_p);

/**
 * Do not comupte dM_du and fill with zeros.
 */
int cm_compute_null_dMdu(const Constitutive_model *m,
                         CM_Ctx &cm_ctx,
                         const double *Grad_op,
                         const int nne,
                         const int ndofn,
                         double *dM_du);

/**
 * Do not cause subdivision (subdiv_param = 0)
 */
int cm_no_subdiv(const Constitutive_model *m,
                 double *subdiv_param,
                 const double dt);

#endif
