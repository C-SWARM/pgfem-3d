/* HEADER */

#pragma once
#ifndef INITIAL_GUESS_H
#define INITIAL_GUESS_H

#include "sig.h"
#include "eps.h"
#include "crpl.h"

  /**
   * Compute an initial guess for the Newton-type methods.
   *
   * Not sure exacly what this function does or how it is used. Seems
   * associated with crystal plasticity.
   */
  long initial_guess (long ii,
              long ip,
              long tim,
              long iter,
              long nne,
              long ndofn,
              long mat,
              long nss,
              CRPL *crpl,
              long STEP,
              long GAMA,
              double *r_r,
              double ****ST,
              double **Fr,
              double **Fn,
              double **FnB,
              double Jr,
              double Tr,
              double L[3][3][3][3],
              SIG *sig,
              EPS *eps,
              double **UU,
              double *HAR,
              double *LAM,
              double dt);

#endif /* #ifndef INITIAL_GUESS_H */
