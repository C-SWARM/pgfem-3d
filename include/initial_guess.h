#ifndef INITIAL_GUESS_H
#define INITIAL_GUESS_H

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif


#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
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

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef INITIAL_GUESS_H */
