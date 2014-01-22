#ifndef TA_GA_H
#define TA_GA_H

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  void TA_GA (long nss,
	      long mat,
	      CRPL *crpl,
	      double Har,
	      double **eFn,
	      double **S,
	      double *TA,
	      double *GA);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef TA_GA_H */
