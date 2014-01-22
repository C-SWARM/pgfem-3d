#ifndef DIFFERENTIATION_OF_UU_H
#define DIFFERENTIATION_OF_UU_H

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  void differentiation_of_UU (long ii,
			      long ip,
			      long mat,
			      CRPL *crpl,
			      double Tr,
			      double Jr,
			      double **Fr,
			      double **FnB,
			      double **eFnB,
			      double **FstB,
			      double **S,
			      double L[3][3][3][3],
			      double **UU,
			      SIG *sig,
			      EPS *eps,
			      double dt,
			      const PGFem3D_opt *opts);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef DIFFERENTIATION_OF_UU_H */
