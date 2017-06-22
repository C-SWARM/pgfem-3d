/* HEADER */
#pragma once
#ifndef DIFFERENTIATION_OF_UU_H
#define DIFFERENTIATION_OF_UU_H

#include "data_structure.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * @todo describe this function.  It has something to do with
   * crystal plasticity.
   */
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
