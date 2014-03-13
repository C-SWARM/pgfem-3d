/* HEADER */
#pragma once
#ifndef PRESS_THETA_H
#define PRESS_THETA_H

#include "element.h"
#include "node.h"
#include "supp.h"
#include "matgeom.h"
#include "hommat.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Update condende variables in the original three field
   * formulation. This is used for qudratic tetrahedra and trininear
   * hexahedron elements. This function will be reimplemented using
   * the Hu_Washizu_element functions in the future.
   */
  void press_theta (long ne,
		    long ndofn,
		    long npres,
		    ELEMENT *elem,
		    NODE *node,
		    double *d_r,
		    double *rr,
		    SUPP sup,
		    MATGEOM matgeom,
		    HOMMAT *hommat,
		    EPS *eps,
		    SIG *sig,
		    long iter,
		    double nor_min,
		    double dt,
		    CRPL *crpl,
		    const PGFem3D_opt *opts);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PRESS_THETA_H */
