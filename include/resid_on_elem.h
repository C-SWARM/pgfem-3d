/* HEADER */
#pragma once
#ifndef RESID_ON_ELEM_H
#define RESID_ON_ELEM_H

#include "data_structure.h"
#include "element.h"
#include "node.h"
#include "matgeom.h"
#include "hommat.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Compute the residual on an element.
   *
   * I (Mosby) believe this is a carry over from the serial code and
   * is unused.
 */
  int resid_on_elem (long ii,
		     long ndofn,
		     long nne,
		     long *nod,
		     ELEMENT *elem,
		     NODE *node,
		     MATGEOM matgeom,
		     HOMMAT *hommat,
		     double *x,
		     double *y,
		     double *z,
		     EPS *eps,
		     SIG *sig,
		     double *r_e,
		     long npres,
		     double nor_min,
		     double *fe,
		     CRPL *crpl,
		     double dt,
		     const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef RESID_ON_ELEM_H */
