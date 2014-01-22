#ifndef RESID_ON_ELEM_H
#define RESID_ON_ELEM_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

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
