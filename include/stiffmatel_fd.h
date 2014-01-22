#ifndef STIFFMATEL_FD_H
#define STIFFMATEL_FD_H

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

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
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
  int stiffmatel_fd (long ii,
		     long ndofn,
		     long nne,
		     long *nod,
		     double *x,
		     double *y,
		     double *z,
		     ELEMENT *elem,
		     MATGEOM matgeom,
		     HOMMAT *hommat,
		     NODE *node,
		     SIG *sig,
		     EPS *eps,
		     double *r_e,
		     long npres,
		     double nor_min,
		     double *Ks,
		     double dt,
		     CRPL *crpl,
		     long FNR,
		     double lm,
		     double *fe,
		     const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef STIFFMATEL_FD_H */
