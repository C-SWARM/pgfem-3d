#ifndef PRESS_THETA_H
#define PRESS_THETA_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef SUPP_H
#include "supp.h"
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

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
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
