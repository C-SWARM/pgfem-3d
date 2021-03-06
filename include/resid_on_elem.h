#ifndef PGFEM3D_RESID_ON_ELEM_H
#define PGFEM3D_RESID_ON_ELEM_H

#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "sig.h"

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
                   Element *elem,
                   Node *node,
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

#endif /* #define PGFEM3D_RESID_ON_ELEM_H */
