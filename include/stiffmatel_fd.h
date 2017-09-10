#ifndef PGFEM3D_STIFFMATEL_FD_H
#define PGFEM3D_STIFFMATEL_FD_H

#include "cohesive_element.h"
#include "crpl.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "sig.h"

/**
 * Compute the element stiffness matrix for the original three-field
 * formulation. To be reimplemented using the Hu_Washizu_element
 * functions.
 */
int stiffmatel_fd (long ii,
                   long ndofn,
                   long nne,
                   long *nod,
                   double *x,
                   double *y,
                   double *z,
                   Element *elem,
                   MATGEOM matgeom,
                   HOMMAT *hommat,
                   Node *node,
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

#endif /* #define PGFEM3D_STIFFMATEL_FD_H */
