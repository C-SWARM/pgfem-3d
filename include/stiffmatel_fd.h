/* HEADER */
#pragma once
#ifndef STIFFMATEL_FD_H
#define STIFFMATEL_FD_H

#include "element.h"
#include "node.h"
#include "matgeom.h"
#include "hommat.h"
#include "cohesive_element.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"

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

#endif /* #ifndef STIFFMATEL_FD_H */
