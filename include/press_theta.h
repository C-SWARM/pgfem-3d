/* HEADER */
#pragma once
#ifndef PGFEM3D_PRESS_THETA_H
#define PGFEM3D_PRESS_THETA_H

#include "data_structure.h"
#include "element.h"
#include "node.h"
#include "supp.h"
#include "matgeom.h"
#include "hommat.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "PGFem3D_options.h"

/**
 * Update condende variables in the original three field
 * formulation. This is used for qudratic tetrahedra and trininear
 * hexahedron elements. This function will be reimplemented using
 * the Hu_Washizu_element functions in the future.
 */
void press_theta (long ne,
                  long ndofn,
                  long npres,
                  Element *elem,
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
                  const PGFem3D_opt *opts,
                  const int mp_id);

#endif /* #ifndef PGFEM3D_PRESS_THETA_H */
