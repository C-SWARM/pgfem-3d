/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 */
#ifndef PGFEM3D_COMPUTE_REACTIONS_H
#define PGFEM3D_COMPUTE_REACTIONS_H

#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "sig.h"
#include "supp.h"

/** After convergence, compute the reaction forces at nodes with
    prescribed deflections. This is achieved by computing the
    residuals on elements which contain node with prescribed
    deflections, filtering the result by deflection ID and summing
    together. NOTE: Only computes the reaction in the direction of
    the prescribed deflection!!! Since dV_u U dV_t = 0, there is no
    need to subtract external forces */

int compute_reactions(long ne,
                      long ndofn,
                      long npres,
                      double *r,
                      Node *node,
                      Element *elem,
                      MATGEOM matgeom,
                      HOMMAT *hommat,
                      SUPP sup,
                      EPS *eps,
                      SIG *sig,
                      double nor_min,
                      CRPL *crpl,
                      double dt,
                      double stab,
		      pgfem3d::CommunicationStructure *com,
                      const int analysis,
                      const int mp_id);

#endif /* #define PGFEM3D_COMPUTE_REACTIONS_H */
