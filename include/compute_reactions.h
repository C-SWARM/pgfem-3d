/* HEADER */
#ifndef COMPUTE_REACTIONS_H
#define COMPUTE_REACTIONS_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#include "PGFEM_mpi.h"

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

#ifndef SUPP_H
#include "supp.h"
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
		      NODE *node,
		      ELEMENT *elem,
		      MATGEOM matgeom,
		      HOMMAT *hommat,
		      SUPP sup,
		      EPS *eps,
		      SIG *sig,
		      double nor_min,
		      CRPL *crpl,
		      double dt,
		      double stab,
		      MPI_Comm mpi_comm,
		      const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef COMPUTE_REACTIONS_H */
