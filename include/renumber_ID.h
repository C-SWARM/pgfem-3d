#ifndef RENUMBER_ID_H
#define RENUMBER_ID_H

#ifndef NODE_H
#include "node.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

  /** This function is used after creating the reorder mapping to
      re-assign the global ID numbers. The variable ndofd is the
      number of dofs on each node, nn is the number of nodes on the
      domain, and g_order is the global renumber mapping. Note that
      the global DOF numbering starts at 1, not zero.  Thus, for the
      mapping (which starts at 0) to work, we must subtract 1 from the
      previous value, then add it back after the new value is
      determined. */
  void renumber_ID (int ndofn,int nn,NODE *node,int *g_order,MPI_Comm comm, const int mp_id);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef RENUMBER_ID_H */
