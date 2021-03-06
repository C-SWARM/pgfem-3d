#ifndef PGFEM3D_RENUMBER_ID_H
#define PGFEM3D_RENUMBER_ID_H

#include "node.h"

/** This function is used after creating the reorder mapping to
    re-assign the global ID numbers. The variable ndofd is the
    number of dofs on each node, nn is the number of nodes on the
    domain, and g_order is the global renumber mapping. Note that
    the global DOF numbering starts at 1, not zero.  Thus, for the
    mapping (which starts at 0) to work, we must subtract 1 from the
    previous value, then add it back after the new value is
    determined. */
void renumber_ID(int ndofn,
                 int nn,
                 Node *node,
                 int *g_order,
		 int myrank,
                 const int mp_id);

#endif /* #define PGFEM3D_RENUMBER_ID_H */
