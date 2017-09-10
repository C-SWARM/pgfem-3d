#ifndef PGFEM3D_GREDIST_NODE_H
#define PGFEM3D_GREDIST_NODE_H

#include "PGFEM_mpi.h"
#include "comm_hints.h"
#include "node.h"

/** Redistributes the degrees of freedom on the boundary nodes and
    returns the number of boundary nodes on the domain. */
long GRedist_node (const int nproc,
                   const int myrank,
                   const long nn,
                   const long ndofn,
                   Node *node,
                   const Comm_hints *hints,
                   const MPI_Comm Comm,
                   const int mp_id);

#endif /* #define PGFEM3D_GREDIST_NODE_H */
