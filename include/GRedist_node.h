#pragma once
#ifndef GREDIST_NODE_H
#define GREDIST_NODE_H

#include "PGFEM_mpi.h"
#include "node.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Redistributes the degrees of freedom on the boundary nodes and
      returns the number of boundary nodes on the domain. */
  long GRedist_node (const int nproc,
		     const int myrank,
		     const long nn,
		     const long ndofn,
		     NODE *node,
		     const MPI_Comm Comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef GREDIST_NODE_H */
