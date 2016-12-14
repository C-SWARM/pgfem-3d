#ifndef SETGLOBALNODENUMBERS_H
#define SETGLOBALNODENUMBERS_H

#include "PGFEM_mpi.h"

#ifndef NODE_H
#include "node.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

/** Assigns each node a global number. Nombering is done contiguously
    from process 0 keeping track of the boundary nodes.*/
void SetGlobalNodeNumbers(int nNodesDom,NODE *node,MPI_Comm comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
