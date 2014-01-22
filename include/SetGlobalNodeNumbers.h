#ifndef SETGLOBALNODENUMBERS_H
#define SETGLOBALNODENUMBERS_H

#include "PGFEM_mpi.h"

#ifndef NODE_H
#include "node.h"
#endif

/** Assigns each node a global number. Nombering is done contiguously
    from process 0 keeping track of the boundary nodes.*/
void SetGlobalNodeNumbers(int nNodesDom,NODE *node,MPI_Comm comm);
#endif
