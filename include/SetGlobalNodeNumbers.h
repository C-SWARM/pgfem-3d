#ifndef PGFEM3D_SET_GLOBAL_NODE_NUMBERS_H
#define PGFEM3D_SET_GLOBAL_NODE_NUMBERS_H

#include "pgfem3d/Communication.hpp"
#include "node.h"

/** Assigns each node a global number. Nombering is done contiguously
    from process 0 keeping track of the boundary nodes.*/
void SetGlobalNodeNumbers(int nNodesDom,
                          Node *node,
                          pgfem3d::CommunicationStructure *com);

#endif // #define PGFEM3D_SET_GLOBAL_NODE_NUMBERS_H
