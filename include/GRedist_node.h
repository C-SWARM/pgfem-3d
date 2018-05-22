#ifndef PGFEM3D_GREDIST_NODE_H
#define PGFEM3D_GREDIST_NODE_H

#include "pgfem3d/Communication.hpp"
#include "node.h"

/** Redistributes the degrees of freedom on the boundary nodes and
    returns the number of boundary nodes on the domain. */
long GRedist_node (const int nproc,
                   const int myrank,
                   const long nn,
                   const long ndofn,
                   Node *node,
		   const pgfem3d::CommunicationStructure *com,
                   const int mp_id);

#endif /* #define PGFEM3D_GREDIST_NODE_H */
