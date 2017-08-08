#ifndef GCODE_NUM_NODE_H
#define GCODE_NUM_NODE_H

#ifndef NODE_H
#include "node.h"
#endif

  /** Assigns the global DOF numbers on each node and returns the global
      number of DOFs on the processor domain. Nodes assigned to another
      processor (boundary nodes) and nodes with prescribed
      deflections/boundary conditions are skipped. */
  long Gcode_num_node (long myrank, long ndofn, long nn, NODE *node);

#endif /* #ifndef GCODE_NUM_NODE_H */
