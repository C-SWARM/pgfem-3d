#ifndef PGFEM3D_COMPUTE_MACRO_F_H
#define PGFEM3D_COMPUTE_MACRO_F_H

#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "node.h"
#include "pgfem3d/Communication.hpp"

/** Compute volume average deformation gradient.  Contains global
    communication. */
double* computeMacroF(Element *elem,
                      long ne,
                      Node *node,
                      long nn,
                      EPS *eps,
                      double oVolume,
                      pgfem3d::CommunicationStructure *com);

#endif /* #define PGFEM3D_COMPUTE_MACRO_F_H */
