#ifndef PGFEM3D_COMPUTE_MACROS_H
#define PGFEM3D_COMPUTE_MACROS_H

#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "node.h"
#include "sig.h"
#include "pgfem3d/Communication.hpp"

/** Compute volume average 2PK stress.  Contains global
    communication. */
double* computeMacroS(Element *elem,
                      long ne,
                      Node *node,
                      long nn,
                      SIG *sig,
                      double oVolume,
                      pgfem3d::CommunicationStructure *com);

/** Compute volume average 1PK stress.  Contains global
    communication. */
double* computeMacroP(Element *elem,
                      long ne,
                      Node *node,
                      long nn,
                      SIG *sig,
                      EPS *eps,
                      double oVolume,
                      pgfem3d::CommunicationStructure *com);

#endif /* #define PGFEM3D_COMPUTE_MACRO_S_H */
