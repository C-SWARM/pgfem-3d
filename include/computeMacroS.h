#ifndef PGFEM3D_COMPUTE_MACROS_H
#define PGFEM3D_COMPUTE_MACROS_H

#include "PGFEM_mpi.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "node.h"
#include "sig.h"

/** Compute volume average 2PK stress.  Contains global
    communication. */
double* computeMacroS(Element *elem,
                      long ne,
                      NODE *node,
                      long nn,
                      SIG *sig,
                      double oVolume,
                      MPI_Comm mpi_comm);

/** Compute volume average 1PK stress.  Contains global
    communication. */
double* computeMacroP(Element *elem,
                      long ne,
                      NODE *node,
                      long nn,
                      SIG *sig,
                      EPS *eps,
                      double oVolume,
                      MPI_Comm mpi_comm);

#endif /* #define PGFEM3D_COMPUTE_MACRO_S_H */
