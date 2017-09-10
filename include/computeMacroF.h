#ifndef COMPUTE_MACRO_F_H
#define COMPUTE_MACRO_F_H

#include "PGFEM_mpi.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "node.h"

/** Compute volume average deformation gradient.  Contains global
    communication. */
double* computeMacroF(Element *elem,
                      long ne,
                      NODE *node,
                      long nn,
                      EPS *eps,
                      double oVolume,
                      MPI_Comm mpi_comm);

#endif /* #define PGFEM3D_COMPUTE_MACRO_F_H */
