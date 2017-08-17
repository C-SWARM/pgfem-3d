#ifndef COMPUTE_MACROS_H
#define COMPUTE_MACROS_H

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

  /** Compute volume average 2PK stress.  Contains global
      communication. */
  double* computeMacroS(ELEMENT *elem,
            long ne,
            NODE *node,
            long nn,
            SIG *sig,
            double oVolume,
            MPI_Comm mpi_comm);

  /** Compute volume average 1PK stress.  Contains global
      communication. */
  double* computeMacroP(ELEMENT *elem,
            long ne,
            NODE *node,
            long nn,
            SIG *sig,
            EPS *eps,
            double oVolume,
            MPI_Comm mpi_comm);

#endif /* #ifndef  */
