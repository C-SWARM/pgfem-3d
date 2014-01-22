#ifndef COMPUTE_MACRO_F_H
#define COMPUTE_MACRO_F_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

  /** Compute volume average deformation gradient.  Contains global
      communication. */
  double* computeMacroF(ELEMENT *elem,
			long ne,
			NODE *node,
			long nn,
			EPS *eps,
			double oVolume,
			MPI_Comm mpi_comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
