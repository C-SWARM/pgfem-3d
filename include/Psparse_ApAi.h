/* HEADER */

/**
 * AUTHORS:
 *    Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef PSPARSE_APAI_H
#define PSPARSE_APAI_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "pgfem_comm.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Create the global sparsity pattern and commincation structure.
   */
  int* Psparse_ApAi (int nproc,
		     int myrank,
		     long ne,
		     long n_be,
		     long nn,
		     long ndofn,
		     long ndofd,
		     ELEMENT *elem,
		     BOUNDING_ELEMENT *b_elems,
		     NODE *node,
		     int *Ap,
		     long nce,
		     COEL *coel,
		     long *DomDof,
		     int *GDof,
		     COMMUN comm,
		     MPI_Comm Comm_Orig,
		     const int cohesive,
		     const int mp_id);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PSPARSE_APAI_H */
