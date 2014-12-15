/* HEADER */
#pragma once
#ifndef RNPSPARSE_APAI_H
#define RNPSPARSE_APAI_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "node.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "pgfem_comm.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Re-build the global distributed global connectivity and
   * communication structure.
   *
   * This function is almost exactly the same as Psparse_ApAi except
   * that it checks to make sure that the LtoG mapping exists. If it
   * does not, then it adds local degrees of freedom and forces the
   * existance of the LtoG mapping. It is very possible that this
   * function will be merged with Psparse_ApAi and be used for the
   * original matrix as well as the renumbered matrix.
   *
   * @todo Remove this function and reimplement @see Psparse_ApAi
   * properly.
   *
  */
  int* RNPsparse_ApAi (int nproc,
		       int myrank,
		       long ne,
		       long n_be,
		       long nn,
		       long ndofn,
		       long *ndofd,
		       ELEMENT *elem,
		       BOUNDING_ELEMENT *b_elems,
		       NODE *node,
		       int *Ap,
		       long nce,
		       COEL *coel,
		       long *DomDof,
		       int *GDof,
		       COMMUN comm,
		       MPI_Comm Comm_RN,
		       const int cohesive);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef RNPSPARSE_APAI_H */
