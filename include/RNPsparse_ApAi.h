#ifndef RNPSPARSE_APAI_H
#define RNPSPARSE_APAI_H

#include "PGFEM_mpi.h"

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef PGFEM_COMM_H
#include "pgfem_comm.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** This function builds the distributed CSR matrix after
      renumbering. This function is almost exactly the same as
      Psparse_ApAi except that it checks to make sure that the LtoG
      mapping exists.  If it does not, then it adds local degrees of
      freedom and forces the existance of teh LtoG mapping.  It is
      very possible that this function will be merged with
      Psparse_ApAi and be used for the original matrix as well as the
      renumbered matrix. */
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
