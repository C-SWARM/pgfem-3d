#ifndef PSPARSE_APAI_H
#define PSPARSE_APAI_H

#include "PGFEM_mpi.h"

#ifndef ELEMENT_H
#include "element.h"
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

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Creates a parallel matrix that is used to solve the system. $name
      populates Ap and returns Ai already allocated and populated. The
      number of DOFs on each domain is not changed.*/
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
		     const int cohesive);

#ifdef __cplusplus
{
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PSPARSE_APAI_H */
