/**
 * AUTHORS:
 *    Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#ifndef PGFEM3D_PSPARSE_APAI_H
#define PGFEM3D_PSPARSE_APAI_H

#include "PGFEM_mpi.h"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "comm_hints.h"
#include "data_structure.h"
#include "element.h"
#include "pgfem_comm.h"

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
                   Element *elem,
                   BoundingElement *b_elems,
                   Node *node,
                   int *Ap,
                   long nce,
                   COEL *coel,
                   long *DomDof,
                   int *GDof,
                   COMMUN comm,
                   MPI_Comm Comm_Orig,
                   const int cohesive,
                   const Comm_hints *hints,
                   const int mp_id);

#endif /* #define PGFEM3D_PSPARSE_APAI_H */
