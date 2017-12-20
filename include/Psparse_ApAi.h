/**
 * AUTHORS:
 *    Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#ifndef PGFEM3D_PSPARSE_APAI_H
#define PGFEM3D_PSPARSE_APAI_H

#include "bounding_element.h"
#include "cohesive_element.h"
#include "comm_hints.h"
#include "data_structure.h"
#include "element.h"
#include "pgfem3d/Communication.hpp"

/**
 * Create the global sparsity pattern and commincation structure.
 */
int* Psparse_ApAi (long ne,
                   long n_be,
                   long nn,
                   long ndofn,
                   long ndofd,
                   Element *elem,
                   BoundingElement *b_elems,
                   Node *node,
                   long nce,
                   COEL *coel,
                   pgfem3d::CommunicationStructure *com,
                   const int cohesive,
                   const int mp_id);

#endif /* #define PGFEM3D_PSPARSE_APAI_H */
