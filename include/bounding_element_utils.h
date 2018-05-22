#ifndef PGFEM3D_BOUNDING_ELEMENT_UTILS_H
#define PGFEM3D_BOUNDING_ELEMENT_UTILS_H

#include "PGFem3D_data_structure.h"
#include "bounding_element.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "node.h"
#include "sig.h"

/** get the element local node numbers, i.e. ordering in volumetric
    element connectivity */
int bounding_element_set_local_ids(const int n_be,
                                   BoundingElement *b_elems,
                                   const Element *elem);

/** get the reverse mapping (elem -> b_elem) */
int bounding_element_reverse_mapping(const int n_be,
                                     const BoundingElement *b_elems,
                                     Element *elem);

/** get the value of damage from the periodic bounding element. If
    the periodic element is on the other domain communicate the
    value. */
int bounding_element_communicate_damage(const int n_be,
                                        BoundingElement *b_elems,
                                        const int ne,
                                        const EPS *eps,
                                        const pgfem3d::CommunicationStructure *com);

int bounding_element_compute_resulting_traction(const int n_be,
                                                const BoundingElement *b_elems,
                                                const Element *elems,
                                                const Node *nodes,
                                                const EPS *eps,
                                                const SIG *sig,
                                                const int ndofd,
                                                const long *DomDof,
                                                const int GDof,
                                                const pgfem3d::CommunicationStructure *com,
                                                const int analysis,
                                                double *res_trac);

/** Compute the permutation to transform the global stiffness into a
    block-ordered matrix for preconditioning */
int compute_block_permutation(const int n_be,
                              const int ndofn,
                              const BoundingElement *b_elems,
                              const Element *elems,
                              const Node *nodes,
                              const long *DomDof,
			      const pgfem3d::CommunicationStructure *com,
                              long *perm,
                              const int mp_id);

#endif /* #define PGFEM3D_BOUNDING_ELEMENT_UTILS_H */
