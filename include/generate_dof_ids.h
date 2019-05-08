/**
 * \file Routines for assigning DOF ids on nodes and elements.
 *
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#ifndef PGFEM3D_GENERATE_DOF_IDS_H
#define PGFEM3D_GENERATE_DOF_IDS_H

#include "pgfem3d/Communication.hpp"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "element.h"
#include "node.h"

/** Generate the local dof id numbers and return the number of local
    dofs */
long generate_local_dof_ids(const long nelem,
                            const long ncoel,
                            const long nnode,
                            const long ndofn,
                            Node *nodes,
                            Element *elems,
                            COEL *coel,
                            BoundingElement *b_elems,
			    const pgfem3d::CommunicationStructure *com,
                            const int mp_id);

/** Generate the global dof id numbers and return the number of dofs
    owned by the domain */
long generate_global_dof_ids(const long nelem,
                             const long ncoel,
                             const long nnode,
                             const long ndofn,
                             Node *nodes,
                             Element *elems,
                             COEL *coel,
                             BoundingElement *b_elems,
			     const pgfem3d::CommunicationStructure *com,
                             const int mp_id);

/** Increment the global dof id numbers based on the number of dofs
    on the other domains. */
void renumber_global_dof_ids(const long nelem,
                             const long ncoel,    /* UNUSED */
                             const int n_belem,
                             const long nnode,
                             const long ndofn,
                             const long *n_G_dof_on_dom,
                             Node *nodes,
                             Element *elems,
                             COEL *coel,
                             BoundingElement *b_elems,
			     const pgfem3d::CommunicationStructure *com,
                             const int mp_id);

/** Redistributes the degrees of freedom on the boundary and returns
    the number of boundary nodes on the domain. */
long distribute_global_dof_ids(const long nelem,   /* UNUSED */
                               const long ncoel,   /* UNUSED */
                               const int n_belem,
                               const long nnode,
                               const long ndofn,
                               const int ndof_be,
                               Node *nodes,
                               Element *elems,
                               COEL *coel,
                               BoundingElement *b_elems,
			       const pgfem3d::CommunicationStructure *com,
                               const int mp_id);

#endif /* #define PGFEM3D_GENERATE_DOF_IDS_H */
