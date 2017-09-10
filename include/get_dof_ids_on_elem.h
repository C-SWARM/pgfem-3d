#ifndef PGFEM3D_GET_DOF_IDS_ON_ELEM_H
#define PGFEM3D_GET_DOF_IDS_ON_ELEM_H

#include "bounding_element.h"
#include "element.h"
#include "node.h"

/** Get the id numbers for all degrees of freedom which need to be
    assembled from the element. */
void get_all_dof_ids_on_elem(const int global,/* this is a boolean flag */
                             const int n_nodes_on_elem,
                             const int ndof_on_elem,
                             const int ndof_per_node,
                             const long *node_ids_on_elem,
                             const NODE *nodes,
                             const BoundingElement *b_elems,
                             const Element *ptr_cur_elem,
                             long *dof_ids,
                             const int mp_id);

/** Get only the id numbers for degrees of freedom associated with
    nodes. This function can replace calls to 'codenum' when the
    element dofs are not needed. */
void get_dof_ids_on_elem_nodes(const int global,
                               const int n_nodes_on_elem,
                               const int n_dof_per_node,
                               const long *node_ids_on_elem,
                               const NODE *nodes,
                               long *dof_ids,
                               const int mp_id);

/** Get only the id numbers for degrees of freedom associated with the element */
void get_dof_ids_on_elem(const int global,
                         const Element *ptr_cur_elem,
                         long *dof_ids);

/** get the dof ids associated with all bounding elements on the volume element */
void get_dof_ids_on_all_elem_be(const int global,
                                const BoundingElement *b_elems,
                                const Element *ptr_cur_elem,
                                long *dof_ids);

/** getr the dof ids associated with a particular bounding element
    and the nodal dofs on the bounded volume element. DOES NOT
    include element dofs or dofs associated with any additional
    bounding elements on the volume element */
void get_dof_ids_on_bnd_elem(const int global,
                             const int n_dof_per_node,
                             const NODE *nodes,
                             const BoundingElement *ptr_cur_be,
                             const Element *elems,
                             long *dof_ids,
                             const int mp_id);

#endif /* #define PGFEM3D_GET_DOF_IDS_ON_ELEM_H */
