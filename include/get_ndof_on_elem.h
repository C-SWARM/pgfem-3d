#ifndef PGFEM3D_GET_NDOF_ON_ELEM_H
#define PGFEM3D_GET_NDOF_ON_ELEM_H

#include "bounding_element.h"
#include "element.h"
#include "node.h"

/** Get the total number of degrees of freedom on an element. This
    includes degrees of freedom from each of the element nodes and
    element-wise degrees of freedom which cannot be condensed
    out as well as from bounding elements. */
int get_total_ndof_on_elem(const int n_nodes_on_elem,
                           const long *node_ids_on_elem,
                           const Node *nodes,
                           const BoundingElement *b_elems,
                           const Element *ptr_cur_elem,
                           const int ndofn);

/** Get the number of dofs associated with the nodes only */
int get_ndof_on_elem_nodes(const int n_nodes_on_elem,
                           const long *node_ids_on_elem,
                           const Node *nodes,
                           const int ndofn);

/** Get the number of degrees of freedom associated with all
    bounding elements on the element. */
int get_ndof_on_all_elem_be(const BoundingElement *b_elems,
                            const Element *ptr_cur_elem);

/** Get the number of degrees of freedom associated with a bounding
    element. This includes the dofs on the bounding element itself
    and dofs associated with the volume element nodes */
int get_ndof_on_bnd_elem(const Node *nodes,
                         const BoundingElement *ptr_cur_be,
                         const Element *elems,
                         const int ndofn);

#endif /* #define PGFEM3D_GET_NDOF_ON_ELEM_H */
