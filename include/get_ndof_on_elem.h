/* HEADER */
#ifndef GET_NDOF_ON_ELEM_H
#define GET_NDOF_ON_ELEM_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

  /** Get the total number of degrees of freedom on an element. This
      includes degrees of freedom from each of the element nodes and
      element-wise degrees of freedom which cannot be condensed
      out as well as from bounding elements. */
  int get_total_ndof_on_elem(const int n_nodes_on_elem,
			     const long *node_ids_on_elem,
			     const NODE *nodes,
			     const BOUNDING_ELEMENT *b_elems,
			     const ELEMENT *ptr_cur_elem);

  /** Get the number of dofs associated with the nodes only */
  int get_ndof_on_elem_nodes(const int n_nodes_on_elem,
			     const long *node_ids_on_elem,
			     const NODE *nodes);

  /** Get the number of degrees of freedom associated with all
      bounding elements on the element. */
  int get_ndof_on_all_elem_be(const BOUNDING_ELEMENT *b_elems,
			      const ELEMENT *ptr_cur_elem);

  /** Get the number of degrees of freedom associated with a bounding
      element. This includes the dofs on the bounding element itself
      and dofs associated with the volume element nodes */
  int get_ndof_on_bnd_elem(const NODE *nodes,
			   const BOUNDING_ELEMENT *ptr_cur_be,
			   const ELEMENT *elems);


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef GET_NDOF_ON_ELEM */
