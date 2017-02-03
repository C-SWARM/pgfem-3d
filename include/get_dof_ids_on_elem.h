/* HEADER */
#ifndef GET_DOF_IDS_ON_ELEM_H
#define GET_DOF_IDS_ON_ELEM_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Get the id numbers for all degrees of freedom which need to be
      assembled from the element. */
  void get_all_dof_ids_on_elem(const int global,/* this is a boolean flag */
			       const int n_nodes_on_elem,
			       const int ndof_on_elem,
			       const int ndof_per_node,
			       const long *node_ids_on_elem,
			       const NODE *nodes,
			       const BOUNDING_ELEMENT *b_elems,
			       const ELEMENT *ptr_cur_elem,
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
			   const ELEMENT *ptr_cur_elem,
			   long *dof_ids);

  /** get the dof ids associated with all bounding elements on the volume element */
  void get_dof_ids_on_all_elem_be(const int global,
				  const BOUNDING_ELEMENT *b_elems,
				  const ELEMENT *ptr_cur_elem,
				  long *dof_ids);

  /** getr the dof ids associated with a particular bounding element
      and the nodal dofs on the bounded volume element. DOES NOT
      include element dofs or dofs associated with any additional
      bounding elements on the volume element */
  void get_dof_ids_on_bnd_elem(const int global,
			       const int n_dof_per_node,
			       const NODE *nodes,
			       const BOUNDING_ELEMENT *ptr_cur_be,
			       const ELEMENT *elems,
			       long *dof_ids,
			       const int mp_id);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef GET_DOF_IDS_ON_ELEM_H */
