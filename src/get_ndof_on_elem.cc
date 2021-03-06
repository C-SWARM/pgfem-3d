#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "get_ndof_on_elem.h"

/* TESTED MM 1/18/2013 */
int get_total_ndof_on_elem(const int n_nodes_on_elem,
                           const long *node_ids_on_elem,
                           const Node *nodes,
                           const BoundingElement *b_elems,
                           const Element *ptr_cur_elem,
                           const int ndofn)
{
  int ndof_on_elem = 0;

  ndof_on_elem += get_ndof_on_elem_nodes(n_nodes_on_elem,node_ids_on_elem,nodes,ndofn);

  ndof_on_elem += ptr_cur_elem->n_dofs;

  ndof_on_elem += get_ndof_on_all_elem_be(b_elems,ptr_cur_elem);

  return ndof_on_elem;
}

/* TESTED MM 1/18/2013 */
int get_ndof_on_elem_nodes(const int n_nodes_on_elem,
                           const long *node_ids_on_elem,
                           const Node *nodes,
                           const int ndofn)
{
  int ndof_on_elem = 0;
  for(int i=0;i<n_nodes_on_elem; i++){
    ndof_on_elem += ndofn;
  }
  return ndof_on_elem;
}

int get_ndof_on_all_elem_be(const BoundingElement *b_elems,
                            const Element *ptr_cur_elem)
{
  int ndof_on_elem = 0;
  for(int i=0; i<ptr_cur_elem->n_be; i++){
    const int be_id = ptr_cur_elem->be_ids[i];
    ndof_on_elem += b_elems[be_id].n_dofs;
  }

  return ndof_on_elem;
}

int get_ndof_on_bnd_elem(const Node *nodes,
                         const BoundingElement *ptr_cur_be,
                         const Element *elems,
                         const int ndofn)
{
  const int elem_id = ptr_cur_be->vol_elem_id;
  const Element *ptr_elem = &elems[elem_id];
  const int nnodes = ptr_elem->toe;
  const long *elem_nodes = ptr_elem->nod;

  int ndof_on_be = 0;
  ndof_on_be += get_ndof_on_elem_nodes(nnodes,elem_nodes,nodes,ndofn);
  ndof_on_be += ptr_cur_be->n_dofs;

  return ndof_on_be;
}
