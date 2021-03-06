#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "get_dof_ids_on_elem.h"
#include <cstring>

/* TESTED MM 1/18/2013 */
void get_all_dof_ids_on_elem(const int global,/* this is a boolean flag */
                             const long n_nodes_on_elem,
                             const int ndof_on_elem,
                             const int ndof_per_node,
                             const long *node_ids_on_elem,
                             const Node *nodes,
                             const BoundingElement *b_elems,
                             const Element *ptr_cur_elem,
                             long *dof_ids,
                             const int mp_id)
{
  const long ndof_on_elem_nodes = n_nodes_on_elem*ndof_per_node;

  /* node dof ids go in the first part of the dof_ids vector */
  get_dof_ids_on_elem_nodes(global,n_nodes_on_elem,ndof_per_node,
                            node_ids_on_elem,nodes,dof_ids, mp_id);

  /* element dofs are appended to the list in dof_ids */
  long *ptr_elem_dof_ids = dof_ids + ndof_on_elem_nodes;
  get_dof_ids_on_elem(global,ptr_cur_elem,ptr_elem_dof_ids);

  long *ptr_be_dof_ids = ptr_elem_dof_ids + ptr_cur_elem->n_dofs;
  get_dof_ids_on_all_elem_be(global,b_elems,ptr_cur_elem,ptr_be_dof_ids);
}

/* TESTED MM 1/18/2013 */
void get_dof_ids_on_elem_nodes(const int global,
                               const long n_nodes_on_elem,
                               const int n_dof_per_node, 
                               const long *node_ids_on_elem,
                               const Node *nodes,
                               long *dof_ids,
                               const int mp_id)
{
  const size_t n_copy = n_dof_per_node * sizeof(long);
  if(global){
    for(long i=0; i<n_nodes_on_elem; i++){
      long *ptr_cur = dof_ids + i*n_dof_per_node;
      memcpy(ptr_cur,nodes[ node_ids_on_elem[i] ].id_map[mp_id].Gid,n_copy);
    }
  } else {
    for(long i=0; i<n_nodes_on_elem; i++){
      long *ptr_cur = dof_ids + i*n_dof_per_node;
      memcpy(ptr_cur,nodes[ node_ids_on_elem[i] ].id_map[mp_id].id,n_copy);
    }
  }
}

/* TESTED MM 1/18/2013 */
void get_dof_ids_on_elem(const int global,
                         const Element *ptr_cur_elem,
                         long *dof_ids)
{
  const size_t n_copy = ptr_cur_elem->n_dofs * sizeof(long);
  if(n_copy > 0){
    if(global){
      memcpy(dof_ids,ptr_cur_elem->G_dof_ids,n_copy);
    } else {
      memcpy(dof_ids,ptr_cur_elem->L_dof_ids,n_copy);
    }
  }
}

void get_dof_ids_on_all_elem_be(const int global,
                                const BoundingElement *b_elems,
                                const Element *ptr_cur_elem,
                                long *dof_ids)
{
  long *ptr_copy = dof_ids;
  for(int i=0; i<ptr_cur_elem->n_be; i++){
    const int be_id = ptr_cur_elem->be_ids[i];
    const BoundingElement *ptr_be = &b_elems[be_id];
    const int n_be_dof = ptr_be->n_dofs;
    const size_t n_copy = n_be_dof*sizeof(long);
    if(global){
      memcpy(ptr_copy,ptr_be->G_dof_ids,n_copy);
    } else {
      memcpy(ptr_copy,ptr_be->L_dof_ids,n_copy);
    }
    ptr_copy += n_be_dof;
  }
}

/* TESTED 2/14/2013 MM */
void get_dof_ids_on_bnd_elem(const int global,
                             const int n_dof_per_node,
                             const Node *nodes,
                             const BoundingElement *ptr_cur_be,
                             const Element *elems,
                             long *dof_ids,
                             const int mp_id)
{
  const int elem_id = ptr_cur_be->vol_elem_id;
  const Element *ptr_elem = &elems[elem_id];
  const int nnodes = ptr_elem->toe;
  const long *elem_nodes = ptr_elem->nod;

  /* get dof ids on elemnodes */
  get_dof_ids_on_elem_nodes(global,nnodes,n_dof_per_node,elem_nodes,nodes,dof_ids,mp_id);

  long *ptr_dof_ids = dof_ids + nnodes*n_dof_per_node;
  const size_t n_copy = ptr_cur_be->n_dofs * sizeof(long);
  if(global){
    memcpy(ptr_dof_ids,ptr_cur_be->G_dof_ids,n_copy);
  } else {
    memcpy(ptr_dof_ids,ptr_cur_be->L_dof_ids,n_copy);
  }
}
