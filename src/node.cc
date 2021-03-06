/** This file defines some routines related to the node structure such
    as allocation, deallocation, reading and writing */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "node.h"
#include "allocation.h"
#include "utils.h"
#include <cassert>

using namespace pgfem3d;
using namespace multiscale::net;

Node* build_node(const long nn,
                 const int ndofn)
{
  return build_node_multi_physics(nn,&ndofn,1);
}

void destroy_node(const long nn,
                  Node* node)
{
  destroy_node_multi_physics(nn, node, 1);
}

/// build node array.
/// Multiphysics needs many ids for nodal variables.
/// To assign local and global ids according to the number physics,
/// an array of id_map object is created as many as number of physics.
/// In id_map, ids is also created based on the number of dofs for
/// the individual physics.
///
/// \param[in] nn number of nodes
/// \param[in] ndofn number of dofs on a node
/// \param[in] physicsno number of physics
/// \return node array
Node* build_node_multi_physics(const long nn,
                               const int *ndofn,
                               const int physicsno)
{
  Node *pom = PGFEM_calloc (Node, nn);
  for (int i = 0 ; i < nn; ++i) {
    pom[i].id_map = new NodeIDMap[physicsno];
    for (int j = 0; j < physicsno; ++j) {
      pom[i].id_map[j].id  = PGFEM_calloc (long, ndofn[j]);
      pom[i].id_map[j].Gid = PGFEM_calloc (long, ndofn[j]);
    }
  }
  return pom;
}

/// destroy node array.
///
/// \param[in] nn number of nodes
/// \param[in] physicsno number of physics
///
/// \return non-zero on internal error
int destroy_node_multi_physics(const long nn,
                               Node* node,
                               const int physicsno)
{
  int err = 0;
  for(int ia=0; ia<nn; ia++)
  {
    for(int ib=0; ib<physicsno; ib++)
    {
      PGFEM_free(node[ia].id_map[ib].id);
      PGFEM_free(node[ia].id_map[ib].Gid);
    }
    delete [] node[ia].id_map;
  }
  PGFEM_free(node);
  return err;
}

long read_nodes (FILE *in,
                 const long nn,
                 Node *node,
                 const int legacy,
                 const CommunicationStructure *com)
/*
  in   - Input file
  nn   - Number of nodes
  node - Structure type of Node

  returns: total number of nodes over all domains, counting nodes
  on boundries only once
  %%%%%%%%%%%%%%%% TESTED 6.12.99 %%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%% MODIFIED 7.20.05 %%%%%%%%%%%%%%%
*/
{
  int myrank = com->rank;
  long tnn = nn;
  for (long i = 0; i < nn; ++i) {
    long id = 0;
    long Gnn = 0;
    long Dom = 0;
    CHECK_SCANF(in,"%ld %ld %ld",&Gnn,&Dom,&id);

    assert(0 <= id and id < nn);
    assert(0 <= Dom and Dom < com->nproc);

    Node *p_node = &node[id];
    p_node->loc_id = id;
    p_node->Gnn = Gnn;
    p_node->Dom = Dom;

    /* Input file error checking */
    if (p_node->Gnn < 0 && p_node->Dom != myrank) {
      PGFEM_printerr("[%d] ERROR: incorrect node domain info (node %ld)!"
                     " %s:%s:%d\n",myrank,i,__func__,__FILE__,__LINE__);
      PGFEM_Abort();
    }

    /* If we get a global node that doesn't live on this domain,
       subtract it from tnn */
    if (p_node->Gnn != -1 && p_node->Dom != myrank) {
      assert(p_node->Gnn >= 0);
      tnn--;
    }

    if (legacy) {
      CHECK_SCANF(in,"%lf %lf %lf %ld",
                  &p_node->x1,
                  &p_node->x2,
                  &p_node->x3,
                  &p_node->pr);
    } else {
      CHECK_SCANF(in,"%lf %lf %lf %d %d %ld",
                  &p_node->x1,
                  &p_node->x2,
                  &p_node->x3,
                  &p_node->model_type,
                  &p_node->model_id,
                  &p_node->pr);
    }

    p_node->x1_fd = p_node->x1;
    p_node->x2_fd = p_node->x2;
    p_node->x3_fd = p_node->x3;

    /* error check read */
    if (ferror(in)) {
      PGFEM_printerr("[%d]ERROR:fscanf returned error"
                     " reading node %ld!\n",myrank,i);
      PGFEM_Abort();
    }

    if (feof(in)) {
      PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
                     myrank);
      PGFEM_Abort();
    }
  }

  /* Gather tnn from all domains */
  long Gtnn = 0;
  com->net->allreduce(&tnn,&Gtnn,1,NET_DT_LONG,NET_OP_SUM,com->comm);
  assert(0 < Gtnn);
  return Gtnn;
}

void write_node_fname(const char *filename,
                      const int nnodes,
                      const Node *nodes,
                      const int ndofn,
                      const int mp_id)
{
  FILE *ofile = fopen(filename,"w");
  if(ofile == NULL){
    PGFEM_printerr("ERROR: cannot open file %s in %s\n",filename,__func__);
    PGFEM_Abort();
  }

  write_node(ofile,nnodes,nodes,mp_id,ndofn);

  fclose(ofile);
}

void write_node(FILE *ofile,
                const int nnodes,
                const Node *nodes,
                const int ndofn,
                const int mp_id)
{
  /* write header describing format */
  PGFEM_fprintf(ofile,"  Gnn DOM   Lnn       X            Y            Z"
                "              X_fd         Y_fd         Z_fd         Lid::Gid ...\n");
  PGFEM_fprintf(ofile,"===================================================="
                "================================================================\n");
  for(int i=0; i<nnodes; i++){
    const Node *p_node = &nodes[i];
    PGFEM_fprintf(ofile,"%5ld %3ld %5d ",p_node->Gnn,p_node->Dom,i);
    PGFEM_fprintf(ofile,"%12.5e %12.5e %12.5e    %12.5e %12.5e %12.5e    ",
                  p_node->x1,p_node->x2,p_node->x3,
                  p_node->x1_fd,p_node->x2_fd,p_node->x3_fd);
    for(int j=0; j<ndofn; j++){
      PGFEM_fprintf(ofile,"%5ld::%-5ld ",p_node->id_map[mp_id].id[j],p_node->id_map[mp_id].Gid[j]);
    }
    PGFEM_fprintf(ofile,"\n");
  }
}

static int node_comp_loc_id(const void *a,
                            const void *b)
{
  return (((Node*) a)->loc_id - ((Node*) b)->loc_id);
}

static int node_comp_own(const void *a,
                         const void *b)
{
  return (((Node*) a)->Dom - ((Node*) b)->Dom);
}

static int node_comp_Gnn(const void *a,
                         const void *b)
{
  return (((Node*) a)->Gnn - ((Node*) b)->Gnn);
}

static int node_comp_own_Gnn(const void *a,
                             const void *b)
{
  int own = node_comp_own(a,b);
  if (own) return own;
  else return node_comp_Gnn(a,b);
}

static int node_comp_Gnn_loc(const void *a,
                             const void *b)
{
  int Gnn = node_comp_Gnn(a,b);
  if (Gnn) return Gnn;
  else return node_comp_loc_id(a,b);
}

static int node_comp_own_Gnn_loc(const void *a,
                                 const void *b)
{
  int own_gnn = node_comp_own_Gnn(a,b);
  if (own_gnn) return own_gnn;
  else return node_comp_loc_id(a,b);
}

void nodes_sort_loc_id(const int nnode,
                       Node *nodes)
{
  qsort(nodes, nnode, sizeof(*nodes), node_comp_loc_id);
}

void nodes_sort_own_Gnn_loc(const int nnode,
                            Node *nodes)
{
  qsort(nodes, nnode, sizeof(*nodes), node_comp_own_Gnn_loc);
}

/**
 * Compute the index range of Global Nodes owned by the specified domain.
 *
 * For valid results, requires the nodes to be sorted by
 * `nodes_sort_own_Gnn`. Index range may include duplicate global node
 * numbersin the case of periodicity.
 *
 * \return non-zero if no shared/global nodes are owned by the
 * specified domain. On success, `range` specifies matches in
 * [range[0], range[1]).
 */
int nodes_get_shared_idx_range(const int nnode,
                               const Node *nodes,
                               const int dom,
                               int range[2])
{
  int err = 0;
  /* create a node for comparison */
  Node comp_node = {0};
  comp_node.Dom = dom;

  /* search for a node with matching ownership */
  const Node *ptr_lb = static_cast<const Node*>(bsearch(&comp_node,
                                                        nodes,
                                                        nnode,
                                                        sizeof(*nodes),
                                                        node_comp_own));

  /* exit early if no match found */
  if (!ptr_lb) return 1;

  /* linearly search for bounds */
  const Node *ptr_ub = ptr_lb;
  while (ptr_ub->Dom == dom) {
    /* limit the search by the length of the array */
    if((++ptr_ub - nodes) == nnode) break;
  }
  if (ptr_lb->Gnn < 0) {
    /* matched node is purely local -> search forward for first owned
       boundary node */
    while (ptr_lb->Gnn < 0 && ptr_lb->Dom == dom) {
      /* limit the search by the length of the array */
      if((++ptr_lb - nodes) == nnode) break;
    }
  } else {
    /* matched node is on the boundary -> search backward for first
       owned boundary node */
    while (ptr_lb->Gnn >= 0 && ptr_lb->Dom == dom) {
      /* limit the search by the length of the array */
      if((--ptr_lb - nodes) < 0) break;
    }

    /* Lower-bound is inclusive -> increment pointer */
    ++ptr_lb;
  }

  assert(ptr_lb != ptr_ub);
  if (ptr_lb == ptr_ub) err++;

  /* perform pointer arithmetic w.r.t nodes to get the index range */
  range[0] = ptr_lb - nodes;
  range[1] = ptr_ub - nodes;

  return err;
}

int nodes_filter_shared_nodes(const int nnode,
                              Node *nodes,
                              int *n_shared,
                              const Node **shared)
{
  int err = 0;

  /* sort by Gnn */
  qsort(nodes, nnode, sizeof(*nodes), node_comp_Gnn_loc);

  /* perform linear search from the end to find the beginning of the
     shared nodes list. We start from the end as typically there are
     fewer boundary nodes than local nodes. */
  Node *ptr = &nodes[nnode-1];
  while (ptr->Gnn >= 0) if ((--ptr - nodes) < 0) break;
  ++ptr; /* lower bound is inclusive, increment pointer */

  /* compute number of shared nodes and starting index */
  *n_shared = nnode - (ptr - nodes);
  *shared = ptr;
  assert(*n_shared >= 0); /* check for implementation error */

  /* re-sort shared nodes by own->Gnn->loc_id */
  qsort(ptr, *n_shared, sizeof(*ptr), node_comp_own_Gnn_loc);

  return err;
}
