#ifndef PGFEM3D_NODE_H
#define PGFEM3D_NODE_H

#include "pgfem3d/Communication.hpp"
#include "PGFEM_io.h"

struct NodeIDMap {
  long *id;
  long *Gid;
};

/** Structure of node properties */
struct Node {
  /** Coordinates of node */
  double x1,x2,x3,x1_fd,x2_fd,x3_fd;

  long loc_id; //Local node number, required for sorting
  NodeIDMap *id_map; // id maps for physics

  long Gnn; /**< Global node number. */
  long Dom; /**< Which processor the node lives on. This is how
               interfacial nodes are identified. */
  long Pr;  /**< Identifies the node as a periodic node. If the node is
               periodic, this variable contains the node number it is periodic
               with/to. */
  //    long ndofn;/**< Number of DOFs on the node */
  long pr;   /**< Property of node. */

  int model_type;
  int model_id;
};

Node* build_node(const long nn,
                 const int ndofn);

void destroy_node(const long nn,
                  Node* node);

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
                               const int physicsno);

/// destroy node array.
///
/// \param[in] nn number of nodes
/// \param[in] physicsno number of physics
///
/// \return non-zero on internal error
int destroy_node_multi_physics(const long nn,
                               Node* node,
                               const int physicsno);


/* Function reads parameters of nodes */
long read_nodes (FILE *in,
                 const long nn,
                 Node *node,
                 const int legacy,
                 const pgfem3d::CommunicationStructure *com);

/** */
void write_node_fname(const char *filename,
                      const int nnodes,
                      const Node *nodes,
                      const int ndofn,
                      const int mp_id);

void write_node(FILE *ofile,
                const int nnodes,
                const Node *nodes,
                const int ndofn,
                const int mp_id);

/**
 * Sort the nodes by their local node numbers (loc_id).
 */
void nodes_sort_loc_id(const int nnode,
                       Node *nodes);

/**
 * Sort the nodes by their ownership (Dom) with ascending global
 * node numbers (Gnn) then ascending local node numbers (loc_id).
 *
 * NOTE: The nodes MUST be re-sorted using `nodes_sort_loc_id` when
 * this ordering is no longer needed. Chaos will ensue otherwise.
 */
void nodes_sort_own_Gnn_loc(const int nnode,
                            Node *nodes);

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
                               int range[2]);

/**
 * Filter the shared nodes and return a pointer to the beginning of
 * the list of shared nodes and the number of shared nodes.
 *
 * Leaves the node list sorted in the following way: local nodes (Gnn
 * < 0) in loc_id order, shared nodes sorted by owning->Gnn->loc_id
 *
 * The node list *MUST* be re-sorted using nodes_sort_loc_id before
 * calling any FE routines.
 *
 * Upon completion `shared` points to the beginning of the shared
 * node list. `n_shared` is the number of shared nodes. Note that
 * there is no allocation or copying involved.
 *
 * \return non-zero on internal error
 */
int nodes_filter_shared_nodes(const int nnode,
                              Node *nodes,
                              int *n_shared,
                              const Node **shared);

#endif /* #define PGFEM3D_NODE_H */
