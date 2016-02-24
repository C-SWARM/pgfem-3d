/* HEADER */

#pragma once
#ifndef NODE_H
#define NODE_H

#include "PGFEM_mpi.h"
#include "PGFEM_io.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of node properties */
  struct NODE{
    /** Coordinates of node */
    double x1,x2,x3,x1_fd,x2_fd,x3_fd;

    long loc_id; /**< Local node number, required for sorting */

    long *id, /**< Local ID of DOFs on node. */

      *Gid, /**< Global ID of DOFs on node. */

      Gnn, /**< Global node number. */

      Dom, /**< Which processor the node lives on. This is how
	      interfacial nodes are identified. */

      Pr, /**< Identifies the node as a periodic node. If the node is
	     periodic, this variable contains the node number it is periodic
	     with/to. */

      ndofn, /**< Number of DOFs on the node */

      pr; /**< Property of node. */

    int model_type;
    int model_id;
  };
  typedef struct NODE NODE;

  NODE* build_node(const long nn,
		   const long ndofn);

  void destroy_node(const long nn,
		    NODE* node);

  /* Function reads parameters of nodes */
  long read_nodes (FILE *in,
		   const long nn,
		   NODE *node,
		   const int legacy,
		   MPI_Comm comm);

  /** */
  void write_node_fname(const char *filename,
			const int nnodes,
			const NODE *nodes);

  void write_node(FILE *ofile,
		  const int nnodes,
		  const NODE *nodes);

  /**
   * Sort the nodes by their local node numbers (loc_id).
   */
  void nodes_sort_loc_id(const int nnode,
                         NODE *nodes);

  /**
   * Sort the nodes by their ownership (Dom) with ascending global
   * node numbers (Gnn) then ascending local node numbers (loc_id).
   *
   * NOTE: The nodes MUST be re-sorted using `nodes_sort_loc_id` when
   * this ordering is no longer needed. Chaos will ensue otherwise.
   */
  void nodes_sort_own_Gnn(const int nnode,
                          NODE *nodes);

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
  int nodes_get_Gnn_idx_range(const int nnode,
                              const NODE *nodes,
                              const int dom,
                              int range[2]);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef NODE_H */
