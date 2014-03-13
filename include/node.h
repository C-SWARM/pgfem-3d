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

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef NODE_H */
