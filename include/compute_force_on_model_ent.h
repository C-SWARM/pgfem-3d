/* HEADER */
#ifndef COMPUTE_FORCE_ON_MODEL_ENT
#define COMPUTE_FORCE_ON_MODEL_ENT

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#include "PGFEM_mpi.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

  typedef struct{
    int n_entities;
    int n_dim;
    int *ent_type;
    int *ent_id;
    int *node_ent;
    int *elem_related;
  } MODEL_ENTITY;

  /** Read the entity list from a file with a given path */
  int read_model_entity_list(const char *filename,
			     MODEL_ENTITY **me,
			     const int nnodes,
			     const int nelems);

  /** read the entity list from an already open stream */
  int read_model_entity_listf(FILE *in,
			      MODEL_ENTITY **me,
			      const int nnodes,
			      const int nelems);

  /** Allocate the Model Entity object */
  int create_model_entity(const int n_entities,
			  const int n_dim,
			  const int nnode,
			  const int nelem,
			  MODEL_ENTITY **me);

  /** Destroy a model entity object */
  int destroy_model_entity(MODEL_ENTITY *me);

  /** Mark nodes and elements as associated with model entities for
      fast search later */
  int model_entity_mark_nodes_elems(const int nnodes,
				    const NODE *nodes,
				    const int nelems,
				    const ELEMENT *elems,
				    MODEL_ENTITY *me);

  /** Compute the incrememental force on the Model Entity object. Note
      that this is called after every iteration and updates the vector
      'force'. The user must pass the total, updated increment, and
      current increment in displacements to the function. */
  int compute_force_increment_on_model_entity(const int nelems,
					      const int ndofn,
					      const double *total_disp_n,
					      const double *disp_k,
					      const double *del_disp_k,
					      const MODEL_ENTITY *me,
					      const ELEMENT *elems,
					      const NODE *nodes,
					      const SUPP sup,
					      const HOMMAT *hommat,
					      const EPS *eps,
					      const SIG *sig,
					      const MPI_Comm mpi_comm,
					      const int analysis,
					      double *forces);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef COMPUTE_FORCE_ON_MODEL_ENT */
