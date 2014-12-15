/* HEADER */
#pragma once
#ifndef BOUNDING_ELEMENT_UTILS_H
#define BOUNDING_ELEMENT_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#include "PGFEM_mpi.h"
#include "element.h"
#include "node.h"
#include "bounding_element.h"
#include "sig.h"
#include "eps.h"
#include "pgfem_comm.h"

  /** get the element local node numbers, i.e. ordering in volumetric
      element connectivity */
  int bounding_element_set_local_ids(const int n_be,
				     BOUNDING_ELEMENT *b_elems,
				     const ELEMENT *elem);

  /** get the reverse mapping (elem -> b_elem) */
  int bounding_element_reverse_mapping(const int n_be,
				       const BOUNDING_ELEMENT *b_elems,
				       ELEMENT *elem);

  /** get the value of damage from the periodic bounding element. If
      the periodic element is on the other domain communicate the
      value. */
  int bounding_element_communicate_damage(const int n_be,
					  BOUNDING_ELEMENT *b_elems,
					  const int ne,
					  const EPS *eps,
					  const MPI_Comm mpi_comm);

  int bounding_element_compute_resulting_traction(const int n_be,
						  const BOUNDING_ELEMENT *b_elems,
						  const ELEMENT *elems,
						  const NODE *nodes,
						  const EPS *eps,
						  const SIG *sig,
						  const int ndofd,
						  const long *DomDof,
						  const int GDof,
						  const COMMUN comm,
						  const MPI_Comm mpi_comm,
						  const int analysis,
						  double *res_trac);

  /** Compute the permutation to transform the global stiffness into a
      block-ordered matrix for preconditioning */
  int compute_block_permutation(const int n_be,
				const int ndofn,
				const BOUNDING_ELEMENT *b_elems,
				const ELEMENT *elems,
				const NODE *nodes,
				const long *DomDof,
				const MPI_Comm mpi_comm,
				long *perm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef BOUNDING_ELEMENT_UTILS_H */
