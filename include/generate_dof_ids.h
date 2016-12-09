/* HEADER */
/**
 * \file Routines for assigning DOF ids on nodes and elements.
 *
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef GENERATE_DOF_IDS_H
#define GENERATE_DOF_IDS_H


#include "PGFEM_mpi.h"
#include "node.h"
#include "element.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "comm_hints.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */



  /** Generate the local dof id numbers and return the number of local
      dofs */
  int generate_local_dof_ids(const int nelem,
			     const int ncoel,
			     const int nnode,
			     const int ndofn,
			     NODE *nodes,
			     ELEMENT *elems,
			     COEL *coel,
			     BOUNDING_ELEMENT *b_elems,
			     MPI_Comm mpi_comm);

  /** Generate the global dof id numbers and return the number of dofs
      owned by the domain */
  int generate_global_dof_ids(const int nelem,
			      const int ncoel,
			      const int nnode,
			      const int ndofn,
			      NODE *nodes,
			      ELEMENT *elems,
			      COEL *coel,
			      BOUNDING_ELEMENT *b_elems,
			      MPI_Comm mpi_comm);

  /** Increment the global dof id numbers based on the number of dofs
      on the other domains. */
  void renumber_global_dof_ids(const int nelem,
			       const int ncoel,
			       const int n_belem,
			       const int nnode,
			       const int ndofn,
			       const long *n_G_dof_on_dom,
			       NODE *nodes,
			       ELEMENT *elems,
			       COEL *coel,
			       BOUNDING_ELEMENT *b_elems,
			       MPI_Comm mpi_comm);

  /** Redistributes the degrees of freedom on the boundary and returns
      the number of boundary nodes on the domain. */
  int distribute_global_dof_ids(const int nelem,
				const int ncoel,
				const int n_belem,
				const int nnode,
				const int ndofn,
				const int ndof_be,
				NODE *nodes,
				ELEMENT *elems,
				COEL *coel,
				BOUNDING_ELEMENT *b_elems,
                                const Comm_hints *hints,
				MPI_Comm mpi_comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef GENERATE_DOF_IDS_H */
