/* HEADER */
#pragma once
#ifndef STIFFMAT_FD_H
#define STIFFMAT_FD_H

#include "PGFEM_mpi.h"
#include "blocksolve_interface.h"
#include "element.h"
#include "node.h"
#include "matgeom.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "pgfem_comm.h"
#include "PGFem3D_options.h"
#include "hypre_global.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Computes element stiffness matrices and assembles local
   * part. Off-process portions of the matrix are communicated via
   * non-blocking point-to-point send/receives using information in
   * COMMUN. Elements with global DOFs are computed first to overlap
   * communication with computation of fully local elements.
   */
  int stiffmat_fd (BSspmat *K,
		   int *Ap,
		   int *Ai,
		   long ne,
		   int n_be,
		   long ndofn,
		   ELEMENT *elem,
		   BOUNDING_ELEMENT *b_elems,
		   long nbndel,
		   long *bndel,
		   NODE *node,
		   HOMMAT *hommat,
		   MATGEOM matgeom,
		   SIG *sig,
		   EPS *eps,
		   double *d_r,
		   double *r,
		   long npres,
		   SUPP sup,
		   long iter,
		   double nor_min,
		   double dt,
		   CRPL *crpl,
		   double stab,
		   long nce,
		   COEL *coel,
		   long FNR,
		   double lm,
		   double *f_u,
		   int myrank,
		   int nproc,
		   long *DomDof,
		   long GDof,
		   COMMUN comm,
		   MPI_Comm mpi_comm,
		     PGFEM_HYPRE_solve_info *PGFEM_hypre,
		   const PGFem3D_opt *opts);

/** Assemble non-local parts as they arrive */
int assemble_nonlocal_stiffmat(const COMMUN pgfem_comm,
			       MPI_Status *sta_r,
			       MPI_Request *req_r,
			       PGFEM_HYPRE_solve_info *PGFEM_hypre,
			       double **recv);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef STIFFMAT_FD_H */
