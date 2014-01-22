#ifndef STIFFMAT_FD_H
#define STIFFMAT_FD_H

#include "PGFEM_mpi.h"
#include "BSprivate.h"

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef PGFEM_COMM_H
#include "pgfem_comm.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifndef HYPRE_GLOBAL_H
#include "hypre_global.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Computes global stiffness matrix.
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
