/* HEADER */
#pragma once
#ifndef STIFFMAT_FD_H
#define STIFFMAT_FD_H

#include "PGFEM_mpi.h"
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
#include "PGFem3D_data_structure.h"
#include "femlib.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Computes an element stiffness matrices Each type of analysis will
   * be excuted in this function Contribution to the residual
   * calculated using this function when prescribed displacement is
   * applied.
  */
  int el_compute_stiffmat(int i,
                          double *lk,
                          long ndofn,
                          long nne,
                          long npres,
                          int nVol,
                          int nsd,
                          ELEMENT *elem,
                          NODE *node,
                          HOMMAT *hommat,
                          MATGEOM matgeom,
                          SIG *sig,
                          EPS *eps,
                          SUPP sup,
                          double dt,
                          double nor_min,
                          double stab,
                          CRPL *crpl,
                          long FNR,
                          double lm,
                          double *x,
                          double *y,
                          double *z,
                          double *fe,
                          long *nod,
                          double *r_n,
                          double *r_e,
                          double alpha,
                          int include_inertia,
                          const int analysis,
                          const int cm);
	                      
  /**
   * Computes element stiffness matrices and assembles local
   * part. Off-process portions of the matrix are communicated via
   * non-blocking point-to-point send/receives using information in
   * COMMUN. Elements with global DOFs are computed first to overlap
   * communication with computation of fully local elements.
   */
  int stiffmat_fd (int *Ap,
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
		   const PGFem3D_opt *opts,double alpha, double *r_n, double *r_n_1,
		   const int mp_id);

int el_compute_stiffmat_MP(FEMLIB *fe,
                           double *lk,
                           GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *fv,
                           SOLVER_OPTIONS *sol,
                           LOADING_STEPS *load,
                           CRPL *crpl,
                           const PGFem3D_opt *opts,
                           MULTIPHYSICS *mp,
                           int mp_id,
                           double dt,
                           double lm,
                           double *be,
                           double *r_e);
		   
/// Compute stiffnes
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com communication object
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dt time step
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int stiffmat_fd_MP(GRID *grid,
                   MATERIAL_PROPERTY *mat,
                   FIELD_VARIABLES *fv,
                   SOLVER_OPTIONS *sol,
                   LOADING_STEPS *load,
                   COMMUNICATION_STRUCTURE *com,
                   CRPL *crpl,
                   MPI_Comm mpi_comm,
                   const PGFem3D_opt *opts,
                   MULTIPHYSICS *mp,
                   int mp_id,
                   double dt,
                   long iter,
                   int myrank);
                   
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
