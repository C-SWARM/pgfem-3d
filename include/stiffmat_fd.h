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
#include "macro_micro_functions.h"

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
                   
/// Multiscale simulation interface to compute stiffness matrix
///
/// \param[in] c structure of macroscale information
/// \param[in,out] s contains the information for the history-dependent solution
/// \param[in] opts structure PGFem3D option
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] nor_min nonlinear convergence tolerance
/// \param[in] FNR if 1: Full Newton-Raphson
///                   0: only compute stiffnes at the 1st iteration
/// \param[in] myrank current process rank
/// \param[in] nproc   number of total process
/// \return non-zero on internal error
int stiffmat_fd_multiscale(COMMON_MACROSCALE *c,
                           MACROSCALE_SOLUTION *s,
                           const PGFem3D_opt *opts,
                           long iter,
                           double nor_min,
                           long FNR,
                           int myrank,
                           int nproc);                   
                   
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
