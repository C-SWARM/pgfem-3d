/* HEADER */
#ifndef SOLVE_SYSTEM_H
#define SOLVE_SYSTEM_H

#include "hypre_global.h"
#include "PGFEM_io.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  const enum {SOLVE_SUCCESS=0,
	      SOLVE_ERROR=1,
	      UNREC_SOLVER,
	      UNREC_SOLVER_TYPE,
	      UNREC_PRECOND,
	      BAD_PRECOND} PFEM_SOLVE_ERR;

  /** define the structure for general access to solver information */
  typedef struct SOLVER_INFO{
    double res_norm;
    int n_iter;
    int err;
  } SOLVER_INFO;

  double solve_system(const PGFem3D_opt *opts,
		      double *loc_rhs,
		      double *loc_sol,
		      const int tim,
		      const int iter,
		      const long *DomDof,
		      SOLVER_INFO *info,
		      PGFEM_HYPRE_solve_info *PGFEM_hypre,
		      MPI_Comm mpi_comm);

  /** Solve the system of equations without calling setup
      routines. Use this function for subsequent solves are made using
      the same matrix. */
  double solve_system_no_setup(const PGFem3D_opt *opts,
			       double *loc_rhs,
			       double *loc_sol,
			       const int tim,
			       const int iter,
			       const long *DomDof,
			       SOLVER_INFO *info,
			       PGFEM_HYPRE_solve_info *PGFEM_hypre,
			       MPI_Comm mpi_comm);

  int solve_system_check_error(FILE *out,
			       const SOLVER_INFO info);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef SOLVE_SYSTEM_H */
