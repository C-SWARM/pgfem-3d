/* HEADER */
#ifndef SOLVE_SYSTEM_H
#define SOLVE_SYSTEM_H

#include "BSprivate.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef HYPRE_GLOBAL_H
#include "hypre_global.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

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
		      BSprocinfo *BSinfo,
		      BSspmat *k,
		      BSpar_mat **pk,
		      BSpar_mat **f_pk,
		      MPI_Comm mpi_comm);

  int solve_system_check_error(FILE *out,
			       const SOLVER_INFO info);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef SOLVE_SYSTEM_H */
