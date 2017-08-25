#pragma once
#ifndef HYPRE_GLOBAL_H
#define HYPRE_GLOBAL_H

#include "PGFEM_mpi.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"

#include "incl.h"
#include "PGFem3D_options.h"

/* HYPRE */

/**
 * @name Structure to contain HYPRE objects and pass them around easily
 */
namespace pgfem3d {
namespace solvers {
struct Hypre {
  HYPRE_IJMatrix       hypre_k = nullptr;       //!< The matrix handle
  HYPRE_ParCSRMatrix  hypre_pk = nullptr;       //!< The matrix data handle
  HYPRE_Solver    hypre_solver = nullptr;       //!< The solver
  HYPRE_Solver        hypre_pc = nullptr;       //!< The preconditioner
  HYPRE_Solver hypre_pc_gotten = nullptr;       //!<
  HYPRE_IJVector     hypre_rhs = nullptr;       //!< The RHS vector handle
  HYPRE_ParVector   hypre_prhs = nullptr;       //!< The RHS vector data handle
  HYPRE_IJVector     hypre_sol = nullptr;       //!< The solution vector handle
  HYPRE_ParVector   hypre_psol = nullptr;       //!< The solution vector data handle

  int  *ncol = nullptr;                         //!<
  int *grows = nullptr;                         //!<
  int ilower = 0;                               //!< Row lower bound
  int iupper = 0;                               //!< Row upper bound
  int jlower = 0;                               //!< Column lower bound
  int jupper = 0;                               //!< Column upper bound

  int precond_type = 0;                         //!< The type of preconditioner
  int  solver_type = 0;                         //!< The type of solver

  ~Hypre();

  void initialize(int *Ap,
                  int *Ai,
                  int size,
                  int maxit,
                  double err,
                  const PGFem3D_opt *options,
                  MPI_Comm mpi_comm);

  void zero(int *Ai, int size) {
    setK(Ai, size, 0);
  }

  void setRowColBounds(long g_n_col, const long* n_row_proc, int rank);

  int createPreconditioner(MPI_Comm comm);
  int destroyPreconditioner();

 private:
  void setK(int *Ai, int size, double val);
}; // struct Hypre
} // namespace solvers
} // namespace pgfem3d

using PGFEM_HYPRE_solve_info = pgfem3d::solvers::Hypre;

/* HYPRE utils */
static inline void
hypre_initialize(int *Ap, int *Ai, int size, int maxit, double err,
                 PGFEM_HYPRE_solve_info *PGFEM_hypre,
                 const PGFem3D_opt *options, MPI_Comm mpi_comm)
{
  PGFEM_hypre->initialize(Ap, Ai, size, maxit, err, options, mpi_comm);
}

static inline void
ZeroHypreK(PGFEM_HYPRE_solve_info *PGFEM_hypre, int *Ai, int size)
{
  PGFEM_hypre->zero(Ai, size);
}

/** Sets the extends of ownership for the process. We partition the
    global matrix row-wise. Each process owns all columns for a
    contiguous set of rows. */
static inline void
set_HYPRE_row_col_bounds(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                         const long g_n_col,
                         const long *n_row_proc,
                         const int myrank)
{
  PGFEM_hypre->setRowColBounds(g_n_col, n_row_proc, myrank);
}

#endif /* #ifndef HYPRE_GLOBAL_H */

