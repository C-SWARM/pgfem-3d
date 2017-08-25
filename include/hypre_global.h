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
struct PGFEM_HYPRE_solve_info {
  HYPRE_IJMatrix hypre_k;        /**< */
  HYPRE_ParCSRMatrix hypre_pk;   /**< */
  HYPRE_Solver hypre_solver;     /**< */
  HYPRE_Solver hypre_pc;         /**< */
  HYPRE_Solver hypre_pc_gotten;  /**< */
  HYPRE_IJVector hypre_rhs;      /**< */
  HYPRE_ParVector hypre_prhs;    /**< */
  HYPRE_IJVector hypre_sol;      /**< */
  HYPRE_ParVector hypre_psol;    /**< */

  int *ncol;                     /**< */
  int *grows;                    /**< */
  int ilower;                    /**< */
  int iupper;                    /**< */
  int jlower;                    /**< */
  int jupper;                    /**< */

  int precond_type;              /**< */
  int solver_type;               /**< */
};

void initialize_PGFEM_HYPRE_solve_info(PGFEM_HYPRE_solve_info **info);
void destroy_PGFEM_HYPRE_solve_info(PGFEM_HYPRE_solve_info *info);

/* HYPRE utils */
void hypre_initialize(int *Ap,
                      int *Ai,
                      int size,
                      int maxit,
                      double err,
                      PGFEM_HYPRE_solve_info *PGFEM_hypre,
                      const PGFem3D_opt *options,
                      MPI_Comm mpi_comm);

void ZeroHypreK(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                int *Ai,
                int size);

/** Sets the extends of ownership for the process. We partition the
    global matrix row-wise. Each process owns all columns for a
    contiguous set of rows. */
void set_HYPRE_row_col_bounds(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                              const long g_n_col,
                              const long *n_row_proc,
                              const int myrank);

int PGFEM_HYPRE_create_preconditioner(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                                      const MPI_Comm mpi_comm);

int destroy_HYPRE_preconditioner(PGFEM_HYPRE_solve_info *PGFEM_hypre);

#endif /* #ifndef HYPRE_GLOBAL_H */

