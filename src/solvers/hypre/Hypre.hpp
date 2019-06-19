// -------------------------------------------------------------------*- C++ -*-
// -----------------------------------------------------------------------------
#ifndef PGFEM3D_SOLVERS_HYPRE_HYPRE_HPP
#define PGFEM3D_SOLVERS_HYPRE_HYPRE_HPP

#include "incl.h"
#include "PGFem3D_options.h"
#include "pgfem3d/SparseSystem.hpp"
#include <HYPRE.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <krylov.h>
#include <string>

#include <mpi.h> // Hypre bound to MPI for now?

struct SOLVER_INFO;

/**
 * @name Structure to contain HYPRE objects and pass them around easily
 */
namespace pgfem3d {
namespace solvers {
namespace hypre {
struct Hypre;

/* create aliases for HYPRE functions so we can set and run with fewer
   switches/code */

typedef HYPRE_t (*ptr_solve_t)(HYPRE_Solver solver,
                               HYPRE_ParCSRMatrix A,
                               HYPRE_ParVector b,
                               HYPRE_ParVector x);

typedef HYPRE_t (*ptr_set_pre_t)(HYPRE_Solver solver,
                                 ptr_solve_t pre_solve,
                                 ptr_solve_t pre_setup,
                                 HYPRE_Solver precond);

typedef HYPRE_t (*ptr_check_pre_t)(HYPRE_Solver solver,
                                   HYPRE_Solver *precond);

typedef HYPRE_t (*ptr_iter_t)(HYPRE_Solver solver,
                              HYPRE_t *num_it);

typedef HYPRE_t (*ptr_norm_t)(HYPRE_Solver solver,
                              double *norm);

class Preconditioner
{
  friend struct Hypre;                          // temporary for conversion

 public:
  Preconditioner(MPI_Comm comm);
  virtual ~Preconditioner();

  virtual void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const;
  virtual void reset();

 protected:
  const MPI_Comm _comm;
  HYPRE_Solver _solver = nullptr;
}; // class Preconditioner

class Solver
{
  friend struct Hypre;                          // temporary for conversion

 public:
  Solver(MPI_Comm comm);
  virtual ~Solver();

  virtual void getFuncs(ptr_solve_t& solve,
                        ptr_solve_t& setup,
                        ptr_set_pre_t& set_precond,
                        ptr_check_pre_t& get_precond,
                        ptr_iter_t& get_num_iter,
                        ptr_norm_t& get_res_norm) const = 0;

 protected:
  const MPI_Comm _comm;
  HYPRE_Solver   _solver = nullptr;
}; // class Solver

struct Hypre : public SparseSystem
{
 public:
  /// Construct a system of equations, solver, and preconditioner.
  ///
  /// @param options     Solver options read from the input data.
  /// @param comm        The communicator on which this object lives.
  /// @param Ap
  /// @param Ai
  /// @param iMin        The minimum row index for this rank.
  /// @param iMax        The maximum row index for this rank (half-open)
  /// @param maxit       The maximum number of iterations.
  /// @param err         The error tolerance used during solving.
  Hypre(const PGFem3D_opt& options, MPI_Comm comm, const int Ap[],
        const Index Ai[], Index iMin, Index iMax, int maxit, double err);

  ~Hypre();

  /// Assemble the matrix.
  void assemble();

  /// Add partial sums to values.
  void add(int nrows, Index ncols[], Index const rids[], const Index cids[],
           const double vals[]);

  /// Reset the prconditioner.
  void resetPreconditioner();

  /// Zero the underlying matrix data.
  void zero();

  /// Print system matrix and vector (LHS and RHS) for debugging.
  ///
  /// This function will generate files as many as number of processes.
  ///
  /// @param basename   The base filename to use as output.
  /// @param rhs        Array of data for the right hand side.
  /// @param ndofd      # of degree freedom of domain
  /// @param myrank     Current process rank.
  void printWithRHS(std::string&& basename, const double rhs[], int ndofd,
                    int rank) const;

  double solveSystem(const PGFem3D_opt *opts,
                     double *loc_rhs,
                     double *loc_sol,
                     const int tim,
                     const int iter,
                     const long *DomDof,
                     SOLVER_INFO *info);

  double solveSystemNoSetup(const PGFem3D_opt *opts,
                            double *loc_rhs,
                            double *loc_sol,
                            const int tim,
                            const int iter,
                            const long *DomDof,
                            SOLVER_INFO *info);

 private:
  /// Encapsulate creation of the solver.
  void createSolver(int type, int maxit, double err, int kdim);

  /// Encapsulate the creation of the preconditioner.
  void createPreconditioner(int type);

  /// A utility function that will set all of the values in the matrix to @p
  /// val.
  ///
  /// @param Ai
  /// @param size
  /// @param val
  void set(double val);

  /** setup the HYPRE solver environment */
  HYPRE_t setupSolverEnv(ptr_solve_t precond_solve,
                     ptr_solve_t precond_setup,
                     ptr_solve_t solver_solve,
                     ptr_solve_t solver_setup,
                     ptr_set_pre_t set_precond,
                     ptr_check_pre_t get_precond,
                     SOLVER_INFO *info,
                     const PGFem3D_opt *opts);

  /** solve using the HYPRE environment */
  HYPRE_t solve(ptr_solve_t solver_solve,
            ptr_iter_t get_num_iter,
            ptr_norm_t get_res_norm,
            SOLVER_INFO *info);

  /// Hypre doesn't want half-open intervals so we adapt the superclass' notion
  // of min and max.
  auto ilower() const {
    return iMin_;
  }

  auto iupper() const {
    return iMax_ - 1;
  }

  // The linear system data (Ax=b)
  HYPRE_IJMatrix        _k = nullptr;           //!< The matrix handle
  HYPRE_IJVector      _rhs = nullptr;           //!< The RHS vector handle
  HYPRE_IJVector _solution = nullptr;           //!< The solution vector handle

  // Sparse matrix data relevant to the distribution of data.
  const MPI_Comm _comm;                         //!<
  const Index * const _Ai;                       //!< reference to column array

  HYPRE_t *_gRows = nullptr;                    //!< Local row indices
  HYPRE_t *_nCols = nullptr;                    //!< Number of columns per row

  Preconditioner* _preconditioner = nullptr;
  Solver*                 _solver = nullptr;
}; // struct Hypre
} // namespace hypre
} // namespace solvers
} // namespace pgfem3d

#endif // #define PGFEM3D_SOLVERS_HYPRE_HYPRE_HPP

