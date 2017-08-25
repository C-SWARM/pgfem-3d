#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Solvers.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Solver;
using pgfem3d::solvers::hypre::Hybrid;

Hybrid::Hybrid(MPI_Comm comm, int maxit, double err, int kdim) : Solver(comm)
{
  if (HYPRE_ParCSRHybridCreate(&_solver) or
      HYPRE_ParCSRHybridSetKDim(_solver, kdim) or
      HYPRE_ParCSRHybridSetTol(_solver, err) or
      HYPRE_ParCSRHybridSetDSCGMaxIter(_solver, maxit) or
      HYPRE_ParCSRHybridSetSolverType(_solver, 2) or
      HYPRE_ParCSRHybridSetConvergenceTol(_solver, 1.0))
  {
    throw std::runtime_error("failed to initialize Hybrid solver");
  }
}

Hybrid::~Hybrid()
{
  // @todo[ld] shouldn't we do something here?
}


void
Hybrid::getFuncs(ptr_solve_t& solve,
               ptr_solve_t& setup,
               ptr_set_pre_t& set_precond,
               ptr_check_pre_t& get_precond,
               ptr_iter_t& get_num_iter,
               ptr_norm_t& get_res_norm) const
{
  set_precond = HYPRE_ParCSRHybridSetPrecond;
  get_precond = nullptr;                        // @todo[ld] is this correct?
  setup = HYPRE_ParCSRHybridSetup;
  solve = HYPRE_ParCSRHybridSolve;
  get_num_iter = HYPRE_ParCSRHybridGetNumIterations;
  get_res_norm = HYPRE_ParCSRHybridGetFinalRelativeResidualNorm;
}
