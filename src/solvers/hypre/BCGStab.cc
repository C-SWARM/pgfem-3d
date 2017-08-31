#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Solvers.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Solver;
using pgfem3d::solvers::hypre::BCGStab;

BCGStab::BCGStab(MPI_Comm comm, int maxit, double err) : Solver(comm)
{
  if (HYPRE_ParCSRBiCGSTABCreate(_comm, &_solver) or
      HYPRE_BiCGSTABSetMaxIter(_solver, maxit) or
      HYPRE_BiCGSTABSetTol(_solver, err) or
      HYPRE_BiCGSTABSetLogging(_solver, 1))
  {
    throw std::runtime_error("failed to initialize BCGStab solver");
  }
}

BCGStab::~BCGStab()
{
  HYPRE_ParCSRBiCGSTABDestroy(_solver);
}

void
BCGStab::getFuncs(ptr_solve_t& solve,
                  ptr_solve_t& setup,
                  ptr_set_pre_t& set_precond,
                  ptr_check_pre_t& get_precond,
                  ptr_iter_t& get_num_iter,
                  ptr_norm_t& get_res_norm) const
{
  set_precond = HYPRE_ParCSRBiCGSTABSetPrecond;
  get_precond = HYPRE_ParCSRBiCGSTABGetPrecond;
  setup = HYPRE_ParCSRBiCGSTABSetup;
  solve = HYPRE_ParCSRBiCGSTABSolve;
  get_num_iter = HYPRE_ParCSRBiCGSTABGetNumIterations;
  get_res_norm = HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm;
}
