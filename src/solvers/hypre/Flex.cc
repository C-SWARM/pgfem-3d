#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Solvers.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Solver;
using pgfem3d::solvers::hypre::Flex;

Flex::Flex(MPI_Comm comm, int maxit, double err, int kdim) : Solver(comm)
{
  if (HYPRE_ParCSRFlexGMRESCreate(_comm, &_solver) or
      HYPRE_ParCSRFlexGMRESSetMaxIter(_solver, maxit) or
      HYPRE_ParCSRFlexGMRESSetTol(_solver, err) or
      HYPRE_ParCSRFlexGMRESSetKDim(_solver, kdim) or
      HYPRE_ParCSRFlexGMRESSetLogging(_solver, 1) or
      HYPRE_ParCSRFlexGMRESSetPrintLevel(_solver, 0))
  {
    throw std::runtime_error("failed to initialize Flex solver");
  }
}

Flex::~Flex()
{
  // @todo[ld] shouldn't we do something here?
}

void
Flex::getFuncs(ptr_solve_t& solve,
               ptr_solve_t& setup,
               ptr_set_pre_t& set_precond,
               ptr_check_pre_t& get_precond,
               ptr_iter_t& get_num_iter,
               ptr_norm_t& get_res_norm) const
{
  set_precond = HYPRE_ParCSRFlexGMRESSetPrecond;
  get_precond = HYPRE_ParCSRFlexGMRESGetPrecond;
  setup = HYPRE_ParCSRFlexGMRESSetup;
  solve = HYPRE_ParCSRFlexGMRESSolve;
  get_num_iter = HYPRE_ParCSRFlexGMRESGetNumIterations;
  get_res_norm = HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm;
}
