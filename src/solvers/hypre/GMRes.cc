#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Solvers.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Solver;
using pgfem3d::solvers::hypre::GMRes;

GMRes::GMRes(MPI_Comm comm, int maxit, double err, int kdim) : Solver(comm)
{
  if (HYPRE_ParCSRGMRESCreate(_comm, &_solver) or
      HYPRE_ParCSRGMRESSetMaxIter(_solver, maxit) or
      HYPRE_ParCSRGMRESSetTol(_solver, err) or
      HYPRE_ParCSRGMRESSetKDim(_solver, kdim) or
      HYPRE_ParCSRGMRESSetLogging(_solver, 1) or
      HYPRE_ParCSRGMRESSetPrintLevel(_solver, 0))
  {
    throw std::runtime_error("failed to initialize GMRes solver");
  }
}

GMRes::~GMRes()
{
  HYPRE_ParCSRGMRESDestroy(_solver);
}

void
GMRes::getFuncs(ptr_solve_t& solve,
                ptr_solve_t& setup,
                ptr_set_pre_t& set_precond,
                ptr_check_pre_t& get_precond,
                ptr_iter_t& get_num_iter,
                ptr_norm_t& get_res_norm) const
{
  set_precond = HYPRE_ParCSRGMRESSetPrecond;
  get_precond = HYPRE_ParCSRGMRESGetPrecond;
  setup = HYPRE_ParCSRGMRESSetup;
  solve = HYPRE_ParCSRGMRESSolve;
  get_num_iter = HYPRE_ParCSRGMRESGetNumIterations;
  get_res_norm = HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm;
}
