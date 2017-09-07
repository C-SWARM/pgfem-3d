#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "boomerAMGInterface.hpp"
#include "Solvers.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Solver;
using pgfem3d::solvers::hypre::AMG;
using pgfem3d::solvers::hypre::boomerAMGOptions;

AMG::AMG(MPI_Comm comm, double err) : Solver(comm)
{
  boomerAMGOptions AMGOptions;
  pgfem3d::solvers::hypre::setBoomerAMGOptions(&AMGOptions);
  if (pgfem3d::solvers::hypre::initializeBoomerAMG(&_solver, &AMGOptions, _comm) or
      HYPRE_BoomerAMGSetTol(_solver, err) or
      HYPRE_BoomerAMGSetMaxIter(_solver, 20))
  {
    throw std::runtime_error("failed to initialize AMG solver");
  }
}

AMG::~AMG()
{
  HYPRE_BoomerAMGDestroy(_solver);
}

void
AMG::getFuncs(ptr_solve_t& solve,
              ptr_solve_t& setup,
              ptr_set_pre_t& set_precond,
              ptr_check_pre_t& get_precond,
              ptr_iter_t& get_num_iter,
              ptr_norm_t& get_res_norm) const
{
  /* no preconditioner */
  set_precond = nullptr;
  get_precond = nullptr;
  setup = HYPRE_BoomerAMGSetup;
  solve = HYPRE_BoomerAMGSolve;
  // *get_num_iter = HYPRE_BoomerAMGGetNumIterations;
  // *get_res_norm = HYPRE_BoomerAMGGetFinalRelativeResidualNorm;
  get_num_iter = HYPRE_ParCSRFlexGMRESGetNumIterations;
  get_res_norm = HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm;
}
