#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"
#include "boomerAMGInterface.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Preconditioner;
using pgfem3d::solvers::hypre::Boomer;
using pgfem3d::solvers::hypre::boomerAMGOptions;

Boomer::Boomer(MPI_Comm comm) : Preconditioner(comm)
{
  boomerAMGOptions AMGOptions;
  setBoomerAMGOptions(&AMGOptions);
  if (initializeBoomerAMG(&_solver, &AMGOptions, _comm)) {
    throw std::runtime_error("failed to initialize Boomer preconditioner");
  }
}

Boomer::~Boomer()
{
  HYPRE_BoomerAMGDestroy(_solver);
}

void
Boomer::getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int) const
{
  solve = HYPRE_BoomerAMGSolve;
  setup = HYPRE_BoomerAMGSetup;
}

void
Boomer::reset()
{
  HYPRE_BoomerAMGDestroy(_solver);
  boomerAMGOptions AMGOptions;
  setBoomerAMGOptions(&AMGOptions);
  if (initializeBoomerAMG(&_solver, &AMGOptions, _comm)) {
    throw std::runtime_error("failed to initialize Boomer preconditioner");
  }
}
