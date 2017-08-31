#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Preconditioner;
using pgfem3d::solvers::hypre::Euclid;

/// Defines the way we configure the Euclid preconditioner.
static const char *config[] = {
  "-bj",
  "-rowScale",
  "-maxNzPerRow",
  "100",
  "-level",
  "1",
  ""
};

Euclid::Euclid(MPI_Comm comm) : Preconditioner(comm)
{
  if (HYPRE_EuclidCreate(_comm, &_solver) or
      HYPRE_EuclidSetParams(_solver, 7, const_cast<char**>(config))) {
    throw std::runtime_error("failed to initialize Euclid preconditioner");
  }
}

Euclid::~Euclid()
{
  HYPRE_EuclidDestroy(_solver);
}

void
Euclid::getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int) const
{
  solve = HYPRE_EuclidSolve;
  setup = HYPRE_EuclidSetup;
}

void
Euclid::reset()
{
  HYPRE_EuclidDestroy(_solver);
  if (HYPRE_EuclidCreate(_comm, &_solver) or
      HYPRE_EuclidSetParams(_solver, 7, const_cast<char**>(config))) {
    throw std::runtime_error("failed to initialize Euclid preconditioner");
  }
}
