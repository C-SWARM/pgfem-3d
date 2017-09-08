#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Preconditioner;
using pgfem3d::solvers::hypre::Pilut;

/// Constants used to initialize the Pilut preconditioner.
/// @{
static constexpr double TOL = 0.1;
static constexpr int    ROW = 20;
static constexpr int  MAXIT = 50;
/// @}

Pilut::Pilut(MPI_Comm comm) : Preconditioner(comm)
{
  if (HYPRE_ParCSRPilutCreate(_comm, &_solver) or
      HYPRE_ParCSRPilutSetMaxIter(_solver, MAXIT) or
      HYPRE_ParCSRPilutSetDropTolerance(_solver, TOL) or
      HYPRE_ParCSRPilutSetFactorRowSize(_solver, ROW)) {
    throw std::runtime_error("failed to initialize Pilut preconditioner");
  }
}

Pilut::~Pilut()
{
  HYPRE_ParCSRPilutDestroy(_solver);
}

void
Pilut::getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int) const
{
  solve = HYPRE_ParCSRPilutSolve;
  setup = HYPRE_ParCSRPilutSetup;
}

void
Pilut::reset()
{
  HYPRE_ParCSRPilutDestroy(_solver);
  if (HYPRE_ParCSRPilutCreate(_comm, &_solver) or
      HYPRE_ParCSRPilutSetMaxIter(_solver, MAXIT) or
      HYPRE_ParCSRPilutSetDropTolerance(_solver, TOL) or
      HYPRE_ParCSRPilutSetFactorRowSize(_solver, ROW)) {
    throw std::runtime_error("failed to initialize Pilut preconditioner");
  }
}
