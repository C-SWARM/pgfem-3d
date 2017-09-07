#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Preconditioner;
using pgfem3d::solvers::hypre::ParaSails;

/// Constants used to initialize the Parasails preconditioner.
/// @{
static constexpr int        SYM = 0;
static constexpr int     NLEVEL = 1;
static constexpr double  THRESH = 0.1;
static constexpr double  FILTER = 0.05;
static constexpr double LOADBAL = 0.9;
/// @}

ParaSails::ParaSails(MPI_Comm comm) : Preconditioner(comm)
{
  if (HYPRE_ParaSailsCreate(_comm, &_solver) or
      HYPRE_ParaSailsSetSym(_solver, SYM) or
      HYPRE_ParaSailsSetParams(_solver, THRESH, NLEVEL) or
      HYPRE_ParaSailsSetFilter(_solver, FILTER) or
      HYPRE_ParaSailsSetLoadbal(_solver, LOADBAL))
  {
    throw std::runtime_error("failed to initialize Parasails preconditioner");
  }
}

ParaSails::~ParaSails()
{
  HYPRE_ParaSailsDestroy(_solver);
}

void
ParaSails::getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const
{
  HYPRE_ParaSailsSetReuse(_solver, (iter <= 1) ? 0 : 1);
  solve = HYPRE_ParCSRParaSailsSolve;
  setup = HYPRE_ParCSRParaSailsSetup;
}

void
ParaSails::reset()
{
  HYPRE_ParaSailsDestroy(_solver);
  if (HYPRE_ParaSailsCreate(_comm, &_solver) or
      HYPRE_ParaSailsSetSym(_solver, SYM) or
      HYPRE_ParaSailsSetParams(_solver, THRESH, NLEVEL) or
      HYPRE_ParaSailsSetFilter(_solver, FILTER) or
      HYPRE_ParaSailsSetLoadbal(_solver, LOADBAL))
  {
    throw std::runtime_error("failed to initialize Parasails preconditioner");
  }
}
