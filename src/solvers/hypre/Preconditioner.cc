#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"

using pgfem3d::solvers::hypre::Preconditioner;

Preconditioner::Preconditioner(MPI_Comm comm) : _comm(comm) {
}

Preconditioner::~Preconditioner() {
}

void
Preconditioner::getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int) const
{
  solve = nullptr;
  setup = nullptr;
}

void
Preconditioner::reset()
{
}
