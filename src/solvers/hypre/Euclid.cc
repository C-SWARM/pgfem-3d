#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"
#include <stdexcept>

using pgfem3d::solvers::hypre::Preconditioner;
using pgfem3d::solvers::hypre::Euclid;

/// Defines the way we configure the Euclid preconditioner.
static const char *config[] = {
 "-level", 
 "1",               //factorization level for ILU(k): 1 (default, we use) varies up to 8
 "-bj", 
 "1",               //parallelization method: 0 (false = PILU, default), 1 (true = BJ-ILU, we use)
 "-rowScale", 
 "0",               //scale values prior to factorization: 0 (false, default), 1 (true, good for large/complex system)
 "-ilut", 
 "0",               //factorization method: 0 (faluse = ILU(k), default, we use), 1 (true = ILUT)
 "-sparseA", 
 "0",               //drop-tolerance for ILU(k) factorizaion: 0 (false, default, we use), 1 (true)
 "-eu_stats", 
 "0"                //print runtime parameters and compute time info: 0 (false, default); 1 (true)
};

Euclid::Euclid(MPI_Comm comm) : Preconditioner(comm)
{
  if (HYPRE_EuclidCreate(_comm, &_solver) or
      HYPRE_EuclidSetParams(_solver, 12, const_cast<char**>(config))) {
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
      HYPRE_EuclidSetParams(_solver, 12, const_cast<char**>(config))) {
    throw std::runtime_error("failed to initialize Euclid preconditioner");
  }
}
