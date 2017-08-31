#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Solvers.hpp"

using pgfem3d::solvers::hypre::Solver;

Solver::Solver(MPI_Comm comm) : _comm(comm) {
}

Solver::~Solver() {
}
