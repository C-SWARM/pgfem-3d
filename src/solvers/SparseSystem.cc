#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Network.hpp"
#include "pgfem3d/SparseSystem.hpp"
#include "hypre/Hypre.hpp"

#ifdef HAVE_TRILINOS
# include "trilinos/Trilinos.hpp"
#endif

using pgfem3d::net::PGFem3D_Comm;
using pgfem3d::solvers::SparseSystem;

SparseSystem::~SparseSystem()
{
}

SparseSystem*
SparseSystem::Create(const PGFem3D_opt& opts, const PGFem3D_Comm comm, const int Ap[],
                     const Ai_t Ai[], const long rowsPerProc[], long maxit,
                     double err)
{
  SparseSystem* system = NULL;

  if(opts.solverpackage == HYPRE)
    system = new hypre::Hypre { opts, (MPI_Comm)comm, Ap, Ai, rowsPerProc, maxit, err };

  #ifdef HAVE_TRILINOS
  else if(opts.solverpackage == TRILINOS)
    system = new trilinos::TrilinosWrap { opts, (MPI_Comm)comm, Ap, Ai, rowsPerProc, maxit, err };
  #endif
  
  else
    std::cerr << "Unsupported solverpackage " << opts.solverpackage << std::endl;
    
  
  return system;
}
