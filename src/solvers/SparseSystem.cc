#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Network.hpp"
#include "pgfem3d/SparseSystem.hpp"
#include "hypre/Hypre.hpp"
#include "trilinos/Trilinos.hpp"

using pgfem3d::net::PGFem3D_Comm;
using pgfem3d::solvers::SparseSystem;

SparseSystem::~SparseSystem()
{
}

SparseSystem*
SparseSystem::Create(const PGFem3D_opt& opts, const PGFem3D_Comm comm, const int Ap[],
                     const int Ai[], const long rowsPerProc[], long maxit,
                     double err)
{
  SparseSystem* system=NULL;

  if(opts.solverpackage == HYPRE)
    system = new hypre::Hypre { opts, (MPI_Comm)comm, Ap, Ai, rowsPerProc, maxit, err };
  
  else if(opts.solverpackage == TRILINOS)
    system = new trilinos::TrilinosWrap { opts, (MPI_Comm)comm, Ap, Ai, rowsPerProc, maxit, err };

  else
    std::cerr << "Unsupported solverpackage " << opts.solverpackage << std::endl;
    
  
  return system;
}
