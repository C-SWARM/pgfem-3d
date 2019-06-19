#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/BlockedRowDistribution.hpp"
#include "pgfem3d/Network.hpp"
#include "pgfem3d/SparseSystem.hpp"
#include "hypre/Hypre.hpp"

#ifdef HAVE_TRILINOS
# include "trilinos/Trilinos.hpp"
#endif

using pgfem3d::net::Network;
using pgfem3d::net::PGFem3D_Comm;
using pgfem3d::solvers::SparseSystem;

SparseSystem::~SparseSystem()
{
}

SparseSystem*
SparseSystem::Create(const PGFem3D_opt& opts,
                     Network& net,
                     PGFem3D_Comm comm,
                     const int Ap[],
                     const Index Ai[],
                     const long rowsPerProc[],
                     long maxit,
                     double err)
{
  auto procs = net.size(comm);
  BlockedRowDistribution<long> rows(procs, rowsPerProc);

  auto rank = net.rank(comm);
  auto  min = rows.min(rank);
  auto  max = rows.max(rank);

  if (opts.solverpackage == HYPRE) {
    return new hypre::Hypre(opts, (MPI_Comm)comm, Ap, Ai, min, max, maxit, err);
  }

#ifdef HAVE_TRILINOS
  if (opts.solverpackage == TRILINOS) {
    return new trilinos::TrilinosWrap(opts, (MPI_Comm)comm, Ap, Ai, min, max, maxit, err);
  }
#endif

  std::cerr << "Unsupported solverpackage " << opts.solverpackage << "\n";
  return nullptr;
}
