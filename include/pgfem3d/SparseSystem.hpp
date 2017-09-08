#ifndef PGFEM3D_SPARSE_SYSTEM_HPP
#define PGFEM3D_SPARSE_SYSTEM_HPP

#include "PGFem3D_options.h"
#include <mpi.h>

/** define the structure for general access to solver information */
struct SOLVER_INFO {
  double res_norm;
  int n_iter;
  int err;
};

namespace pgfem3d {
namespace solvers {
class SparseSystem
{
 public:

  static SparseSystem* Create(const PGFem3D_opt& opts, MPI_Comm comm,
                              const int Ap[], const int Ai[],
                              const long rowsPerProc[], long maxit, double err);

  virtual ~SparseSystem();

  /// Assemble the matrix.
  virtual void assemble() = 0;

  /// Add partial sums to values.
  virtual void add(int nrows, int ncols[], int const rids[], const int cids[],
                  const double vals[]) = 0;

  /// Reset the preconditioner.
  virtual void resetPreconditioner() = 0;

  /// Check to see if the row is owned by the local rank.
  ///
  /// @param i          The global row index to check.
  /// @returns          TRUE if the row is local, FALSE otherwise.
  virtual bool isLocalRow(int i) const = 0;

  /// Zero the underlying matrix data.
  virtual void zero() = 0;

  virtual double solveSystem(const PGFem3D_opt *opts,
                             double *loc_rhs,
                             double *loc_sol,
                             const int tim,
                             const int iter,
                             const long *DomDof,
                             SOLVER_INFO *info) = 0;

  virtual double solveSystemNoSetup(const PGFem3D_opt *opts,
                                    double *loc_rhs,
                                    double *loc_sol,
                                    const int tim,
                                    const int iter,
                                    const long *DomDof,
                                    SOLVER_INFO *info) = 0;
};
}
}

#endif // #define PGFEM3D_SPARSE_SYSTEM_HPP
