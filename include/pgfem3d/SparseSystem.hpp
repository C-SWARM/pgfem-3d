#ifndef PGFEM3D_SPARSE_SYSTEM_HPP
#define PGFEM3D_SPARSE_SYSTEM_HPP

#include "msnet/Network.hpp"
#include "PGFem3D_options.h"
#include "datatype.hpp"

/** define the structure for general access to solver information */
struct SOLVER_INFO {
  double res_norm;
  Ai_t n_iter;
  int err;
};

namespace pgfem3d {
namespace solvers {
class SparseSystem
{
 public:
  using Index = Ai_t;

  static SparseSystem* Create(const PGFem3D_opt& opts,
                              multiscale::net::Network& net,
                              multiscale::net::MSNET_Comm comm,
                              const int Ap[],
                              const Index Ai[],
                              const long rowsPerProc[],
                              int maxit,
                              double err);

  virtual ~SparseSystem();

  /// Assemble the matrix.
  virtual void assemble() = 0;

  /// Add partial sums to values.
  virtual void add(int nrows, Index ncols[], Index const rids[], const Index cids[],
                   const double vals[]) = 0;

  /// Reset the preconditioner.
  virtual void resetPreconditioner() = 0;

  /// Check to see if the row is owned by the local rank.
  ///
  /// @param i          The global row index to check.
  /// @returns          TRUE if the row is local, FALSE otherwise.
  bool isLocal(Index i) const {
    return (iMin_ <= i and i < iMax_);
  }

  auto nRows() const {
    return iMax_ - iMin_;
  }

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

 protected:
  SparseSystem(Index iMin, Index iMax) : iMin_(iMin), iMax_(iMax) {
  }

  Index iMin_ = 0;
  Index iMax_ = 0;
};
}
}

#endif // #define PGFEM3D_SPARSE_SYSTEM_HPP
