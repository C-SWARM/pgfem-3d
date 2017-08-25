// -------------------------------------------------------------------*- C++ -*-
// -----------------------------------------------------------------------------
#ifndef PGFEM3D_SOLVERS_HYPRE_PRECONDITIONERS_H
#define PGFEM3D_SOLVERS_HYPRE_PRECONDITIONERS_H

#include "Hypre.hpp"

namespace pgfem3d {
namespace solvers {
namespace hypre {
class Euclid : public Preconditioner {
 public:
  Euclid(MPI_Comm comm);
  ~Euclid();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const;
  void reset();
}; // class Euclid

class ParaSails : public Preconditioner {
 public:
  ParaSails(MPI_Comm comm);
  ~ParaSails();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const;
  void reset();
}; // class ParaSails

class Pilut : public Preconditioner {
 public:
  Pilut(MPI_Comm comm);
  ~Pilut();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const;
  void reset();
}; // class Pilut

class Boomer : public Preconditioner {
 public:
  Boomer(MPI_Comm comm);
  ~Boomer();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const;
  void reset();
}; // class Boomer

class Jacobi : public Preconditioner {
 public:
  Jacobi(MPI_Comm comm);
  ~Jacobi();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const;
  void reset();
}; // class Jacobi

class ScaleDiag : public Preconditioner {
 public:
  ScaleDiag(MPI_Comm comm);
  ~ScaleDiag();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int iter) const;
  void reset();
}; // class ScaleDiag
} // namespace hypre
} // namespace solvers
} // namespace pgfem3d

#endif // #define PGFEM3D_SOLVERS_HYPRE_PRECONDITIONERS_H
