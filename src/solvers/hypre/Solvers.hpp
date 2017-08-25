// -------------------------------------------------------------------*- C++ -*-
// -----------------------------------------------------------------------------
#ifndef PGFEM3D_SOLVERS_HYPRE_SOLVERS_H
#define PGFEM3D_SOLVERS_HYPRE_SOLVERS_H

#include "Hypre.hpp"

namespace pgfem3d {
namespace solvers {
namespace hypre {
class GMRes : public Solver {
 public:
  GMRes(MPI_Comm comm, int maxit, double err, int kdim);
  ~GMRes();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup,
                ptr_set_pre_t& set_precond, ptr_check_pre_t& get_precond,
                ptr_iter_t& get_num_iter, ptr_norm_t& get_res_norm) const;

  void setupEnvironment(const Preconditioner& pc);
}; // class GMRes

class BCGStab : public Solver {
 public:
  BCGStab(MPI_Comm comm, int maxit, double err);
  ~BCGStab();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup,
                ptr_set_pre_t& set_precond, ptr_check_pre_t& get_precond,
                ptr_iter_t& get_num_iter, ptr_norm_t& get_res_norm) const;

}; // class BCGStab

class AMG : public Solver {
 public:
  AMG(MPI_Comm comm, double err);
  ~AMG();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup,
                ptr_set_pre_t& set_precond, ptr_check_pre_t& get_precond,
                ptr_iter_t& get_num_iter, ptr_norm_t& get_res_norm) const;

  void setupEnvironment(const Preconditioner& pc);
}; // class AMG

class Flex : public Solver {
 public:
  Flex(MPI_Comm comm, int maxit, double err, int kdim);
  ~Flex();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup,
                ptr_set_pre_t& set_precond, ptr_check_pre_t& get_precond,
                ptr_iter_t& get_num_iter, ptr_norm_t& get_res_norm) const;

  void setupEnvironment(const Preconditioner& pc);
}; // class Flex

class Hybrid : public Solver {
 public:
  Hybrid(MPI_Comm comm, int maxit, double err, int kdim);
  ~Hybrid();

  void getFuncs(ptr_solve_t& solve, ptr_solve_t& setup,
                ptr_set_pre_t& set_precond, ptr_check_pre_t& get_precond,
                ptr_iter_t& get_num_iter, ptr_norm_t& get_res_norm) const;

  void setupEnvironment(const Preconditioner& pc);
}; // class Hybrid
} // namespace hypre
} // namespace solvers
} // namespace pgfem3d

#endif // #define PGFEM3D_SOLVERS_HYPRE_SOLVERS_H
