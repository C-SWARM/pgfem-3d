// -------------------------------------------------------------------*- C++ -*-
// -----------------------------------------------------------------------------
#ifndef PGFEM3D_SOLVER_H
#define PGFEM3D_SOLVER_H

/// @brief This file defines the sparse linear system solver API used in PGFEM.

// Forward declare the types we need.
struct ARC_LENGTH_VARIABLES;
struct PGFEM_HYPRE_solve_info;
namespace pgfem3d { struct CommandLineOptions; }

namespace pgfem3d {
struct Solver {
  virtual ~Solver();

  int n_step = 0; //!< the number of nonlinear steps taken to solve the give increment
  double nor_min = 0.0;  //!< nonlinear convergence tolerance for Newton Raphson
  long iter_max = 0;     //!< maximum number of iterations for Newton Raphson
  double alpha = 0.0;    //!< midpoint rule alpha
  void *microscale = nullptr;           //!< Container of microscale information
  PGFEM_HYPRE_solve_info *PGFEM_hypre = nullptr; //!< HYPRE solver object
  long FNR = 0; //!< "Full Newton-Raphson" == 0, only compute tangent on 1st iteration
  double gama = 0.0;   //!< related to linesearch, but is modified internally...
  double err = 0.0;    //!< linear solve tolerance
  long iter_max_sol = 0;     //!< maximum number of iterations for linear solver
  double computer_zero = 0.0;                   //!< computer zero
  int run_integration_algorithm = 0; //!< if yes, run integration algorithm when compute residuals
  int max_NR_staggering = 0; //!< maximum number of Newton Raphson staggering when physic is coulpled, default = 5
  int max_subdivision = 0; //!< maximum number of subdivision, if -1 (default): take maximum
  int is_subdivided = 0;   //!< if yes, has been subdivided
  double last_residual = 0.0; //!< last residual achieved during Newton-Raphson iterations
  int set_initial_residual = 0; //!< if yes, compute residual before the first NR iteration
  double du = 0.0;      //!< perturbation value for computing the first residual
  ARC_LENGTH_VARIABLES *arc = nullptr; //!< Container of Arc length related variables
}; // struct Solver
} // namespace pgfem3d

#endif // #define PGFEM3D_SOLVER_H
