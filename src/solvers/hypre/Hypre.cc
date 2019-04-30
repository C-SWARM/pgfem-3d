#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Hypre.hpp"
#include "Preconditioners.hpp"
#include "Solvers.hpp"
#include "enumerations.h"
#include "pgfem3d/Solver.hpp"
#include "pgfem3d/Communication.hpp"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <assert.h>

#ifndef PFEM_HYPRE_DEBUG
#define PFEM_HYPRE_DEBUG 0
#endif

using pgfem3d::PGFEM_Abort;
using pgfem3d::solvers::SparseSystem;
using namespace pgfem3d::solvers::hypre;

Hypre::Hypre(const PGFem3D_opt& opts, MPI_Comm comm, const int Ap[],
             const Ai_t Ai[], const long rowsPerProc[], long maxit, double err)
    : SparseSystem(),
      _comm(comm),
      _Ai(Ai)
{
  /* Rank and nproc */
  int rank, nproc;
  if (MPI_Comm_size(_comm, &nproc) or
      MPI_Comm_rank(_comm, &rank)) {
    PGFEM_Abort();
  }

  if (rank == 0) {
    PGFEM_printf("Iterative solver info:\n"
                 "Kdim     = %d\n"
                 "Max. It. = %d\n", opts.kdim, opts.maxit);
  }

  const auto rows = rowsPerProc[rank];
  _gRows = new HYPRE_t[rows];
  _nCols = new HYPRE_t[rows];

  // Prefix sum to find the index of my first row, and then compute the end of
  // my closed interval range.
  _ilower = std::accumulate(rowsPerProc, rowsPerProc + rank, 0);
  _iupper = _ilower + rows - 1;

  // Generate the rows that we own (dense integer range starting at _ilower)
  std::iota(_gRows, _gRows + rows, _ilower);    // ksaha

  // Record the number of non-zero columns in each of the rows.
  for (long i = 0, e = rows; i < e; ++i) {
    _nCols[i] = Ap[i + 1] - Ap[i];
  }

  // Create the matrix.
  HYPRE_IJMatrixCreate(_comm, _ilower, _iupper, _ilower, _iupper, &_k);
  HYPRE_IJMatrixSetObjectType(_k, HYPRE_PARCSR);
  HYPRE_IJMatrixSetRowSizes(_k, _nCols);

  // HYPRE would like to know how the columns are partitioned in each row, i.e.,
  // how many of the columns are locally owned (diag) and how many columns are
  // not (offd). We compute this by looping through every element we own, and
  // counting the number of diagonal elements (the offd elements are the
  // remainder).
  auto diag = new HYPRE_t[rows]{};
  auto offd = new HYPRE_t[rows];

  // NB: n tracks the linear offset into Ai where we find the column id and is
  //     incremented in the inner loop.
  for (long i = 0, n = 0, e = rows; i < e; ++i) {
    for (HYPRE_t j = 0, e = _nCols[i]; j < e; ++j, ++n) {
      diag[i] += (isLocalRow(Ai[n])) ? 1 : 0;
    }
    offd[i] = _nCols[i] - diag[i];
  }

  HYPRE_IJMatrixSetDiagOffdSizes(_k, diag, offd);
  //HYPRE_IJMatrixSetDiagOffdSizes(_k, const_cast<long long *> (diag), (HYPRE_t *) offd);

  delete [] diag;
  delete [] offd;

  HYPRE_IJMatrixInitialize(_k);

  if (PFEM_HYPRE_DEBUG) {
    /* print pattern */
    set(1.0);
    HYPRE_IJMatrixAssemble(_k);
    HYPRE_IJMatrixPrint(_k, "stiff_pattern.out");
    HYPRE_IJMatrixInitialize(_k);
  }

  HYPRE_IJVectorCreate(_comm, _ilower, _iupper, &_rhs);
  HYPRE_IJVectorSetObjectType(_rhs, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(_rhs);

  HYPRE_IJVectorCreate(_comm, _ilower, _iupper, &_solution);
  HYPRE_IJVectorSetObjectType(_solution, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(_solution);

  createPreconditioner(opts.precond);
  createSolver(opts.solver, opts.maxit, err, opts.kdim);
}

Hypre::~Hypre()
{
  delete _preconditioner;
  delete _solver;

  delete [] _nCols;
  delete [] _gRows;

  /* destroy vectors/matrix */
  HYPRE_IJVectorDestroy(_rhs);
  HYPRE_IJVectorDestroy(_solution);
  HYPRE_IJMatrixDestroy(_k);
}

void
Hypre::assemble()
{
  HYPRE_IJMatrixAssemble(_k);
}

void
Hypre::add(int nrows, sp_idx ncols[], sp_idx const rids[], const sp_idx cids[],
           const double vals[])
{
  HYPRE_IJMatrixAddToValues(_k, nrows, (HYPRE_t *) ncols, (HYPRE_t *) rids, (HYPRE_t *) cids, vals);
}

void
Hypre::resetPreconditioner()
{
  assert(_preconditioner);
  _preconditioner->reset();
}

bool
Hypre::isLocalRow(Ai_t i) const
{
  return (_ilower <= i and i < _iupper + 1);
}

void
Hypre::createSolver(int type, int maxit, double err, int kdim)
{
  switch (type) {
   default:
   case SOLVER_GMRES:
    _solver = new GMRes(_comm, maxit, err, kdim);
    break;
   case SOLVER_BCG_STAB:
    _solver = new BCGStab(_comm, maxit, err);
    break;
   case SOLVER_AMG:
    _solver = new AMG(_comm, err);
    break;
   case SOLVER_FLEX:
    _solver = new Flex(_comm, maxit, err, kdim);
    break;
   case SOLVER_HYBRID:
    _solver = new Hybrid(_comm, maxit, err, kdim);
    break;
  }
}

void
Hypre::createPreconditioner(int type)
{
  switch (type) {
   default:
    _preconditioner = new Preconditioner(_comm);
    break;
   case PRECOND_PARA_SAILS:
    _preconditioner = new ParaSails(_comm);
    break;
   case PRECOND_PILUT:
    _preconditioner = new Pilut(_comm);
    break;
   case PRECOND_EUCLID:
    _preconditioner = new Euclid(_comm);
    break;
   case PRECOND_BOOMER:
    _preconditioner = new Boomer(_comm);
    break;
   case PRECOND_DIAG_SCALE:
    _preconditioner = new ScaleDiag(_comm);
    break;
   case PRECOND_JACOBI:
    _preconditioner = new Jacobi(_comm);
    break;
  }
}

void
Hypre::zero()
{
  set(0);
}

void
Hypre::set(double val)
{
  HYPRE_Int cols[1] = {1};

  for(int i = 0, e = (_iupper - _ilower + 1), n = 0; i < e; ++i) {
    for(HYPRE_t j = 0, e = _nCols[i]; j < e; ++j, ++n) {
      HYPRE_IJMatrixSetValues(_k, 1, cols, &_gRows[i], (HYPRE_t *) &_Ai[n], &val);
    }
  }
}

void
Hypre::printWithRHS(std::string&& basename, const double rhs[], int ndofd,
                    int rank) const
{
  std::stringstream filename;
  filename.str(basename);
  filename << "k.txt";
  HYPRE_IJMatrixPrint(_k, filename.str().c_str());

  filename.str(basename);
  filename << "f.txt." << std::setfill('0') << std::setw(5) << rank;

  std::ofstream os;
  os.open(filename.str());
  os << _ilower << " " << _iupper << "\n";
  for (int i = 0; i < ndofd; ++i) {
    os << _ilower + i << " " << rhs[i] << "\n";
  }
  os.close();
}

double
Hypre::solveSystem(const PGFem3D_opt *opts,
                   double *loc_rhs,
                   double *loc_sol,
                   const int tim,
                   const int iter,
                   const long *DomDof,
                   SOLVER_INFO *info)
{
  assert(opts->solverpackage == HYPRE);
  double func_time = -MPI_Wtime();
  int rank=0;
  MPI_Comm_rank(_comm, &rank);

  /* Assemble the rhs and solution vector */
  nulld(loc_sol, DomDof[rank]);
  HYPRE_IJVectorSetValues(_rhs, DomDof[rank], _gRows, loc_rhs);
  HYPRE_IJVectorSetValues(_solution, DomDof[rank], _gRows, loc_sol);
  HYPRE_IJVectorAssemble(_rhs);
  HYPRE_IJVectorAssemble(_solution);

  /* get preconditioner functions */
  ptr_solve_t precond_solve = NULL;
  ptr_solve_t precond_setup = NULL;
  _preconditioner->getFuncs(precond_solve, precond_setup, iter);

  /* Get solver functions */
  ptr_solve_t    solver_solve = NULL;
  ptr_solve_t    solver_setup = NULL;
  ptr_set_pre_t   set_precond = NULL;
  ptr_check_pre_t get_precond = NULL;
  ptr_iter_t     get_num_iter = NULL;
  ptr_norm_t     get_res_norm = NULL;
  _solver->getFuncs(solver_solve, solver_setup, set_precond, get_precond,
                    get_num_iter, get_res_norm);

  /* Setup the solver */
  if (setupSolverEnv(precond_solve, precond_setup, solver_solve, solver_setup,
                     set_precond, get_precond, info, opts)) {
    return 0.0;
  }

  /* Solve */
  solve(solver_solve, get_num_iter, get_res_norm, info);

  /* Get the solution */
  HYPRE_IJVectorGetValues(_solution, DomDof[rank], _gRows, loc_sol);

  /* reset hypre error flag so we can restart cleanly based on OUR
     error tolerances etc. */
  hypre__global_error = 0;

  /* update timer and return */
  func_time += MPI_Wtime();
  return func_time;
}

double
Hypre::solveSystemNoSetup(const PGFem3D_opt *opts,
                          double *loc_rhs,
                          double *loc_sol,
                          const int tim,
                          const int iter,
                          const long *DomDof,
                          SOLVER_INFO *info)
{
  assert(opts->solverpackage == HYPRE);
  double func_time = -MPI_Wtime();
  int rank=0;
  MPI_Comm_rank(_comm, &rank);

  /* Assemble the rhs and solution vector */
  nulld(loc_sol,DomDof[rank]);
  HYPRE_IJVectorSetValues(_rhs, DomDof[rank], _gRows, loc_rhs);
  HYPRE_IJVectorSetValues(_solution, DomDof[rank], _gRows, loc_sol);
  HYPRE_IJVectorAssemble(_rhs);
  HYPRE_IJVectorAssemble(_solution);

  /* Get solver functions */
  ptr_solve_t solver_solve = NULL;
  ptr_solve_t solver_setup = NULL;

  ptr_set_pre_t set_precond = NULL;
  ptr_check_pre_t get_precond = NULL;
  ptr_iter_t get_num_iter = NULL;
  ptr_norm_t get_res_norm = NULL;
  _solver->getFuncs(solver_solve, solver_setup, set_precond, get_precond,
                    get_num_iter, get_res_norm);

  /* solver should already be set up, just call solve */
  solve(solver_solve, get_num_iter, get_res_norm, info);

  /* Get the solution */
  HYPRE_IJVectorGetValues(_solution, DomDof[rank], _gRows, loc_sol);

  /* update timer and return */
  func_time += MPI_Wtime();
  return func_time;
}

HYPRE_t
Hypre::setupSolverEnv(ptr_solve_t precond_solve,
                      ptr_solve_t precond_setup,
                      ptr_solve_t solver_solve,
                      ptr_solve_t solver_setup,
                      ptr_set_pre_t set_precond,
                      ptr_check_pre_t get_precond,
                      SOLVER_INFO *info,
                      const PGFem3D_opt *opts)
{
  int err = 0;
  /* set preconditioner */
  if (precond_solve and precond_setup and set_precond) {
    auto pc = _preconditioner->_solver;
    set_precond(_solver->_solver, precond_solve, precond_setup, pc);
    if (get_precond) {
      HYPRE_Solver read = nullptr;
      get_precond(_solver->_solver, &read);
      if (read != pc) {
        info->err = BAD_PRECOND;
        err = 1;
      }
    }
  }

  /* call setup on the solver environment */
  HYPRE_ParCSRMatrix     k = nullptr;
  HYPRE_ParVector      rhs = nullptr;
  HYPRE_ParVector solution = nullptr;
  HYPRE_IJMatrixGetObject(_k, reinterpret_cast<void**>(&k));
  HYPRE_IJVectorGetObject(_rhs, reinterpret_cast<void**>(&rhs));
  HYPRE_IJVectorGetObject(_solution, reinterpret_cast<void**>(&solution));
  solver_setup(_solver->_solver, k, rhs, solution);
  return err;
}

HYPRE_t
Hypre::solve(ptr_solve_t solver_solve,
             ptr_iter_t get_num_iter,
             ptr_norm_t get_res_norm,
             SOLVER_INFO *info)
{
  int err = 0;
  HYPRE_ParCSRMatrix     k = nullptr;
  HYPRE_ParVector      rhs = nullptr;
  HYPRE_ParVector solution = nullptr;
  HYPRE_IJMatrixGetObject(_k, reinterpret_cast<void**>(&k));
  HYPRE_IJVectorGetObject(_rhs, reinterpret_cast<void**>(&rhs));
  HYPRE_IJVectorGetObject(_solution, reinterpret_cast<void**>(&solution));
  auto solver = _solver->_solver;
  info->err = solver_solve(solver, k, rhs, solution);
  get_num_iter(solver, (HYPRE_t *) &(info->n_iter));    //ksaha
  get_res_norm(solver, &(info->res_norm));

  /* reset hypre error flag so we can restart cleanly based on OUR
     error tolerances etc. */
  hypre__global_error = 0;

  return err;
}

