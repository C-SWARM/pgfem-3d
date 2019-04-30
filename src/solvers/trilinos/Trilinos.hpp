// -------------------------------------------------------------------*- C++ -*-
// -----------------------------------------------------------------------------
#ifndef PGFEM3D_SOLVERS_TRILINOS_TRILINOS_HPP
#define PGFEM3D_SOLVERS_TRILINOS_TRILINOS_HPP


#include "incl.h"
#include "PGFem3D_options.h"
#include "pgfem3d/SparseSystem.hpp"
#include <string>

#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosTpetraAdapter.hpp>

#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>

#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>



struct SOLVER_INFO;

/**
 * @name Structure to contain TRILINOS objects and pass them around easily
 */
namespace pgfem3d {
namespace solvers {
namespace trilinos {

using ST =                                                 double;
using LO =                                                    int;
using GO =                                                   Ai_t;
using Node =          KokkosClassic::DefaultNode::DefaultNodeType;
using TCrsMatrix =               Tpetra::CrsMatrix<ST,LO,GO,Node>;
using TRowMatrix =               Tpetra::RowMatrix<ST,LO,GO,Node>;
using Map =                               Tpetra::Map<LO,GO,Node>;
using MV =                     Tpetra::MultiVector<ST,LO,GO,Node>;
using OP =                        Tpetra::Operator<ST,LO,GO,Node>;
using prec_type =          Ifpack2::Preconditioner<ST,LO,GO,Node>;
using OT =          Teuchos::OrdinalTraits<Tpetra::global_size_t>;


struct TrilinosWrap;
struct TrilinosWrap : public SparseSystem
{
 public:
  /// Construct a system of equations, solver, and preconditioner.
  ///
  /// @param options     Solver options read from the input data.
  /// @param comm        The communicator on which this object lives.
  /// @param Ap
  /// @param Ai
  /// @param rowsPerProc The row partitioning for this solver.
  /// @param maxit       The maximum number of iterations.
  /// @param err         The error tolerance used during solving.
  TrilinosWrap(const PGFem3D_opt& options, MPI_Comm comm, const int Ap[],
        const GO Ai[], const long rowsPerProc[], long maxit, double err);

  ~TrilinosWrap();

  /// Assemble the matrix.
  void assemble();

  /// Add partial sums to values.
  void add(LO nrows, GO ncols[], GO const rids[], const GO cids[],
           const ST vals[]);

  /// Reset the prconditioner.
  void resetPreconditioner();

  /// Check to see if the row is owned by the local rank.
  ///
  /// @param i          The global row index to check.
  /// @returns          TRUE if the row is local, FALSE otherwise.
  bool isLocalRow(GO i) const;
  
  /// Zero the underlying matrix data.
  void zero();

  /// Print system matrix and vector (LHS and RHS) for debugging.
  ///
  /// This function will generate files as many as number of processes.
  ///
  /// @param basename   The base filename to use as output.
  /// @param rhs        Array of data for the right hand side.
  /// @param ndofd      # of degree freedom of domain
  /// @param myrank     Current process rank.
  void printWithRHS(std::string&& basename, const ST rhs[], int ndofd,
                    int rank) const;

  double solveSystem(const PGFem3D_opt *opts,
                     ST *loc_rhs,
                     ST *loc_sol,
                     const int tim,
                     const int iter,
                     const long *DomDof,
                     SOLVER_INFO *info);

  double solveSystemNoSetup(const PGFem3D_opt *opts,
                            ST *loc_rhs,
                            ST *loc_sol,
                            const int tim,
                            const int iter,
                            const long *DomDof,
                            SOLVER_INFO *info);

 private:
  /// Encapsulate creation of the solver.
  void createSolver(int type, int maxit, ST err, int kdim);

  /// Encapsulate the creation of the preconditioner.
  void createPreconditioner();

  /// A utility function that will set all of the values in the matrix to val.
  void set(ST val);

  /** setup the TRILINOS solver environment */
  int setupSolverEnv();

  /** solve using the TRILINOS environment */
  int solve(SOLVER_INFO *info);

  /** set preconditioner options */
  void setTeuchosParameters(const char *prec_xmlpath, const char *sol_xmlpath);
  
  // Teuchos::RCP is the reference counting smart pointer used throughout
  // Trilinos.
  
  // The linear system data (Ax=b)
  Teuchos::RCP<TCrsMatrix> _localk;               //!< The local matrix handle
  Teuchos::RCP<TCrsMatrix>      _k;               //!< The matrix handle
  Teuchos::RCP<MV>            _rhs;               //!< The RHS vector handle
  Teuchos::RCP<MV>       _solution;               //!< The solution vector handle

  // Belos linear problem and solver manager
  Teuchos::RCP<Belos::LinearProblem<ST,MV,OP> > _problem;
  Teuchos::RCP<Belos::SolverManager<ST,MV,OP> >  _solver;
  
  // Trilinos parameter lists
  Teuchos::RCP<Teuchos::ParameterList> _belosList;
  Teuchos::RCP<Teuchos::ParameterList> _precList;
  //Teuchos::RCP<Teuchos::ParameterList>  _mueluList;

  // Ifpack2 Preconditioner
  Teuchos::RCP<prec_type> _ifpackprec;

  // MueLu Preconditioner
  Teuchos::RCP<MueLu::TpetraOperator<ST,LO,GO> > _mueluprec;
  
  // Preconditioner operator for either Ifpack2 or MueLu
  Teuchos::RCP<OP> _preconditioner;

  // Sparse matrix data relevant to the distribution of data.
  const MPI_Comm _comm;                         //!<
  const GO * const _Ai;                        //!< reference to column array

  GO *_gRows = nullptr;                        //!< Local row indices
  size_t *_nCols = nullptr;                    //!< Number of columns per row

  GO _ilower = 0;                              //!< Row lower bound
  GO _iupper = 0;                              //!< Row upper bound

  int _prectype;
  std::string _ifpacktype;
}; // struct Trilinos
} // namespace trilinos
} // namespace solvers
} // namespace pgfem3d

#endif // #define PGFEM3D_SOLVERS_TRILINOS_TRILINOS_HPP

