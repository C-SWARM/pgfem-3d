#ifdef HAVE_CONFIG_H
# include "config.h"
#endif


#include "Trilinos.hpp"
#include "enumerations.h"
#include "pgfem3d/Solver.hpp"
#include "pgfem3d/Communication.hpp"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <assert.h>


#ifndef PGFEM_TRILINOS_DEBUG
#define PGFEM_TRILINOS_DEBUG 0
#endif

using pgfem3d::PGFEM_Abort;
using pgfem3d::solvers::SparseSystem;
using namespace pgfem3d::solvers::trilinos;


TrilinosWrap::TrilinosWrap(const PGFem3D_opt& opts, MPI_Comm comm, const int Ap[],
             const int Ai[], const long rowsPerProc[], long maxit, double err)
    : SparseSystem(),
      _comm(comm),
      _Ai(Ai),
      _prectype(opts.precond)
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

  Teuchos::RCP<Node> node = Teuchos::null;
  
  const auto localrows = rowsPerProc[rank];

  _gRows = new int[localrows];
  _nCols = new size_t[localrows];
  
  // Prefix sum to find the index of my first row, and then compute the end of
  // my closed interval range.
  _ilower = std::accumulate(rowsPerProc, rowsPerProc + rank, 0);
  _iupper = _ilower + localrows - 1;

  // Generate the rows that we own (dense integer range starting at _ilower)
  std::iota(_gRows, _gRows + localrows, _ilower);

  int ourcols = 0;
  // Record the number of non-zero columns in each of the rows.
  for (int i = 0, e = localrows; i < e; ++i) {
    _nCols[i] = Ap[i + 1] - Ap[i];
    ourcols += _nCols[i];
  }

  Teuchos::RCP<const Teuchos::Comm<int> > mycomm =
    Teuchos::rcp(new Teuchos::MpiComm<int>(_comm));

  Teuchos::RCP<const Map> map =
    Teuchos::RCP<const Map>(new Map(OT::invalid(), localrows, 0, mycomm, node));

  const Teuchos::ArrayView<const int> colids(&Ai[Ap[0]], ourcols);
  Teuchos::RCP<const Map> colmap =
    Teuchos::RCP<const Map>(new Map(OT::invalid(), colids, 0, mycomm, node));
  

  const Teuchos::ArrayRCP< size_t > nCols(_nCols, 0, localrows, false); 
  _localk =
    Teuchos::rcp(new TCrsMatrix(map, colmap, nCols, Tpetra::DynamicProfile));
  
  int numrhs = 1;
  _rhs = Teuchos::rcp(new MV(map, numrhs));
  _solution = Teuchos::rcp(new MV(map, numrhs));

  createPreconditioner(opts.trilxml);
  createSolver(opts.solver, opts.maxit, err, opts.kdim);
}

TrilinosWrap::~TrilinosWrap()
{
  delete [] _nCols;
  delete [] _gRows;
}

void
TrilinosWrap::assemble()
{
  _localk->fillComplete();

  Tpetra::Export<int,int,Node> rowexporter(_localk->getRowMap(),
					   _localk->getRowMap());
  
  _k = Tpetra::exportAndFillCompleteCrsMatrix(_localk.getConst(),
					      rowexporter,
					      _localk->getRowMap(),
					      _localk->getRowMap());
}

void
TrilinosWrap::add(int nrows, int ncols[], int const rids[], const int cids[],
           const double vals[])
{
  int idx=0;
  for (int i = 0; i < nrows; ++i){
    const Teuchos::ArrayView<const int> _cids(cids + idx, ncols[i]);
    const Teuchos::ArrayView<const ST> _vals(vals + idx, ncols[i]);
    
    _localk->sumIntoGlobalValues(rids[i], _cids, _vals);
    idx += ncols[i];
  }
}

void
TrilinosWrap::resetPreconditioner()
{
}

bool
TrilinosWrap::isLocalRow(int i) const
{
  return (_ilower <= i and i < _iupper + 1);
}

void
TrilinosWrap::createSolver(int type, int maxit, double err, int kdim)
{
  // The belos options should probably have their own xml read
  _belosList=Teuchos::rcp(new Teuchos::ParameterList());
  _belosList->set("Maximum Iterations", maxit);
  _belosList->set( "Convergence Tolerance", err );
  _belosList->set("Num Blocks", kdim);
  _belosList->set("Block Size", 1);

  int verbLevel = 0;
  // verbLevel+=Belos::Errors + Belos::Warnings + Belos::TimingDetails
  // + Belos::FinalSummary + Belos::StatusTestDetails + Belos::Debug;

  _belosList->set( "Verbosity", verbLevel );
  _belosList->set( "Output Frequency", -1 );
  _belosList->set( "Maximum Restarts", 20);
  _problem = Teuchos::rcp(new Belos::LinearProblem<ST,MV,OP>());
}

void
TrilinosWrap::createPreconditioner(const char *xmlpath)
{
  switch(_prectype) {
  case PRECOND_MUELU: {
    _mueluList = Teuchos::getParametersFromXmlFile(std::string(xmlpath));
  }
    break;
    
  case PRECOND_IFPACK2: {
    _ifpackList = Teuchos::rcp(new Teuchos::ParameterList());
    _ifpackList = Teuchos::getParametersFromXmlFile(std::string(xmlpath));
    _ifpacktype = _ifpackList->get<std::string>("Preconditioner");

    // Ifpack::Relaxation will complain if this is in the parameter
    // list so remove it.
    _ifpackList->remove("Preconditioner");

    Ifpack2::Factory factory;
    _ifpackprec = factory.create(_ifpacktype,_k.getConst());
    _ifpackprec->setParameters(*_ifpackList);
  }
    break;
  }
}

void
TrilinosWrap::zero()
{
  set(0);
}

void
TrilinosWrap::set(double val)
{
  _localk->resumeFill();

  ST *vals = new ST[_iupper - _ilower +1];
  for (int i = 0; i < (_iupper-_ilower+1) ; ++i) {
    vals[i] = val;
  }
  
  int idx = 0;
  for (int i = 0, e = (_iupper - _ilower + 1); i < e; ++i) {
    const Teuchos::ArrayView<const int> _cids(&_Ai[idx], _nCols[i]);
    const Teuchos::ArrayView<const ST> _val(&vals[0], _nCols[i]);
    
    if (_localk->getProfileType() == Tpetra::DynamicProfile) {
    _localk->insertGlobalValues(_gRows[i], _cids, _val);
    }
    else{
    _localk->replaceGlobalValues(_gRows[i], _cids, _val);
    }
    
    idx += _nCols[i];
  }
  
  delete [] vals;
}

void
TrilinosWrap::printWithRHS(std::string&& basename, const double rhs[], int ndofd,
                    int rank) const
{
  std::stringstream filename;
  filename.str(basename);
  filename << "k.txt";
  
  std::ofstream os;
  os.open(filename.str());
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (os));
  _k->describe(*out,Teuchos::VERB_EXTREME);
  filename.str(basename);
  filename << "f.txt." << std::setfill('0') << std::setw(5) << rank;

  
  os << _ilower << " " << _iupper << "\n";
  for (int i = 0; i < ndofd; ++i) {
    os << _ilower + i << " " << rhs[i] << "\n";
  }
  os.close();
}

double
TrilinosWrap::solveSystem(const PGFem3D_opt *opts,
                   double *loc_rhs,
                   double *loc_sol,
                   const int tim,
                   const int iter,
                   const long *DomDof,
                   SOLVER_INFO *info)
{
  assert(opts->solverpackage == TRILINOS);
  double func_time = -MPI_Wtime();
  int rank=0;
  MPI_Comm_rank(_comm, &rank);

  /* Assemble the rhs and solution vector */
  nulld(loc_sol, DomDof[rank]);

  _rhs->sync<Kokkos::HostSpace>();
  auto rhs_2d = _rhs->getLocalView<Kokkos::HostSpace>();
  auto rhs_1d = Kokkos::subview(rhs_2d, Kokkos::ALL(), 0);
  _rhs->modify<Kokkos::HostSpace>();
  const size_t localLength = _rhs->getLocalLength();

  if (DomDof[rank] != (long)localLength) {
    std::cerr << "DomDof[rank] does not match rhs vector length." << std::endl;
  }
  
  for (size_t k = 0; k < localLength; ++k) {
    rhs_1d(k) = loc_rhs[k];
  }
  
  using memory_space = MV::device_type::memory_space;
  _rhs->sync<memory_space> ();

  _solution->sync<Kokkos::HostSpace>();
  auto solution_2d = _solution->getLocalView<Kokkos::HostSpace>();
  auto solution_1d = Kokkos::subview(solution_2d, Kokkos::ALL(), 0);
  _solution->modify<Kokkos::HostSpace>();
  
  for (size_t k = 0; k < localLength; ++k) {
    solution_1d(k) = loc_sol[k];
  }
  
  _solution->sync<memory_space>();

  setupSolverEnv();

  solve(info);

  for (size_t k = 0; k < localLength; ++k) {
    loc_sol[k] = solution_1d(k);
  }
  
  /* update timer and return */
  func_time += MPI_Wtime();
  
  return func_time;
}

double
TrilinosWrap::solveSystemNoSetup(const PGFem3D_opt *opts,
                          double *loc_rhs,
                          double *loc_sol,
                          const int tim,
                          const int iter,
                          const long *DomDof,
                          SOLVER_INFO *info)
{
  assert(opts->solverpackage == TRILINOS);
  double func_time = -MPI_Wtime();
  int rank=0;
  MPI_Comm_rank(_comm, &rank);

  /* Assemble the rhs and solution vector */
  nulld(loc_sol, DomDof[rank]);

  _rhs->sync<Kokkos::HostSpace>();
  auto rhs_2d = _rhs->getLocalView<Kokkos::HostSpace>();
  auto rhs_1d = Kokkos::subview (rhs_2d, Kokkos::ALL(), 0);
  _rhs->modify<Kokkos::HostSpace>();
  const size_t localLength = _rhs->getLocalLength();

  if (DomDof[rank] != (long)localLength) {
    std::cerr << "DomDof[rank] does not match rhs vector length." << std::endl;
  }

  for (size_t k = 0; k < localLength; ++k) {
    rhs_1d(k) = loc_rhs[k];
  }
  
  using memory_space = MV::device_type::memory_space;
  _rhs->sync<memory_space>();

  _solution->sync<Kokkos::HostSpace>();
  auto solution_2d = _solution->getLocalView<Kokkos::HostSpace>();
  auto solution_1d = Kokkos::subview (solution_2d, Kokkos::ALL(), 0);
  _solution->modify<Kokkos::HostSpace>();
  
  for (size_t k = 0; k < localLength; ++k) {
    solution_1d(k) = loc_sol[k];
  }
  
  _solution->sync<memory_space>();

  solve(info);

  for (size_t k = 0; k < localLength; ++k) {
    loc_sol[k] = solution_1d(k);
  }
  
  /* update timer and return */
  func_time += MPI_Wtime();
  return func_time;
}

int
TrilinosWrap::setupSolverEnv()
{
  int err = 0;
  assert(_k->isFillComplete());

  switch (_prectype) {
  case PRECOND_MUELU: {
    if (_mueluprec == Teuchos::null)
      _mueluprec = MueLu::CreateTpetraPreconditioner(_k, *_mueluList);
    else
      MueLu::ReuseTpetraPreconditioner(_k,*_mueluprec);
    
    _preconditioner = _mueluprec;
  }
    break;
  case PRECOND_IFPACK2: {
    // This dynamic casting shouldn't be necessary but the
    // Ifpack2::Preconditioner base class does not have member setMatrix().
    // The derived classes inherit this from Ifpack2::Details::CanChangeMatrix.
    // If other Ifpack2 preconditioners are used, they need to be added here
    // until a better fix is found (writing our own preconditioner wrapper?).

    if (_ifpacktype == "ILUT") {
      typedef Ifpack2::ILUT<TRowMatrix> prec_cast;
      (Teuchos::rcp_dynamic_cast<prec_cast, prec_type> (_ifpackprec))
	->setMatrix(_k.getConst());
    }
    else if (_ifpacktype == "RILUK") {
      typedef Ifpack2::RILUK<TRowMatrix> prec_cast;
      (Teuchos::rcp_dynamic_cast<prec_cast, prec_type> (_ifpackprec))
	->setMatrix(_k.getConst());
    }
    else if (_ifpacktype == "SCHWARZ") {
      typedef Ifpack2::AdditiveSchwarz<TRowMatrix> prec_cast;
      (Teuchos::rcp_dynamic_cast<prec_cast, prec_type> (_ifpackprec))
	->setMatrix(_k.getConst());
    }
    else if (_ifpacktype == "RELAXATION") {
      typedef Ifpack2::Relaxation<TRowMatrix> prec_cast;
      (Teuchos::rcp_dynamic_cast<prec_cast, prec_type> (_ifpackprec))
	->setMatrix(_k.getConst());
    }
    
    _ifpackprec->initialize();
    _ifpackprec->compute();
    _preconditioner = _ifpackprec;
  }
    break;
  }

  _problem->setOperator(_k);
  
  bool set = _problem->setProblem(_solution,_rhs);
  if (set == false) {
    std::cerr << "ERROR: Belos::LinearProblem failed to set up correctly" << std::endl;
  }
  
  if (_preconditioner != Teuchos::null) {
    _problem->setLeftPrec(_preconditioner);
  }
  
  _solver =
    Teuchos::rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>( _problem, _belosList ));

  return err;
}

int
TrilinosWrap::solve(SOLVER_INFO *info)
{
  int err=0;
  
  Belos::ReturnType ret = _solver->solve();
  info->res_norm=_solver->achievedTol();
  info->n_iter=_solver->getNumIters();
  if (ret == Belos::Unconverged) {
    std::cerr << "Trilinos did not converge with residual "
	      << info->res_norm << std::endl;
  }
  
  return err;
}

