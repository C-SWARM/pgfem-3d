/* HEADER */
#include "hypre_global.h"
#include "boomerAMGInterface.h"
#include "enumerations.h"
#include "allocation.h"
#include "PGFEM_HYPRE_preconditioners.h"

#ifndef PFEM_HYPRE_DEBUG
#define PFEM_HYPRE_DEBUG 0
#endif

using pgfem3d::solvers::Hypre;

/*==== STATIC LOCAL FUNCTIONS ====*/
static int create_EUCLID(MPI_Comm comm, HYPRE_Solver& pre) {
  /* Euclid parameters */
  static const char *euclid_argv[] = {
    "-bj",
    "-rowScale",
    "-maxNzPerRow",
    "100",
    "-level",
    "1",
    ""
  };

  /* int euclid_argc = 5; */
  /* char *euclid_argv[] = {"-maxNzPerRow","100", */
  /*             "-level","1",""}; */

  /* {"-level","1","-bj","-eu_stats","1","-rowScale",
     "-sparseA","0.0","-maxNzPerRow","100"};*/

  /* if(myrank == 0){ */
  /*   PGFEM_printf ("Preconditioner : Euclid\n"); */
  /* } */


  return (HYPRE_EuclidCreate(comm, &pre) or
          HYPRE_EuclidSetParams(pre, 7, const_cast<char**>(euclid_argv)));
}

static int create_PARASAILS(MPI_Comm comm, HYPRE_Solver& pre) {
  /* ParaSails */
  static constexpr int        sym = 0;
  static constexpr int     nlevel = 1;
  static constexpr double  thresh = 0.1;
  static constexpr double  filter = 0.05;
  static constexpr double loadbal = 0.9;

  /* if(myrank == 0){ */
  /*   PGFEM_printf ("Preconditioner : ParaSails\n"); */
  /*   PGFEM_printf ("nlevel  = %d\n",nlevel); */
  /*   PGFEM_printf ("thresh  = %f\n",thresh); */
  /*   PGFEM_printf ("filter  = %f\n",filter); */
  /*   PGFEM_printf ("loadbal = %f\n",loadbal); */
  /* } */
  return (HYPRE_ParaSailsCreate(comm, &pre) or
          HYPRE_ParaSailsSetSym(pre, sym) or
          HYPRE_ParaSailsSetParams(pre, thresh, nlevel) or
          HYPRE_ParaSailsSetFilter(pre, filter) or
          HYPRE_ParaSailsSetLoadbal(pre, loadbal)
          /* or HYPRE_ParaSailsSetLogging(hypre_pc,1) */);
}

static int create_PILUT(MPI_Comm comm, HYPRE_Solver& pre) {
  static constexpr double pilut_tol = 0.1;
  static constexpr int    pilut_row = 20;
  static constexpr int  pilut_maxit = 50;

  /* if(myrank == 0){ */
  /*   PGFEM_printf ("Preconditioner : Pilut\n"); */
  /*   PGFEM_printf ("maximum iterations= %d\n",pilut_maxit); */
  /*   PGFEM_printf ("drop tolerance    = %f\n",pilut_tol); */
  /*   PGFEM_printf ("nonzeros per row  = %f\n",pilut_row); */
  /* } */

  return (HYPRE_ParCSRPilutCreate(comm, &pre) or
          HYPRE_ParCSRPilutSetMaxIter(pre, pilut_maxit) or
          HYPRE_ParCSRPilutSetDropTolerance(pre, pilut_tol) or
          HYPRE_ParCSRPilutSetFactorRowSize(pre, pilut_row));
}

static int create_BOOMER(MPI_Comm comm, HYPRE_Solver& pre) {
  boomerAMGOptions AMGOptions;
  setBoomerAMGOptions(&AMGOptions);
  return initializeBoomerAMG(&pre, &AMGOptions, comm);
}

static void destroy_preconditioner(int type, HYPRE_Solver& pre) {
  if (!pre) return;

  switch (type) {
   default:
    PGFEM_printerr("Unexpected preconditioner type (%d)\n", type);
    break;
   case PARA_SAILS:
    HYPRE_ParaSailsDestroy(pre);
    break;
   case PILUT:
    HYPRE_ParCSRPilutDestroy(pre);
    break;
   case EUCLID:
    HYPRE_EuclidDestroy(pre);
    break;
   case BOOMER:
    HYPRE_BoomerAMGDestroy(pre);
    break;
   case DIAG_SCALE:
    PGFEM_HYPRE_ScaleDiagDestroy(pre);
    break;
   case JACOBI:
    PGFEM_HYPRE_JacobiDestroy(pre);
    break;
  }

  pre = nullptr;
}

static HYPRE_Solver create_preconditioner(MPI_Comm comm, int type) {
  int error = 0;
  HYPRE_Solver pre = nullptr;
  switch (type) {
   default:
    /* if(myrank == 0) PGFEM_printf("No Preconditioner\n"); */
    break;
   case PARA_SAILS:
    error = create_PARASAILS(comm, pre);
    break;
   case PILUT:
    error = create_PILUT(comm, pre);
    break;
   case EUCLID:
    error = create_EUCLID(comm, pre);
    break;
   case BOOMER:
    error = create_BOOMER(comm, pre);
    break;
   case DIAG_SCALE:
    /* if(myrank == 0) PGFEM_printf ("Preconditioner : Diagonal Scale\n"); */
    error = PGFEM_HYPRE_ScaleDiagCreate(&pre);
    break;
   case JACOBI:
    /* if(myrank == 0) PGFEM_printf ("Preconditioner : Jacobi\n"); */
    error = PGFEM_HYPRE_JacobiCreate(&pre);
    break;
  }

  if (error) {
    PGFEM_printerr("Failed to initialize preconditioner (type %d)\n", type);
    destroy_preconditioner(type, pre);
    pre = nullptr;
  }

  return pre;
}

static void destroy_solver(int type, HYPRE_Solver& solver) {
  if (!solver) return;

  switch (type) {
   default:
    PGFEM_printerr("Unexpected solver type (%d)\n", type);
    break;
   case HYPRE_GMRES:
    HYPRE_ParCSRGMRESDestroy(solver);
    break;
   case HYPRE_BCG_STAB:
    HYPRE_ParCSRBiCGSTABDestroy(solver);
    break;
   case HYPRE_AMG:
    HYPRE_BoomerAMGDestroy(solver);
    break;
  }

  solver = nullptr;
}

static void ap_to_ncols(int *Ap, int *ncols, int size) {
  for (int n = 0; n < size; n++) {
    ncols[n] = Ap[n + 1] - Ap[n];
  }
}

static void get_diag_offd_sizes(int *ncols, int *Ai, int *dsize, int *odsize,
                                const int size, const int il, const int iu)
{
  int j = 0;

  for (int i = 0; i < size; i++) {
    dsize[i] = 0;
    for(int c = 0, e = ncols[i]; c < e; c++) {
      int col_id = Ai[j++];
      if (col_id >= il && col_id <= iu) {
        dsize[i]++;
      }
    }
    odsize[i] = ncols[i] - dsize[i];
  }
}

/*===== API FUNCTION DEFINITIONS ====*/
Hypre::~Hypre()
{
  destroy_preconditioner(precond_type, hypre_pc);
  destroy_solver(solver_type, hypre_solver);

  delete [] ncol;
  delete [] grows;

  /* destroy vectors/matrix */
  HYPRE_IJVectorDestroy(hypre_rhs);
  HYPRE_IJVectorDestroy(hypre_sol);
  HYPRE_IJMatrixDestroy(hypre_k);
}

void
Hypre::initialize(int *Ap, int *Ai, int size, int maxit, double err,
                  const PGFem3D_opt *options, MPI_Comm comm)
{
  /* override maxit with command line options */
  maxit = options->maxit;

  /* Rank and nproc */
  int myrank, nproc;
  if (MPI_Comm_size(comm, &nproc) or
      MPI_Comm_rank(comm, &myrank)) {
    PGFEM_Abort();
  }

  if (myrank == 0) {
    PGFEM_printf("Iterative solver info:\n"
                 "Kdim     = %d\n"
                 "Max. It. = %d\n", options->kdim, maxit);
  }

  /* Allocations */
  ncol = new int[size];
  grows = new int[size];

  for (int i = 0; i < size; ++i) {
    grows[i] = ilower + i;
  }

  /* Compute ncol */
  ap_to_ncols(Ap, ncol, size);

  /* Matrix */
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &hypre_k);
  HYPRE_IJMatrixSetObjectType(hypre_k, HYPRE_PARCSR);

  int *diag_sizes = new int[size];
  int *offd_sizes = new int[size];

  get_diag_offd_sizes(ncol, Ai, diag_sizes, offd_sizes, size, ilower, iupper);
  HYPRE_IJMatrixSetRowSizes(hypre_k, ncol);
  HYPRE_IJMatrixSetDiagOffdSizes(hypre_k, diag_sizes, offd_sizes);

  delete[] diag_sizes;
  delete[] offd_sizes;

  HYPRE_IJMatrixInitialize(hypre_k);
  HYPRE_IJMatrixGetObject(hypre_k, reinterpret_cast<void**>(&hypre_pk));

  if (PFEM_HYPRE_DEBUG) {
    /* print pattern */
    setK(Ai, size, 1.0);
    HYPRE_IJMatrixAssemble(hypre_k);
    HYPRE_IJMatrixPrint(hypre_k, "stiff_pattern.out");
    HYPRE_IJMatrixInitialize(hypre_k);
  }

  /* RHS */
  HYPRE_IJVectorCreate(comm, ilower, iupper, &hypre_rhs);
  HYPRE_IJVectorSetObjectType(hypre_rhs, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(hypre_rhs);
  HYPRE_IJVectorGetObject(hypre_rhs, reinterpret_cast<void**>(&hypre_prhs));

  /* Solution vector */
  HYPRE_IJVectorCreate(comm, ilower, iupper, &hypre_sol);
  HYPRE_IJVectorSetObjectType(hypre_sol, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(hypre_sol);
  HYPRE_IJVectorGetObject(hypre_sol, reinterpret_cast<void**>(&hypre_psol));

  /* set up the preconditioner */
  precond_type = options->precond;
  createPreconditioner(comm);

  /* set up the solver */
  solver_type = options->solver;

  switch (options->solver) {
   default:
   case HYPRE_GMRES:
    HYPRE_ParCSRGMRESCreate(comm, &hypre_solver);
    HYPRE_ParCSRGMRESSetMaxIter(hypre_solver, maxit);
    HYPRE_ParCSRGMRESSetTol(hypre_solver, err);
    HYPRE_ParCSRGMRESSetKDim(hypre_solver, options->kdim);
    HYPRE_ParCSRGMRESSetLogging(hypre_solver, 1);
    HYPRE_ParCSRGMRESSetPrintLevel(hypre_solver, 0);
    break;

   case HYPRE_BCG_STAB:
    HYPRE_ParCSRBiCGSTABCreate(comm, &hypre_solver);
    HYPRE_BiCGSTABSetMaxIter(hypre_solver, maxit);
    HYPRE_BiCGSTABSetTol(hypre_solver, err);
    HYPRE_BiCGSTABSetLogging(hypre_solver, 1);
    break;

   case HYPRE_AMG:
     {
       boomerAMGOptions AMGOptions;
       setBoomerAMGOptions(&AMGOptions);
       initializeBoomerAMG(&hypre_solver, &AMGOptions, comm);
       HYPRE_BoomerAMGSetTol(hypre_solver, err);
       HYPRE_BoomerAMGSetMaxIter(hypre_solver, 20);
     }
     break;

   case HYPRE_FLEX:
    HYPRE_ParCSRFlexGMRESCreate(comm, &hypre_solver);
    HYPRE_ParCSRFlexGMRESSetMaxIter(hypre_solver, maxit);
    HYPRE_ParCSRFlexGMRESSetTol(hypre_solver, err);
    HYPRE_ParCSRFlexGMRESSetKDim(hypre_solver, options->kdim);
    HYPRE_ParCSRFlexGMRESSetLogging(hypre_solver, 1);
    HYPRE_ParCSRFlexGMRESSetPrintLevel(hypre_solver, 0);
    break;

   case HYPRE_HYBRID:
    HYPRE_ParCSRHybridCreate(&hypre_solver);
    HYPRE_ParCSRHybridSetKDim(hypre_solver, options->kdim);
    HYPRE_ParCSRHybridSetTol(hypre_solver, err);
    HYPRE_ParCSRHybridSetDSCGMaxIter(hypre_solver, maxit);
    HYPRE_ParCSRHybridSetSolverType(hypre_solver, 2);
    HYPRE_ParCSRHybridSetConvergenceTol(hypre_solver, 1.0);
    break;
  }
}

void
Hypre::setRowColBounds(long g_n_col, const long *n_row_proc, int rank)
{
  jlower = 0;
  jupper = g_n_col - 1;
  for (int i = 0; i < rank; ++i) {
    ilower += n_row_proc[i];
  }
  iupper = ilower + n_row_proc[rank] - 1;
}

int
Hypre::createPreconditioner(MPI_Comm comm)
{
  hypre_pc = create_preconditioner(comm, precond_type);
  return (hypre_pc != nullptr);
}

int
Hypre::destroyPreconditioner()
{
  destroy_preconditioner(precond_type, hypre_pc);
  return 0;
}

void
Hypre::setK(int *Ai, int size, double val)
{
  HYPRE_Int cols[1] = {1};
  int curpos = 0;

  for(int i = 0, e = size; i < e; ++i) {
    for(int j = 0, e = ncol[i]; j < e; ++j) {
      HYPRE_IJMatrixSetValues(hypre_k, 1, cols, &grows[i], &Ai[curpos++], &val);
    }
  }
}
