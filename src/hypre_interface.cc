/* HEADER */
#include "hypre_global.h"
#include "boomerAMGInterface.h"
#include "enumerations.h"
#include "allocation.h"
#include "PGFEM_HYPRE_preconditioners.h"

#ifndef PFEM_HYPRE_DEBUG
#define PFEM_HYPRE_DEBUG 0
#endif

/*==== STATIC FUNCTION PROTOTYPES ====*/
static int create_precond_EUCLID(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                 const MPI_Comm mpi_comm);

static int create_precond_PARASAILS(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                    const MPI_Comm mpi_comm);

static int create_precond_PILUT(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                const MPI_Comm mpi_comm);

static int create_precond_BOOMER(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                 const MPI_Comm mpi_comm);

/*===== API FUNCTION DEFINITIONS ====*/
void initialize_PGFEM_HYPRE_solve_info(PGFEM_HYPRE_solve_info **info)
{
  (*info) = PGFEM_calloc(PGFEM_HYPRE_solve_info, 1);
  (*info)->hypre_k = NULL;
  (*info)->hypre_pk = NULL;
  (*info)->hypre_solver = NULL;
  (*info)->hypre_pc = NULL;
  (*info)->hypre_pc_gotten = NULL;
  (*info)->hypre_rhs = NULL;
  (*info)->hypre_prhs = NULL;
  (*info)->hypre_sol = NULL;
  (*info)->hypre_psol = NULL;

  (*info)->ncol = NULL;
  (*info)->grows = NULL;
  (*info)->ilower = 0;
  (*info)->iupper = 0;
  (*info)->jlower = 0;
  (*info)->jupper = 0;

  (*info)->precond_type = 0;
  (*info)->solver_type = 0;
}

void destroy_PGFEM_HYPRE_solve_info(PGFEM_HYPRE_solve_info *info)
{
  if(info != NULL){
    destroy_HYPRE_preconditioner(info);

    switch(info->solver_type){
    case HYPRE_GMRES:
      HYPRE_ParCSRGMRESDestroy (info->hypre_solver);
      break;
    case HYPRE_BCG_STAB:
      HYPRE_ParCSRBiCGSTABDestroy (info->hypre_solver);
      break;
    case HYPRE_AMG:
      HYPRE_BoomerAMGDestroy (info->hypre_solver);
      break;
    }

    free(info->ncol);
    free(info->grows);

    /* destroy vectors/matrix */
    HYPRE_IJVectorDestroy(info->hypre_rhs);
    HYPRE_IJVectorDestroy(info->hypre_sol);
    HYPRE_IJMatrixDestroy(info->hypre_k);
  }
  free(info);
}

void hypre_initialize(int *Ap,
              int *Ai,
              int size,
              int maxit,
              double err,
              PGFEM_HYPRE_solve_info *PGFEM_hypre,
              const PGFem3D_opt *options,
              MPI_Comm mpi_comm)
{
  int kdim = 0;
  int *diag_sizes = NULL;
  int *offd_sizes = NULL;
  int i = 0;
  int func_err = 0;

  /* override maxit with command line options */
  maxit = options->maxit;
  kdim = options->kdim;

  /* Rank and nproc */
  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  if(myrank == 0){
    PGFEM_printf("Iterative solver info:\n"
         "Kdim     = %d\n"
         "Max. It. = %d\n",kdim,maxit);
  }


  /* Allocations */
  PGFEM_hypre->ncol  = PGFEM_calloc(int, size);
  diag_sizes         = PGFEM_calloc(int, size);
  offd_sizes         = PGFEM_calloc(int, size);
  PGFEM_hypre->grows = PGFEM_calloc(int, size);

  i = 0;
  while(i < size){
    PGFEM_hypre->grows[i] = PGFEM_hypre->ilower+i;
    i++;
  }

  /* Compute ncol */
  Ap2ncols(Ap,PGFEM_hypre->ncol,size);

  /* Matrix */
  HYPRE_IJMatrixCreate(mpi_comm,PGFEM_hypre->ilower,PGFEM_hypre->iupper,
               PGFEM_hypre->ilower,PGFEM_hypre->iupper,
               &PGFEM_hypre->hypre_k);
  HYPRE_IJMatrixSetObjectType(PGFEM_hypre->hypre_k,HYPRE_PARCSR);
  GetDiagOffdSizes(PGFEM_hypre->ncol,Ai,diag_sizes,offd_sizes,size,
           PGFEM_hypre->ilower,PGFEM_hypre->iupper);
  HYPRE_IJMatrixSetRowSizes(PGFEM_hypre->hypre_k,PGFEM_hypre->ncol);
  HYPRE_IJMatrixSetDiagOffdSizes(PGFEM_hypre->hypre_k,
                 diag_sizes,offd_sizes);
  HYPRE_IJMatrixInitialize(PGFEM_hypre->hypre_k);
  HYPRE_IJMatrixGetObject(PGFEM_hypre->hypre_k,
              (void **) &PGFEM_hypre->hypre_pk);

  if(PFEM_HYPRE_DEBUG){
    /* print pattern */
    SetHypreK(PGFEM_hypre,Ai,size,1.0);
    HYPRE_IJMatrixAssemble(PGFEM_hypre->hypre_k);
    HYPRE_IJMatrixPrint(PGFEM_hypre->hypre_k,"stiff_pattern.out");
    HYPRE_IJMatrixInitialize(PGFEM_hypre->hypre_k);
  }


  /* RHS */
  HYPRE_IJVectorCreate(mpi_comm,PGFEM_hypre->ilower,PGFEM_hypre->iupper,
               &PGFEM_hypre->hypre_rhs);
  HYPRE_IJVectorSetObjectType(PGFEM_hypre->hypre_rhs,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(PGFEM_hypre->hypre_rhs);
  HYPRE_IJVectorGetObject(PGFEM_hypre->hypre_rhs,
              (void **)&PGFEM_hypre->hypre_prhs);

  /* Solution vector */
  HYPRE_IJVectorCreate(mpi_comm,PGFEM_hypre->ilower,PGFEM_hypre->iupper,
               &PGFEM_hypre->hypre_sol);
  HYPRE_IJVectorSetObjectType(PGFEM_hypre->hypre_sol,HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(PGFEM_hypre->hypre_sol);
  HYPRE_IJVectorGetObject(PGFEM_hypre->hypre_sol,
              (void **)&PGFEM_hypre->hypre_psol);


  /* set up the preconditioner */
  PGFEM_hypre->precond_type = options->precond;
  func_err += PGFEM_HYPRE_create_preconditioner(PGFEM_hypre,
                       mpi_comm);

  /* set up the solver */
  PGFEM_hypre->solver_type = options->solver;
  switch(options->solver){
  default:
  case HYPRE_GMRES:
    HYPRE_ParCSRGMRESCreate(mpi_comm,&PGFEM_hypre->hypre_solver);
    HYPRE_ParCSRGMRESSetMaxIter(PGFEM_hypre->hypre_solver,maxit);
    HYPRE_ParCSRGMRESSetTol(PGFEM_hypre->hypre_solver,err);
    HYPRE_ParCSRGMRESSetKDim(PGFEM_hypre->hypre_solver,kdim);
    HYPRE_ParCSRGMRESSetLogging(PGFEM_hypre->hypre_solver,1);
    HYPRE_ParCSRGMRESSetPrintLevel(PGFEM_hypre->hypre_solver,0);
    break;

  case HYPRE_BCG_STAB:
    HYPRE_ParCSRBiCGSTABCreate(mpi_comm, &PGFEM_hypre->hypre_solver);
    HYPRE_BiCGSTABSetMaxIter(PGFEM_hypre->hypre_solver,maxit);
    HYPRE_BiCGSTABSetTol(PGFEM_hypre->hypre_solver,err);
    HYPRE_BiCGSTABSetLogging(PGFEM_hypre->hypre_solver,1);
    break;

  case HYPRE_AMG:
    {
      boomerAMGOptions AMGOptions;
      setBoomerAMGOptions(&AMGOptions);
      initializeBoomerAMG(&PGFEM_hypre->hypre_solver,&AMGOptions,mpi_comm);
      HYPRE_BoomerAMGSetTol(PGFEM_hypre->hypre_solver,err);
      HYPRE_BoomerAMGSetMaxIter(PGFEM_hypre->hypre_solver,20);
    }
    break;

  case HYPRE_FLEX:
    HYPRE_ParCSRFlexGMRESCreate(mpi_comm,&PGFEM_hypre->hypre_solver);
    HYPRE_ParCSRFlexGMRESSetMaxIter(PGFEM_hypre->hypre_solver,maxit);
    HYPRE_ParCSRFlexGMRESSetTol(PGFEM_hypre->hypre_solver,err);
    HYPRE_ParCSRFlexGMRESSetKDim(PGFEM_hypre->hypre_solver,kdim);
    HYPRE_ParCSRFlexGMRESSetLogging(PGFEM_hypre->hypre_solver,1);
    HYPRE_ParCSRFlexGMRESSetPrintLevel(PGFEM_hypre->hypre_solver,0);
    break;

  case HYPRE_HYBRID:
    HYPRE_ParCSRHybridCreate(&PGFEM_hypre->hypre_solver);
    HYPRE_ParCSRHybridSetKDim(PGFEM_hypre->hypre_solver,kdim);
    HYPRE_ParCSRHybridSetTol(PGFEM_hypre->hypre_solver,err);
    HYPRE_ParCSRHybridSetDSCGMaxIter(PGFEM_hypre->hypre_solver,maxit);
    HYPRE_ParCSRHybridSetSolverType(PGFEM_hypre->hypre_solver,2);
    HYPRE_ParCSRHybridSetConvergenceTol(PGFEM_hypre->hypre_solver,1.0);
    break;
  }

  free(diag_sizes); free(offd_sizes);
}

void Ap2ncols(int *Ap,
          int *ncols,
          int size)
{
  int n = 0;

  for(n = 0;n < size;n++)
    ncols[n] = Ap[n+1] - Ap[n];
}

void GetDiagOffdSizes(int *ncols,
              int *Ai,
              int *dsize,
              int *odsize,
              int size,
              int il,
              int iu)
{
  int i = 0;
  int c = 0;
  int j = 0;
  int col_id = 0;
  int diag_count = 0;

  for(i = 0;i < size;i++){
    diag_count = 0;
    for(c = 0;c < ncols[i];c++){
      col_id = Ai[j++];
      if(col_id >= il && col_id <= iu)
    diag_count++;
    }
    dsize[i] = diag_count;
    odsize[i] = ncols[i] - diag_count;
  }
}

void ZeroHypreK(PGFEM_HYPRE_solve_info *PGFEM_hypre,
        int *Ai,
        int size)
{
  int i = 0;
  int j = 0;
  int curpos = 0;
  double values[1] = {0};
  int cols[1] = {1};

  for(i = 0;i < size;i++)
    for(j = 0;j < PGFEM_hypre->ncol[i];j++)
      HYPRE_IJMatrixSetValues(PGFEM_hypre->hypre_k,1,cols,
                  &PGFEM_hypre->grows[i],
                  &Ai[curpos++],values);
}

void SetHypreK(PGFEM_HYPRE_solve_info *PGFEM_hypre,
           int *Ai,
           int size,
           double val)
{
  int i = 0;
  int j = 0;
  int curpos = 0;
  double *values = &val;
  int cols[1] = {1};

  for(i = 0;i < size;i++)
    for(j = 0;j < PGFEM_hypre->ncol[i];j++)
      HYPRE_IJMatrixSetValues(PGFEM_hypre->hypre_k,1,cols,
                  &PGFEM_hypre->grows[i],
                  &Ai[curpos++],values);
}

double* hypreGetDiagonal(HYPRE_IJMatrix A,
             int lower,
             int upper)
{
  int nrows = upper - lower + 1;
  int *ncols, *rows;
  int i;
  ncols = PGFEM_calloc(int, nrows);
  rows = PGFEM_calloc(int, nrows);

  for(i=0; i<nrows;i++){
    ncols[i] = 1;
    rows[i] = i+lower;
  }

  double *values = PGFEM_calloc(double, nrows);

  HYPRE_IJMatrixGetValues(A,nrows,ncols,rows,rows,values);

  return values;

}

void set_HYPRE_row_col_bounds(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                  const long g_n_col,
                  const long *n_row_proc,
                  const int myrank)
{
  int idx = 0;
  PGFEM_hypre->jlower = 0;
  PGFEM_hypre->jupper = g_n_col-1;
  while(idx < myrank){
    PGFEM_hypre->ilower += n_row_proc[idx];
    idx++;
  }
  PGFEM_hypre->iupper = PGFEM_hypre->ilower + n_row_proc[idx] - 1;

}

int PGFEM_HYPRE_create_preconditioner(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                      const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);

  switch(PGFEM_hypre->precond_type) {
  case PARA_SAILS:
    err += create_precond_PARASAILS(PGFEM_hypre,mpi_comm);
    break;

  case PILUT:
    err += create_precond_PILUT(PGFEM_hypre,mpi_comm);
    break;

  case EUCLID:
    err += create_precond_EUCLID(PGFEM_hypre,mpi_comm);
    break;

  case BOOMER:
    err += create_precond_BOOMER(PGFEM_hypre,mpi_comm);
    break;

  case DIAG_SCALE:
    /* if(myrank == 0) PGFEM_printf ("Preconditioner : Diagonal Scale\n"); */
    err += PGFEM_HYPRE_ScaleDiagCreate(&PGFEM_hypre->hypre_pc);
    break;

  case JACOBI:
    /* if(myrank == 0) PGFEM_printf ("Preconditioner : Jacobi\n"); */
    err += PGFEM_HYPRE_JacobiCreate(&PGFEM_hypre->hypre_pc);
    break;

  default:
    /* if(myrank == 0) PGFEM_printf("No Preconditioner\n"); */
    break;
  }

  return err;
}

int destroy_HYPRE_preconditioner(PGFEM_HYPRE_solve_info *info)
{
  int err = 0;
  if(info->hypre_pc != NULL){
    switch(info->precond_type){
    case PARA_SAILS:
      HYPRE_ParaSailsDestroy (info->hypre_pc);
      break;
    case PILUT:
      HYPRE_ParCSRPilutDestroy (info->hypre_pc);
      break;
    case EUCLID:
      HYPRE_EuclidDestroy (info->hypre_pc);
      break;
    case BOOMER:
      HYPRE_BoomerAMGDestroy (info->hypre_pc);
      break;
    case DIAG_SCALE:
      PGFEM_HYPRE_ScaleDiagDestroy(info->hypre_pc);
      break;
    case JACOBI:
      PGFEM_HYPRE_JacobiDestroy(info->hypre_pc);
      break;
    }
  }

  /* reset to NULL in case this function is called again from full
     destructor */
  info->hypre_pc = NULL;
  return err;
}

/*==== STATIC FUNCTION DEFINITIONS ====*/
static int create_precond_EUCLID(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                 const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);
  /* Euclid parameters */
  int euclid_argc = 7;
  const char *euclid_argv[] = {"-bj",
             "-rowScale",
             "-maxNzPerRow","100",
             "-level","1",""};
  /* int euclid_argc = 5; */
  /* char *euclid_argv[] = {"-maxNzPerRow","100", */
  /*             "-level","1",""}; */

  /* {"-level","1","-bj","-eu_stats","1","-rowScale",
     "-sparseA","0.0","-maxNzPerRow","100"};*/

  /* if(myrank == 0){ */
  /*   PGFEM_printf ("Preconditioner : Euclid\n"); */
  /* } */

  err += HYPRE_EuclidCreate(mpi_comm, &PGFEM_hypre->hypre_pc);
  err += HYPRE_EuclidSetParams (PGFEM_hypre->hypre_pc,
      euclid_argc,const_cast<char**>(euclid_argv));

  return err;
}

static int create_precond_PARASAILS(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                    const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);

  /* ParaSails */
  int sym     = 0;
  int nlevel  = 1;
  double thresh  = 0.1;
  double filter  = 0.05;
  double loadbal = 0.9;

  /* if(myrank == 0){ */
  /*   PGFEM_printf ("Preconditioner : ParaSails\n"); */
  /*   PGFEM_printf ("nlevel  = %d\n",nlevel); */
  /*   PGFEM_printf ("thresh  = %f\n",thresh); */
  /*   PGFEM_printf ("filter  = %f\n",filter); */
  /*   PGFEM_printf ("loadbal = %f\n",loadbal); */
  /* } */

  err += HYPRE_ParaSailsCreate(mpi_comm,&PGFEM_hypre->hypre_pc);
  err += HYPRE_ParaSailsSetSym(PGFEM_hypre->hypre_pc,sym);
  err += HYPRE_ParaSailsSetParams(PGFEM_hypre->hypre_pc,thresh,nlevel);
  err += HYPRE_ParaSailsSetFilter(PGFEM_hypre->hypre_pc,filter);
  err += HYPRE_ParaSailsSetLoadbal(PGFEM_hypre->hypre_pc,loadbal);
  /* err += HYPRE_ParaSailsSetLogging(hypre_pc,1); */

  return err;
}

static int create_precond_PILUT(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);

  double pilut_tol = 0.1;
  int pilut_row = 20;
  int pilut_maxit = 50;

  /* if(myrank == 0){ */
  /*   PGFEM_printf ("Preconditioner : Pilut\n"); */
  /*   PGFEM_printf ("maximum iterations= %d\n",pilut_maxit); */
  /*   PGFEM_printf ("drop tolerance    = %f\n",pilut_tol); */
  /*   PGFEM_printf ("nonzeros per row  = %f\n",pilut_row); */
  /* } */

  err += HYPRE_ParCSRPilutCreate(mpi_comm,&PGFEM_hypre->hypre_pc);
  err += HYPRE_ParCSRPilutSetMaxIter(PGFEM_hypre->hypre_pc,pilut_maxit);
  err += HYPRE_ParCSRPilutSetDropTolerance(PGFEM_hypre->hypre_pc,pilut_tol);
  err += HYPRE_ParCSRPilutSetFactorRowSize(PGFEM_hypre->hypre_pc,pilut_row);

  return err;
}

static int create_precond_BOOMER(PGFEM_HYPRE_solve_info *PGFEM_hypre,
                 const MPI_Comm mpi_comm)
{
  int err = 0;
  boomerAMGOptions AMGOptions;
  setBoomerAMGOptions(&AMGOptions);
  err += initializeBoomerAMG(&PGFEM_hypre->hypre_pc,&AMGOptions,mpi_comm);
  return err;
}
