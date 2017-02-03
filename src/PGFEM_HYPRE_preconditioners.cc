/* HEADER */

/*
 * This file contains definitions of special preconditioners developed
 * for use with the HYPRE package
 */

#include "PGFEM_HYPRE_preconditioners.h"
#include <math.h>
#include "_hypre_parcsr_mv.h"
#include "PGFEM_mpi.h"

#include "allocation.h"

typedef struct{
  double *data;
  int size;
  int setup;
  double thresh;
  double replacement;
} PGFEM_diag_pc;

int PGFEM_HYPRE_ScaleDiagCreate(HYPRE_Solver *vdiag_pc)
{
  PGFEM_diag_pc **diag_pc = (PGFEM_diag_pc**) vdiag_pc;
  /* if(diag_pc != NULL){ */
  /*   fprintf(stderr,"WARNING: Preconditioner already exists!\n"); */
  /*   /\* PGFEM_Abort(); *\/ */
  /* } */
  (*diag_pc) = PGFEM_calloc(1,sizeof(PGFEM_diag_pc));
  (*diag_pc)->data = NULL;
  (*diag_pc)->thresh = 1.0e-5;
  (*diag_pc)->replacement = 1.1;
  (*diag_pc)->setup = 0;
  return 0;
}

int PGFEM_HYPRE_ScaleDiagSetup(HYPRE_Solver vdiag_pc,
			       HYPRE_ParCSRMatrix vA,
			       HYPRE_ParVector vb,
			       HYPRE_ParVector vx)
{
  /* setup scale vector accounting for zeros on the diagonal */
  PGFEM_diag_pc *diag_pc = (PGFEM_diag_pc*) vdiag_pc;
  hypre_ParCSRMatrix *A = (hypre_ParCSRMatrix *) vA;
  hypre_ParVector    *x = (hypre_ParVector *) vx;
  const double *A_data = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(A));
  const int *A_i = hypre_CSRMatrixI(hypre_ParCSRMatrixDiag(A));

  const int local_size = hypre_VectorSize(hypre_ParVectorLocalVector(x));
  if(diag_pc->setup){
    free(diag_pc->data);
    diag_pc->setup = 0;
  }
  diag_pc->data = PGFEM_calloc(local_size,sizeof(double));
  diag_pc->size = local_size;

  double val;
  for(int i=0; i<local_size; i++){
    /* diagonal value is always stored first on the local row */
    val = A_data[A_i[i]];
    if (fabs(val) > diag_pc->thresh){
      diag_pc->data[i] = 1./val;
    } else {
      diag_pc->data[i] = diag_pc->replacement;
    }
  }
  diag_pc->setup = 1;
  return diag_pc->setup;
}

int PGFEM_HYPRE_ScaleDiagSolve(HYPRE_Solver vdiag_pc,
			       HYPRE_ParCSRMatrix vA,
			       HYPRE_ParVector vb,
			       HYPRE_ParVector vx)
{
  PGFEM_diag_pc *diag_pc = (PGFEM_diag_pc*) vdiag_pc;
  hypre_ParVector    *x = (hypre_ParVector *) vx;
  hypre_ParVector    *b = (hypre_ParVector *) vb;
  double *x_data = hypre_VectorData(hypre_ParVectorLocalVector(x));
  const double *b_data = hypre_VectorData(hypre_ParVectorLocalVector(b));
  for(int i=0; i<diag_pc->size; i++){
    x_data[i] = b_data[i]*diag_pc->data[i];
  }
  return 0;
}

int PGFEM_HYPRE_ScaleDiagDestroy(HYPRE_Solver vdiag_pc)
{
  PGFEM_diag_pc *diag_pc = (PGFEM_diag_pc*) vdiag_pc;
  free(diag_pc->data);
  free(diag_pc);
  return 0;
}

/*===================================================================
 *                     JACOBI PRECONDITIONER
 *=================================================================*/


typedef struct{
  double scale;
  double thresh;
  double replace;
  int global;
  int size;
  MPI_Comm comm;
} PGFEM_jacobi_pc;

int PGFEM_HYPRE_JacobiCreate(HYPRE_Solver *vjacobi_pc)
{
  PGFEM_jacobi_pc **jacobi_pc = (PGFEM_jacobi_pc **) vjacobi_pc;
  (*jacobi_pc) = PGFEM_calloc(1,sizeof(PGFEM_jacobi_pc));
  (*jacobi_pc)->thresh = 1.e-5;
  (*jacobi_pc)->replace = 1.0;
  (*jacobi_pc)->scale = 1.0;
  (*jacobi_pc)->global = 0;
  return 0;
}

int PGFEM_HYPRE_JacobiSetup(HYPRE_Solver vjacobi_pc,
			    HYPRE_ParCSRMatrix vA,
			    HYPRE_ParVector vb,
			    HYPRE_ParVector vx)
{
  PGFEM_jacobi_pc *jacobi_pc = (PGFEM_jacobi_pc *) vjacobi_pc;
  hypre_ParCSRMatrix *A = (hypre_ParCSRMatrix *) vA;
  hypre_ParVector    *x = (hypre_ParVector *) vx;

  const double *A_data = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(A));
  const int *A_i = hypre_CSRMatrixI(hypre_ParCSRMatrixDiag(A));
  const int local_size = hypre_VectorSize(hypre_ParVectorLocalVector(x));

  jacobi_pc->scale = 0.0;
  jacobi_pc->comm = hypre_ParCSRMatrixComm(A);
  jacobi_pc->size = local_size;
  for(int i=0; i<local_size; i++){
    jacobi_pc->scale += A_data[A_i[i]];
  }

  if(jacobi_pc->global){
    MPI_Allreduce(MPI_IN_PLACE,&(jacobi_pc->scale),1,
		  MPI_DOUBLE,MPI_SUM,jacobi_pc->comm);
  }

  if(fabs(jacobi_pc->scale) < jacobi_pc->thresh){
    jacobi_pc->scale = jacobi_pc->replace;
  } else {
    jacobi_pc->scale = 1.0/jacobi_pc->scale;
  }

  return 0;
}

int PGFEM_HYPRE_JacobiSolve(HYPRE_Solver vjacobi_pc,
			    HYPRE_ParCSRMatrix vA,
			    HYPRE_ParVector vb,
			    HYPRE_ParVector vx)
{
  PGFEM_jacobi_pc *jacobi_pc = (PGFEM_jacobi_pc *) vjacobi_pc;
  hypre_ParVector    *x = (hypre_ParVector *) vx;
  hypre_ParVector    *b = (hypre_ParVector *) vb;
  double *x_data = hypre_VectorData(hypre_ParVectorLocalVector(x));
  const double *b_data = hypre_VectorData(hypre_ParVectorLocalVector(b));
  for(int i=0; i<jacobi_pc->size; i++){
    x_data[i] = b_data[i]*jacobi_pc->scale;
  }
  return 0;
}

int PGFEM_HYPRE_JacobiDestroy(HYPRE_Solver vjacobi_pc)
{
  PGFEM_jacobi_pc *jacobi_pc = (PGFEM_jacobi_pc *) vjacobi_pc;
  free(jacobi_pc);
  return 0;
}
