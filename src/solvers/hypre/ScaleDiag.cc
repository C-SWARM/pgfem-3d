#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"
#include <_hypre_parcsr_mv.h>
#include <stdexcept>

using pgfem3d::solvers::hypre::Preconditioner;
using pgfem3d::solvers::hypre::ScaleDiag;

namespace {
struct PGFEM_diag_pc {
  double *data;
  int size;
  int setup;
  double thresh;
  double replacement;
};

HYPRE_t PGFEM_HYPRE_ScaleDiagCreate(HYPRE_Solver *vdiag_pc)
{
  PGFEM_diag_pc **diag_pc = (PGFEM_diag_pc**) vdiag_pc;
  /* if(diag_pc != NULL){ */
  /*   fprintf(stderr,"WARNING: Preconditioner already exists!\n"); */
  /*   /\* PGFEM_Abort(); *\/ */
  /* } */
  (*diag_pc) = PGFEM_calloc(PGFEM_diag_pc, 1);
  (*diag_pc)->data = NULL;
  (*diag_pc)->thresh = 1.0e-5;
  (*diag_pc)->replacement = 1.1;
  (*diag_pc)->setup = 0;
  return 0;
}

HYPRE_t PGFEM_HYPRE_ScaleDiagSetup(HYPRE_Solver vdiag_pc,
                                   HYPRE_ParCSRMatrix vA,
                                   HYPRE_ParVector vb,
                                   HYPRE_ParVector vx)
{
  /* setup scale vector accounting for zeros on the diagonal */
  PGFEM_diag_pc *diag_pc = (PGFEM_diag_pc*) vdiag_pc;
  hypre_ParCSRMatrix *A = (hypre_ParCSRMatrix *) vA;
  hypre_ParVector    *x = (hypre_ParVector *) vx;
  const double *A_data = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(A));
  const HYPRE_t *A_i = hypre_CSRMatrixI(hypre_ParCSRMatrixDiag(A));      // ksaha ???

  const int local_size = hypre_VectorSize(hypre_ParVectorLocalVector(x));
  if(diag_pc->setup){
    free(diag_pc->data);
    diag_pc->setup = 0;
  }
  diag_pc->data = PGFEM_calloc(double, local_size);
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

HYPRE_t PGFEM_HYPRE_ScaleDiagSolve(HYPRE_Solver vdiag_pc,
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

HYPRE_t PGFEM_HYPRE_ScaleDiagDestroy(HYPRE_Solver vdiag_pc)
{
  PGFEM_diag_pc *diag_pc = (PGFEM_diag_pc*) vdiag_pc;
  free(diag_pc->data);
  free(diag_pc);
  return 0;
}
}

ScaleDiag::ScaleDiag(MPI_Comm comm) : Preconditioner(comm)
{
  if (PGFEM_HYPRE_ScaleDiagCreate(&_solver)) {
    throw std::runtime_error("failed to initialize ScaleDiag preconditioner");
  }
}

ScaleDiag::~ScaleDiag()
{
  PGFEM_HYPRE_ScaleDiagDestroy(_solver);
}

void
ScaleDiag::getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int) const
{
  solve = PGFEM_HYPRE_ScaleDiagSolve;
  setup = PGFEM_HYPRE_ScaleDiagSetup;
}

void
ScaleDiag::reset()
{
  PGFEM_HYPRE_ScaleDiagDestroy(_solver);
  PGFEM_HYPRE_ScaleDiagCreate(&_solver);
}
