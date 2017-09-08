#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Preconditioners.hpp"
#include <_hypre_parcsr_mv.h>
#include <stdexcept>

using pgfem3d::solvers::hypre::Preconditioner;
using pgfem3d::solvers::hypre::Jacobi;

/*===================================================================
 *                     JACOBI PRECONDITIONER
 *=================================================================*/

namespace {
struct PGFEM_jacobi_pc {
  double scale;
  double thresh;
  double replace;
  int global;
  int size;
  MPI_Comm comm;
};

int PGFEM_HYPRE_JacobiCreate(HYPRE_Solver *vjacobi_pc)
{
  PGFEM_jacobi_pc **jacobi_pc = (PGFEM_jacobi_pc **) vjacobi_pc;
  (*jacobi_pc) = PGFEM_calloc(PGFEM_jacobi_pc, 1);
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
}

Jacobi::Jacobi(MPI_Comm comm) : Preconditioner(comm)
{
  if (PGFEM_HYPRE_JacobiCreate(&_solver)) {
    throw std::runtime_error("failed to initialize Jacobi preconditioner");
  }
}

Jacobi::~Jacobi()
{
  PGFEM_HYPRE_JacobiDestroy(_solver);
}

void
Jacobi::getFuncs(ptr_solve_t& solve, ptr_solve_t& setup, int) const
{
  solve = PGFEM_HYPRE_JacobiSolve;
  setup = PGFEM_HYPRE_JacobiSetup;
}

void
Jacobi::reset()
{
  PGFEM_HYPRE_JacobiDestroy(_solver);
  PGFEM_HYPRE_JacobiCreate(&_solver);
}
