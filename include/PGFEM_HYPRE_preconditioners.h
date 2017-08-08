/* HEADER */
#pragma once
#ifndef PGFEM_HYPRE_PRECOND_H
#define PGFEM_HYPRE_PRECOND_H

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

  /** Diagonal scale preconditioner */
  int PGFEM_HYPRE_ScaleDiagCreate(HYPRE_Solver *vdiag_pc);

  int PGFEM_HYPRE_ScaleDiagSetup(HYPRE_Solver vdiag_pc,
                 HYPRE_ParCSRMatrix vA,
                 HYPRE_ParVector vb,
                 HYPRE_ParVector vx);

  int PGFEM_HYPRE_ScaleDiagSolve(HYPRE_Solver vdiag_pc,
                 HYPRE_ParCSRMatrix vA,
                 HYPRE_ParVector vb,
                 HYPRE_ParVector vx);

  int PGFEM_HYPRE_ScaleDiagDestroy(HYPRE_Solver vdiag_pc);

  /** Jacobi preconditioner */
  int PGFEM_HYPRE_JacobiCreate(HYPRE_Solver *vjacobi_pc);

  int PGFEM_HYPRE_JacobiSetup(HYPRE_Solver vjacobi_pc,
                  HYPRE_ParCSRMatrix vA,
                  HYPRE_ParVector vb,
                  HYPRE_ParVector vx);

  int PGFEM_HYPRE_JacobiSolve(HYPRE_Solver vjacobi_pc,
                  HYPRE_ParCSRMatrix vA,
                  HYPRE_ParVector vb,
                  HYPRE_ParVector vx);

  int PGFEM_HYPRE_JacobiDestroy(HYPRE_Solver vjacobi_pc);

#endif /* #ifndef PGFEM_HYPRE_PRECOND_H */
