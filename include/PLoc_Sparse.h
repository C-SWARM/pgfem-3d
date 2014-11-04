/* HEADER */
#pragma once
#ifndef PLOC_SPARSE_H
#define PLOC_SPARSE_H

#include "pgfem_comm.h"
#include "hypre_global.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Sparse nonsymmetric row storage format.
   *
   * A->rows[i]->nz = (FLOAT *) calloc (A->rows[i]->length,sizeof(FLOAT));
   *
   *       (j)
   *      xxxxxx
   *  (i) xxxxxx
   *      xxxxxx
   */

  void PLoc_Sparse (double **Lk,
		    double *lk,
		    int *Ai, /**< UNUSED */
		    int *Ap, /**< UNUSED */
		    long *cnL, /**< UNUSED */
		    long *cnG,
		    long ndofe,
		    int *Ddof, /**< UNUSED */
		    long GDof, /**< UNUSED */
		    int myrank,
		    int nproc,
		    COMMUN comm,
		    int interior,
		    PGFEM_HYPRE_solve_info *PGFEM_hypre,
		    const int analysis);

/** PLoc_Sparse_rec for general (non-square) matrices. */
void PLoc_Sparse_rec (double **Lk,
		      double *lk,
		      int *Ai,
		      int *Ap,
		      long *cnG_row,
		      long *cnG_col,
		      long nrow,
		      long ncol,
		      int *Ddof,
		      long GDof,
		      int myrank,
		      int nproc,
		      COMMUN comm,
		      int interior,
		      PGFEM_HYPRE_solve_info *PGFEM_hypre);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PLOC_SPARSE_H */
