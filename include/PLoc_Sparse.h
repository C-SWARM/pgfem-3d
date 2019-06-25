/* HEADER */
#pragma once
#ifndef PLOC_SPARSE_H
#define PLOC_SPARSE_H

#include "data_structure.h"
#include "pgfem3d/Solver.hpp"
#include "pgfem3d/Communication.hpp"

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
                  Ai_t *Ai, /**< UNUSED */
                  int *Ap, /**< UNUSED */
                  long *cnL, /**< UNUSED */
                  long *cnG,
                  long ndofe,
                  long *Ddof, /**< UNUSED */
                  long GDof, /**< UNUSED */
		  int myrank,
		  int nproc,
                  pgfem3d::SparseComm *comm,
                  int interior,
                  pgfem3d::solvers::SparseSystem *system,
                  const int analysis);

#endif /* #ifndef PLOC_SPARSE_H */
