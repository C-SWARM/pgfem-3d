/* HEADER */

/**
 * AUTHORS:
 * Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef PGFEM_PAR_MATVEC_H
#define PGFEM_PAR_MATVEC_H

#include "PGFEM_mpi.h"

  typedef void* PGFEM_par_matrix_comm;

  /** Parallel dense matrix object */
  typedef struct PGFEM_par_matrix{
    MPI_Comm mpi_comm;
    int assembled;
    int add_values;
    int n_rows;
    int n_cols;
    int *n_own_rows; /* [nproc] */
    int *idx_starts; /* [nproc + 1] */
    double *data; /* n_loc_row[myrank] x n_col */
    void *send_info;
    void *recv_info;
    void *s_rows;
    void *r_rows;
  } PGFEM_par_matrix;

  /** allocate a PGFEM_par_matrix object and initialize the
      communication structure. Contains collective (blocking)
      communication. */
  int initialize_PGFEM_par_matrix(const int n_rows,
                  const int n_cols,
                  const int n_own_rows,
                  const int n_entries,
                  const int *row_idx,
                  const int *col_idx,
                  MPI_Comm mpi_comm,
                  PGFEM_par_matrix **mat);

  /** Destroy the matrix object. No communication */
  void destroy_PGFEM_par_matrix(PGFEM_par_matrix *mat);

  /** Set values in matrix. Note that subseqent calls to *set_values
      and *assemble may overwrite previous information. Overwritting
      data at assembly can be averted by calling *add_to_values
      first. Duplicate entries will be overwritten with last in
      row-sorted list. No communication */
  int PGFEM_par_matrix_set_values(const int n_entries,
                  const int *row_idx,
                  const int *col_idx,
                  const double *values,
                  PGFEM_par_matrix *mat);

  /** Add to values in matrix. Assembly will add values from other
      processors. No communication*/
  int PGFEM_par_matrix_add_to_values(const int n_entries,
                     const int *row_idx,
                     const int *col_idx,
                     const double *values,
                     PGFEM_par_matrix *mat);

  /** Zero values in the matrix. No communication */
  int PGFEM_par_matrix_zero_values(PGFEM_par_matrix *mat);

  /** Starts a non-blocking communication sequence to assemble the
      matrix. NOTE: the user must not modify the matrix after calling
      this function until after calling *end_assembly */
  int PGFEM_par_matrix_start_assembly(PGFEM_par_matrix *mat,
                      PGFEM_par_matrix_comm *comm);

  /** Finishes a non-blocking communication sequence to assemble the
      matrix */
  int PGFEM_par_matrix_end_assembly(PGFEM_par_matrix *mat,
                    PGFEM_par_matrix_comm comm);

  /** Get the local part of the column. No communication */
  int PGFEM_par_matrix_get_column(const PGFEM_par_matrix *mat,
                  const int col_idx,
                  double *col);

  /** Compute the dot product of 2 distributed vectors. The number of
      elements in each vector must match on the domain. The result is
      distributed to all members in the communicator. Collective
      blocking communication */
  int PGFEM_par_vec_dot(const int len,
            const double *vec_a,
            const double *vec_b,
            MPI_Comm comm,
            double *result);

#endif /* #ifndef  */
