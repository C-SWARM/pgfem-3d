#ifndef MATRIX_PRINTING_H
#define MATRIX_PRINTING_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Prints the matrix indicies by columns. Note that this function
      cannot print the renumbered columns until after the matrix has
      been rebuilt. */
  void print_columns(char *name, int ncol, int *Ap, int *Ai, int start);

  /** This function prints the sparse CSR matrix in a format readable
   by MATLAB. Discription of arguments:\n
   name -- filename to printto\n
   ncol -- the number of columns\n
   Ap/Ai -- the matrix indicies in CSR format\n
   start -- the index of the first column of the matrix\n
   The important thing to note is that the I index is actually
   the column.  Thus when importing into MATLAB or the like, the
   second column should be the row variable, not the first. Also,
   since we only care about the pattern, we just set the values of the
   non-zero entries to 1. */
  void print_matrix(char *name, int ncol, int *Ap, int *Ai, int start);

  /** This function is the same as print_matrix except that it prints
      the matrix using the renumber mapping. This function can print
      the renumbered pattern before the matrix has been rebuilt. */
  void print_rn_matrix(char *name, int ncol, int *Ap, int *Ai, int *order, int start);

  /** Prints the matrix to a file in CSR format */
  void print_CSR(char *name, int ncol, int *Ap, int *Ai);

  /** prints a vector to a file */
  void printVector(double *vector, int length,char *filename);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef MATRIX_PRINTING_H */
