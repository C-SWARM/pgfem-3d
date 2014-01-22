#ifndef REMOVE_DIAGONAL_H
#define REMOVE_DIAGONAL_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Removes the diagonal elements from the graph to be
      renumbered. */
void remove_diagonal(int ncol, int start, int *Ap,
		     int *Ai, int *Apd, int *Aid);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef REMOVE_DIAGONAL_H */
