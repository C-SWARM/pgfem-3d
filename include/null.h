#ifndef NULL_H
#define NULL_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Set a 4d tensor (3^4 elements) equal to zero. */
  void null_4d (double A[3][3][3][3]);

  /** set a 2d tensor (3^2 elements) equal to zero. */
  void null_2d (double A[3][3]);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef NULL_H */
