#ifndef PRINT_DIST_H
#define PRINT_DIST_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** This function is intended to print the dist vector to a file to
      be read by another program, but it is general enough to print
      any vector really. */
  void print_dist(char *name,int length, int *dist);

#ifdef __cplusplus
{
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PRINT_DIST_H */
