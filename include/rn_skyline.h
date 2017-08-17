#ifndef RN_SKYLINE_H
#define RN_SKYLINE_H

  /** This function returns the skyline profile (including the
      diagonal) of a matrix renumbered with ParMETIS */
  long rn_skyline(int ncol, int *Ap, int *Ai, int *order, int start);

#endif /* #ifndef RN_SKYLINE_H */
