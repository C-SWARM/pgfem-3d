#ifndef SKYLINE_H
#define SKYLINE_H

#include "datatype.hpp"

  /** This function calculates and returns the symmetric skyline
      storage requirement (including the diagonal) of sparse
      matrix. For the total skyline of a distributed matrix, gather
      the result of this function. Description of arguments:\n
      ncol -- the number of columns on the domain\n
      Ap/Ai -- the sparse matrix IMPORTANT: The indicies (values)
      of Ai MUST be in numerical order!\n
      start -- the value of the first index on the domain */
  long skyline(long ncol, int *Ap, Ai_t *Ai, long start);

#endif /* #ifndef SKYLINE_H */
