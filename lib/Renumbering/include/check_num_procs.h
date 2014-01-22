#ifndef CHECK_NUM_PROCS_H
#define CHECK_NUM_PROCS_H

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Determines if the number of processors renumbering is to be run
      on is a power of 2. */
void check_num_procs(MPI_Comm comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef CHECK_NUM_PROCS_H */
