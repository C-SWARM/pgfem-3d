#ifndef SWITCH_DOMAINS_H
#define SWITCH_DOMAINS_H

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Creates a new communicator with ranks determined to minimize communication. */
  void switch_domains(int ncol, int *l_order, int *dist, MPI_Comm Old, MPI_Comm *New);

#ifdef __cplusplus
{
#endif /* #ifdef __cplusplus */

#endif /* #ifndef SWITCH_DOMAINS_H */
