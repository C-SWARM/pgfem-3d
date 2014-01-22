#ifndef COMM_PROFILER_H
#define COMM_PROFILER_H

#include "sort_container.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Writes a report to file which outlines the amount of
      communication required. CAUTION: this function re-sorts the sort
      container to its original state. */
  void comm_profiler(int ndom, int *dist, int *new_ranks, sort_container *info);

#ifdef __cplusplus
{
#endif /* #ifdef __cplusplus */

#endif /* #ifndef COMM_PROFILER_H */
