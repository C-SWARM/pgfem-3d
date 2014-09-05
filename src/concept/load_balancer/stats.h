/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#pragma once
#ifndef STATS_H
#define STATS_H
#include <stdio.h>

struct STATS{
  double avg;
  double std;
  double min;
  double max;
  double total;
};
#ifndef TYPEDEF_STATS
#define TYPEDEF_STATS
typedef struct STATS STATS;
#endif

/**
 * Compute min, mav, avg, and std. of arr. Variable arr is sorted in
 * ascending order on exit
 */
void STATS_compute(STATS *stats,
		   double *arr,
		   const size_t len);

/** Reset stats object to zeros. */
void STATS_reset(STATS *stats);

/** print statistics info to file */
void STATS_print(FILE *out,
		 const STATS *stats);

#endif
