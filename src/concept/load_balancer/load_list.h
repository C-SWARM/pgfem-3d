/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#pragma once
#ifndef LOAD_LIST_H
#define LOAD_LIST_H

#include <stdlib.h>
#include <stdio.h>

typedef struct LOAD{
  double time;
  size_t client_proc;
  size_t server_proc;
}LOAD;

typedef struct STATS{
  double avg;
  double std;
  double min;
  double max;
}STATS;

/**
 * Compute min, mav, avg, and std. of arr. Variable arr is sorted in
 * ascending order on exit
 */
void stats_compute(STATS *stats,
		   double *arr,
		   const size_t len);

/** Reset stats object to zeros. */
void stats_reset(STATS *stats);

typedef struct LOAD_LIST{
  LOAD *loads;
  size_t n_loads;
  size_t n_servers;
  STATS time_stats;
}LOAD_LIST;

/**
 * Allocate the load list.
 *
 * Allocates internal loads and times buffers and sets all members to
 * 0.
 */
void build_LOAD_LIST(LOAD_LIST *list,
		     const size_t n_loads,
		     const size_t n_servers);

/** Free memory associated with a LOAD LIST object */
void destroy_LOAD_LIST(LOAD_LIST *list);

/** Get pointer to the load at idx */
LOAD* load_list_get_load(LOAD_LIST *list,
			 const size_t idx);

/**
 * Sort loads by their time in descending order.
 */
void load_list_sort_load_time(LOAD_LIST *list);

/**
 * Sort loads by their server_proc.
 */
void load_list_sort_load_server(LOAD_LIST *list);

/**
 * Extract the times from the list of loads.
 *
 * Returns buffer of times with the same internal ordering of loads.
 */
double* load_list_extract_load_times(const LOAD_LIST *list);

/** Print out a nicely formatted listing of the loads */
void load_list_print(FILE *out,
		     const LOAD_LIST *list);

/**
 * Compute time statistics from load.
 *
 * Computes and stores average, standard deviation, min and max of
 * times from the list of loads.
 */
void load_list_compute_stats(LOAD_LIST *list);

/**
 * Reset the loads and avg/std times to 0.
 */
void load_list_reset(LOAD_LIST *list);

#endif
