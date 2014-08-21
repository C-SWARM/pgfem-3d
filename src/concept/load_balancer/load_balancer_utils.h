/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef LOAD_BALANCER_UTILS_H
#define LOAD_BALANCER_UTILS_H

#include <stdlib.h>

/** compare two doubles (as in for qsort) */
int double_comp(const void *lhs,
		const void *rhs);

/** compare two integers (as in for qsort) */
int int_comp(const void *lhs,
	     const void *rhs);

/** compare two *unsigned* integers (as in for qsort) */
int size_t_comp(const void *lhs,
		const void *rhs);

/** compute average value of arr */
double compute_avg(double *restrict arr,
		   const size_t len);

/** compute standard deviation of arr given average of arr */
double compute_std(double *restrict arr,
		   const size_t len,
		   const double avg);

/** return the index of the minimum value */
size_t min_arr_idx(const double *restrict arr,
		   const size_t len);

/**
 * Get pointer to the first entry in arr that compares >= *val using
 * the compare function. Returns NULL if a lower bound does not exist.
 */
void *lower_bound(const void *val,
		  void *arr,
		  const size_t nmemb,
		  const size_t size,
		  int (*compare)(const void*,const void*));
#endif
