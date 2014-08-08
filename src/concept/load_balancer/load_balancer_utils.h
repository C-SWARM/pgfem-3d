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

#endif
