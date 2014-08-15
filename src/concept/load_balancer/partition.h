/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#pragma once
#ifndef PARTITION_H
#define PARTITION_H

#include "load.h"
#include "stats.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct PARTITION{
  LOAD *loads;
  STATS stats;
  size_t size;
  size_t max_size;
} PARTITION;

void PARTITION_build(PARTITION *P,
		     const size_t max_size);

void PARTITION_destroy(PARTITION *P);

void PARTITION_copy(PARTITION *dest,
		    const PARTITION *src);

/**
 * Append the partition B onto A, i.e., A += B.
 *
 * Returns non-zero if A is full and operation is not performed.
 */
int PARTITION_push_partition(PARTITION *A,
			     const PARTITION *B);

/**
 * Append a LOAD to the PARTITION object.
 *
 * Returns 0 on success. Returns non-zero if PARTITION is full and
 * load is not appended.
 */
int PARTITION_push_load(PARTITION *P,
			const LOAD *load);

/**
 * Remove the last LOAD in the PARTITION object.
 *
 * Returns 0 on success. Returns non-zero if PARTITION is empty.
 */
int PARTITION_pop_load(PARTITION *P);

/**
 * Give _const_ access to last LOAD in the PARTITION object.
 *
 * Returns NULL if PARTITION is empty.
 */
const LOAD *PARTITION_top_load(const PARTITION *P);

/**
 * Remove the LOAD at idx in PARTITION.
 *
 * Behaves as if the element was removed from the list, but does not
 * free memory. Returns 0 on success. Returns non-zero if idx is
 * invalid.
 */
int PARTITION_remove_load(PARTITION *P,
			  const size_t idx);

/**
 * Insert LOAD at idx in PARTITION.
 *
 * Behaves as if the element was inserted before the element currently
 * at idx. Returns 0 on success. Returns non-zero (and performs no
 * action) if insertion cases size > max_size.
 */
int PARTITION_insert_load(PARTITION *P,
			  const LOAD *load,
			  const size_t idx);

/**
 * T/F (1/0) if PARTITION object is full.
 */
int PARTITION_is_full(const PARTITION *P);

/**
 * T/F (1/0) if PARTITION object is empty.
 */
int PARTITION_is_empty(const PARTITION *P);

/**
 * Sort PARTITION loads by LOAD id.
 */
void PARTITION_sort_load_id(PARTITION *P);

/**
 * Sort PARTITION loads by size/weight.
 */ 
void PARTITION_sort_load(PARTITION *P);

/**
 * Compute the indices of intersection between two partitions based on
 * LOAD_ID.
 *
 * Requires that A and B are sorted by PARTITION_sort_load_id before
 * calling this function. Returns two allocated arrays containing the
 * corresponding indices of elements A and B in the intersection well
 * as the size of the intersection. Returns non-zero on error.
 */
int PARTITION_compute_intersection(const PARTITION *A,
				   const PARTITION *B,
				   size_t **A_idx,
				   size_t **B_idx,
				   size_t *size);

/**
 * Extract the list of loads from the partition.
 */
double *PARTITION_extract_loads(const PARTITION *P);

/**
 * Compute and store statistics information for the PARTITION object.
 */
void PARTITION_compute_stats(PARTITION *P);

/**
 * Retutrn the total load/time on the PARTITION object.
 *
 * The user must ensure that PARTITION_compute_stats has been called
 * on the current object for valid return value.
 */
double PARTITION_stats_total(const PARTITION *P);

/**
 * Print a PARTITION object to a file.
 */
void PARTITION_print(FILE *out,
		     const PARTITION *P);

/** Container of PARTITION objects */
typedef struct PARTITION_LIST{
  PARTITION *partitions;
  size_t n_parts;
  size_t n_loads;
} PARTITION_LIST;

void PARTITION_LIST_build(PARTITION_LIST *PL,
			  const size_t n_parts,
			  const size_t part_max_size);

void PARTITION_LIST_destroy(PARTITION_LIST *PL);

void PARTITION_LIST_print(FILE *out,
			  const PARTITION_LIST *PL);

size_t PARTITION_LIST_total_size(const PARTITION_LIST *PL);

void PARTITION_LIST_copy(PARTITION_LIST *dest,
			 const PARTITION_LIST *src);
#endif
