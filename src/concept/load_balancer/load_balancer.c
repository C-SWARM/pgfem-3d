/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer.h"
#include "partition.h"
#include "load.h"
#include "load_balancer_utils.h"
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>

static void greedy_add_to_partition(PARTITION *S,
				    PARTITION *P,
				    double *sum)
{
  const LOAD *ld = PARTITION_top_load(S);
  if( ! PARTITION_is_full(P) ){
    PARTITION_push_load(P, ld);
    *sum += ld->load;
    PARTITION_pop_load(S); /* invalidates ld */
  } else {
    *sum = DBL_MAX;
  }
}

int load_balancer_greedy(PARTITION_LIST *PL)
{
  int err = 0;
  const size_t n_parts = PL->n_parts;
  const size_t set_size = PARTITION_LIST_total_size(PL);
  PARTITION *parts = PL->partitions; /* alias */
  assert(set_size >= n_parts);

  /* return early if only one (or none) partitions */
  if(n_parts <= 1) return err;

  /* Allocate the set */
  PARTITION *S = calloc(1,sizeof(*S));
  PARTITION_build(S,set_size);

  /* Fill the set from the partitions */
  for(size_t i=0; i<n_parts; i++){
    /* append partitions to the set */
    PARTITION_push_partition(S,parts + i);

    /* reset the partitions to empty */
    parts[i].size = 0;
  }

  /* sort the set by load */
  PARTITION_sort_load(S);

  /* set initial partitions */
  double *sums = calloc(n_parts,sizeof(*sums));
  for(size_t i=0; i<n_parts; i++){
    greedy_add_to_partition(S,parts + i,sums + i);    
  }

  /* proceed with greedy partitioning */
  while ( ! PARTITION_is_empty(S) ){
    size_t idx = min_arr_idx(sums,n_parts);
    greedy_add_to_partition(S, parts + idx, sums + idx);
  }

  PARTITION_destroy(S);
  free(S);
  free(sums);
  return err;
}


/**
 * Search list of loads in A previously belonging to B (at time n) for
 * equivalent loads in B previously belonging to A (at time n). If
 * found, swap jobs.
 */
static size_t reduce_comm_swap(PARTITION *restrict A,
			       PARTITION *restrict B,
			       const size_t A_id,
			       const size_t B_id,
			       const double tol)
{
  /* get pointers to min/max loads */
  LOAD *AL = lower_bound(&B_id,A->loads,A->size,sizeof(*AL),LOAD_compare_part_id);
  LOAD *AL_end = lower_bound(&(B_id+1),A->loads,A->size,sizeof(*AL),LOAD_compare_part_id);
  if(AL_end == NULL) AL_end = A->loads + A->size;

  LOAD *BL = lower_bound(&A_id,B->loads,B->size,sizeof(*BL),LOAD_compare_part_id);
  LOAD *BL_end = lower_bound(&(B_id+1),B->loads,B->size,sizeof(*BL),LOAD_compare_part_id);
  if(BL_end == NULL) BL_end = B->loads + B->size;

  for(; AL != AL_end; AL++){
    for(; BL != BL_end; BL++){

    }
    if(BL == BL_end) break;
  }
}			       

int load_balancer_reduce_comm(PARTITION_LIST *PL,
			      const double tol)
{
  int err = 0;
  const size_t n_part = PL->n_parts;

  /* exit early if no comm required */
  if(n_part <= 1) return err;

  /* alias to partitions */
  PARTITION *parts = PL->partitions;

  /* sort each PARTITION's loads by part_id (from previous step) */
  for(size_t i=0; i<n_part; i++){
    PARTITION_sort_load_part_id(parts + i);
  }

  return err;
}
