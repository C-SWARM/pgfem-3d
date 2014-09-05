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
#include <math.h>
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
static size_t reduce_comm_swap(PARTITION *A,
			       PARTITION *B,
			       const size_t A_id,
			       const size_t B_id,
			       const double tol)
{
  size_t err = 0;
  /* sort partition by owner then load */
  PARTITION_sort_load_part_id(A);
  PARTITION_sort_load_part_id(B);

  LOAD match;
  match.load = 0.0;
  match.part_id = B_id;
  const size_t a_lower = lower_bound_idx(&match,A->loads,A->size,sizeof(*A),LOAD_compare_part_id);
  match.part_id ++;
  const size_t a_upper = lower_bound_idx(&match,A->loads,A->size,sizeof(*A),LOAD_compare_part_id);
  /* return early if none in range */
  if(a_lower == a_upper || a_lower == A->size) return err;

  match.part_id = A_id;
  const size_t b_lower = lower_bound_idx(&match,B->loads,B->size,sizeof(*B),LOAD_compare_part_id);
  match.part_id ++;
  const size_t b_upper = lower_bound_idx(&match,B->loads,B->size,sizeof(*B),LOAD_compare_part_id);
  /* return if none in range */
  if(b_lower == b_upper || b_lower == B->size) return err;

  /* aliases */
  LOAD *A_load = A->loads;
  LOAD *B_load = B->loads;

  /* skip array */
  size_t b_len = b_upper - b_lower;
  size_t *restrict skip_arr = calloc(b_len,sizeof(*skip_arr));
  size_t n_skip = 0;

  /* loop through range in A */
  for(size_t i = a_lower; i < a_upper; i++){
    /* break if no possible matches */
    if(n_skip >= b_len) break;
    int swap_idx = -1;
    double min = tol;
    /* loop through range in B */
    for(size_t j = b_lower; j < b_upper; j++){
      /* search for minimum index */
      if(!skip_arr[j - b_lower]){
	double diff = fabs(A_load[i].load - B_load[j].load);
	if(diff < min){
	  min = diff;
	  swap_idx = j;
	}
      }
    }

    if(swap_idx >= 0){
      skip_arr[swap_idx - b_lower] ++; /* mark index to skip */
      n_skip ++;
      LOAD_swap(A_load + i, B_load + swap_idx);
    }

  }

  free(skip_arr);
  return err;
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
  for(size_t i=0; i<n_part; i++){
    for(size_t j=0; j<n_part; j++){
      if(i == j) continue;
      reduce_comm_swap(parts + i,parts + j,i,j,tol);
    }
  }
  return err;
}
