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
