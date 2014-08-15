/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "partition.h"
#include "load.h"
#include "load_balancer_utils.h"
#include <stdlib.h>
#include <string.h>

void PARTITION_build(PARTITION *P,
		     const size_t max_size)
{
  P->max_size = max_size;
  P->loads = calloc(max_size,sizeof(*(P->loads)));
  P->size = 0;
}

void PARTITION_destroy(PARTITION *P)
{
  free(P->loads); P->loads = NULL;
  P->max_size = P->size = 0;
}

void PARTITION_copy(PARTITION *dest,
		    const PARTITION *src)
{
  if(src == dest) return;
  PARTITION_destroy(dest);
  PARTITION_build(dest,src->max_size);
  dest->size = src->size;
  /* only copy valid data [0,size) */
  memcpy(dest->loads,src->loads,
	 (dest->size)*sizeof(*(dest->loads)));
}

void PARTITION_set_load_part_id(PARTITION *P,
				const size_t idx)
{
  for(size_t i=0, end = P->size; i<end; i++){
    P->loads[i].part_id = idx;
  }
}

int PARTITION_push_partition(PARTITION *A,
			     const PARTITION *B)
{
  if(A->size + B->size <= A->max_size){
    memcpy(A->loads + A->size,B->loads,(B->size)*sizeof(*(A->loads)));
    A->size += B->size;
    return 0;
  } else return 1;
}

int PARTITION_push_load(PARTITION *P,
			const LOAD *load)
{
  if(P->size >= P->max_size) return -1;
  LOAD_copy(P->loads + P->size,load);
  P->size++;
  return 0;
}


int PARTITION_pop_load(PARTITION *P)
{
  if(P->size == 0) return -1;
  P->size--;
  return 0;
}

const LOAD *PARTITION_top_load(const PARTITION *P)
{
  if(P->size == 0) return NULL;
  else return (P->loads + P->size - 1);
}

int PARTITION_remove_load(PARTITION *P,
			  const size_t idx)
{
  if(idx >= P->size) return -1;

  /* only move the currently valid data, (P->size - idx -1) */
  memmove(P->loads + idx,P->loads + idx + 1,
	  (P->size - idx -1)*sizeof(*(P->loads)));
  P->size--;
  return 0;
}

int PARTITION_insert_load(PARTITION *P,
			  const LOAD *load,
			  const size_t idx)
{
  if(P->size >= P->max_size) return -1;
  if(idx >= P->size) return PARTITION_push_load(P,load);

  memmove(P->loads + idx + 1, P->loads + idx,
	  (P->size - idx)*sizeof(*(P->loads)));
  LOAD_copy(P->loads + idx,load);
  P->size++;
  return 0;
}

int PARTITION_is_full(const PARTITION *P)
{
  if(P->size >= P->max_size) return 1;
  else return 0;
}

int PARTITION_is_empty(const PARTITION *P)
{
  if(P->size <= 0) return 1;
  else return 0;
}

void PARTITION_sort_load_id(PARTITION *P)
{
  qsort(P->loads,P->size,sizeof(*(P->loads)),LOAD_compare_id);
}

void PARTITION_sort_load(PARTITION *P)
{
  qsort(P->loads,P->size,sizeof(*(P->loads)),LOAD_compare_load);
}

static int private_PARTITION_search_intersection(const PARTITION *A,
						  const PARTITION *B,
						  const size_t max_int_size,
						  size_t *A_idx,
						  size_t *B_idx,
						  size_t *size)
{
  for(size_t i=0, end = A->size; i<end; i++){
    LOAD *match = bsearch(A->loads + i,B->loads,B->size,
			  sizeof(*(B->loads)),LOAD_compare_id);
    if(match != NULL){
      if(*size < max_int_size){
	A_idx[*size] = i;
	B_idx[*size] = (size_t) (match - B->loads);
	(*size)++;
      } else return 1; /* ERROR: too many elements in intersection */
    }
  }
  return 0;
}

int PARTITION_compute_intersection(const PARTITION *A,
				   const PARTITION *B,
				   size_t **A_idx,
				   size_t **B_idx,
				   size_t *size)
{
  int err = 0;

  /* allocate intersections to minimum partition size */
  const size_t max_size = (A->size <= B->size) ? A->size : B->size;
  *A_idx = calloc(max_size,sizeof(**A_idx));
  *B_idx = calloc(max_size,sizeof(**B_idx));

  /* initialize intersection to empty */
  *size = 0;

  /* search for intersection from smallest partition */
  if(A->size <= B->size){
   err = private_PARTITION_search_intersection(A,B,max_size,*A_idx,*B_idx,size);
  } else {
   err = private_PARTITION_search_intersection(B,A,max_size,*B_idx,*A_idx,size);
  }

  return err;
}

double *PARTITION_extract_loads(const PARTITION *P)
{
  const LOAD *restrict loads = P->loads;
  double *restrict result = calloc(P->size,sizeof(*result));
  for(size_t i=0, end = P->size; i<end; i++){
    result[i] = loads[i].load;
  }
  return result;
}

void PARTITION_compute_stats(PARTITION *P)
{
  STATS_reset(&P->stats);
  double *load = PARTITION_extract_loads(P);
  STATS_compute(&P->stats,load,P->size);
  free(load);
}

inline double PARTITION_stats_total(const PARTITION *P)
{
  return P->stats.total;
}

void PARTITION_print(FILE *out,
		     const PARTITION *P)
{
  fprintf(out,"size: %ld max_size: %ld\n",P->size,P->max_size);
  for(size_t i=0, end = P->size; i<end; i++){
    LOAD_print(out,P->loads + i);
    fprintf(out,"\n");
  }
  STATS_print(out,&P->stats);
  fprintf(out,"\n\n");
}

/* === PARTITION_LIST === */

void PARTITION_LIST_build(PARTITION_LIST *PL,
			  const size_t n_parts,
			  const size_t part_max_size)
{
  PL->partitions = calloc(n_parts,sizeof(*(PL->partitions)));
  PL->n_parts = n_parts;
  for(size_t i=0; i<n_parts; i++){
    PARTITION_build(PL->partitions + i,part_max_size);
  }
}

void PARTITION_LIST_destroy(PARTITION_LIST *PL)
{
  for(size_t i=0, n_parts = PL->n_parts; i<n_parts; i++){
    PARTITION_destroy(PL->partitions + i);
  }
  free(PL->partitions);
  PL->partitions = NULL;
  PL->n_parts = 0;
}

void PARTITION_LIST_print(FILE *out,
			  const PARTITION_LIST *PL)
{
  size_t  n_parts = PL->n_parts; 
  fprintf(out,"No. Partitions: %ld\n",n_parts);
  PARTITION *parts = PL->partitions; /* alias */
  double *restrict totals = calloc(n_parts,sizeof(*totals));
  for(size_t i=0; i<n_parts; i++){
    PARTITION_compute_stats(parts + i);
    PARTITION_print(out,parts + i);
    totals[i] = PARTITION_stats_total(parts + i);
  }
  qsort(totals,n_parts,sizeof(*totals),double_comp);
  fprintf(out,"Max. Diff.: %11.3e\n",totals[n_parts-1] - totals[0]);
  fprintf(out,"==============================\n\n");
  free(totals);
}

size_t PARTITION_LIST_total_size(const PARTITION_LIST *PL)
{
  size_t result = 0;
  const PARTITION *restrict parts = PL->partitions;
  for(size_t i=0, n_parts = PL->n_parts; i<n_parts; i++){
    result += parts[i].size;
  }
  return result;
}

void PARTITION_LIST_copy(PARTITION_LIST *dest,
			 const PARTITION_LIST *src)
{
  if(src == dest) return;
  PARTITION_LIST_destroy(dest);
  dest->n_parts = src->n_parts;
  dest->partitions = calloc(dest->n_parts,sizeof(*(dest->partitions)));
  for(size_t i=0, end = dest->n_parts; i<end; i++){
    dest->partitions[i].loads = NULL; /* allows free without error */
    PARTITION_copy(dest->partitions + i, src->partitions + i);
  }
}
