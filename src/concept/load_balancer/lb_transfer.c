/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "lb_transfer.h"
#include "partition.h"
#include <string.h>
#include <assert.h>

inline void TRANSFER_copy(TRANSFER *dest,
			  const TRANSFER *src)
{
  if(src == dest) return;
  memcpy(dest,src,sizeof(*dest));
}

inline void TRANSFER_print(FILE *out,
			   const TRANSFER *T)
{
  fprintf(out,"%3ld --> %3ld \t",T->from,T->to);
  LOAD_ID_print(out,&(T->id));
  fprintf(out,"\n");
}

/* === TRANSFER_LIST === */
void TRANSFER_LIST_build(TRANSFER_LIST *TL,
			 const size_t max_size)
{
  TL->transfers = calloc(max_size,sizeof(*(TL->transfers)));
  TL->max_size = max_size;
  TL->size = 0;
}

void TRANSFER_LIST_destroy(TRANSFER_LIST *TL)
{
  free(TL->transfers); TL->transfers = NULL;
  TL->size = TL->max_size = 0;
}

static int private_TRANSFER_LIST_search_partition(const size_t part_id,
						  const PARTITION *P,
						  TRANSFER_LIST *TL)
{
  const LOAD *restrict load = P->loads; /* alias */
  TRANSFER T;
  T.to = part_id;
  for(size_t i=0, n_load=P->size; i<n_load; i++){
    if(load[i].part_id != part_id){
      T.from = load[i].part_id;
      LOAD_ID_copy(&(T.id),&(load[i].id));
      if(TRANSFER_LIST_push(TL,&T) != 0) return 1;
    }
  }
  return 0;
}

void TRANSFER_LIST_compute(const PARTITION_LIST *PL,
			   TRANSFER_LIST *TL)
{
  /* destroy TL and re-build for worst-case scenario -> all loads need
     to be transferred. */
  TRANSFER_LIST_destroy(TL);
  const size_t total_loads = PARTITION_LIST_total_size(PL);
  TRANSFER_LIST_build(TL,total_loads);
  for(size_t i=0, n_part = PL->n_parts; i<n_part; i++){
    int err = private_TRANSFER_LIST_search_partition(i,(PL->partitions) + i,TL);
    assert(err == 0);
  }
}

int TRANSFER_LIST_push(TRANSFER_LIST *TL,
		       const TRANSFER *T)
{
  if(TRANSFER_LIST_is_full(TL)) return 1;

  TRANSFER_copy((TL->transfers) + TL->size,T);
  TL->size ++;
  return 0;
}

int TRANSFER_LIST_is_full(const TRANSFER_LIST *TL)
{
  return (TL->size >= TL->max_size);
}

int TRANSFER_LIST_is_empty(const TRANSFER_LIST *TL)
{
  return (TL->size == 0);
}

void TRANSFER_LIST_print(FILE *out,
			 const TRANSFER_LIST *TL)
{
  fprintf(out,"No. Transfers: %ld\n",TL->size);
  for(size_t i=0, nt=TL->size; i<nt; i++){
    TRANSFER_print(out,TL->transfers + i);
  }
  fprintf(out,"\n");
}
