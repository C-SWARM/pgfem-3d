/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load.h"
#include "load_balancer_utils.h"
#include <string.h>
#include <math.h>
#include <assert.h>

void LOAD_ID_copy(LOAD_ID *dest,
		  const LOAD_ID *src)
{
  if(src == dest) return;
  memcpy(dest,src,sizeof(*dest));
}

int LOAD_ID_compare(const void *a,
		    const void *b)
{
  const LOAD_ID *restrict A = (const LOAD_ID *)a;
  const LOAD_ID *restrict B = (const LOAD_ID *)b;
  int result = size_t_comp((const void*) &(A->proc),
			   (const void*) &(B->proc));
  if(!result){
    result = size_t_comp((const void*) &(A->elem),
			 (const void*) &(B->elem));
  }
  return result;
}

void LOAD_ID_print(FILE *out,
		   const LOAD_ID *id)
{
  fprintf(out,"id: %3ld::%ld",id->proc,id->elem);
}

/* LOAD FUNCTIONS */
void LOAD_copy(LOAD *dest,
	       const LOAD *src)
{
  if(src != dest)
    memcpy(dest,src,sizeof(*dest));
}

void LOAD_swap(LOAD *A,
	       LOAD *B)
{
  LOAD buffer;
  LOAD_copy(&buffer,A);
  LOAD_copy(A,B);
  LOAD_copy(B,&buffer);
}

int LOAD_compare_id(const void *a,
		    const void *b)
{
  return LOAD_ID_compare((void*) &(((const LOAD*)a)->id),
			 (void*) &(((const LOAD*)b)->id));
}

inline int LOAD_compare_part_id(const void *a,
				const void *b)
{
  return size_t_comp((void*) &(((const LOAD*)a)->part_id),
		     (void*) &(((const LOAD*)b)->part_id));
}

int LOAD_compare_part_id_load(const void *a,
			      const void *b)
{
  int part_id = LOAD_compare_part_id(a,b);
  if(part_id == 0){
    return LOAD_compare_load(a,b);
  } else {
    return part_id;
  }
}

inline int LOAD_compare_load(const void *a,
			     const void *b)
{
  return double_comp((void*) &(((LOAD *) a)->load),
		     (void*) &(((LOAD *) b)->load));
}

int LOAD_compare_r_load(const void *a,
			const void *b)
{
  return LOAD_compare_load(b,a);
}

int LOAD_approx_equal(const LOAD *a,
		      const LOAD *b,
		      const double tol)
{
  assert(tol >= 0.0);
  return (fabs(a->load - b->load) >= tol);
}

void LOAD_print(FILE *out,
		const LOAD *load)
{
  fprintf(out,"%11.3e p: %3ld ",load->load,load->part_id);
  LOAD_ID_print(out,&(load->id));
}
