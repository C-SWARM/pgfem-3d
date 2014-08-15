/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load.h"
#include "load_balancer_utils.h"
#include <string.h>

int LOAD_ID_compare(const void *a,
		    const void *b)
{
  const LOAD_ID *restrict A = (const LOAD_ID *)a;
  const LOAD_ID *restrict B = (const LOAD_ID *)b;
  int result = size_t_comp((void*) &(A->proc),(void*) &(B->proc));
  if(!result){
    result = size_t_comp((void*) &(A->elem),(void*) &(B->elem));
  }
  return result;
}

void LOAD_ID_print(FILE *out,
		   const LOAD_ID *id)
{
  fprintf(out,"id: %ld::%ld",id->proc,id->elem);
}

/* LOAD FUNCTIONS */
void LOAD_copy(LOAD *dest,
	       const LOAD *src)
{
  if(src != dest)
    memcpy(dest,src,sizeof(*dest));
}

int LOAD_compare_id(const void *a,
		    const void *b)
{
  return LOAD_ID_compare((void*) &(((const LOAD*)a)->id),
			 (void*) &(((const LOAD*)b)->id));
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

void LOAD_print(FILE *out,
		const LOAD *load)
{
  fprintf(out,"%11.3e ",load->load);
  LOAD_ID_print(out,&load->id);
}
