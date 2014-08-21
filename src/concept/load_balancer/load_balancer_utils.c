/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer_utils.h"
#include <math.h>
#include <assert.h>
#include <float.h>

int double_comp(const void *lhs,
		const void *rhs)
{
  double diff = *((double *) lhs) - *((double *) rhs);
  if(diff < 0.0) return -1;
  else if(diff > 0.0) return 1;
  else return 0;
}

int int_comp(const void *lhs,
	     const void *rhs)
{
  return *((int*)lhs) - *((int*)rhs);
}

int size_t_comp(const void *lhs,
		const void *rhs)
{
  size_t a = *(size_t*)lhs;
  size_t b = *(size_t*)rhs;
  if(a<b) return -1;
  else if(a>b) return 1;
  else return 0;
}

double compute_avg(double *restrict arr,
		   const size_t len)
{
  double avg = 0.0;
  for(size_t i=0; i<len; i++){
    avg += arr[i];
  }
  avg /= len;
  return avg;
}

double compute_std(double *restrict arr,
		   const size_t len,
		   const double avg)
{
  assert(len > 0);
  double std = 0;
  for(size_t i=0; i<len; i++){
    std += (arr[i] - avg)*(arr[i]-avg);
  }
  std = sqrt(std/len);
  return std;
}


size_t min_arr_idx(const double *restrict arr,
		   const size_t len)
{
  double min = DBL_MAX;
  size_t idx = -1;
  for(size_t i=0; i<len; i++){
    if(min > arr[i]){
      min = arr[i];
      idx = i;
    }
  }
  return idx;
}

void *lower_bound(const void *val,
		  void *arr,
		  const size_t nmemb,
		  const size_t size,
		  int (*compare)(const void*,const void*))
{
  char *c = arr;
  for(size_t i=0; i<nmemb; i++){
    c += i*size;
    if(compare(c,val) >= 0) return c;
  }
  return NULL;
}
		  
