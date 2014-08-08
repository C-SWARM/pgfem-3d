/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer.h"
#include <math.h>

static double compute_avg(double *restrict arr,
			  const size_t len)
{
  double avg = 0.0;
  for(size_t i=0; i<len; i++){
    avg += arr[i];
  }
  avg /= len;
}

static double compute_std(double *restrict arr,
			  const size_t len,
			  const double avg)
{
  double std = 0;
  for(size_t i=0; i<len; i++){
    std += (arr[i] - avg)*(arr[i]-avg);
  }
  std = sqrt(std/len);
}


/**
 * Compute the average and standard deviation of values in arr of
 * length len. Return in avg and std.
 */
static void compute_avg_std(double *restrict arr,
			    const size_t len,
			    double *restrict avg,
			    double *restrict std)
{
  *avg = 0.0;
  for(size_t i=0; i<len; i++){
    *avg += arr[i];
  }
  *avg /= len;

  *std = 0;
  for(size_t i=0; i<len; i++){
    *std += arr[i];
  }

}

int load_balancer()
{
  return 0;
}
