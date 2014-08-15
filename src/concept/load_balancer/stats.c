/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "stats.h"
#include "load_balancer_utils.h"
#include <stdlib.h>
#include <string.h>

/* STATS FUNCTIONS */
void STATS_compute(STATS *stats,
		   double *arr,
		   const size_t len)
{
  qsort(arr,len,sizeof(*arr),double_comp);
  stats->avg = compute_avg(arr,len);
  stats->total = stats->avg*len;
  stats->std = compute_std(arr,len,stats->avg);
  stats->min = arr[0];
  stats->max = arr[len-1];
}

void STATS_reset(STATS *stats)
{
  memset(stats,0,sizeof(*stats));
}

void STATS_print(FILE *out,
		 const STATS *stats)
{
  printf("Minimum time:   %f\n",stats->min);
  printf("Maximum time:   %f\n",stats->max);
  printf("Average time:   %f\n",stats->avg);
  printf("Std. Dev. time: %f\n",stats->std);
  printf("Total time:     %f\n",stats->total);
}
