/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer.h"
#include <stdlib.h>

int load_balancer(LOAD_LIST *list)
{
  load_list_sort_load_time(list);
  load_list_compute_stats(list);
  return 0;
}
