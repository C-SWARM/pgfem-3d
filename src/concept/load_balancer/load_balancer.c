/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer.h"

int load_balancer(LOAD_LIST *servers,
		  const size_t n_servers)
{
  for(size_t i=0; i<n_servers; i++){
    load_list_sort_load_time(servers + i);
    load_list_compute_stats(servers + i);
  }
  return 0;
}
