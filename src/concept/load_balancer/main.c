/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer.h"
#include "load_list.h"

#include <stdlib.h>

int main(int argc, char **argv)
{
  LOAD_LIST list;
  {
    const size_t n_loads = 10;
    const size_t n_servers = 3;
    build_LOAD_LIST(&list,n_loads,n_servers);
  }

  /* set the loads */
  srand(1);
  for(size_t i=0; i<list.n_loads; i++){
    LOAD *restrict load = load_list_get_load(&list,i);
    load->time = 1.5*i*i;
    load->client_proc = 0;
    load->server_proc = rand()%2;
  }

  load_balancer(&list);
  load_list_print(stdout,&list);
  return 0;
}
