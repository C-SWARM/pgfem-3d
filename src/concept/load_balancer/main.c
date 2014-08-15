/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer.h"

#include <stdlib.h>

static const size_t n_servers = 3;
int main(int argc, char **argv)
{
  LOAD_LIST server[n_servers];
  for(size_t i=0; i<n_servers; i++){
    const size_t n_loads = 10;
    build_LOAD_LIST(server+i,n_loads);
  }

  /* set the loads */
  srand(1);
  for(size_t s = 0; s<n_servers; s++){
    for(size_t i=0, nl = server[s].n_loads; i<nl; i++){
      LOAD *restrict load = load_list_get_load(server+s,i);
      load->time = (1.5 + s)*i*i;
      load->client_proc = 0;
      load->server_proc = s;
    }
  }

  load_balancer(server,n_servers);
  for(size_t i=0; i<n_servers; i++){
    load_list_print(stdout,server +i);
    printf("\n");
  }


  /* destroy server loads */
  for(size_t i=0; i<n_servers; i++){
    destroy_LOAD_LIST(server+i);
  }
  return 0;
}
