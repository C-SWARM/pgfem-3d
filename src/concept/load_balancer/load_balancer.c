/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_balancer.h"
#include "load_balancer_utils.h"
#include <string.h>
#include <stdio.h>
#include <float.h>

static int get_loads_from_servers(LOAD_LIST *global,
				  const LOAD_LIST *servers,
				  const size_t n_servers)
{
  size_t n_loads = 0;
  for(size_t i=0; i<n_servers; i++){
    n_loads += servers[i].n_loads;
  }
  build_LOAD_LIST(global,n_loads);

  size_t idx = 0;
  for(size_t i=0; i<n_servers; i++){
    n_loads = servers[i].n_loads;
    memcpy(global->loads+idx,servers[i].loads,
	   n_loads*sizeof(*(global->loads)));
    idx += n_loads;
  }
  return 0;
}

/** reassign loads by round robin assuming decreasing order in global
    list. */
static void round_robin_reassignment(const LOAD_LIST *global,
				     const size_t n_servers,
				     LOAD **dest,
				     LOAD **ends)
{
  for(size_t i=0, end = global->n_loads; i<end; i++){
    for(size_t j = 0; j<n_servers; j++){
      if(dest[j] != ends[j]){
	copy_LOAD(dest[j],global->loads+i);
	dest[j]++; /* increment pointer */
	i++; /* increment global counter */
      }
    }
  }
}


/** greedy reassignment. THis is not working quite properly... */
static void greedy_reassignment(const LOAD_LIST *global,
				const size_t n_servers,
				const size_t K_greedy,
				LOAD **dest,
				LOAD **ends)
{
  double *sum = calloc(n_servers,sizeof(*sum));
  const size_t n_loads = global->n_loads;

  if(K_greedy == 0){
    round_robin_reassignment(global,n_servers,dest,ends);
  } else {
    if(n_loads < n_servers*K_greedy) abort();

    size_t load_idx = 0;
    for(size_t i=0; i<n_servers; i++){
      for(size_t j=0; j<K_greedy; j++){
	copy_LOAD(dest[i],global->loads + load_idx);
	sum[i] += dest[i]->time;
	dest[i]++;
	load_idx++;
      }
    }

    for(; load_idx < n_loads; load_idx++){
      size_t server_idx = min_arr(sum,n_servers);
      if(dest[server_idx] != ends[server_idx]){
	copy_LOAD(dest[server_idx],global->loads + load_idx);
	sum[server_idx] += dest[server_idx]->time;
	dest[server_idx]++;
      } else {
	sum[server_idx] = DBL_MAX;
      }
    }

    free(sum);
  }
}

/* From G. Horton, "A multi-level diffusion method for dynamic load
   balancing", Parallel Computing, 19(2)L1993, 209--218 

   !! requires pow(2) servers !!
*/
/* Assume that global is sorted */
/* static void multi_level_balancer(const LOAD_LIST *global, */
/* 				 const size_t n_servers, */
/* 				 const size_t n_levels, */
/* 				 const size_t level, */
/* 				 LOAD **dest, */
/* 				 LOAD **ends) */
/* { */
/*   if(global->n_loads == 1) return; */

/*   /\* recursively balance *\/ */
/*   if(level < n_levels){ */
/*     /\* */
/*       multi_level_balancer(lhs ... level + 1); */
/*       multi_level_balancer(rhs ... level + 1); */
/*     *\/ */
/*   } */
/* } */

				

int load_balancer(LOAD_LIST *servers,
		  const size_t n_servers)
{
  /* compile loads into single list */
  LOAD_LIST global;
  get_loads_from_servers(&global,servers,n_servers);
  load_list_sort_load_time(&global);
  load_list_compute_stats(&global);
  load_list_print(stdout,&global);
  printf("\n\n\n");

  /* array of pointers to current server load */
  LOAD **dest = calloc(n_servers,sizeof(*dest));
  LOAD **ends = calloc(n_servers,sizeof(*ends));
  for(size_t i=0; i<n_servers; i++){
    ends[i] = servers[i].loads + servers[i].n_loads;
    dest[i] = servers[i].loads;
  }

  /* round_robin_reassignment(&global,n_servers,dest,ends); */
  greedy_reassignment(&global,n_servers,1,dest,ends);

  for(size_t i=0; i<n_servers; i++){
    load_list_compute_stats(servers+i);
  }

  destroy_LOAD_LIST(&global);
  free(dest);
  free(ends);
  return 0;
}

