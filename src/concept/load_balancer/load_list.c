/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_list.h"
#include "load_balancer_utils.h"
#include <string.h>

/* LOAD FUNCTIONS */
static int load_time_comp(const void *lhs,
			  const void *rhs)
{
  return double_comp((void*) &(((LOAD *) lhs)->time),
		     (void*) &(((LOAD *) rhs)->time));
}

static int load_server_comp(const void *lhs,
			    const void *rhs)
{
  return size_t_comp((void*) &(((LOAD *) lhs)->server_proc),
		     (void*) &(((LOAD *) rhs)->server_proc));
}

static int load_time_comp_reverse(const void *lhs,
				  const void *rhs)
{
  return load_time_comp(rhs,lhs);
}

static void load_print(FILE *out,
		       const LOAD *load)
{
  fprintf(out,"%11.3e c:%ld s:%ld",load->time,
	  load->client_proc,load->server_proc);
}

/* STATS FUNCTIONS */
void stats_compute(STATS *stats,
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

void stats_reset(STATS *stats)
{
  memset(stats,0,sizeof(*stats));
}

void stats_print(FILE *out,
		 const STATS *stats)
{
  printf("Minimum time:   %f\n",stats->min);
  printf("Maximum time:   %f\n",stats->max);
  printf("Average time:   %f\n",stats->avg);
  printf("Std. Dev. time: %f\n",stats->std);
  printf("Total time:     %f\n",stats->total);
}

/* LOAD_LIST FUNCTIONS */
void build_LOAD_LIST(LOAD_LIST *list,
		     const size_t n_loads)
{
  list->loads = calloc(n_loads,sizeof(*(list->loads)));
  list->n_loads = n_loads;
}

void destroy_LOAD_LIST(LOAD_LIST *list)
{
  free(list->loads); list->loads = NULL;
}

LOAD* load_list_get_load(LOAD_LIST *list,
			 const size_t idx)
{
  return list->loads + idx;
}

void load_list_sort_load_time(LOAD_LIST *list)
{
  qsort(list->loads,list->n_loads,
	sizeof(*(list->loads)),load_time_comp_reverse);
}

void load_list_sort_load_server(LOAD_LIST *list)
{
  qsort(list->loads,list->n_loads,
	sizeof(*(list->loads)),load_server_comp);
}

double* load_list_extract_load_times(const LOAD_LIST *list)
{
  LOAD *restrict loads = list->loads;
  double *restrict times = calloc(list->n_loads,sizeof(*times));
  for(int i=0, end=list->n_loads; i<end; i++){
    times[i] = loads[i].time;
  }
  return times;
}

void load_list_print(FILE *out,
		     const LOAD_LIST *list)
{
  const LOAD *loads = list->loads;
  fprintf(out,"Number of loads:   %ld\n",list->n_loads);
  stats_print(out,&list->time_stats);
  fprintf(out,"LOADS:\n");
  for(int i=0, end=list->n_loads; i<end; i++){
    load_print(out,loads + i);
    fprintf(out,"\n");
  }
}

void load_list_compute_stats(LOAD_LIST *list)
{
  double *times = load_list_extract_load_times(list);
  stats_compute(&list->time_stats,times,list->n_loads);
  free(times);
}

void load_list_reset(LOAD_LIST *list)
{
  memset(list->loads,0,list->n_loads*sizeof(*(list->loads)));
  stats_reset(&list->time_stats);
}
