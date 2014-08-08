/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "load_list.h"
#include <math.h>
#include <assert.h>
#include <string.h>

/* LOAD FUNCTIONS */
static int load_time_comp(const void *lhs,
			  const void *rhs)
{
  double diff = ((LOAD *) lhs)->time - ((LOAD *) rhs)->time;
  if(diff < 0.0) return -1;
  else if(diff > 0.0) return 1;
  else return 0;
}

static int load_server_comp(const void *lhs,
			    const void *rhs)
{
  return (int) (((LOAD *) lhs)->server_proc
		- ((LOAD *) rhs)->server_proc);
}

/* static int load_server_comp_reverse(const void *lhs, */
/* 				  const void *rhs) */
/* { */
/*   return load_server_comp(rhs,lhs); */
/* } */

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
static int double_comp(const void *lhs,
		       const void *rhs)
{
  double diff = *((double *) lhs) - *((double *) rhs);
  if(diff < 0.0) return -1;
  else if(diff > 0.0) return 1;
  else return 0;
}

static double compute_avg(double *restrict arr,
			  const size_t len)
{
  double avg = 0.0;
  for(size_t i=0; i<len; i++){
    avg += arr[i];
  }
  avg /= len;
  return avg;
}

static double compute_std(double *restrict arr,
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

void stats_compute(STATS *stats,
		   double *arr,
		   const size_t len)
{
  qsort(arr,len,sizeof(*arr),double_comp);
  stats->avg = compute_avg(arr,len);
  stats->std = compute_std(arr,len,stats->avg);
  stats->min = arr[0];
  stats->max = arr[len-1];
}

void stats_reset(STATS *stats)
{
  memset(stats,0,sizeof(*stats));
}

/* LOAD_LIST FUNCTIONS */
void build_LOAD_LIST(LOAD_LIST *list,
		     const size_t n_loads,
		     const size_t n_servers)
{
  list->loads = calloc(n_loads,sizeof(*(list->loads)));
  list->n_loads = n_loads;
  list->n_servers = n_servers;
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
  fprintf(out,"Number of servers: %ld\n",list->n_servers);
  printf("Minimum time:   %f\n",list->time_stats.min);
  printf("Maximum time:   %f\n",list->time_stats.max);
  printf("Average time:   %f\n",list->time_stats.avg);
  printf("Std. Dev. time: %f\n",list->time_stats.std);
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
