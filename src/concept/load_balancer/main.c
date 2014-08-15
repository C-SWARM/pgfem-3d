/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#include "partition.h"
#include "load.h"
#include "lb_transfer.h"
#include "load_balancer.h"
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

static const char *opts = "hp:l:f:";
static const char *usage = 
  "USAGE: load_bal [OPTIONS]\n"
  "OPTIONS:\n"
  "\t-p (int > 0) = 1\t number of partitions\n"
  "\t-l (int >= 1)   \t number of loads > number of partitions\n"
  "\t-f (flt >= 1.0) \t factor for maximum number of loads per partition.\n"
  "\t                  \t Gives max_size = factor*n_load/n_part\n\n";

static double get_rand_in_range(const double min,
				const double max)
{
  return ((max-min) * rand())/RAND_MAX - min;
}

int main(int argc, char **argv)
{
  size_t n_part = 1;
  size_t n_load = 0;
  double factor = 1.0;

  /* parse options from command line */
  {
    int opt;
    while((opt = getopt(argc,argv,opts)) != -1){
      switch(opt){
      case 'p': n_part = atoi(optarg); break;
      case 'l': n_load = atoi(optarg); break;
      case 'f': factor = atof(optarg); break;
      default:
	fprintf(stderr,"%s",usage);
	exit(0);
      }
    }
    /* error check options */
    assert(factor >= 1.0);
    assert(n_part > 0);
    if(n_load <= 0) n_load = n_part;
    assert(n_load >= n_part);
  }

  /* create the initial partitioning */
  PARTITION_LIST *PL = calloc(1,sizeof(*PL));
  PARTITION_LIST_build(PL,n_part,factor*n_load/n_part);

  /* set initial partitions */
  {
    srand(1); /* seed random number generator */
    LOAD tl;
    const double min = 0.0;
    const double max = 1000.0;
    PARTITION *parts = PL->partitions;

    /* quotient */
    for(size_t j=0; j<n_part; j++){
      tl.id.proc = j;
      tl.part_id = j;
      for(size_t i=0, end=n_load/n_part; i<end; i++){
	tl.id.elem = i;
	tl.load = get_rand_in_range(min,max);
	PARTITION_push_load(parts + j,&tl);
      }
    }

    /* remainder */
    for(size_t i=0, end=n_load%n_part; i<end; i++){
      tl.id.proc = i;
      tl.part_id = i;
      tl.id.elem = i + n_load/n_part;
      tl.load = get_rand_in_range(min,max);
      PARTITION_push_load(parts + i,&tl);
    }
  }

  /* print the partition list */
  PARTITION_LIST_print(stdout,PL);

  /* load balance */
  load_balancer_greedy(PL);

  /* print the new partitions */
  PARTITION_LIST_print(stdout,PL);

  TRANSFER_LIST *TL = calloc(1,sizeof(*TL));
  TRANSFER_LIST_compute(PL,TL);
  TRANSFER_LIST_print(stdout,TL);

  /* clean up and exit */
  PARTITION_LIST_destroy(PL);
  free(PL);
  TRANSFER_LIST_destroy(TL);
  free(TL);
  return 0;
}
