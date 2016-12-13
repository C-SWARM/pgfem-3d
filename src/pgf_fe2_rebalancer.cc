/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */
#include "pgf_fe2_rebalancer.h"

#include "pgf_fe2_job.h"
#include "pgf_fe2_server_rebalance.h"
#include "pgf_fe2_micro_server.h"
#include "PGFEM_mpi.h"
#include "PGFEM_io.h"
#include <string.h>
#include <assert.h>
#include <limits.h>

/* turn on/off adaptive rebalancing print statements */
static const int LOGGING = 1;

/* a helper structure */
enum new_partition_idx{
  NEW_PART_JOB_ID,
  NEW_PART_TIME,
  NEW_PART_SR,
  NEW_PART_NUMEL
};

static const size_t NEW_PART_JOB_SIZE = NEW_PART_NUMEL*sizeof(size_t);

typedef struct new_partition{
  size_t max_size;
  size_t n_job;
  size_t n_keep;
  size_t n_send;
  size_t n_recv;
  size_t total_time;
  size_t *keep;
  size_t *send;
  size_t *recv;
} new_partition;

/*** static function declarations ***/
/**
 * compare two objects of helper structure new_partition by
 * time. (decreasing order)
 */
static int new_partition_buf_compare_time(const void *a,
					  const void *b);

/**
 * Build a new_partition object with space for max_n_jobs.
 */
static void new_partition_build(new_partition *np,
				const size_t max_n_jobs);

/**
 * Destroy a new_partition object.
 */
static void new_partition_destroy(new_partition *np);

/**
 * Sort a new_partition object keep buffer by time.
 */
static void new_partition_sort_keep_time(new_partition *np);

/**
 * Set a new_partition object to behave as if it is empty.
 */
static void new_partition_set_empty(new_partition *np);

/**
 * Get an array offset for appending to the "keep" buffer of a
 * new_partition object.
 */
static size_t new_partition_get_offset_keep(const new_partition *np);

/**
 * Get an array offset for appending to the "send" buffer of a
 * new_partition object.
 */
static size_t new_partition_get_offset_send(const new_partition *np);

/**
 * Get an array offset for appending to the "recv" buffer of a
 * new_partition object.
 */
static size_t new_partition_get_offset_recv(const new_partition *np);

/**
 * Get a pointer to the idx-th job in the "keep" buffer of the
 * new_partition object.
 */
static size_t* new_partition_get_ptr_to_job_keep(const new_partition *np,
						 const size_t idx);

/**
 * Format conversion from the helper struct to one more communication
 * friendly.
 */
static void new_partition_to_pgf_FE2_server_rebalance(const new_partition *np,
						      pgf_FE2_server_rebalance **rb);

/**
 * Assign a job to a given partition. Updates send buffers on other
 * partitions as needed.
 */
static void rebalance_partitions_assign_job(const size_t *job,
					    const size_t part_id,
					    const size_t n_parts,
					    new_partition *parts);

/**
 * Get the index to the partition with the smallest total time that is
 * not full.
 *
 * Returns a poisoned value if there are no valid partitions, i.e.,
 * max_size specified at allocation is insufficient for the number of
 * jobs/partitions.
 */
static size_t get_smallest_part_idx(const size_t n_parts,
				    const new_partition *parts);


/**
 * Keeps the same partitioning as the previous step.
 */
static void rebalance_partitions_none(const size_t n_parts,
				      new_partition *all_parts,
				      new_partition *parts);

/**
 * Chooses between rebalancing via greedy algorithm and original
 * partitioning. This allows the system to decide whether it is worth
 * it to migrate data.
 */
static void rebalance_partitions_adaptive(const size_t n_parts,
					  const double *orig_totals,
					  new_partition *all_parts,
					  new_partition *parts);

/**
 * Format conversion server->new_partition (reduces footprint).
 */
static void new_partition_extract_server_as_keep(new_partition *all_parts,
						 const size_t server_id,
						 const pgf_FE2_micro_server *server);

/*** API ***/
/**
 * Very ugly initialization of data structure I don't want to keep around...
 */
void new_partition_build_set_keep(void **np,
				  const size_t max_n_job,
				  const size_t n_job,
				  const int *job_id,
				  const int *job_time,
				  const int *proc_id)
{
  *np = malloc(sizeof(new_partition));
  new_partition *NP = *np;
  new_partition_build(NP,max_n_job);

  NP->n_keep = n_job;
  for(size_t i=0; i<n_job; i++){
    size_t *ptr = new_partition_get_ptr_to_job_keep(NP,i);
    ptr[NEW_PART_JOB_ID] = job_id[i];
    ptr[NEW_PART_TIME] = job_time[i];
    ptr[NEW_PART_SR] = proc_id[i];
  }
}

/**
 * destroy via opaque handle
 */
void new_partition_destroy_void(void *np)
{
  new_partition_destroy(np);
  free(np);
  np = NULL;
}

void new_partitions_void(void **parts,
			 const size_t n_parts,
			 const size_t n_max_job)
{
  *parts = malloc(n_parts*sizeof(new_partition));
  new_partition *NP = *parts; /* alias */
  for(size_t i=0; i<n_parts; i++){
    new_partition_build(NP + i,n_max_job);
  }
}

void new_partitions_destroy_void(void *np,
				 const size_t n_parts)
{
  new_partition *NP = np;
  for(size_t i=0; i<n_parts; i++){
    new_partition_destroy(NP+i);
  }
  free(np);
  np = NULL;
}

void new_partitions_void_to_pgf_FE2_server_rebalance(const int n_parts,
						     const void *np,
						     pgf_FE2_server_rebalance **rb)
{
  auto NP = (new_partition *) np;
  for(int i=0; i<n_parts; i++){
    new_partition_to_pgf_FE2_server_rebalance(NP + i, rb + i);
  }
}

void rebalance_partitions_greedy(const size_t n_parts,
				 void *All_parts,
				 void *Parts)
{
  /* cast pointers */
  new_partition *all_parts = All_parts;
  new_partition *parts = Parts;
  /* sort the full job list by time */
  new_partition_sort_keep_time(all_parts);

  /* reset the partitions */
  for(size_t i=0; i<n_parts; i++){
    new_partition_set_empty(parts + i);
  }

  /* give each partition one job */
  for(size_t i=0; i<n_parts; i++){
    const size_t *job = new_partition_get_ptr_to_job_keep(all_parts,i);
    rebalance_partitions_assign_job(job,i,n_parts,parts);
  }

  /* assign remaining jobs iteratively to smallest partition */
  for(size_t i=n_parts, total_n_jobs = all_parts->n_keep;
      i<total_n_jobs; i++){
    const size_t *job = new_partition_get_ptr_to_job_keep(all_parts,i);
    const size_t part_id = get_smallest_part_idx(n_parts,parts);
    rebalance_partitions_assign_job(job,part_id,n_parts,parts);
  }

  /* compute stats on rebalanced partitions */
  /* ... */
}

pgf_FE2_server_rebalance** pgf_FE2_rebalancer(const PGFEM_mpi_comm *mpi_comm,
					     const size_t total_n_jobs,
					     const size_t max_n_jobs,
					     const int heuristic)
{
  /* get rank and number of macro and micro porcs. */
  int n_macro_proc = 0;
  int n_micro_proc = 0;
  MPI_Comm_size(mpi_comm->mm_inter,&n_micro_proc);
  MPI_Comm_size(mpi_comm->macro,&n_macro_proc);
  n_micro_proc -= n_macro_proc;

  char **buf = malloc(n_micro_proc*sizeof(*buf));
  MPI_Status stat;
  MPI_Request *req = malloc(n_micro_proc*sizeof(*req));
  /* probe for communication, allocate and post matching receives. */
  /* Can improve asynchrony here with some work. */
  for(int i=0; i<n_micro_proc; i++){
    int src = i + n_macro_proc;
    MPI_Probe(src,FE2_MICRO_SERVER_REBALANCE,
	      mpi_comm->mm_inter,&stat);
    int count = 0;
    MPI_Get_count(&stat,MPI_CHAR,&count);
    buf[i] = malloc(count);
    MPI_Irecv(buf[i],count,MPI_CHAR,src,FE2_MICRO_SERVER_REBALANCE,
	      mpi_comm->mm_inter,req + i);
  }

  /* build some helper data structures */
  pgf_FE2_micro_server *server = NULL;
  pgf_FE2_micro_server_stats *all_stats = malloc(n_micro_proc*sizeof(*all_stats));
  new_partition *all_parts = malloc(sizeof(*all_parts));
  new_partition_build(all_parts,total_n_jobs);
  new_partition *parts = malloc(n_micro_proc*sizeof(*parts));
  for(int i=0; i<n_micro_proc; i++){
    new_partition_build(parts + i,max_n_jobs);
  }

  /* receive buffers and unpack into data structure */
  int n_finished = 0;
  while(n_finished < n_micro_proc){
    int idx = 0;
    MPI_Waitany(n_micro_proc,req,&idx,MPI_STATUS_IGNORE);
    pgf_FE2_micro_server_unpack_summary(&server,buf[idx]);

    /* copy the server stats */
    memcpy(all_stats+idx,server->stats,sizeof(*all_stats));

    /* push the job information into all_parts */
    new_partition_extract_server_as_keep(all_parts,idx,server);

    /* destroy the server object, but keep handle. */
    pgf_FE2_micro_server_destroy(server);
    n_finished++;
  }

  /* rebalance partitions according to the heuristic */
  switch(heuristic){
  case FE2_REBALANCE_NONE:
    rebalance_partitions_none(n_micro_proc,all_parts,parts);
    break;
  case FE2_REBALANCE_GREEDY:
    rebalance_partitions_greedy(n_micro_proc,all_parts,parts);
    break;
  case FE2_REBALANCE_ADAPTIVE:
    {
      /* push unbalanced server totals into array */
      double *orig_totals = malloc(n_micro_proc*sizeof(*orig_totals));
      for(int i=0; i<n_micro_proc; i++){
	orig_totals[i] = all_stats[i].total;
      }

      /* adaptively rebalance */
      rebalance_partitions_adaptive(n_micro_proc,orig_totals,all_parts,parts);
      free(orig_totals);
    }
    break;
  default:
    /* Should not get here */
    assert(0);
    break;
  }

  /* push new_partitions to rebalance data structure for communication */
  pgf_FE2_server_rebalance **rb = calloc(n_micro_proc,sizeof(rb));
  for(int i=0; i<n_micro_proc; i++){
    new_partition_to_pgf_FE2_server_rebalance(parts + i,rb + i);
  }


  for(int i=0; i<n_micro_proc; i++){
    free(buf[i]);
    new_partition_destroy(parts+i);
    /* pgf_FE2_server_rebalance_destroy(rb+i); */
  }
  free(buf);
  free(req);
  free(parts);
  /* free(rb); */
  new_partition_destroy(all_parts);
  free(all_parts);
  free(all_stats);
  return rb;
}

/*** STATIC FUNCTIONS ****/
static int new_partition_buf_compare_time(const void *a,
					  const void *b)
{
  auto A = (size_t *) a;
  auto B = (size_t *) b;
  return ((A[NEW_PART_TIME] < B[NEW_PART_TIME])
	  - (A[NEW_PART_TIME] > B[NEW_PART_TIME]));
}

static void new_partition_build(new_partition *np,
				const size_t max_n_jobs)
{
  np->max_size = max_n_jobs;
  np->n_job = 0;
  np->n_keep = 0;
  np->n_send = 0;
  np->n_recv = 0;
  np->total_time = 0;
  const size_t len = max_n_jobs*NEW_PART_JOB_SIZE;
  np->keep = malloc(len);
  np->send = malloc(len);
  np->recv = malloc(len);
}

static void new_partition_destroy(new_partition *np)
{
  free(np->keep);
  free(np->send);
  free(np->recv);
  np->max_size = np->n_keep = np->n_recv = np->n_send = 0;
}

static void new_partition_sort_keep_time(new_partition *np)
{
  qsort(np->keep,np->n_keep,NEW_PART_JOB_SIZE,
	new_partition_buf_compare_time);
}

static void new_partition_set_empty(new_partition *np)
{
  np->n_job = np->n_keep = np->n_send = np->n_recv = np->total_time = 0;
}

static size_t new_partition_get_offset_keep(const new_partition *np)
{
  return np->n_keep * NEW_PART_NUMEL;
}

static size_t new_partition_get_offset_send(const new_partition *np)
{
  return np->n_send * NEW_PART_NUMEL;
}

static size_t new_partition_get_offset_recv(const new_partition *np)
{
  return np->n_recv * NEW_PART_NUMEL;
}


static size_t* new_partition_get_ptr_to_job_keep(const new_partition *np,
						 const size_t idx)
{
  const size_t i = idx * NEW_PART_NUMEL;
  return (np->keep + i);
}


static void new_partition_to_pgf_FE2_server_rebalance(const new_partition *np,
						      pgf_FE2_server_rebalance **rb)
{
  pgf_FE2_server_rebalance_build(rb,np->n_keep,np->n_send,np->n_recv);
  int *keep = pgf_FE2_server_rebalance_keep_buf(*rb);
  int *send = pgf_FE2_server_rebalance_send_buf(*rb);
  int *send_to = pgf_FE2_server_rebalance_send_dest(*rb);
  int *recv = pgf_FE2_server_rebalance_recv_buf(*rb);
  int *recv_from = pgf_FE2_server_rebalance_recv_src(*rb);

  for(size_t i=0, n_keep = np->n_keep; i<n_keep; i++){
    keep[i] = np->keep[i*NEW_PART_NUMEL + NEW_PART_JOB_ID];
  }

  for(size_t i=0, n_send = np->n_send; i<n_send; i++){
    send[i] = np->send[i*NEW_PART_NUMEL + NEW_PART_JOB_ID];
    send_to[i] = np->send[i*NEW_PART_NUMEL + NEW_PART_SR];
  }

  for(size_t i=0, n_recv = np->n_recv; i<n_recv; i++){
    recv[i] = np->recv[i*NEW_PART_NUMEL + NEW_PART_JOB_ID];
    recv_from[i] = np->recv[i*NEW_PART_NUMEL + NEW_PART_SR];
  }
}


static void rebalance_partitions_assign_job(const size_t *job,
					    const size_t part_id,
					    const size_t n_parts,
					    new_partition *parts)
{
  /* job is a pointer to some job information to copy */
  const size_t job_sr = job[NEW_PART_SR];
  if(job_sr == part_id){ /* job is kept */
    size_t keep_off = new_partition_get_offset_keep(parts+part_id);
    memcpy(parts[part_id].keep + keep_off,job,NEW_PART_JOB_SIZE);
    parts[part_id].n_keep++;
  } else { /* job is communicated */
    /* set receive */
    {
      const size_t recv_off = new_partition_get_offset_recv(parts+part_id);
      memcpy(parts[part_id].recv + recv_off,job,NEW_PART_JOB_SIZE);
      parts[part_id].n_recv++;
    }

    /* set send on other partition */
    {
      const size_t send_off = new_partition_get_offset_send(parts+job_sr);
      memcpy(parts[job_sr].send + send_off,job,NEW_PART_JOB_SIZE);
      /* reset destination */
      (parts[job_sr].send + send_off)[NEW_PART_SR] = part_id;
      parts[job_sr].n_send++;
    }
  }

  /* update number of jobs on partition and total time on partition */
  parts[part_id].n_job++;
  parts[part_id].total_time += job[NEW_PART_TIME];
}


static size_t get_smallest_part_idx(const size_t n_parts,
				    const new_partition *parts)
{
  size_t idx = -1; /* poisoned */
  size_t min = -1; /* overflow --> max */
  for(size_t i=0; i<n_parts; i++){
    size_t t = parts[i].total_time;
    if(t < min && parts[i].n_job < parts[i].max_size){
      min = t;
      idx = i;
    }
  }
  return idx;
}

static void rebalance_partitions_none(const size_t n_parts,
				      new_partition *all_parts,
				      new_partition *parts)
{

  /* reset the partitions */
  for(int i=0; i<n_parts; i++){
    new_partition_set_empty(parts + i);
  }

  for(size_t i=0, total_n_jobs = all_parts->max_size;
      i<total_n_jobs; i++){
    const size_t *job = new_partition_get_ptr_to_job_keep(all_parts,i);
    const size_t part_id = job[NEW_PART_SR];
    rebalance_partitions_assign_job(job,part_id,n_parts,parts);
  }

  /* compute stats on rebalanced partitions */
  /* ... */
}

static void new_partition_extract_server_as_keep(new_partition *all_parts,
						 const size_t server_id,
						 const pgf_FE2_micro_server *server)
{
  const pgf_FE2_job *j = server->jobs;
  size_t *keep = all_parts->keep;
  for(size_t i=0, e=server->n_jobs; i<e; i++){
    size_t off = new_partition_get_offset_keep(all_parts);
    keep[off + NEW_PART_JOB_ID] = j[i].id;
    keep[off + NEW_PART_TIME] = j[i].time;
    keep[off + NEW_PART_SR] = server_id;
    all_parts->n_keep ++;
  }
}


static void rebalance_partitions_adaptive(const size_t n_parts,
					  const double *orig_totals,
					  new_partition *all_parts,
					  new_partition *parts)
{
  int rank = -1;
  PGFEM_Error_rank(&rank);
  /* compute the greedy partitioning. */
  rebalance_partitions_greedy(n_parts,all_parts,parts);

  /* compute statistics on original partitioning and new
     partitioning. Note that total and average are unchanged between
     the old and new partitionings. Most important statistic is max
     time and wait time */
  double orig_min = INT_MAX; /* should be big enough */
  double orig_max = -1.0; /* actual values are (+) */
  double orig_wait = 0.0;
  double new_min = INT_MAX; /* should be big enough */
  double new_max = -1.0; /* actual values are (+) */
  double new_wait = 0.0;

  for(size_t i=0; i<n_parts; i++){
    /* min */
    if(orig_min > orig_totals[i]) orig_min = orig_totals[i];
    if(new_min > parts[i].total_time) new_min = (double) parts[i].total_time;

    /* max */
    if(orig_max < orig_totals[i]) orig_max = orig_totals[i];
    if(new_max < parts[i].total_time) new_max = (double) parts[i].total_time;
  }

  orig_wait = orig_max - orig_min;
  new_wait = new_max - new_min;

  /*
    Only rebalance if the estimated max server time is reduced by some
    factor. Wait time, while a good metric for resource utilization,
    does not matter in terms of scalability or overall
    efficiency. Faster simulation = lower CPU hours. period.
 */
  static const double max_factor = 0.95; /* save ~1hr on a 24hr simulation */
  if ( new_max < max_factor * orig_max ){ /* reduced total time */
    if ( rank == 0 && LOGGING){
      PGFEM_printerr("REBALANCING:\n"
		     "       max       min       wait\n"
		     "OLD: %8.3e %8.3e %8.3e\n"
		     "NEW: %8.3e %8.3e %8.3e\n"
		     "N_JOB:\n",
		     orig_max,orig_min,orig_wait,
		     new_max,new_min,new_wait);
      for(size_t i = 0; i < n_parts; i++){
	PGFEM_printerr("%ld ",parts[i].n_job);
	if( !((i + 1) % 20) ) PGFEM_printerr("\n");
      }
      PGFEM_printerr("\n");
    }
  } else { /* do not rebalance */
    /* reset the parts */
    for(size_t i=0; i<n_parts; i++){
      new_partition_set_empty(parts + i);
    }

    /* extract original partition */
    rebalance_partitions_none(n_parts,all_parts,parts);
  }

}
