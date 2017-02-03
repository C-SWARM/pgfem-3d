/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_micro_server.h"

#include "utils.h"
#include "PGFEM_io.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

void pgf_FE2_micro_server_stats_print(const pgf_FE2_micro_server *server)
{

  const pgf_FE2_micro_server_stats *stat = server->stats;
  PGFEM_printf("SERVER STATS: njob total avg std min max\n");
  PGFEM_printf("%ld %0.5e %0.5e %0.5e %0.5e %0.5e\n\n",
	       server->n_jobs,
	       stat->total,stat->avg,
	       stat->std,stat->min,stat->max);
}

void pgf_FE2_micro_server_init(pgf_FE2_micro_server **server)
{
  *server = malloc(sizeof(**server));

  (*server)->n_jobs = 0;
  (*server)->jobs = NULL;
  (*server)->stats = NULL;
}

void pgf_FE2_micro_server_build(pgf_FE2_micro_server *server,
				const pgf_FE2_server_rebalance *rebal)
{
  const size_t keep = pgf_FE2_server_rebalance_n_keep(rebal);
  const size_t recv = pgf_FE2_server_rebalance_n_recv(rebal);
  server->n_jobs = keep + recv;
  
  server->jobs = malloc(server->n_jobs*sizeof(*(server->jobs)));
  server->stats = malloc(sizeof(*(server->stats)));
  
  /* initialize the jobs that do not need rebalancing */
  const int *job_ids_keep = pgf_FE2_server_rebalance_keep_buf(rebal);
  for(size_t i=0; i<keep; i++){
    pgf_FE2_job_init((server->jobs) + i,job_ids_keep[i],FE2_STATE_NEED_INFO);
  }

  /* initialize the jobs that do not need rebalancing */
  const int *job_ids_recv = pgf_FE2_server_rebalance_recv_buf(rebal);
  for(size_t i=0; i<recv; i++){
    pgf_FE2_job_init((server->jobs) + keep + i,job_ids_recv[i],
		     FE2_STATE_NEED_INFO_REBALANCE);
  }
}

void pgf_FE2_micro_server_destroy(pgf_FE2_micro_server *server)
{
  for(size_t i=0,e=server->n_jobs; i<e; i++){
    pgf_FE2_job_destroy((server->jobs) + i);
  }
  free(server->jobs);
  free(server->stats);
  free(server);
  server = NULL;
}

/**
 * Check to see if all jobs are finished. Returns true or false.
 */
static int pgf_FE2_micro_server_done(pgf_FE2_micro_server *server)
{
  pgf_FE2_job *restrict jobs = server->jobs; /* alias */
  for(size_t i=0, n=server->n_jobs; i<n; i++){
    if(pgf_FE2_job_complete(jobs + i) != FE2_STATE_DONE) return 0;
  }
  return 1;
}

/**
 * Attempt to get information from the microscale.
 */
static void pgf_FE2_micro_server_get_info(pgf_FE2_micro_server *server,
					  const PGFEM_mpi_comm *mpi_comm)
{
  pgf_FE2_job *restrict jobs = server->jobs; /* alias */
  for(size_t i=0, n=server->n_jobs; i<n; i++){
    pgf_FE2_job_get_info(jobs + i,mpi_comm);
  }
}

/**
 * Attempt to compute ready jobs.
 */
static void pgf_FE2_micro_server_compute_ready(pgf_FE2_micro_server *server,
					       MICROSCALE *micro,
					       const PGFEM_mpi_comm *mpi_comm,
					       const int mp_id)
{
  pgf_FE2_job *restrict jobs = server->jobs; /* alias */
  for(size_t i=0, n=server->n_jobs; i<n; i++){
    pgf_FE2_job_compute(jobs + i,micro,mpi_comm,mp_id);
  }
}


/**
 * Busy loop looking for message to start the server cycle.
 */
static void pgf_FE2_micro_server_probe_start(const PGFEM_mpi_comm *mpi_comm,
					     MPI_Status *stat)
{
  int msg_waiting = 0;
  while (1){

    /* exit */
    MPI_Iprobe(MPI_ANY_SOURCE,FE2_MICRO_SERVER_EXIT,
	       mpi_comm->mm_inter,&msg_waiting,stat);
    if(msg_waiting) break;

    /* rebalance */
    MPI_Iprobe(MPI_ANY_SOURCE,FE2_MICRO_SERVER_REBALANCE,
	       mpi_comm->mm_inter,&msg_waiting,stat);
    if(msg_waiting) break;

    /* other ... */
  }

  /* this assert will fail if we accidentally overlap with
     MPI_ANY_TAG. FYI: often MPI_ANY_TAG = -1.*/
  assert(stat->MPI_TAG == FE2_MICRO_SERVER_EXIT 
	 || stat->MPI_TAG == FE2_MICRO_SERVER_REBALANCE);
}


/**
 * On the MASTER server process look for info from the macroscale and
 * propogate to the workers.
 */
static void pgf_FE2_micro_server_start_cycle(const PGFEM_mpi_comm *mpi_comm,
					     pgf_FE2_micro_server *server,
					     pgf_FE2_server_rebalance **rebal,
					     int *exit_server)
{
  static const size_t n_info = 2;
  long info[n_info]; /* = {0 /\* tag *\/, */
		    /*    0 /\* buffer_len *\/}; */
  *exit_server = 0;
  MPI_Status stat;

  /* probe for incomming message on comm->mm_inter */
  pgf_FE2_micro_server_probe_start(mpi_comm,&stat);

  /* analyze source of the message and its content */
#ifndef NDEBUG
  int n_micro_proc = 0;
  int n_macro_proc = 0;
  MPI_Comm_size(mpi_comm->worker_inter,&n_micro_proc);
  MPI_Comm_size(mpi_comm->mm_inter,&n_macro_proc);
  n_macro_proc -= n_micro_proc;
  assert(stat.MPI_SOURCE < n_macro_proc);
#endif
  int buf_len = 0;
  MPI_Get_count(&stat,MPI_CHAR,&buf_len);
  info[0] = stat.MPI_TAG;
  info[1] = buf_len;

  /* allocate buffer for receive */
  char *buf = malloc(buf_len);

  /* post non-blocking receive matching the probed message */
  MPI_Request req;
  MPI_Irecv(buf,buf_len,MPI_CHAR,
	    stat.MPI_SOURCE,stat.MPI_TAG,
	    mpi_comm->mm_inter,&req);

  /* broadcast information to workers !!Signature must match that in
     worker busy loop!!*/
  MPI_Bcast(info,n_info,MPI_LONG,0,mpi_comm->micro);

  /* complete communication w/ macroscale */
  MPI_Wait(&req,MPI_STATUS_IGNORE);

  switch(stat.MPI_TAG){
  case FE2_MICRO_SERVER_EXIT:
    /* set flag and move on to exit function */
    *exit_server = 1;

    /* free the buffer */
    free(buf);
    break;
 
 case FE2_MICRO_SERVER_REBALANCE:
    /* broadcast information to microscale (could change to
       non-blocking?) */
    MPI_Bcast(buf,buf_len,MPI_CHAR,0,mpi_comm->micro);

    /* build rebalance */
    pgf_FE2_server_rebalance_build_from_buffer(rebal,(void**) &buf);

    /* build server list */
    pgf_FE2_micro_server_build(server,*rebal);

    /* do not free the buffer. Maintain access through rebal. Buffer
       is free'd when rebal is destroyed */
    break;
  default: assert(0); /* should never get here */ break;
  }
}

static void pgf_FE2_micro_server_compute_stats(pgf_FE2_micro_server *server)
{
  /* aliases */
  pgf_FE2_micro_server_stats *s = server->stats;
  pgf_FE2_job *j = server->jobs;

  /* initialization */
  const size_t njob = server->n_jobs;
  s->total = s->avg = s->std = 0.0;

  /* total, min, max, and average */
  double total = 0.0;
  double min = INT_MAX;
  double max = 0.0;
  for(size_t i=0; i<njob; i++){
    double cur = (double) j[i].time;
    total += cur;
    if(min > cur) min = cur;
    if(max < cur) max = cur;
  }
  const double avg = total / njob;

  /* standard deviation */
  double std = 0.0;
  for(size_t i=0; i<njob; i++){
    double cur = (double) j[i].time;
    std += (cur-avg)*(cur-avg);
  }
  s->total = total;
  s->avg = avg;
  s->std = sqrt(1./njob*std);
  s->min = min;
  s->max = max;
}					     

/**
 * Make a fairly shallow copy of the job data into a communication
 * buffer.
 */
static void pgf_FE2_micro_server_pack_summary(pgf_FE2_micro_server *server,
					      size_t *buf_len,
					      char **buf)
{
  /* compute stats */
  pgf_FE2_micro_server_compute_stats(server);

  /* compute buffer length and allocate */
  const size_t n_jobs = server->n_jobs;
  static const size_t size_job = sizeof(*(server->jobs));
  static const size_t size_stats = sizeof(*(server->stats));
  *buf_len = sizeof(server->n_jobs) + n_jobs*size_job + size_stats;
  (*buf) = malloc(*buf_len);

  /* pack the server job summary */
  size_t pos = 0;
  pack_data(&n_jobs,*buf,&pos,1,sizeof(n_jobs));
  pack_data(server->stats,*buf,&pos,1,size_stats);
  pack_data(server->jobs,*buf,&pos,n_jobs,size_job);
}

static void pgf_FE2_micro_server_finish_cycle(const PGFEM_mpi_comm *mpi_comm,
					      pgf_FE2_micro_server *server)
{
  int n_micro_proc = 0;
  int n_macro_proc = 0;
  MPI_Comm_size(mpi_comm->worker_inter,&n_micro_proc);
  MPI_Comm_size(mpi_comm->mm_inter,&n_macro_proc);
  n_macro_proc -= n_micro_proc;
  MPI_Request *req = malloc(n_macro_proc*sizeof(*req));
  char *buf = NULL;
  size_t buf_len = 0;
  pgf_FE2_micro_server_pack_summary(server,&buf_len,&buf);

  for(int i=0; i<n_macro_proc; i++){
    MPI_Isend(buf,buf_len,MPI_CHAR,i,
	      FE2_MICRO_SERVER_REBALANCE,
	      mpi_comm->mm_inter,req+i);
  }

  MPI_Waitall(n_macro_proc,req,MPI_STATUS_IGNORE);
  free(req);
  free(buf);
}

void pgf_FE2_micro_server_unpack_summary(pgf_FE2_micro_server **Server,
					 const char *buf)
{
  pgf_FE2_micro_server_init(Server);
  pgf_FE2_micro_server *server = *Server; /* alias */
  static const size_t size_jobs = sizeof(*(server->jobs));
  static const size_t size_stats = sizeof(*(server->stats));

  server->stats = malloc(size_stats);
  size_t pos = 0;
  unpack_data(buf,&(server->n_jobs),&pos,1,sizeof(server->n_jobs));
  server->jobs = malloc(server->n_jobs*size_jobs);
  unpack_data(buf,server->stats,&pos,1,size_stats);
  unpack_data(buf,server->jobs,&pos,server->n_jobs,size_jobs);
  for(int i=0,e=server->n_jobs; i<e; i++){
    server->jobs[i].comm_buf = NULL;
  }
}

static int pgf_FE2_micro_server_master(const PGFEM_mpi_comm *mpi_comm,
				       MICROSCALE *micro,
				       const int mp_id)
{
  int err = 0;
  int exit_server = 0;

  while(1){
    /* begin server loop */
    pgf_FE2_micro_server *server = NULL;
    pgf_FE2_server_rebalance *rebal = NULL;

    pgf_FE2_micro_server_init(&server);
    pgf_FE2_micro_server_start_cycle(mpi_comm,server,&rebal,&exit_server);
    if(exit_server){
      pgf_FE2_micro_server_destroy(server);
      pgf_FE2_server_rebalance_destroy(rebal);
      break;
    }

    /* swap microscale information according to rebal */
    /* Future implementation/improvements will overlay this with
       communication of job information from macroscale and computation
       of microstructures that are not to be rebalanced. */
    pgf_FE2_server_rebalance_post_exchange(rebal,mpi_comm,micro);
    /* blocks until all exchanges are done. Could improve
       async. comm here. */
    pgf_FE2_server_rebalance_finalize_exchange(rebal,mpi_comm);

    /* for now, go ahead and mark all jobs as needing information */
    for(size_t i=0,e=server->n_jobs; i<e; i++){
      pgf_FE2_job_set_state(server->jobs + i,FE2_STATE_NEED_INFO);
    }

    while(!pgf_FE2_micro_server_done(server)){
      /* Get info from macroscale */
      pgf_FE2_micro_server_get_info(server,mpi_comm);

      /* in improved implementation, advance rebalancing here. I.e.,
	 check for rebalancing jobs where the communication has
	 completed and update their state appropriately. */

      /* compute the jobs that are ready. Posts sends */
      pgf_FE2_micro_server_compute_ready(server,micro,mpi_comm,mp_id);

    }

    /* send server/job statistics to macroscale for rebalancing */
    pgf_FE2_micro_server_finish_cycle(mpi_comm,server);
    pgf_FE2_micro_server_stats_print(server);

    pgf_FE2_micro_server_destroy(server);
    pgf_FE2_server_rebalance_destroy(rebal);
  }

  return err;
}

/**
 * Worker busy loop. Some logic as to which jobs to initiate.
 */
static int pgf_FE2_micro_server_worker(const PGFEM_mpi_comm *mpi_comm,
				       MICROSCALE *micro,
				       const int mp_id)
{
  int err = 0;
  int exit_server = 0;
  static const size_t n_info = 2;
  long info[n_info]; /* = {0 /\*id*\/, */
		     /*   0 /\*msg_size in bytes *\/}; */
  while( !exit_server ){
    /* get information from master on server */
    MPI_Bcast(info,n_info,MPI_LONG,0,micro->common->mpi_comm);

    switch(info[0]){
    case FE2_MICRO_SERVER_EXIT:
      exit_server = 1;
      break;

    case FE2_MICRO_SERVER_REBALANCE:
      /* rebalance on workers */
      {
	/* allocate buffer for receive */
	int buf_len = info[1];
	char *buf = malloc(buf_len);
	pgf_FE2_server_rebalance *rebal = NULL;

	/* get information from MASTER */
	MPI_Bcast(buf,buf_len,MPI_CHAR,0,mpi_comm->micro);

	/* build rebalancing data structure */
	pgf_FE2_server_rebalance_build_from_buffer(&rebal,(void**) &buf);

	/* perform rebalancing */
	pgf_FE2_server_rebalance_post_exchange(rebal,mpi_comm,micro);
	/* blocks until all exchanges are done. Could improve
	   async. comm here. */
	pgf_FE2_server_rebalance_finalize_exchange(rebal,mpi_comm);

	/* implicitly free's buf */
	pgf_FE2_server_rebalance_destroy(rebal);
      }
      break;

    default: /* valid job, compute work */
      pgf_FE2_job_compute_worker(info[0],info[1],micro,mp_id);
      break;
    }
  }

  return err;
}

int pgf_FE2_micro_server_START(const PGFEM_mpi_comm *mpi_comm,
			       MICROSCALE *micro,
			       const int mp_id)
{
  int err = 0;
  assert(mpi_comm->valid_micro);
  if(mpi_comm->valid_mm_inter){
   err += pgf_FE2_micro_server_master(mpi_comm,micro,mp_id);
  } else {
    err += pgf_FE2_micro_server_worker(mpi_comm,micro,mp_id);
  }
  return err;
}
