/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_micro_server.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

void pgf_FE2_micro_server_init(pgf_FE2_micro_server *server)
{
  server->n_jobs = 0;
  server->jobs = NULL;
  server->stats = NULL;
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

  /* initialize the jobs that do no need rebalancing */
  const int *job_ids_recv = pgf_FE2_server_rebalance_recv_buf(rebal);
  for(size_t i=0; i<recv; i++){
    pgf_FE2_job_init((server->jobs) + keep + i,job_ids_keep[i],
		     FE2_STATE_NEED_INFO_REBALANCE);
  }
}

void pgf_FE2_micro_server_destroy(pgf_FE2_micro_server *server)
{
  for(size_t i=0,e=server->n_jobs; i<e; i++){
    pgf_FE2_job_destroy((server->jobs) + i);
  }
  server->n_jobs = 0;
  free(server->jobs);
  free(server->stats);
  server->jobs = NULL;
  server->stats = NULL;
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
					       const PGFEM_mpi_comm *mpi_comm)
{
  pgf_FE2_job *restrict jobs = server->jobs; /* alias */
  for(size_t i=0, n=server->n_jobs; i<n; i++){
    pgf_FE2_job_compute(jobs + i,micro,mpi_comm);
  }
}


/**
 * On the MASTER server process look for info from the macroscale and
 * propogate to the workers.
 */
static void pgf_FE2_micro_server_start_cycle(const PGFEM_mpi_comm *mpi_comm,
					     pgf_FE2_micro_server *server,
					     pgf_FE2_server_rebalance *rebal,
					     int *exit_server)
{
  static const size_t n_info = 2;
  long info[n_info]; /* = {0 /\* tag *\/, */
		    /*    0 /\* buffer_len *\/}; */
  *exit_server = 0;
  MPI_Status stat;

  /* probe for incomming message on comm->mm_inter */
  MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,mpi_comm->mm_inter,&stat);

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
    break;
 
 case FE2_MICRO_SERVER_REBALANCE:
    /* broadcast information to microscale (could change to
       non-blocking?) */
    MPI_Bcast(buf,buf_len,MPI_CHAR,0,mpi_comm->micro);

    /* build rebalance */
    pgf_FE2_server_rebalance_build_from_buffer(rebal,buf);

    /* build server list */
    pgf_FE2_micro_server_build(server,rebal);

    break;
  default: assert(0); /* should never get here */ break;
  }

  free(buf);
}


static int pgf_FE2_micro_server_master(const PGFEM_mpi_comm *mpi_comm,
				       MICROSCALE *micro)
{
  int err = 0;
  pgf_FE2_micro_server *server = malloc(sizeof(*server));
  pgf_FE2_micro_server_init(server);
  pgf_FE2_server_rebalance *rebal = malloc(sizeof(*rebal));
  int exit_server = 0;

  while(1){
    /* begin server loop */
    pgf_FE2_micro_server_start_cycle(mpi_comm,server,rebal,&exit_server);
    if(exit_server) break;

    /* swap microscale information according to rebal */
    /* Future implementation/improvements will overlay this with
       communication of job information from macroscale and computation
       of microstructures that are not to be rebalanced. */
    pgf_FE2_server_rebalance_post_exchange(rebal,mpi_comm,micro);
    /* blocks until all exchanges are done. Could improve
       async. comm here. */
    pgf_FE2_server_rebalance_finalize_exchange(mpi_comm);

    /* for now, go ahead and mark all jobs as needing information */
    for(size_t i=0,e=server->n_jobs; i<e; i++){
      server->jobs->state = FE2_STATE_NEED_INFO;
    }

    while(!pgf_FE2_micro_server_done(server)){
      /* Get info from macroscale */
      pgf_FE2_micro_server_get_info(server,mpi_comm);

      /* in improved implementation, advance rebalancing here. I.e.,
	 check for rebalancing jobs where the communication has
	 completed and update their state appropriately. */

      /* compute the jobs that are ready. Posts sends */
      pgf_FE2_micro_server_compute_ready(server,micro,mpi_comm);

    }
  }

  pgf_FE2_micro_server_destroy(server);
  free(server);

  pgf_FE2_server_rebalance_destroy(rebal);
  free(rebal);

  return err;
}

/**
 * Worker busy loop. Some logic as to which jobs to initiate.
 */
static int pgf_FE2_micro_server_worker(const PGFEM_mpi_comm *mpi_comm,
				       MICROSCALE *micro)
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
	pgf_FE2_server_rebalance *rebal = malloc(sizeof(*rebal));
	/* get information from MASTER */
	MPI_Bcast(buf,buf_len,MPI_CHAR,0,mpi_comm->micro);

	/* build rebalancing data structure */
	pgf_FE2_server_rebalance_build_from_buffer(rebal,buf);

	/* perform rebalancing */
	pgf_FE2_server_rebalance_post_exchange(rebal,mpi_comm,micro);
	/* blocks until all exchanges are done. Could improve
	   async. comm here. */
	pgf_FE2_server_rebalance_finalize_exchange(mpi_comm);

	free(buf);
	pgf_FE2_server_rebalance_destroy(rebal);
	free(rebal);
      }
      break;

    default: /* valid job, compute work */
      pgf_FE2_job_compute_worker(info[0],info[1],micro);
      break;
    }
  }

  return err;
}

/**
 * This is the main function that starts the master/worker server
 * processes.
 */
int pgf_FE2_micro_server_START(const PGFEM_mpi_comm *mpi_comm,
			       MICROSCALE *micro)
{
  int err = 0;
  assert(mpi_comm->valid_micro);
  if(mpi_comm->valid_mm_inter){
   err += pgf_FE2_micro_server_master(mpi_comm,micro);
  } else {
    err += pgf_FE2_micro_server_worker(mpi_comm,micro);
  }
  return err;
}
