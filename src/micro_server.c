/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "micro_server.h"
#include "microscale_information.h"

#include <stdlib.h>
#include <string.h>

void pgf_FE2_micro_server_init(pgf_FE2_micro_server *server)
{
  server->n_jobs = 0;
  server->jobs = NULL;
}

void pgf_FE2_micro_server_build(pgf_FE2_micro_server *server,
				const pgf_FE2_server_rebalance *rebal)
{
  const size_t keep = pgf_FE2_server_rebalance_n_keep(rebal);
  const size_t recv = pgf_FE2_server_rebalance_n_recv(rebal);
  server->n_jobs = keep + recv;
  
  server->jobs = malloc(server->n_jobs*sizeof(*(server->jobs)));
  
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
  server->jobs = NULL;
}

static void pgf_FE2_micro_server_get_job_list(pgf_FE2_micro_server *server,
					      pgf_FE2_server_rebalance *rebal,
					      const PGFEM_mpi_comm *mpi_comm)
{
  pgf_FE2_server_rebalance_build_from_message(rebal,mpi_comm);
  pgf_FE2_micro_server_build(server,rebal);
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

int pgf_new_micro_server_master(const PGFEM_mpi_comm *mpi_comm,
				MICROSCALE *micro)
{
  int err = 0;
  pgf_FE2_micro_server *server = malloc(sizeof(*server));
  pgf_FE2_micro_server_init(server);
  pgf_FE2_server_rebalance *rebal = malloc(sizeof(*rebal));

  int exit_server = 0;
  while(!exit_server){
    /* begin server loop */
    pgf_FE2_micro_server_get_job_list(server,rebal,mpi_comm);

    /* swap microscale information according to rebal */
    /* Future implementation/improvements will overlay this with
       communication of job information from macroscale and computation
       of microstructures that are not to be rebalanced. */
    pgf_FE2_server_rebalance_post_exchange(rebal,mpi_comm,micro);
    pgf_FE2_server_rebalance_finalize_exchange();

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

