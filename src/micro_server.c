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

int pgf_new_micro_server_master(const PGFEM_mpi_comm *comm)
{
  int err = 0;

  /* look for partition/server info on pgf_mpi_comm->mm_inter */
  {
    int len = 0;
    char *buffer = NULL;
    {
      int flag = 0;
      MPI_Status rebal_stat;
      while(!flag){
	err += MPI_Iprobe(MPI_ANY_SOURCE,-1,comm->mm_inter,&flag,&rebal_stat);
      }

      /* post appropriate receive for partition information from macroscale */
      err += MPI_Get_count(&rebal_stat,MPI_CHAR,&len);
      buffer = calloc(len,sizeof(*buffer));
      MPI_Request rebal_req;
      err += MPI_Irecv(buffer,len,MPI_CHAR,rebal_stat.MPI_SOURCE,
		       rebal_stat.MPI_TAG,comm->mm_inter,&rebal_req);
      err += MPI_Wait(&rebal_req,MPI_STATUS_IGNORE);
    }

    /* build work_queue for current step */

    /* spawn rebalancing on pgf_mpi_comm->micro */

    /* rebalance master(s) on pgf_mpi_comm->worker_inter */

    /* Complete rebalancing -- mark as NEED_INFO */

    free(buffer);
  }

  /* post receives for jobs to compute (mm_inter)*/

  while(1 /* !work_queue_is_empty(queue) */){

    /* MPI_Waitsome for information. */

    /* Mark jobs in work_queue as READY */

    /* Compute READY jobs in work_queue */

    /* Post response and mark REPLY */

    /* MPI_Testsome for finished replies. Mark DONE */

  }

  return err;
}

