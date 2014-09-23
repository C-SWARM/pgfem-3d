/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_macro_client.h"
#include "PGFEM_mpi.h"
#include "pgf_fe2_job.h"
#include "pgf_fe2_server_rebalance.h"
#include "ms_cohe_job_info.h"
#include "ms_cohe_job_list.h"
#include "macro_micro_functions.h"

#include <stdlib.h>
#include <assert.h>

/* fully define the macro client structure */
struct pgf_FE2_macro_client{
  size_t n_jobs_loc;
  size_t n_jobs_glob;
  size_t n_jobs_max;
  size_t n_server;
  MS_COHE_JOB_INFO *jobs;  /* macroscale job(s) */

  /* communication buffers */
  PGFEM_server_ctx *send;
  PGFEM_server_ctx *recv;
};

void pgf_FE2_macro_client_init(pgf_FE2_macro_client **client)
{
  *client = malloc(sizeof(**client));
  (*client)->n_jobs_loc = 0;
  (*client)->n_jobs_glob = 0;
  (*client)->n_jobs_max = 0;
  (*client)->n_server = 0;
  (*client)->jobs = NULL;

  /* should be done by PGFEM_server_ctx init function... */
  (*client)->send = malloc(sizeof(*((*client)->send)));
  (*client)->recv = malloc(sizeof(*((*client)->recv)));
  initialize_PGFEM_server_ctx((*client)->send);
  initialize_PGFEM_server_ctx((*client)->recv);

  /* other initialization stuff */
}

void pgf_FE2_macro_client_destroy(pgf_FE2_macro_client *client)
{
  destroy_PGFEM_server_ctx(client->send);
  destroy_PGFEM_server_ctx(client->recv);
  free(client->send);
  free(client->recv);

  for(int i=0,e=client->n_jobs_loc; i<e; i++){
    destroy_MS_COHE_JOB_INFO((client->jobs) + i);
  }
  free(client->jobs);
  /* destroy internal objects */

  /* destroy the handle */
  free(client);
  client = NULL;
}

void pgf_FE2_macro_client_create_job_list(pgf_FE2_macro_client *client,
					  const int n_jobs_max,
					  const MACROSCALE *macro,
					  const PGFEM_mpi_comm *mpi_comm)
{
  /* compute and store number of servers */
  {
    int n_proc_macro = 0;
    int n_proc_inter = 0;
    MPI_Comm_size(mpi_comm->macro,&n_proc_macro);
    MPI_Comm_size(mpi_comm->mm_inter,&n_proc_inter);
    client->n_server = n_proc_inter - n_proc_macro;
  }

  const COMMON_MACROSCALE *c = macro->common;

  /* get number of jobs and buffer sizes. */
  int n_jobs = 0;
  int *job_buf_sizes = NULL;
  compute_n_job_and_job_sizes(c,&n_jobs,&job_buf_sizes);
  client->n_jobs_loc = n_jobs;

  long Gn_jobs = 0;
  long *n_job_dom = NULL;
  create_group_ms_cohe_job_list(c->nce,c->coel,c->node,
				mpi_comm->macro,MPI_COMM_SELF,
				0,&Gn_jobs,&n_job_dom,
				&(client->jobs));
  /* error check */
  assert(Gn_jobs == n_jobs);

  /* free unused memory */
  free(n_job_dom);

  /* compute total number of jobs and set maximum number of jobs per
     server */
  MPI_Allreduce(MPI_IN_PLACE,&n_jobs,1,MPI_INT,MPI_SUM,mpi_comm->macro);
  client->n_jobs_glob = n_jobs;
  client->n_jobs_max = n_jobs_max;

  /* create server contexts for send/recv of job information. */
  build_PGFEM_server_ctx(client->send,client->n_jobs_loc,job_buf_sizes);
  build_PGFEM_server_ctx(client->recv,client->n_jobs_loc,job_buf_sizes);

  /* do initial server assignment here? */

  /* set tags for jobs */
  {
    const MS_COHE_JOB_INFO *j = client->jobs;
    for(int i=0; i<n_jobs; i++){
      const int tag = pgf_FE2_job_compute_encoded_id(j[i].proc_id,
						     j[i].elem_id,
						     j[i].int_pt);

      PGFEM_server_ctx_set_tag_at_idx(client->send,tag,i);
      PGFEM_server_ctx_set_tag_at_idx(client->recv,tag,i);
    }
  }
  /* final cleanup */
  free(job_buf_sizes);
}

void pgf_FE2_macro_client_assign_initial_servers(pgf_FE2_macro_client *client,
						 const PGFEM_mpi_comm *mpi_comm)
{
  /* create initial partition and send to microscale */
  int nproc_macro = 0;
  const int rank = mpi_comm->rank_macro;
  MPI_Comm_size(mpi_comm->macro,&nproc_macro);
  int *restrict n_jobs = malloc(nproc_macro*sizeof(*n_jobs));
  int *restrict displ = malloc(nproc_macro*sizeof(*n_jobs));
  int *restrict id = malloc(client->n_jobs_glob*sizeof(*id));
  int *restrict proc = malloc(client->n_jobs_glob*sizeof(*proc));
  int *restrict time = malloc(client->n_jobs_glob*sizeof(*time));

  n_jobs[rank] = client->n_jobs_loc;
  MPI_Allgather(MPI_IN_PLACE,1,MPI_INT,n_jobs,1,MPI_INT,mpi_comm->macro);

  /* compute partial sum for displ and set proc_id and time */
  displ[0] = 0;
  for(int i=0; i<n_jobs[0]; i++){
    proc[i] = 0;
    time[i] = 1;
  }
  for(int i=1; i<nproc_macro; i++){
    displ[i] = displ[i-1] + n_jobs[i-1];
    const int start = displ[i];
    for(int j=0; j<n_jobs[i]; j++){
      proc[start + j] = i;
      time[start + j] = 1;
    }
  }

  /* set job ids */
  const int start = displ[rank];
  memcpy(id+start,client->send->tags,sizeof(*id));

  /* get job ids from other procs */
  MPI_Allgatherv(MPI_IN_PLACE,n_jobs[rank],MPI_INT,
		 id,n_jobs,displ,MPI_INT,mpi_comm->macro);

  void *all_part = NULL;
  void *parts = NULL;

  new_partition_build_set_keep(&all_part,client->n_jobs_max,
			       client->n_jobs_glob,id,time,proc);
  new_partitions_void(&parts,client->n_server,client->n_jobs_max);

  /* partition the jobs using the greedy algorithm. Since everything
     has the same weight, amounts to round robin assignment. */
  rebalance_partitions_greedy(client->n_server,all_part,parts);

  /* push partitions into structure I will send */
  pgf_FE2_server_rebalance *rb = malloc(client->n_server*sizeof(*rb));
  new_partitions_void_to_pgf_FE2_server_rebalance(client->n_server,parts,rb);

  /* cleanup */
  new_partition_destroy_void(all_part);
  new_partitions_destroy_void(parts,client->n_server);

  /* final cleanup */
  for(size_t i=0, e=client->n_server; i<e; i++){
    pgf_FE2_server_rebalance_destroy(rb + i);
  }
  free(rb);
}

void pgf_FE2_macro_client_rebalance_servers(pgf_FE2_macro_client *client
					    /* TBD */)
{
  /* receive message and rebalance according to flag (NONE or GREEDY) */
}

void pgf_FE2_macro_client_send_jobs(pgf_FE2_macro_client *client
				    /* TBD */)
{

}

void pgf_FE2_macro_client_recv_jobs(pgf_FE2_macro_client *client
				    /* TBD */)
{

}

void pgf_FE2_macro_client_send_exit(pgf_FE2_macro_client *client
				     /* TBD */)
{

}

