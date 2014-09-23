/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_macro_client.h"
#include "PGFEM_mpi.h"
#include "pgf_fe2_job.h"
#include "ms_cohe_job_list.h"
#include "macro_micro_functions.h"

#include <stdlib.h>
#include <assert.h>

/* fully define the macro client structure */
struct pgf_FE2_macro_client{
  size_t n_jobs;
  MS_COHE_JOB_INFO *jobs;  /* macroscale job(s) */

  /* communication buffers */
  PGFEM_server_ctx *send;
  PGFEM_server_ctx *recv;
};

void pgf_FE2_macro_client_init(pgf_FE2_macro_client **client)
{
  *client = malloc(sizeof(**client));
  (*client)->n_jobs = 0;
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

  for(int i=0,e=client->n_jobs; i<e; i++){
    destroy_MS_COHE_JOB_INFO((client->jobs) + i);
  }
  free(client->jobs);
  /* destroy internal objects */

  /* destroy the handle */
  free(client);
  client = NULL;
}

void pgf_FE2_macro_client_create_job_list(pgf_FE2_macro_client *client,
					  const MACROSCALE *macro,
					  const PGFEM_mpi_comm *mpi_comm)
{
  const COMMON_MACROSCALE *c = macro->common;

  /* get number of jobs and buffer sizes. */
  int n_jobs = 0;
  int *job_buf_sizes = NULL;
  compute_n_job_and_job_sizes(c,&n_jobs,&job_buf_sizes);
  client->n_jobs = n_jobs;

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

  /* create server contexts for send/recv of job information. */
  build_PGFEM_server_ctx(client->send,client->n_jobs,job_buf_sizes);
  build_PGFEM_server_ctx(client->recv,client->n_jobs,job_buf_sizes);

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

void pgf_FE2_macro_client_assign_initial_servers(pgf_FE2_macro_client *client
						 /* TBD */)
{

}

void pgf_FE2_macro_client_rebalance_servers(pgf_FE2_macro_client *client
					    /* TBD */)
{

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

