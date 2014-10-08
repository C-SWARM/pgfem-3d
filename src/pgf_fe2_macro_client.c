/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_macro_client.h"
#include "pgf_fe2_micro_server.h"
#include "PGFEM_mpi.h"
#include "pgf_fe2_job.h"
#include "pgf_fe2_server_rebalance.h"
#include "ms_cohe_job_info.h"
#include "ms_cohe_job_list.h"
#include "macro_micro_functions.h"

#include <stdlib.h>
#include <assert.h>

/**
 * Comparator function for array of ints as object. Compares n-th in
 * array object.
 */
static int compare_nth_int(const void *n,
			   const void *a,
			   const void *b)
{
  const size_t N = *(const size_t*) n;
  return ((int*) a)[N] - ((int*) b)[N];
}


/**
 * Comparator function for array of ints as object. Compares first in
 * array object.
 */
static int compare_first_int(const void *a,
			     const void *b)
{
  static const size_t n = 0;
  return compare_nth_int(&n,a,b);
}

/**
 * Comparator function for array of ints as object. Compares second in
 * array object.
 */
static int compare_second_int(const void *a,
			      const void *b)
{
  static const size_t n = 1;
  return compare_nth_int(&n,a,b);
}

/**
 * Comparator function for array of ints as object. Compares third in
 * array object.
 */
static int compare_third_int(const void *a,
			     const void *b)
{
  static const int n = 2;
  return compare_nth_int(&n,a,b);
}

/* fully define the macro client structure */
struct pgf_FE2_macro_client{
  size_t n_jobs_loc;      /* how many jobs are on THIS macro-domain */
  size_t n_jobs_glob;     /* how many jobs are on ALL macro-domains */
  size_t n_jobs_max;      /* maximum number of jobs on each server */
  size_t n_server;        /* number of servers */
  MS_COHE_JOB_INFO *jobs; /* macroscale job(s) */

  /* communication buffers */
  PGFEM_server_ctx *send;
  PGFEM_server_ctx *recv;
};

static int pgf_FE2_macro_client_update_send_recv(pgf_FE2_macro_client *client,
						 const pgf_FE2_server_rebalance *rb_list,
						 const size_t nproc_macro)
{
  int err = 0;

  /* Get list of tags from ctx. Put into tuple for {tag, idx, server} */
  const size_t len = 3;
  const size_t s_tags_len = len*client->send->n_comms;
  int *s_tags = malloc(s_tags_len*sizeof(*s_tags));
  for(size_t i = 0; i < s_tags_len; i += len){
    s_tags[i] = client->send->tags[i%len];
    s_tags[i+1] = i%len;
    s_tags[i+2] = -1;
  }

  /* sort tags by {tag} */
  qsort(s_tags,client->send->n_comms,len*sizeof(*s_tags),compare_first_int);

  /* loop through each rebalance and extract list of jobs on the
     server. */
  for(size_t i = 0, e = client->n_server; i<e; i++){
    const size_t n_keep = pgf_FE2_server_rebalance_n_keep(&(rb_list[i]));
    const size_t n_recv = pgf_FE2_server_rebalance_n_recv(&(rb_list[i]));
    const size_t serv_tags_len = (n_keep + n_recv);
    int *serv_tags = malloc(serv_tags_len * sizeof(*serv_tags));
    const int *k = pgf_FE2_server_rebalance_keep_buf(&(rb_list[i]));
    const int *r = pgf_FE2_server_rebalance_recv_buf(&(rb_list[i]));
    memcpy(serv_tags,k,n_keep*sizeof(*serv_tags));
    memcpy(serv_tags + n_keep,r,n_recv*sizeof(*serv_tags));

    /* sort by job id */
    qsort(serv_tags,serv_tags_len,sizeof(*serv_tags),compare_first_int);

    for(size_t j = 0; j < s_tags_len; j += len){
      /* do not search for match if already found */
      if(s_tags[j+2] >= nproc_macro) continue;

      void *ptr = bsearch(s_tags + j,serv_tags,serv_tags_len,
			  sizeof(*serv_tags),compare_first_int);

      /* Check for match, set src/dest proc id */
      if(ptr != NULL){
	s_tags[j+2] = i + nproc_macro;
      }
    }
    free(serv_tags);
  }

#ifndef NDEBUG
  /* sanity check/debugging. ensure that the src/dest is valid */
  /* sort s_tags by {server} */
  qsort(s_tags,client->send->n_comms,len*sizeof(*s_tags),compare_third_int);

  /* ensure that the smallest server id is valid, i.e., we assigned all
     jobs on this domain. */
  assert(s_tags[2] >= nproc_macro);
  assert(s_tags[s_tags_len - 1] < nproc_macro + client->n_server);
#endif

  /* sort s_tags by {idx} */
  qsort(s_tags,client->send->n_comms,len*sizeof(*s_tags),compare_second_int);

  /* assign proc on send and recv */
  for(size_t i = 0; i < s_tags_len; i += len){
    const int tag = s_tags[i];
    const int idx = s_tags[i+1];
    const int proc = s_tags[i+2];

    /* error checking. Compiled out with NDEBUG */
    assert(client->send->tags[idx] == tag);
    assert(client->recv->tags[idx] == tag);

    /* set procs for send and recv */
    PGFEM_server_ctx_set_proc_at_idx(client->send,proc,idx);
    PGFEM_server_ctx_set_proc_at_idx(client->recv,proc,idx);
  }

  free(s_tags);

  return err;
}

/**
 * Determine what servers this domain is responsible for broadcasting
 * information to.
 *
 * Returns number of messages and allocated list of process ids on
 * mpi_comm->mm_inter.
 */
static int pgf_FE2_macro_client_bcast_list(const PGFEM_mpi_comm *mpi_comm,
					   size_t *n_comm,
					   int **ranks)
{
  int err = 0;
  const size_t rank = mpi_comm->rank_macro;
  int nproc_macro = 0;
  int nproc_inter = 0;
  int n_server = 0;
  err += MPI_Comm_size(mpi_comm->macro,&nproc_macro);
  err += MPI_Comm_size(mpi_comm->mm_inter,&nproc_inter);
  n_server = nproc_inter - nproc_macro;

  if(rank >= n_server){
    /* this domain is not responsible for broadcasting any information */
    *n_comm = 0;
    *ranks = NULL;
  } else {
    /* determine what ranks (in mm_inter) this domain is responsible
       for broadcasting info to. */
    if(nproc_macro >= n_server) *n_comm = 1;
    else {
      *n_comm = n_server/nproc_macro;
      int rem = n_server % nproc_macro;
      if(rem > 0 && rem > rank) (*n_comm)++;
    }

    /* allocate/populate ranks  */
    *ranks = malloc((*n_comm)*sizeof(**ranks));
    {
      int * restrict r = *ranks; /* restrict alias */
      r[0] = rank;
      for(size_t i=1, e=*n_comm; i<e; i++){
	r[i] = r[i-1] + nproc_macro;
      }
    }

    /* check max rank for validity */
    assert((*ranks)[*n_comm - 1] < nproc_inter);
  }

  return err;
}
					   

/**
 * Send the rebalancing information to the servers.
 *
 * Rebalancing information sent from servers in round-robin. Not all
 * macro-clients may perform communication.
 */
static int pgf_FE2_macro_client_bcast_rebal_to_servers(pgf_FE2_macro_client *client,
						       const pgf_FE2_server_rebalance *rb_list,
						       const PGFEM_mpi_comm *mpi_comm)
{
  int err = 0;

  size_t n_send = 0;
  int *ranks = NULL;
  err += pgf_FE2_macro_client_bcast_list(mpi_comm,&n_send,&ranks);

  /* should put request et al. in local buffer for overlay of
     comp/comm but will just waitall for now. */
  int nproc_macro = 0;
  const int rank = mpi_comm->rank_macro;
  err += MPI_Comm_size(mpi_comm->macro,&nproc_macro);
  MPI_Request *req = malloc(n_send*sizeof(*req));
  for(size_t i = 0; i < n_send; i++){
    size_t idx = rank + i*nproc_macro;
    assert(idx < client->n_server);
    size_t len = pgf_FE2_server_rebalance_n_bytes(rb_list + idx);
    err += MPI_Isend(rb_list + idx,len,MPI_CHAR,ranks[i],
		     FE2_MICRO_SERVER_REBALANCE,
		     mpi_comm->mm_inter,req+i);
  }

  err += MPI_Waitall(n_send,req,MPI_STATUS_IGNORE);

  free(req);
  free(ranks);
  return err;
}

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

  /* Update server context (who I am sending to)*/
  pgf_FE2_macro_client_update_send_recv(client,rb,nproc_macro);

  /* Broadcast the rebalancing information to the servers */
  pgf_FE2_macro_client_bcast_rebal_to_servers(client,rb,mpi_comm);

  /* final cleanup */
  for(size_t i=0, e=client->n_server; i<e; i++){
    pgf_FE2_server_rebalance_destroy(rb + i);
  }
  free(rb);
}

void pgf_FE2_macro_client_rebalance_servers(pgf_FE2_macro_client *client,
					    const PGFEM_mpi_comm *mpi_comm,
					    const int heuristic)
{
  int nproc_macro = 0;
  MPI_Comm_rank(mpi_comm->macro,&nproc_macro);

  /* receive message and rebalance according to heuristic */
  pgf_FE2_server_rebalance *rb_list = pgf_FE2_rebalancer(mpi_comm,
							 client->n_jobs_glob,
							 client->n_jobs_max,
							 heuristic);

  /* update the server context (send/recv) */
  pgf_FE2_macro_client_update_send_recv(client,rb_list,nproc_macro);

  /* Broadcast rebalancing info to servers */
  pgf_FE2_macro_client_bcast_rebal_to_servers(client,rb_list,mpi_comm);

  /* cleanup */
  for(size_t i=0, e=client->n_server; i<e; i++){
    pgf_FE2_server_rebalance_destroy(rb_list + i);
  }
  free(rb_list);
}

void pgf_FE2_macro_client_send_jobs(pgf_FE2_macro_client *client
				    /* TBD */)
{
  /* see start_macroscale_compute_jobs */
}

void pgf_FE2_macro_client_recv_jobs(pgf_FE2_macro_client *client
				    /* TBD */)
{
  /* see finish_macroscale_compute_jobs */
}

void pgf_FE2_macro_client_send_exit(pgf_FE2_macro_client *client,
				    const PGFEM_mpi_comm *mpi_comm)
{
  size_t n_send = 0;
  int *ranks = NULL;
  void *empty = NULL;
  pgf_FE2_macro_client_bcast_list(mpi_comm,&n_send,&ranks);

  MPI_Request *req = malloc(n_send*sizeof(*req));
  for(size_t i=0; i<n_send; i++){
    MPI_Isend(empty,0,MPI_CHAR,ranks[i],FE2_MICRO_SERVER_EXIT,
	      mpi_comm->mm_inter,req + i);
  }

  /* wait for exit code to be received */
  MPI_Waitall(n_send,req,MPI_STATUS_IGNORE);
  free(req);
  free(ranks);
}

