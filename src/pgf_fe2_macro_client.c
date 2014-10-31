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
#include "PLoc_Sparse.h"
#include "stiffmat_fd.h"

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
  static const size_t n = 2;
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

  /* for broadcasting */
  /* When using this datastructure to communicate with servers, set
     active to true before 1st communication and set to 0 after
     communication compled by MPI_Test/Wait functions. */
  struct BCAST{
    size_t active;
    size_t n_comm;
    int *ranks;
    MPI_Request *req;
  } bcast;
};

static int pgf_FE2_macro_client_update_send_recv(pgf_FE2_macro_client *client,
						 pgf_FE2_server_rebalance **rb_list,
						 const size_t nproc_macro)
{
  int err = 0;
  assert(nproc_macro > 0);

  /* Get list of tags from ctx. Put into tuple for {tag, idx, server} */
  const size_t len = 3;
  const size_t s_tags_nel = client->send->n_comms;

  /* exit early if this domain does not contain jobs */
  if(s_tags_nel == 0) return err;

  int *s_tags = malloc(s_tags_nel*len*sizeof(*s_tags));
  for(size_t i = 0; i < s_tags_nel; i++){
    s_tags[i*len] = client->send->tags[i];
    s_tags[i*len + 1] = i;
    s_tags[i*len + 2] = 0; /* nproc_macro always > 0 */
  }

  /* sort s_tags by {tag} */
  qsort(s_tags,s_tags_nel,len*sizeof(*s_tags),compare_first_int);

  /* loop through each rebalance and extract list of jobs on the
     server. */
  for(size_t i = 0, e = client->n_server; i<e; i++){
    const size_t n_keep = pgf_FE2_server_rebalance_n_keep((rb_list[i]));
    const size_t n_recv = pgf_FE2_server_rebalance_n_recv((rb_list[i]));
    const size_t serv_tags_len = (n_keep + n_recv);
    int *serv_tags = malloc(serv_tags_len * sizeof(*serv_tags));
    const int *k = pgf_FE2_server_rebalance_keep_buf((rb_list[i]));
    const int *r = pgf_FE2_server_rebalance_recv_buf((rb_list[i]));
    memcpy(serv_tags,k,n_keep*sizeof(*serv_tags));
    memcpy(serv_tags + n_keep,r,n_recv*sizeof(*serv_tags));

    /* sort by job id */
    qsort(serv_tags,serv_tags_len,sizeof(*serv_tags),compare_first_int);

    for(size_t j = 0; j < s_tags_nel; j++){
      /* do not search for match if already found */
      if(s_tags[j*len + 2] >= nproc_macro) continue;

      void *ptr = bsearch(s_tags + j*len,serv_tags,serv_tags_len,
			  sizeof(*serv_tags),compare_first_int);

      /* Check for match, set src/dest proc id */
      if(ptr != NULL){
	s_tags[j*len + 2] = i + nproc_macro;
      }
    }
    free(serv_tags);
  }

#ifndef NDEBUG
  /* sanity check/debugging. ensure that the src/dest is valid */
  /* sort s_tags by {server} */
  qsort(s_tags,s_tags_nel,len*sizeof(*s_tags),compare_third_int);

  /* ensure that the smallest server id is valid, i.e., we assigned all
     jobs on this domain. */
  assert(s_tags[2] >= nproc_macro);
  assert(s_tags[s_tags_nel*len - 1] < (nproc_macro + client->n_server));
#endif

  /* sort s_tags by {idx} */
  qsort(s_tags,s_tags_nel,len*sizeof(*s_tags),compare_second_int);

  /* assign proc on send and recv */
  for(size_t i = 0; i < s_tags_nel; i++){
    const int tag = s_tags[i*len];
    const int idx = s_tags[i*len + 1];
    const int proc = s_tags[i*len + 2];

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
 * Initializes client->bcast.
 */
static int pgf_FE2_macro_client_bcast_list(pgf_FE2_macro_client *client,
					   const PGFEM_mpi_comm *mpi_comm)
{
  int err = 0;
  const size_t rank = mpi_comm->rank_macro;
  const size_t n_server = client->n_server;

  /* exit early if not responsible for signaling servers */
  if(rank >= n_server) return err;

  int n_comm = 0;
  int nproc_macro = 0;
  int nproc_inter = 0;
  err += MPI_Comm_size(mpi_comm->macro,&nproc_macro);
  err += MPI_Comm_size(mpi_comm->mm_inter,&nproc_inter);

  /* determine what ranks (in mm_inter) this domain is responsible
     for broadcasting info to. */
  if(nproc_macro >= n_server) n_comm = 1;
  else {
    n_comm = n_server/nproc_macro;
    int rem = n_server % nproc_macro;
    if(rem > 0 && rem > rank) n_comm++;
  }

  /* allocate/populate client->bcast */
  client->bcast.n_comm = n_comm;
  client->bcast.ranks = malloc(n_comm*sizeof(*(client->bcast.ranks)));
  client->bcast.req = calloc(n_comm,sizeof(*(client->bcast.req)));
  {
    int * restrict r = client->bcast.ranks; /* restrict alias */
    r[0] = rank + nproc_macro;
    for(size_t i=1; i<n_comm; i++){
      r[i] = r[i-1] + nproc_macro;
    }
  }

  /* check max rank for validity */
  assert(client->bcast.ranks[n_comm - 1] < nproc_inter);

  return err;
}
					   

/**
 * Send the rebalancing information to the servers.
 *
 * Rebalancing information sent from servers in round-robin. Not all
 * macro-clients may perform communication.
 */
static int pgf_FE2_macro_client_bcast_rebal_to_servers(pgf_FE2_macro_client *client,
						       pgf_FE2_server_rebalance **rb_list,
						       const PGFEM_mpi_comm *mpi_comm)
{
  int err = 0;

  /* get aliases */
  const size_t n_send = client->bcast.n_comm;
  const int *ranks = client->bcast.ranks;
  MPI_Request *req = client->bcast.req;

  /* should put request et al. in local buffer for overlay of
     comp/comm but will just waitall for now. */
  int nproc_macro = 0;
  const int rank = mpi_comm->rank_macro;
  err += MPI_Comm_size(mpi_comm->macro,&nproc_macro);

  client->bcast.active = 1;

  /* hold addresses of rebal buffers to send */
  void **buffs = calloc(n_send,sizeof(*buffs));

  for(size_t i = 0; i < n_send; i++){
    size_t idx = rank + i*nproc_macro;
    assert(idx < client->n_server);
    size_t len = pgf_FE2_server_rebalance_n_bytes(rb_list[idx]);
    buffs[i] = pgf_FE2_server_rebalance_buff(rb_list[idx]);
    err += MPI_Isend(buffs[i],len,MPI_CHAR,ranks[i],
		     FE2_MICRO_SERVER_REBALANCE,
		     mpi_comm->mm_inter,req+i);
  }

  err += MPI_Waitall(n_send,req,MPI_STATUS_IGNORE);
  free(buffs);
  client->bcast.active = 0;

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

  /* bcast */
  (*client)->bcast.active = 0;
  (*client)->bcast.n_comm = 0;
  (*client)->bcast.ranks = NULL;
  (*client)->bcast.req = NULL;

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
  if(client->bcast.active){
    /* complete communication before destruction */
    MPI_Waitall(client->bcast.n_comm,client->bcast.req,MPI_STATUS_IGNORE);
  }
  free(client->bcast.ranks);
  free(client->bcast.req);

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
  MPI_Allreduce(MPI_IN_PLACE,&Gn_jobs,1,MPI_INT,MPI_SUM,mpi_comm->macro);
  client->n_jobs_glob = Gn_jobs;
  client->n_jobs_max = n_jobs_max;

  /* make sure sufficient resources to compute solution */
  if(client->n_server * client->n_jobs_max < Gn_jobs){
    PGFEM_printerr("ERROR: insufficent microscale resources!\n"
		   "Increase n_jobs_max or number of servers and try again!\n");
    PGFEM_Abort();
  }

  /* create server contexts for send/recv of job information. */
  build_PGFEM_server_ctx(client->send,client->n_jobs_loc,job_buf_sizes);
  build_PGFEM_server_ctx(client->recv,client->n_jobs_loc,job_buf_sizes);

  /* setup the bcast list */
  pgf_FE2_macro_client_bcast_list(client,mpi_comm);

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
  const size_t max_proc = client->n_server;
  size_t proc_id = 0;
  for(int i=0; i<n_jobs[0]; i++){
    proc[i] = proc_id++;
    if(proc_id >= max_proc) proc_id = 0;
    time[i] = 1;
  }
  for(int i=1; i<nproc_macro; i++){
    displ[i] = displ[i-1] + n_jobs[i-1];
    const int start = displ[i];
    for(int j=0; j<n_jobs[i]; j++){
      proc[start + j] = proc_id++;
      if(proc_id >= max_proc) proc_id = 0;
      time[start + j] = 1;
    }
  }

  /* set job ids */
  const int start = displ[rank];
  memcpy(id+start,client->send->tags,client->send->n_comms*sizeof(*id));

  /* get job ids from other procs */
  MPI_Allgatherv(MPI_IN_PLACE,n_jobs[rank],MPI_INT,
		 id,n_jobs,displ,MPI_INT,mpi_comm->macro);

  void *all_part = NULL;
  void *parts = NULL;

  new_partition_build_set_keep(&all_part,client->n_jobs_glob,
			       client->n_jobs_glob,id,time,proc);
  new_partitions_void(&parts,client->n_server,client->n_jobs_max);

  /* partition the jobs using the greedy algorithm. Since everything
     has the same weight, amounts to round robin assignment. */
  rebalance_partitions_greedy(client->n_server,all_part,parts);

  /* push partitions into structure I will send */
  pgf_FE2_server_rebalance **rb = calloc(client->n_server,sizeof(*rb));
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
    pgf_FE2_server_rebalance_destroy(rb[i]);
  }
  free(rb);
}

void pgf_FE2_macro_client_rebalance_servers(pgf_FE2_macro_client *client,
					    const PGFEM_mpi_comm *mpi_comm,
					    const int heuristic)
{
  int nproc_macro = 0;
  MPI_Comm_size(mpi_comm->macro,&nproc_macro);

  /* receive message and rebalance according to heuristic */
  pgf_FE2_server_rebalance **rb_list = pgf_FE2_rebalancer(mpi_comm,
							 client->n_jobs_glob,
							 client->n_jobs_max,
							 heuristic);

  /* update the server context (send/recv) */
  pgf_FE2_macro_client_update_send_recv(client,rb_list,nproc_macro);

  /* Broadcast rebalancing info to servers */
  pgf_FE2_macro_client_bcast_rebal_to_servers(client,rb_list,mpi_comm);

  /* cleanup */
  for(size_t i=0, e=client->n_server; i<e; i++){
    pgf_FE2_server_rebalance_destroy(rb_list[i]);
  }
  free(rb_list);
}

void pgf_FE2_macro_client_send_jobs(pgf_FE2_macro_client *client,
				    const PGFEM_mpi_comm *mpi_comm,
				    const MACROSCALE *macro,
				    const int job_type)
{
  int err = 0;
  /* see start_macroscale_compute_jobs */

  /* Get aliases from client object etc. */
  PGFEM_server_ctx *send = client->send;
  PGFEM_server_ctx *recv = client->recv;
  MS_COHE_JOB_INFO *job_list = client->jobs;

  MPI_Comm comm = mpi_comm->mm_inter;
  /* Wait to complete any pending communcication. The error flag is
     incremented as this should only occur by a logical/programming
     error. Note that buffers may be overwritten. If the communication
     is completed but the flags have not been reset to 0, the
     MPI_Wait* commands should return immediately */
  if(send->in_process || recv->in_process){
    PGFEM_printerr("WARNING: communication in progress!(%s:%s:%d)\n",
		   __func__,__FILE__,__LINE__);
    err++;
    err += MPI_Waitall(recv->n_comms,recv->req,recv->stat);
    err += MPI_Waitall(send->n_comms,send->req,send->stat);
    recv->in_process = 0;
    send->in_process = 0;
  }


  /* post recieves (from the running server) */
  for(int i=0; i<recv->n_comms; i++){
    err += MPI_Irecv(recv->buffer[i],recv->sizes[i],MPI_CHAR,
		     recv->procs[i],recv->tags[i],comm,
		     &(recv->req[i]));
  }

  for(int i=0; i<send->n_comms; i++){
    /* update the job information according to job_type */
    err += macroscale_update_job_info(macro,job_type,macro->sol->f,job_list+i);

    /* pack the job info into the buffer to send */
    err += pack_MS_COHE_JOB_INFO(job_list + i,send->sizes[i],
				 send->buffer[i]);

    /* post the send (to the running server) */
    err += MPI_Isend(send->buffer[i],send->sizes[i],MPI_CHAR,
		     send->procs[i],send->tags[i],comm,
		     &(send->req[i]));
  }

  /* set in_process flags to true (1) */
  recv->in_process = 1;
  send->in_process = 1;

}

void pgf_FE2_macro_client_recv_jobs(pgf_FE2_macro_client *client,
				    MACROSCALE *macro,
				    int *max_micro_sub_step)
{
  int err = 0;

  /* Get aliases from client object etc. */
  PGFEM_server_ctx *send = client->send;
  PGFEM_server_ctx *recv = client->recv;
  MS_COHE_JOB_INFO *job_list = client->jobs;

  /* reset the max number of steps */
  *max_micro_sub_step = 0;

  /* see finish_macroscale_compute_jobs */
  COMMON_MACROSCALE *c = macro->common;
  MACROSCALE_SOLUTION *s = macro->sol;
  int rank_macro = 0;
  int nproc_macro = 0;

  /* exit early if !*->in_process */
  if(!recv->in_process && !send->in_process) return;

  err += MPI_Comm_rank(c->mpi_comm,&rank_macro);
  err += MPI_Comm_size(c->mpi_comm,&nproc_macro);

  /* if expecting to receive buffers */
  if(recv->in_process){
    /* set up the stiffness matrix communication */
    double **Lk = NULL;
    double **receive = NULL;
    MPI_Request *req_r = NULL;
    MPI_Status *sta_r = NULL;
    err += init_and_post_stiffmat_comm(&Lk,&receive,&req_r,&sta_r,
				       c->mpi_comm,c->pgfem_comm);

    /* assemble jobs as they are received */
    for(int i=0; i<recv->n_comms; i++){
      int idx = 0;
      err += MPI_Waitany(recv->n_comms,recv->req,&idx,recv->stat);
      MS_COHE_JOB_INFO *job = job_list + idx;
      err += unpack_MS_COHE_JOB_INFO(job,recv->sizes[idx],
				     recv->buffer[idx]);

      /* compute the number of micro-sub-steps */
      if(job->n_step > *max_micro_sub_step) *max_micro_sub_step = job->n_step;

      /* finish jobs based on job_type */
      switch(job->job_type){
      case JOB_COMPUTE_EQUILIBRIUM:
	/* assemble residual (local) */
	for(int j=0; j<job->ndofe; j++){
	  int dof_id = job->loc_dof_ids[j] - 1;
	  if(dof_id < 0) continue; /* boundary condition */
	  s->f_u[dof_id] += job->traction_res[j];
	}

	/*** Deliberate drop through ***/
      case JOB_NO_COMPUTE_EQUILIBRIUM:
	/* assemble tangent to local and off-proc buffers */
	PLoc_Sparse(NULL,Lk,job->K_00_contrib,NULL,NULL,NULL,job->g_dof_ids,
		    job->ndofe,NULL,c->GDof,rank_macro,nproc_macro,
		    c->pgfem_comm,0,c->SOLVER,macro->opts->analysis_type);
	break;
      case JOB_UPDATE:
	/* update cohesive elements */
	err += macroscale_update_coel(job,macro);
	break;
      default: /* do nothing */ break;
      }
    }

    /* send/finalize communication of the stiffness matrix */
    MPI_Status *sta_s = NULL;
    MPI_Request *req_s = NULL;
    err += send_stiffmat_comm(&sta_s,&req_s,Lk,c->mpi_comm,c->pgfem_comm);

    /* get maximum number of steps from all macro processors */
    err += MPI_Allreduce(MPI_IN_PLACE,max_micro_sub_step,
			 1,MPI_INT,MPI_MAX,c->mpi_comm);

    err += assemble_nonlocal_stiffmat(c->pgfem_comm,sta_r,req_r,
				      c->SOLVER,receive);

    err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,c->pgfem_comm);

    /* re-initialize preconditioner ? */
    /* err += PGFEM_HYPRE_create_preconditioner(c->SOLVER,c->mpi_comm); */

    /* clean up memory */
    for(int i=0; i<nproc_macro; i++){
      if(Lk != NULL) free(Lk[i]);
      if(receive != NULL) free(receive[i]);
    }
    free(Lk);
    free(receive);
    free(sta_r);
    free(req_r);
    free(sta_s);
    free(req_s);

    /* set in_process to false (0) */
    recv->in_process = 0;
  }/* end if(recv->in_process) */

  /* wait for any pending send communication to finish */
  if(send->in_process){
    MPI_Waitall(send->n_comms,send->req,send->stat);
    send->in_process = 0;
  }

}

void pgf_FE2_macro_client_send_exit(pgf_FE2_macro_client *client,
				    const PGFEM_mpi_comm *mpi_comm)
{
  const size_t n_send = client->bcast.n_comm;
  const int *ranks = client->bcast.ranks;
  MPI_Request *req = client->bcast.req;

  client->bcast.active = 1;
  void *empty = NULL;
  for(size_t i=0; i<n_send; i++){
    MPI_Isend(empty,0,MPI_CHAR,ranks[i],FE2_MICRO_SERVER_EXIT,
	      mpi_comm->mm_inter,req + i);
  }

  /* wait for exit code to be received */
  MPI_Waitall(n_send,req,MPI_STATUS_IGNORE);
  client->bcast.active = 0;
}

