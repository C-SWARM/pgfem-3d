/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/MultiscaleCommon.hpp"
#include "allocation.h"
#include "ms_cohe_job_info.h"
#include "ms_cohe_job_list.h"
#include "macro_micro_functions.h"
#include "pgf_fe2_job.h"
#include "pgf_fe2_macro_client.h"
#include "pgf_fe2_micro_server.h"
#include "pgf_fe2_server_rebalance.h"
#include "PLoc_Sparse.h"
#include "stiffmat_fd.h"

#include <stdlib.h>
#include <assert.h>

using namespace pgfem3d;
using namespace pgfem3d::net;

/**
 * Comparator function for array of ints as object. Compares n-th in
 * array object.
 */
static int
compare_nth_int(const void *n, const void *a, const void *b)
{
  const size_t N = *(const size_t*) n;
  return ((int*) a)[N] - ((int*) b)[N];
}


/**
 * Comparator function for array of ints as object. Compares first in
 * array object.
 */
static int
compare_first_int(const void *a, const void *b)
{
  static const size_t n = 0;
  return compare_nth_int(&n,a,b);
}

/**
 * Comparator function for array of ints as object. Compares second in
 * array object.
 */
static int
compare_second_int(const void *a, const void *b)
{
  static const size_t n = 1;
  return compare_nth_int(&n,a,b);
}

/* fully define the macro client structure */
struct pgf_FE2_macro_client{
  size_t n_jobs_loc;      /* how many jobs are on THIS macro-domain */
  size_t n_jobs_loc_ROM;
  size_t n_jobs_glob;     /* how many jobs are on ALL macro-domains */
  size_t n_jobs_glob_ROM;
  size_t n_jobs_max;      /* maximum number of jobs on each server */
  size_t n_server;        /* number of servers */
  size_t n_server_ROM;      //number of micro 2 servers
  MS_COHE_JOB_INFO *jobs; /* macroscale job(s) (pde)*/
  MS_COHE_JOB_INFO *jobs_ROM; /* macroscale job(s) (rom)*/

  /* communication buffers */
  MultiscaleServerContext *send;
  MultiscaleServerContext *recv;
  MultiscaleServerContext *send_ROM;
  MultiscaleServerContext *recv_ROM;
  /* for broadcasting */
  /* When using this datastructure to communicate with servers, set
     active to true before 1st communication and set to 0 after
     communication completed by Test/Wait functions. */
  struct BCAST{
    size_t active;
    size_t n_comm;
    int *ranks;
    Request *req;
  } bcast;

  struct BCAST_ROM : BCAST{ //there is probably a better way to do this
    size_t active;
    size_t n_comm;
    int *ranks;
    Request *req;
  } bcast_ROM;
 
  Network *net; // handle to the network
};

/**
 * Determine the communication pattern between macroscale clients and
 * microscale servers.
 *
 * Only the domains with multiscale elements actually exchange
 * information w/ microscale servers (other than broadcasting which
 * rebalancing information). This function determines for the
 * macroscale the location of cells on the servers after rebalancing
 * and thus where to send/receive data to/from. The servers
 * dynamically compute the src/dest of communications from the cell
 * ID, chich encodes the macroscopic domain.
 *
 * \param[in,out] client Macroscopic domain
 * \param[in] rb_list Rebalancing information for each server
 * \param[in] nproc_macro Number of macroscale domains
 * \return non-zero on error.
 */
static int pgf_FE2_macro_client_update_send_recv(pgf_FE2_macro_client *client,
                         pgf_FE2_server_rebalance **rb_list,
                         const size_t nproc_macro, int micro_model) //1 = pde, 2 = taylor
{
  int err = 0;
  assert(nproc_macro > 0);
    auto recv = client->recv; //default do nothing
    auto send = client->send; //if ROM, then overwrite
    auto n_server = client->n_server;

  if (micro_model == 2) {//will eventually be replaced by an array
    recv = client->recv_ROM;
    send = client->send_ROM;
    n_server = client->n_server_ROM;
  }
  /*
   * Get number of communications (number of multiscale elements on
   * domain). Can exit early if there are no multiscale elements,
   * i.e., no required communication.
   */
  const size_t s_tags_nel = send->n_comms;
  if(s_tags_nel == 0) return err;

  /* push information into tuple workspace {tag, idx, server} */
  const size_t len = 3;
  int *s_tags = PGFEM_malloc<int>(s_tags_nel*len);
  for(size_t i = 0; i < s_tags_nel; i++){
    s_tags[i*len] = send->tags[i];
    s_tags[i*len + 1] = i;
    s_tags[i*len + 2] = 0; /* nproc_macro always > 0 */
  }

  /* sort s_tags by {tag} */
  qsort(s_tags,s_tags_nel,len*sizeof(*s_tags),compare_first_int);

  /*
   * loop through each server and extract list of jobs on the server
   * after rebalancing (kept + received).
   */
  unsigned n_client_matched = 0;
  for(size_t i = 0, e = n_server; i<e; i++){
    const size_t n_keep = pgf_FE2_server_rebalance_n_keep((rb_list[i]));
    const size_t n_recv = pgf_FE2_server_rebalance_n_recv((rb_list[i]));
    const int *k = pgf_FE2_server_rebalance_keep_buf((rb_list[i]));
    const int *r = pgf_FE2_server_rebalance_recv_buf((rb_list[i]));

    /*
     * serv_tags contains the tags (cell ids) on server i after
     * rebalancing
     */
    const size_t serv_tags_len = (n_keep + n_recv);
    int *serv_tags = PGFEM_malloc<int>(serv_tags_len);
    memcpy(serv_tags,k,n_keep*sizeof(*serv_tags));
    memcpy(serv_tags + n_keep,r,n_recv*sizeof(*serv_tags));

    /* sort by job id */
    qsort(serv_tags,serv_tags_len,sizeof(*serv_tags),compare_first_int);

    /*
     * Determine if server i contains tag (cell) j. If so, mark it
     * with the server rank (in mpi_comm->mm_inter, rank = i +
     * nproc_macro) to communicate with.
     */
    unsigned n_server_matched = 0;
    for(size_t j = 0; j < s_tags_nel; j++){
      assert(s_tags[j*len + 2] >= 0);
      if ((unsigned)s_tags[j*len + 2] >= nproc_macro) {
    /* already matched */
    continue;
      } else if (n_server_matched >= serv_tags_len
         || n_client_matched >= s_tags_nel) {
    /* no more possible matches on server or client */
    break;
      }

      /* search server for match on client */
      void *ptr = bsearch(s_tags + j*len,serv_tags,serv_tags_len,
              sizeof(*serv_tags),compare_first_int);

      /* Check for match, set src/dest proc id */
      if (ptr != NULL) {
	s_tags[j*len + 2] = i + nproc_macro;
	n_server_matched++;
	n_client_matched++;
      }
    }
    free(serv_tags);

    /* break loop if all matchess with client are already found */
    if (n_client_matched >= s_tags_nel) break;
  }

  /* sort s_tags by {idx} */
  qsort(s_tags,s_tags_nel,len*sizeof(*s_tags),compare_second_int);

  /* assign proc on send and recv */
  for(size_t i = 0; i < s_tags_nel; i++){
    const int tag = s_tags[i*len];
    const int idx = s_tags[i*len + 1];
    const int proc = s_tags[i*len + 2];

    /* error checking. Compiled out with NDEBUG */
    assert(send->tags[idx] == tag);
    assert(recv->tags[idx] == tag);


    /* set procs for send and recv */
  if (micro_model == 1) {//will eventually be replaced by an array
    client->send->set_proc_at_idx(proc,idx);
    client->recv->set_proc_at_idx(proc,idx);
  } else { //micro model == ROM
    client->send_ROM->set_proc_at_idx(proc,idx);
    client->recv_ROM->set_proc_at_idx(proc,idx);
  }


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
					   const MultiscaleComm *mscom)
{
  int err = 0;
  const size_t rank = mscom->rank_macro;
  const size_t n_server = client->n_server;
  const size_t n_server_ROM = client->n_server_ROM;
  /* exit early if not responsible for signaling servers */
  if(rank >= n_server) return err;

  int n_comm = 0;
  int n_comm_ROM = 0;
  int nproc_macro = 0;
  int nproc_inter = 0;
  int nproc_inter_ROM = 0;
  client->net->comm_size(mscom->macro, &nproc_macro);
  client->net->comm_size(mscom->mm_inter, &nproc_inter);
  client->net->comm_size(mscom->mm_inter_ROM, &nproc_inter_ROM);

  assert(nproc_macro > 0);
  assert(nproc_inter > 0);
  assert(nproc_inter_ROM > 0); //technically doesnt have to be > 0

  /* determine what ranks (in mm_inter) this domain is responsible
     for broadcasting info to. */
  if((unsigned)nproc_macro >= n_server) n_comm = 1;
  else {
    n_comm = n_server/nproc_macro;
    unsigned rem = n_server % nproc_macro;
    if(rem > 0 && rem > rank) n_comm++;

    n_comm_ROM = n_server_ROM/nproc_macro;
    rem = n_server_ROM % nproc_macro;
    if(rem > 0 && rem > rank) n_comm_ROM++;
  }

  /* allocate/populate client->bcast */
  client->bcast.n_comm = n_comm;
  client->bcast.ranks = PGFEM_malloc<int>(n_comm);
  client->net->allocRequestArray(n_comm, &client->bcast.req);
  {
    int * __restrict r = client->bcast.ranks; /* __restrict alias */
    r[0] = rank + nproc_macro;
    for(int i=1; i<n_comm; i++){
      r[i] = r[i-1] + nproc_macro;
    }
  }

  /*allocate second micro bcast */
  client->bcast_ROM.n_comm = n_comm_ROM;
  client->bcast_ROM.ranks = PGFEM_malloc<int>(n_comm_ROM);
  client->net->allocRequestArray(n_comm_ROM, &client->bcast_ROM.req);
  {
    int * __restrict r = client->bcast_ROM.ranks; /* __restrict alias */
    r[0] = rank + nproc_macro;
    for(int i=1; i<n_comm; i++){
      r[i] = r[i-1] + nproc_macro;
    }
  }



  /* check max rank for validity */
  assert(client->bcast.ranks[n_comm - 1] < nproc_inter);

  /* check max rank for validity */
  assert(client->bcast_ROM.ranks[n_comm_ROM - 1] < nproc_inter_ROM);//may not be necessary if no 2nd micro

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
			       const MultiscaleComm *mscom,int micro_model)
{
  int err = 0;
   
   auto n_server = client->n_server;
   auto bcast = client->bcast;
   auto mm_inter = mscom->mm_inter;

  if (micro_model == 2) {//will eventually be replaced by an array
    n_server = client->n_server_ROM;
    bcast = client->bcast_ROM;
    mm_inter = mscom->mm_inter_ROM;
  }

  /* get aliases */
  const size_t n_send = bcast.n_comm;
  const int *ranks = bcast.ranks;
  Request *req = bcast.req;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(client->net);
  
  /* should put request et al. in local buffer for overlay of
     comp/comm but will just waitall for now. */
  int nproc_macro = 0;
  const int rank = mscom->rank_macro;
  net->comm_size(mscom->macro, &nproc_macro);

  client->bcast.active = 1;
  client->bcast_ROM.active = 1;
  /* hold addresses of rebal buffers to send */
  void **buffs = PGFEM_calloc(void*, n_send);

  for(size_t i = 0; i < n_send; i++){
    size_t idx = rank + i*nproc_macro;
    assert(idx < n_server);
    size_t len = pgf_FE2_server_rebalance_n_bytes(rb_list[idx]);
    buffs[i] = pgf_FE2_server_rebalance_buff(rb_list[idx]);
    net->isend(buffs[i],len,NET_DT_CHAR,ranks[i],
	       FE2_MICRO_SERVER_REBALANCE,
	       mm_inter,req+i);
  }

  net->waitall(n_send,req,NET_STATUS_IGNORE);
  free(buffs);
  client->bcast.active = 0;
  client->bcast_ROM.active = 0;

  return err;
}

void pgf_FE2_macro_client_init(pgf_FE2_macro_client **client,
			       Network *n)
{
  *client = PGFEM_malloc<pgf_FE2_macro_client>();
  (*client)->n_jobs_loc = 0;
  (*client)->n_jobs_glob = 0;
  (*client)->n_jobs_max = 0;
  (*client)->n_server = 0;
  (*client)->n_server_ROM = 0;
  (*client)->jobs = NULL;
  (*client)->jobs_ROM = NULL;

  /* should be done by PGFEM_server_ctx init function... */
  (*client)->send = new MultiscaleServerContext(n);
  (*client)->recv = new MultiscaleServerContext(n);

  /* bcast */
  (*client)->bcast.active = 0;
  (*client)->bcast.n_comm = 0;
  (*client)->bcast.ranks = NULL;
  (*client)->bcast.req = NULL;

  /* bcast for micro 2 */
  (*client)->bcast_ROM.active = 0;
  (*client)->bcast_ROM.n_comm = 0;
  (*client)->bcast_ROM.ranks = NULL;
  (*client)->bcast_ROM.req = NULL;


  /* save the network handle */
  (*client)->net = n;
  
  /* other initialization stuff */
}

void pgf_FE2_macro_client_destroy(pgf_FE2_macro_client *client)
{
  delete client->send;
  delete client->recv;

  ISIRNetwork *net = static_cast<ISIRNetwork*>(client->net);
  
  for(int i=0,e=client->n_jobs_loc; i<e; i++){
    destroy_MS_COHE_JOB_INFO((client->jobs) + i);
  }
  free(client->jobs);

  /* destroy internal objects */
  if(client->bcast.active){
    /* complete communication before destruction */
    net->waitall(client->bcast.n_comm,client->bcast.req,
		 NET_STATUS_IGNORE);
  }
  free(client->bcast.ranks);
  delete [] client->bcast.req;

  /* destroy the handle */
  free(client);
  client = NULL;
}

void find_job_sizes_from_map(const Macroscale *macro,int *jobs_ROM, int *pde_jobs) {
  int i;
  for (i = 0; i < macro->nce ; i++ ) { 
    if (macro->opts->methods[i] == 1) { pde_jobs++;}
    else {jobs_ROM++;}
  }

}

void pgf_FE2_macro_client_create_job_list(pgf_FE2_macro_client *client,
                      const int n_jobs_max,
                      const Macroscale *macro,
                      const MultiscaleComm *mscom,
                      const int mp_id)
{
  /* compute and store number of servers */
  {
    int n_proc_macro = 0;
    int n_proc_inter = 0;
    int n_proc_inter_ROM = 0;
    client->net->comm_size(mscom->macro, &n_proc_macro);
    client->net->comm_size(mscom->mm_inter, &n_proc_inter);
    client->net->comm_size(mscom->mm_inter_ROM, &n_proc_inter_ROM);
    client->n_server = n_proc_inter - n_proc_macro;
    client->n_server_ROM= n_proc_inter_ROM - n_proc_macro;
  }

  const Macroscale *c = macro;

  /* get number of jobs and buffer sizes. */
  int n_jobs = 0;
  int jobs_ROM = 0;
  int pde_jobs = 0;
  int *job_buf_sizes = NULL;
  compute_n_job_and_job_sizes(c,&n_jobs,&job_buf_sizes);
  find_job_sizes_from_map(macro,&jobs_ROM, &pde_jobs); 

  client->n_jobs_loc = pde_jobs;
  client->n_jobs_loc_ROM = jobs_ROM;

  long Gn_jobs = 0;
  long *n_job_dom = NULL;//why is this created outside and not inside, if its destroyed right after?
  create_group_ms_cohe_job_list(pde_jobs,jobs_ROM,c->coel,c->node,
				mscom->macro,NET_COMM_SELF,
<<<<<<< HEAD
				0,&n_job_dom,
				&(client->jobs),&(client->jobs_ROM),c->net,mp_id);
=======
				0,&n_job_dom, &n_job_dom_ROM,
				&(client->jobs),&(client->jobs_ROM),c->net,mp_id,macro->opts);
>>>>>>> more 2-comm changes, getting ready for rebase
  /* free unused memory */
  free(n_job_dom);

  /* compute total number of jobs and set maximum number of jobs per
     server */
//  client->net->allreduce(NET_IN_PLACE,&Gn_jobs,1,NET_DT_LONG,
//			 NET_OP_SUM,mscom->macro);
  client->n_jobs_glob = pde_jobs;
  client->n_jobs_max = n_jobs_max;
  client->n_jobs_glob_ROM = jobs_ROM;//just one macro server for now

  /* make sure sufficient resources to compute solution */
  if(client->n_server * client->n_jobs_max < (size_t)Gn_jobs){
    PGFEM_printerr("ERROR: insufficent microscale resources!\n"
           "Increase n_jobs_max or number of servers and try again!\n");
    PGFEM_Abort();
  }

  /* create server contexts for send/recv of job information. */
  client->send->initialize(client->n_jobs_loc,job_buf_sizes);
  client->recv->initialize(client->n_jobs_loc,job_buf_sizes);
  client->send_ROM->initialize(client->n_jobs_loc_ROM,job_buf_sizes);
  client->recv_ROM->initialize(client->n_jobs_loc_ROM,job_buf_sizes);
 
  /* setup the bcast lists */
  pgf_FE2_macro_client_bcast_list(client,mscom);

  /* set tags for jobs */
  {
    const MS_COHE_JOB_INFO *j = client->jobs;
    for(int i=0; i<n_jobs; i++){
      const int tag = pgf_FE2_job_compute_encoded_id(j[i].proc_id,
						     j[i].elem_id,
						     j[i].int_pt);

      client->send->set_tag_at_idx(tag, i);
      client->recv->set_tag_at_idx(tag, i);
    }
    j = client->jobs_ROM;
    for(int i=0; i<((int)client->n_jobs_loc_ROM); i++){
      const int tag = pgf_FE2_job_compute_encoded_id(j[i].proc_id,
                 j[i].elem_id,
                 j[i].int_pt);

      client->send_ROM->set_tag_at_idx(tag, i);
      client->recv_ROM->set_tag_at_idx(tag, i);
    }
  }
  /* final cleanup */
  free(job_buf_sizes);
}

void pgf_FE2_macro_client_assign_initial_servers(pgf_FE2_macro_client *client,
						 const MultiscaleComm *mscom)
{
  /* create initial partition and send to microscale */
  int nproc_macro = 0;
  const int rank = mscom->rank_macro;
  client->net->comm_size(mscom->macro, &nproc_macro);
  int *__restrict n_jobs = PGFEM_malloc<int>(nproc_macro);//maybe nproc_macro is the maximum number of jobs there can be?
  int *__restrict n_jobs_ROM = PGFEM_malloc<int>(nproc_macro);
  int *__restrict displ = PGFEM_malloc<int>(nproc_macro);
  int *__restrict id = PGFEM_malloc<int>(client->n_jobs_glob);
  int *__restrict proc = PGFEM_malloc<int>(client->n_jobs_glob);
  int *__restrict time = PGFEM_malloc<int>(client->n_jobs_glob);

  int *__restrict id_ROM = PGFEM_malloc<int>(client->n_jobs_glob_ROM);
  int *__restrict proc_ROM = PGFEM_malloc<int>(client->n_jobs_glob_ROM);
  int *__restrict time_ROM = PGFEM_malloc<int>(client->n_jobs_glob_ROM);


  /*gather number of jobs given out by each macro*/
  n_jobs[rank] = client->n_jobs_loc;
  n_jobs_ROM[rank] = client->n_jobs_loc_ROM;
  client->net->allgather(NET_IN_PLACE,1,NET_DT_INT,n_jobs,1,NET_DT_INT,
			 mscom->macro);
  client->net->allgather(NET_IN_PLACE,1,NET_DT_INT,n_jobs_ROM,1,NET_DT_INT,
       mscom->macro);


  /* compute partial sum for displ and set proc_id and time*/
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

  /* compute partial sum for displ and set proc_id and time for rom*/
  displ[0] = 0;
  const size_t max_proc_ROM = client->n_server_ROM;
  proc_id = 0;
  for(int i=0; i<n_jobs_ROM[0]; i++){
    proc_ROM[i] = proc_id++;
    if(proc_id >= max_proc_ROM) proc_id = 0;
    time_ROM[i] = 1;
  }
  for(int i=1; i<nproc_macro; i++){
    displ[i] = displ[i-1] + n_jobs[i-1];
    const int start = displ[i];
    for(int j=0; j<n_jobs_ROM[i]; j++){
      proc_ROM[start + j] = proc_id++;
      if(proc_id >= max_proc_ROM) proc_id = 0;
      time_ROM[start + j] = 1;
    }
  }



  /* set job ids */
  const int start = displ[rank];
  memcpy(id+start,client->send->tags,client->send->n_comms*sizeof(*id));
  memcpy(id_ROM+start,client->send_ROM->tags,client->send_ROM->n_comms*sizeof(*id_ROM));

  /* get job ids from other procs */
  client->net->allgatherv(NET_IN_PLACE,n_jobs[rank],NET_DT_INT,
			  id,n_jobs,displ,NET_DT_INT,mscom->macro);
  client->net->allgatherv(NET_IN_PLACE,n_jobs_ROM[rank],NET_DT_INT,
        id_ROM,n_jobs_ROM,displ,NET_DT_INT,mscom->macro);

  void *all_part = NULL;
  void *parts = NULL;
  void *all_part_ROM = NULL;
  void *parts_ROM = NULL;
  /*fill job list (all_parts) with job info*/
  new_partition_build_set_keep(&all_part,client->n_jobs_glob,
                   client->n_jobs_glob,id,time,proc);
  new_partitions_void(&parts,client->n_server,client->n_jobs_max);

  /*fill job info for ROM */
  new_partition_build_set_keep(&all_part_ROM,client->n_jobs_glob_ROM,
                   client->n_jobs_glob_ROM,id_ROM,time_ROM,proc_ROM);
  new_partitions_void(&parts_ROM,client->n_server_ROM,client->n_jobs_max);


  /* partition the jobs using the greedy algorithm. Since everything
     has the same weight, amounts to round robin assignment. */
  rebalance_partitions_greedy(client->n_server,all_part,parts);

  rebalance_partitions_greedy(client->n_server_ROM,all_part_ROM,parts_ROM);


  /* push partitions into structure I will send */
  auto rb = PGFEM_calloc(pgf_FE2_server_rebalance *, client->n_server);
  auto rb_ROM = PGFEM_calloc(pgf_FE2_server_rebalance *, client->n_server_ROM);

  new_partitions_void_to_pgf_FE2_server_rebalance(client->n_server,parts,rb);
  new_partitions_void_to_pgf_FE2_server_rebalance(client->n_server_ROM,parts_ROM,rb_ROM);

  /* cleanup */
  new_partition_destroy_void(all_part);
  new_partitions_destroy_void(parts,client->n_server);

  new_partition_destroy_void(all_part_ROM);
  new_partitions_destroy_void(parts_ROM,client->n_server_ROM);

  /* Update server context (who I am sending to)*/
  pgf_FE2_macro_client_update_send_recv(client,rb,nproc_macro,1);//full pde micros
  pgf_FE2_macro_client_update_send_recv(client,rb_ROM,nproc_macro,2);//taylor model micros


  /* Broadcast the rebalancing information to the servers */
  pgf_FE2_macro_client_bcast_rebal_to_servers(client,rb,mscom,1);//full pde micros
  pgf_FE2_macro_client_bcast_rebal_to_servers(client,rb_ROM,mscom,2);//taylor model micros

  /* final cleanup */
  for(size_t i=0, e=client->n_server; i<e; i++){
    pgf_FE2_server_rebalance_destroy(rb[i]);
  }
  free(rb);
}

void pgf_FE2_macro_client_rebalance_servers(pgf_FE2_macro_client *client,
					    const MultiscaleComm *mscom,
					    const int heuristic)
{
  int nproc_macro = 0;
  client->net->comm_size(mscom->macro, &nproc_macro);

  /* receive message and rebalance according to heuristic */
  pgf_FE2_server_rebalance **rb_list = pgf_FE2_rebalancer(client->net,
							  mscom,
							  client->n_jobs_glob,
							  client->n_jobs_max,
							  heuristic,1);
  pgf_FE2_server_rebalance **rb_list_ROM = pgf_FE2_rebalancer(client->net,
                mscom,
                client->n_jobs_glob_ROM,
                client->n_jobs_max,
                heuristic,2);


  /* update the server context (send/recv) */
  pgf_FE2_macro_client_update_send_recv(client,rb_list,nproc_macro,1);
  pgf_FE2_macro_client_update_send_recv(client,rb_list_ROM,nproc_macro,2);

  /* Broadcast rebalancing info to servers */
  pgf_FE2_macro_client_bcast_rebal_to_servers(client,rb_list,mscom,1);
  pgf_FE2_macro_client_bcast_rebal_to_servers(client,rb_list_ROM,mscom,2);

  /* cleanup */
  for(size_t i=0, e=client->n_server; i<e; i++){
    pgf_FE2_server_rebalance_destroy(rb_list[i]);
    pgf_FE2_server_rebalance_destroy(rb_list_ROM[i]);
  }
  free(rb_list);
  free(rb_list_ROM);
}

void pgf_FE2_macro_client_send_jobs(pgf_FE2_macro_client *client,
                    const MultiscaleComm *mscom,
                    const Macroscale *macro,
                    const int job_type,
                    const int micro_model)
{
  int err = 0;
  /* see start_macroscale_compute_jobs */
    MultiscaleServerContext *recv;
    MultiscaleServerContext *send;
    MS_COHE_JOB_INFO *job_list;
    PGFem3D_Comm comm;;

  if (micro_model == 2) {//will eventually be replaced by an array

    recv = client->recv_ROM;
    send = client->send_ROM;
    job_list = client->jobs_ROM;
    comm = mscom->mm_inter_ROM;
  } else {
    recv = client->recv;
    send = client->send;
    job_list = client->jobs;
    comm = mscom->mm_inter;
  }


  ISIRNetwork *net = static_cast<ISIRNetwork*>(client->net);
  /* Wait to complete any pending communcication. The error flag is
     incremented as this should only occur by a logical/programming
     error. Note that buffers may be overwritten. If the communication
     is completed but the flags have not been reset to 0, the
     Wait* commands should return immediately */
  if(send->in_process || recv->in_process){
    PGFEM_printerr("WARNING: communication in progress!(%s:%s:%d)\n",
           __func__,__FILE__,__LINE__);
    err++;
    net->waitall(recv->n_comms,recv->req,recv->stat);
    net->waitall(send->n_comms,send->req,send->stat);
    recv->in_process = 0;
    send->in_process = 0;
  }


  /* post receives (from the running server) */
  for(int i=0; i<recv->n_comms; i++){
    net->irecv(recv->buffer[i],recv->sizes[i],NET_DT_CHAR,
	       recv->procs[i],recv->tags[i],comm,
	       &(recv->req[i]));
  }

  for(int i=0; i<send->n_comms; i++){
    /* update the job information according to job_type. Note that we
       send the *increment* of the solution (f = d_r + rr) since the
       cohesive elements are implemented assuming updated formulation
       of the displacement field (u_n computed from coordinates) */
    if (micro_model == 2) {
      err += macroscale_update_job_info(macro,job_type,macro->sol->f,job_list+i);
    } else {
      err += macroscale_update_job_info(macro,job_type,macro->sol->f,job_list+i);
    }
    /* pack the job info into the buffer to send */
    err += pack_MS_COHE_JOB_INFO(job_list + i,send->sizes[i],
                 send->buffer[i]);

    /* post the send (to the running server) */
    net->isend(send->buffer[i],send->sizes[i],NET_DT_CHAR,
	       send->procs[i],send->tags[i],comm,
	       &(send->req[i]));
  }

  /* set in_process flags to true (1) */
  recv->in_process = 1;
  send->in_process = 1;
}

void pgf_FE2_macro_client_recv_jobs(pgf_FE2_macro_client *client,
				    Macroscale *macro,
				    int *max_micro_sub_step,
            int micro_model)
{
  /* Get aliases from client object etc. */
  MultiscaleServerContext *recv;
  MultiscaleServerContext *send;
  MS_COHE_JOB_INFO *job_list;
  if (micro_model == 2) {//will eventually be replaced by an array
    recv = client->recv_ROM;
    send = client->send_ROM;
    job_list = client->jobs_ROM;
  } else {
    recv = client->recv;
    send = client->send;
    job_list = client->jobs;
  }

  /* reset the max number of steps */
  *max_micro_sub_step = 0;

  /* see finish_macroscale_compute_jobs */
  Macroscale *c = macro;
  MULTISCALE_SOLUTION *s = macro->sol;
  int rank_macro = 0;
  int nproc_macro = 0;

  ISIRNetwork *net = static_cast<ISIRNetwork*>(client->net);
  
  /* exit early if !*->in_process */
  if(!recv->in_process && !send->in_process) return;

  net->comm_size(c->comm, &rank_macro);
  net->comm_size(c->comm, &nproc_macro);

  /* if expecting to receive buffers */
  if(recv->in_process){
    /* set up the stiffness matrix communication */
    double **Lk = NULL;
    double **receive = NULL;
    c->spc->post_stiffmat(&Lk, &receive);

    /* assemble jobs as they are received */
    for(int i=0; i<recv->n_comms; i++){
      int idx = 0;
      net->waitany(recv->n_comms,recv->req,&idx,recv->stat);
      MS_COHE_JOB_INFO *job = job_list + idx;
      unpack_MS_COHE_JOB_INFO(job,recv->sizes[idx],
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
	      PLoc_Sparse(Lk,job->K_00_contrib,NULL,NULL,NULL,job->g_dof_ids,
		    job->ndofe,NULL,c->GDof,rank_macro,nproc_macro,
		    c->spc,0,c->SOLVER,c->opts->analysis_type);
	      break;
      case JOB_UPDATE:
	/* update cohesive elements */

        if (micro_model == 2) {//ROM
  	      macroscale_update_coel(job,macro);
        } else {  //pde
          macroscale_update_coel(job,macro);
        }
    	break;
      default: /* do nothing */ break;
      }
    }

    /* send/finalize communication of the stiffness matrix */
    c->spc->send_stiffmat();
//    c->spc->print_stiffmat();
    /* get maximum number of steps from all macro processors */
    net->allreduce(NET_IN_PLACE,max_micro_sub_step,
		   1,NET_DT_INT,NET_OP_MAX,c->comm);

    c->spc->assemble_nonlocal_stiffmat(c->SOLVER);
    c->spc->finalize_stiffmat();

    /* set in_process to false (0) */
    recv->in_process = 0;
  }/* end if(recv->in_process) */

  /* wait for any pending send communication to finish */
  if(send->in_process){
    net->waitall(send->n_comms,send->req,send->stat);
    send->in_process = 0;
  }

}

void pgf_FE2_macro_client_send_exit(pgf_FE2_macro_client *client,
				    const MultiscaleComm *mscom)
{
  const size_t n_send = client->bcast.n_comm;
  const int *ranks = client->bcast.ranks;
  Request *req = client->bcast.req;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(client->net);
  
  client->bcast.active = 1;
  void *empty = NULL;
  for(size_t i=0; i<n_send; i++){
    net->isend(empty,0,NET_DT_CHAR,ranks[i],FE2_MICRO_SERVER_EXIT,
	       mscom->mm_inter,req + i);
  }
  
  /* wait for exit code to be received */
  net->waitall(n_send,req,NET_STATUS_IGNORE);
  client->bcast.active = 0;

  const size_t n_send_ROM = client->bcast_ROM.n_comm;
  const int *ranks_ROM = client->bcast_ROM.ranks;
  Request *req_ROM = client->bcast_ROM.req;

  client->bcast_ROM.active = 1;
  void *empty_ROM = NULL;
  for(size_t i=0; i<n_send_ROM; i++){
    net->isend(empty_ROM,0,NET_DT_CHAR,ranks_ROM[i],FE2_MICRO_SERVER_EXIT,
         mscom->mm_inter_ROM,req_ROM + i);
  }

  /* wait for exit code to be received */
  net->waitall(n_send_ROM,req,NET_STATUS_IGNORE);
  client->bcast_ROM.active = 0;

}

