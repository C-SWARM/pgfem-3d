/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "allocation.h"
#include "pgf_fe2_server_rebalance.h"
#include "pgf_fe2_restart.h"
#include <assert.h>

using namespace pgfem3d;
using namespace pgfem3d::net;

/**
 * Internal data structure
 */
namespace {
struct Workspace {
  size_t count = 0;
  char **buffers = nullptr;
  Request *reqs = nullptr;
};
}

/**
 * main encapsulated data structure
 */
struct pgf_FE2_server_rebalance
{
  int *t;
  Workspace *send_wkspc;
  Workspace *recv_wkspc;
};

static void allocate_workspace(Workspace **wkspc, const size_t n_comm,
			       const Microscale *micro)
{
  assert(wkspc != NULL && "wkspc can't be null");
  *wkspc = static_cast<Workspace*>(malloc(sizeof(Workspace)));
  (*wkspc)->count = n_comm;
  (*wkspc)->buffers = static_cast<char**>(malloc(n_comm*sizeof(*((*wkspc)->buffers))));
  micro->net->allocRequestArray(n_comm, &(*wkspc)->reqs);
}

static void destroy_workspace(Workspace *wkspc)
{
  if(wkspc == NULL) return;
  free(wkspc->buffers);
  delete [] wkspc->reqs;
}

static void allocate_recv_workspace(pgf_FE2_server_rebalance *t,
				    const Microscale *micro,
				    const size_t n_recv)
{
  allocate_workspace(&(t->recv_wkspc),n_recv,micro);
}

static void destroy_recv_workspace(pgf_FE2_server_rebalance *t)
{
  destroy_workspace(t->recv_wkspc);
  free(t->recv_wkspc);
  t->recv_wkspc = NULL;
}

static void allocate_send_workspace(pgf_FE2_server_rebalance *t,
				    const Microscale *micro,
				    const size_t n_send,
				    const size_t buff_size)
{
  allocate_workspace(&(t->send_wkspc),n_send,micro);
  for(size_t i=0; i<n_send; i++){
    t->send_wkspc->buffers[i] = static_cast<char*>(malloc(buff_size));
  }
}

static void destroy_send_workspace(pgf_FE2_server_rebalance *t)
{
  if(t->send_wkspc == NULL) return;
  for(size_t i=0,e=(t->send_wkspc)->count; i<e; i++){
    free((t->send_wkspc)->buffers[i]);
  }
  destroy_workspace(t->send_wkspc);
  free(t->send_wkspc);
  t->send_wkspc = NULL;
}

enum {REBAL_N_KEEP=0,
      REBAL_N_SEND,
      REBAL_N_RECV,
      REBAL_KEEP_OFF,
      REBAL_SEND_OFF,
      REBAL_RECV_OFF,
      REBAL_N_META};

static size_t get_size(const size_t n_keep,
               const size_t n_send,
               const size_t n_recv)
{
  return ((size_t) REBAL_N_META) + n_keep + 2*(n_send + n_recv);
}

size_t
pgf_FE2_server_rebalance_n_bytes(const pgf_FE2_server_rebalance *R)
{
  const int *t = R->t; /* alias */
  return get_size((t)[REBAL_N_KEEP],
          (t)[REBAL_N_SEND],
          (t)[REBAL_N_RECV])*sizeof(*t);
}

void*
pgf_FE2_server_rebalance_buff(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return t->t;
  else return NULL;
}

void
pgf_FE2_server_rebalance_build(pgf_FE2_server_rebalance **R,
                               const size_t n_keep,
                               const size_t n_send,
                               const size_t n_recv)
{
  /* currently, pgf_FE2_server_rebalance is simply a typedef for a
     pointer to an integer, but we will have a special layout of the
     data described in the code below. */

  /* allocate buffer */
  *R = static_cast<pgf_FE2_server_rebalance*>(malloc(sizeof(**R)));
  (*R)->t = PGFEM_calloc(int, get_size(n_keep,n_send,n_recv));
  int *t = (*R)->t; /* alias */

  /* encode separate buffer sizes */
  t[REBAL_N_KEEP] = n_keep;
  t[REBAL_N_SEND] = n_send;
  t[REBAL_N_RECV] = n_recv;

  /* encode offsets */
  t[REBAL_KEEP_OFF] = REBAL_N_META; /* keep_offset */
  t[REBAL_SEND_OFF] = t[3] + n_keep; /* send_offset */
  t[REBAL_RECV_OFF] = t[4] + 2*n_send; /* recv_offset */

  /* set workspaces to NULL */
  (*R)->send_wkspc = NULL;
  (*R)->recv_wkspc = NULL;
}


void
pgf_FE2_server_rebalance_build_from_buffer(pgf_FE2_server_rebalance **t, void **buffer)
{
  /* avoid memory leaks */
  pgf_FE2_server_rebalance_destroy(*t);
  *t = static_cast<pgf_FE2_server_rebalance*>(malloc(sizeof(pgf_FE2_server_rebalance)));
  (*t)->t = static_cast<int*>(*buffer);

  /* invalidate pointer to buffer */
  *buffer = NULL;

  /* set workspaces to NULL */
  (*t)->send_wkspc = NULL;
  (*t)->recv_wkspc = NULL;
}

void
pgf_FE2_server_rebalance_destroy(pgf_FE2_server_rebalance *t)
{
  if(t != NULL){
    free(t->t);
    if(t->send_wkspc != NULL) destroy_send_workspace(t);
    if(t->recv_wkspc != NULL) destroy_recv_workspace(t);
  }
  free(t);
}

int
pgf_FE2_server_rebalance_n_keep(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return t->t[REBAL_N_KEEP];
  else return -1;
}

int*
pgf_FE2_server_rebalance_keep_buf(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return (t->t + (t->t)[REBAL_KEEP_OFF]);
  else return NULL;
}

size_t
pgf_FE2_server_rebalance_n_send(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return (t->t)[REBAL_N_SEND];
  else return 0;
}

int*
pgf_FE2_server_rebalance_send_buf(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return (t->t + (t->t)[REBAL_SEND_OFF]);
  else return NULL;
}

int*
pgf_FE2_server_rebalance_send_dest(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return (t->t + (t->t)[REBAL_SEND_OFF] + (t->t)[REBAL_N_SEND]);
  else return NULL;
}

size_t
pgf_FE2_server_rebalance_n_recv(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) {
    return (t->t)[REBAL_N_RECV];
  }
  else {
    return 0;
  }
}

int*
pgf_FE2_server_rebalance_recv_buf(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return (t->t + (t->t)[REBAL_RECV_OFF]);
  else return NULL;
}

int*
pgf_FE2_server_rebalance_recv_src(const pgf_FE2_server_rebalance *t)
{
  if(t != NULL) return (t->t + (t->t)[REBAL_RECV_OFF] + (t->t)[REBAL_N_RECV]);
  else return NULL;
}


static int dbg_cmp_int(const void *a, const void *b)
{
  return *((const int*) a) - *((const int*) b);
}

/** Verify that the id's on the microscale servers are unique */
static int debug_keep_id(const int n_keep,
                         const int *keep_id,
			 Network *n,
                         const MultiscaleComm *mscom)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(n);
  int err = 0;
  int n_servers = -1;
  net->comm_size(mscom->worker_inter,&n_servers);
  
  /* get n_keep from each proc on worker_inter */
  int *counts = PGFEM_calloc(int, n_servers);
  int *displ = PGFEM_calloc(int, n_servers + 1);
  net->allgather(&n_keep,1,NET_DT_INT,counts,1,NET_DT_INT,
		 mscom->worker_inter);

  /* allocate buffer for all keep_id */
  int total_n_keep = counts[0];
  displ[0] = 0;
  for(int i = 1; i < n_servers; i++){
    displ[i] = displ[i-1] + counts[i-1];
    total_n_keep += counts[i];
  }
  int *all_keep_id = PGFEM_calloc(int, total_n_keep);

  /* Gather keep_id from all procs */
  net->allgatherv(keep_id,n_keep,NET_DT_INT,all_keep_id,counts,
		  displ,NET_DT_INT,mscom->worker_inter);
  
  /* sort array */
  qsort(all_keep_id,total_n_keep,sizeof(*all_keep_id),dbg_cmp_int);

  /* linear search for duplicates on master of server 0 */
  if(mscom->server_id == 0 && mscom->rank_micro == 0){
    for(int i = 0; i < n_keep - 1; i++){
      if(all_keep_id[i] == all_keep_id[i+1]){
	PGFEM_printerr("Found duplicate: %d\n",all_keep_id[i]);
	assert(0); /* Cause abort on duplicate. Requires investigation */
	err ++;
      }
    }
  }
  
  /* brodcast error value to slaves */
  net->bcast(&err,1,NET_DT_INT,0,mscom->worker_inter); /* between eq. procs */
  net->bcast(&err,1,NET_DT_INT,0,mscom->micro); /* master to slaves */

  free(counts);
  free(displ);
  free(all_keep_id);
  return err;
}

int
pgf_FE2_server_rebalance_post_exchange(pgf_FE2_server_rebalance *t,
                                       const MultiscaleComm *mscom,
                                       Microscale *micro)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(micro->net);
  static int sol_id_initialized = 0;
  int err = 0;

  /* alias to solution index map */
  sol_idx_map *pmap = &(micro->idx_map);

  /* initialize ids of solutions from keep buffer. This is done once at
     startup and assumes that we are not sending/receiving info from
     other servers. */
  if(!sol_id_initialized){
    int n_keep = pgf_FE2_server_rebalance_n_keep(t);
    const int *keep_id = pgf_FE2_server_rebalance_keep_buf(t);

#ifndef NDEBUG
    err += debug_keep_id(n_keep,keep_id,micro->net,mscom);
#endif

    /* set first n_keep ids to match tags */
    for(int i=0; i<n_keep; i++){
      sol_idx_map_idx_set_id(pmap,i,keep_id[i]);
    }

    /* ensure all others are -1 (available) */
    for(int i=n_keep, e=pmap->size; i<e; i++){
      sol_idx_map_idx_set_id(pmap,i,-1);
    }

    /* read restart information if needed */
    if(micro->opts->restart >= 0){
      const size_t step = micro->opts->restart;
      for(int i=0; i<n_keep; i++){
	err += pgf_FE2_restart_read_micro(micro,step,keep_id[i]);
      }
      /* turn off restart at the microscale */
      micro->opts->restart = -1;
    }

    /* increment flag so we do not do this initialization again */
    sol_id_initialized++;
  }

  const size_t buff_size = micro->sol[0].packed_state_var_len;
  const int n_send = pgf_FE2_server_rebalance_n_send(t);
  const int *send_id = pgf_FE2_server_rebalance_send_buf(t);
  const int *send_to = pgf_FE2_server_rebalance_send_dest(t);

  /* allocate file local buffer for send */
  allocate_send_workspace(t,micro,n_send,buff_size);

  /* copy and send information */
  for(int i=0; i<n_send; i++){
    /* search micro solutions list for matching id and mark as available. */
    int sol_idx = sol_idx_map_get_idx_reset_id(pmap,send_id[i],-1);

    /* copy state at n to buffer */
    void *src = (micro->sol[sol_idx]).packed_state_var_n;
    void *dest = t->send_wkspc->buffers[i];
    memcpy(dest,src,buff_size);

    /* post non-blocking send of the data */
    net->isend(t->send_wkspc->buffers[i],
	       buff_size, NET_DT_CHAR,
	       send_to[i], send_id[i],
	       mscom->worker_inter,
	       &(t->send_wkspc->reqs[i]));
  }
  
  const int n_recv = pgf_FE2_server_rebalance_n_recv(t);
  const int *recv_id = pgf_FE2_server_rebalance_recv_buf(t);
  const int *recv_from = pgf_FE2_server_rebalance_recv_src(t);

  /* allocate file local buffer for recv */
  allocate_recv_workspace(t,micro,n_recv);

  /* set workspace buffers to point to solution space and post receives */
  for(int i=0; i<n_recv; i++){
    /* search micro solutions list for available solution space and assign new id */
    int sol_idx = sol_idx_map_get_idx_reset_id(pmap,-1,recv_id[i]);

    /* set workspace buffer to _point_ to the state at n buffer */
    t->recv_wkspc->buffers[i] = micro->sol[sol_idx].packed_state_var_n;
    
    /* post non-blocking receive of the data */
    net->irecv(t->recv_wkspc->buffers[i],
	       buff_size,NET_DT_CHAR,
	       recv_from[i], recv_id[i],
	       mscom->worker_inter,
	       &(t->recv_wkspc->reqs[i]));
  }
  
  return err;
}

int
pgf_FE2_server_rebalance_finalize_exchange(pgf_FE2_server_rebalance *t,
                                           const MultiscaleComm *mscom,
					   const Microscale *micro)
{
  int err = 0;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(micro->net);
  
  /* complete sends. Use waitsome to aid in splitting up later */
  {
    int outcount = 0;
    int *idx = PGFEM_calloc(int, t->send_wkspc->count);
    while(outcount != NET_UNDEFINED){
      net->waitsome(t->send_wkspc->count,t->send_wkspc->reqs,
		    &outcount,idx,NET_STATUS_IGNORE);
    }
    free(idx);
  }
  destroy_send_workspace(t);

  /* complete receives. Use Waitsome to aid in splitting later. */
  {
    int outcount = 0;
    int *idx = PGFEM_calloc(int, t->recv_wkspc->count);
    while(outcount != NET_UNDEFINED){
      net->waitsome(t->recv_wkspc->count,t->recv_wkspc->reqs,
		    &outcount,idx,NET_STATUS_IGNORE);
    }
    free(idx);
  }

  destroy_recv_workspace(t);

  /* barrier to ensure that transfer comm is complete on all procs in
     the server. */
  net->barrier(mscom->micro);
  return err;
}
