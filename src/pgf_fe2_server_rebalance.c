/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_server_rebalance.h"

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

size_t pgf_FE2_server_rebalance_n_bytes(const pgf_FE2_server_rebalance *t)
{
  return get_size((*t)[REBAL_N_KEEP],
		  (*t)[REBAL_N_SEND],
		  (*t)[REBAL_N_RECV]);
}

void pgf_FE2_server_rebalance_build(pgf_FE2_server_rebalance *t,
				    const size_t n_keep,
				    const size_t n_send,
				    const size_t n_recv)
{
  /* currently, pgf_FE2_server_rebalance is simply a typedef for a
     pointer to an integer, but we will have a special layout of the
     data described in the code below. */

  /* allocate buffer */
  *t = calloc(get_size(n_keep,n_send,n_recv),sizeof(**t));

  /* encode separate buffer sizes */
  (*t)[REBAL_N_KEEP] = n_keep;
  (*t)[REBAL_N_SEND] = n_send;
  (*t)[REBAL_N_RECV] = n_recv;

  /* encode offsets */
  (*t)[REBAL_KEEP_OFF] = REBAL_N_META; /* keep_offset */
  (*t)[REBAL_SEND_OFF] = (*t)[3] + n_keep; /* send_offset */
  (*t)[REBAL_RECV_OFF] = (*t)[4] + 2*n_send; /* recv_offset */
}


void pgf_FE2_server_rebalance_build_from_buffer(pgf_FE2_server_rebalance *t,
						 void *buffer)
{
  /* avoid memory leaks */
  pgf_FE2_server_rebalance_destroy(t);
  *t = buffer;

  /* invalidate buffer */
  buffer = NULL;
}

void pgf_FE2_server_rebalance_destroy(pgf_FE2_server_rebalance *t)
{
  free(*t);
}

int pgf_FE2_server_rebalance_n_keep(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t)[REBAL_N_KEEP];
  else return -1;
}

int* pgf_FE2_server_rebalance_keep_buf(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return *t + (*t)[REBAL_KEEP_OFF];
  else return NULL;
}

int pgf_FE2_server_rebalance_n_send(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t)[REBAL_N_SEND];
  else return -1;
}

int* pgf_FE2_server_rebalance_send_buf(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_SEND_OFF];
  else return NULL;
}

int* pgf_FE2_server_rebalance_send_dest(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_SEND_OFF] + (*t)[REBAL_N_SEND];
  else return NULL;
}

int pgf_FE2_server_rebalance_n_recv(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t)[REBAL_N_RECV];
  else return -1;
}

int* pgf_FE2_server_rebalance_recv_buf(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_RECV_OFF];
  else return NULL;
}

int* pgf_FE2_server_rebalance_recv_src(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_RECV_OFF] + (*t)[REBAL_N_RECV];
  else return NULL;
}


typedef struct{
  size_t count;
  char **buffers;
  MPI_Request *reqs;
} pgf_FE2_server_rebalance_workspace;

pgf_FE2_server_rebalance_workspace *send_wkspc = NULL, *recv_wkspc = NULL;

static void allocate_workspace(pgf_FE2_server_rebalance_workspace **wkspc,
			       const size_t n_comm)
{
  (*wkspc)->count = n_comm;
  *wkspc = malloc(sizeof(**wkspc));
  (*wkspc)->reqs = malloc(n_comm*sizeof(*((*wkspc)->reqs)));
  (*wkspc)->buffers = malloc(n_comm*sizeof(*((*wkspc)->buffers)));
}

static void destroy_workspace(pgf_FE2_server_rebalance_workspace *wkspc)
{
  if(wkspc == NULL) return;
  free(wkspc->buffers);
  free(wkspc->reqs);
}

static void allocate_recv_workspace(const size_t n_recv)
{
  allocate_workspace(&recv_wkspc,n_recv);
}

static void destroy_recv_workspace()
{
  destroy_workspace(recv_wkspc);
  free(recv_wkspc);
  recv_wkspc = NULL;
}

static void allocate_send_workspace(const size_t n_send,
				    const size_t buff_size)
{
  allocate_workspace(&send_wkspc,n_send);
  for(size_t i=0; i<n_send; i++){
    send_wkspc->buffers[i] = malloc(n_send*sizeof(*(send_wkspc->buffers[i])));
  }
}

static void destroy_send_workspace()
{
  if(send_wkspc == NULL) return;
  for(size_t i=0,e=send_wkspc->count; i<e; i++){
    free(send_wkspc->buffers[i]);
  }
  destroy_workspace(send_wkspc);
  free(send_wkspc);
  send_wkspc = NULL;
}

int pgf_FE2_server_rebalance_post_exchange(const pgf_FE2_server_rebalance *t,
					   const PGFEM_mpi_comm *mpi_comm,
					   MICROSCALE *micro)
{
  int err = 0;

  /* alias to solution index map */
  sol_idx_map *pmap = &(micro->idx_map);

  const size_t buff_size = micro->sol[0].packed_state_var_len;
  const int n_send = pgf_FE2_server_rebalance_n_send(t);
  const int *send_id = pgf_FE2_server_rebalance_send_buf(t);
  const int *send_to = pgf_FE2_server_rebalance_send_dest(t);

  /* allocate file local buffer for send */
  allocate_send_workspace(n_send,buff_size);

  /* copy and send information */
  for(int i=0; i<n_send; i++){
    /* search micro solutions list for matching id and mark as available. */
    int sol_idx = sol_idx_map_get_idx_reset_id(pmap,send_id[i],-1);

    /* copy state at n to buffer */
    memcpy(send_wkspc->buffers[i],micro->sol[sol_idx].packed_state_var_n,buff_size);

    /* post non-blocking send of the data */
    err += MPI_Isend(send_wkspc->buffers[i],
		     buff_size, MPI_CHAR,
		     send_to[i], send_id[i],
		     mpi_comm->worker_inter,
		     &(send_wkspc->reqs[i]));
  }

  const int n_recv = pgf_FE2_server_rebalance_n_recv(t);
  const int *recv_id = pgf_FE2_server_rebalance_recv_buf(t);
  const int *recv_from = pgf_FE2_server_rebalance_recv_src(t);

  /* allocate file local buffer for recv */
  allocate_recv_workspace(n_recv);

  /* set workspace buffers to point to solution space and post receives */
  for(int i=0; i<n_recv; i++){
    /* search micro solutions list for available solution space and assign new id */
    int sol_idx = sol_idx_map_get_idx_reset_id(pmap,-1,recv_id[i]);

    /* set workspace buffer to _point_ to the state at n buffer */
    recv_wkspc->buffers[i] = micro->sol[sol_idx].packed_state_var_n;

    /* post non-blocking receive of the data */
    err += MPI_Irecv(recv_wkspc->buffers[i],
		     buff_size,MPI_CHAR,
		     recv_from[i], recv_id[i],
		     mpi_comm->worker_inter,
		     &(recv_wkspc->reqs[i]));
  }

  return err;
}

int pgf_FE2_server_rebalance_finalize_exchange(const PGFEM_mpi_comm *mpi_comm)
{
  int err = 0;

  /* complete sends. Use waitsome to aid in splitting up later */
  {
    int outcount = 0;
    int *idx = calloc(send_wkspc->count,sizeof(*idx));
    while(outcount != MPI_UNDEFINED){
      err += MPI_Waitsome(send_wkspc->count,send_wkspc->reqs,
			  &outcount,idx,MPI_STATUS_IGNORE);
    }
    free(idx);
  }
  destroy_send_workspace();

  /* complete receives. Use Waitsome to aid in splitting later. */
  {
    int outcount = 0;
    int *idx = calloc(recv_wkspc->count,sizeof(*idx));
    while(outcount != MPI_UNDEFINED){
      err += MPI_Waitsome(recv_wkspc->count,recv_wkspc->reqs,
			  &outcount,idx,MPI_STATUS_IGNORE);
    }
    free(idx);
  }
  destroy_recv_workspace();

  /* barrier to ensure that transfer comm is complete on all procs in
     the server. */
  MPI_Barrier(mpi_comm->micro);
  return err;
}
