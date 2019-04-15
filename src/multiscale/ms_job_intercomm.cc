/* HEADER */
#include "ms_job_intercomm.h"
#include <string.h>
#include "allocation.h"
#include "PGFEM_io.h"
#include "utils.h"

/*==== STATIC FUNCTION PROTOTYPES ====*/
static int create_ms_job_intercomm_micro(const int nproc_macro,
                     const int n_jobs,
                     const int *job_buff_sizes,
                     PGFEM_ms_job_intercomm *ic);

static int create_ms_job_intercomm_macro(const int n_jobs,
                     const int *job_buff_sizes,
                     PGFEM_ms_job_intercomm *ic,
                     MSNET_Comm macro_comm);

/*==== API FUNCTION DEFINITIONS ====*/
int create_PGFEM_ms_job_intercomm(const int nproc_macro, //deprecated
                  const MultiscaleCommunicator*mscom,
                  const int n_jobs,
                  const int *job_buff_sizes,
                  PGFEM_ms_job_intercomm **intercomm)
{
  int err = 0;
  *intercomm = NULL;
  if(!comm->valid_mm_inter){
    return err;
  }
  /* working within the intercommunicator only */
  else {
    /* allocate the communication object */
    *intercomm = PGFEM_calloc(PGFEM_ms_job_intercomm, 1);
    PGFEM_ms_job_intercomm *ic = *intercomm; /* alias */

    /* duplicate the MPI_Comm communicator */
    /* err += MPI_Comm_dup(comm->mm_inter,&(ic->comm)); */
    ic->comm = comm->mm_inter;


    if(comm->valid_macro){/*=== MACROSCALE ===*/
      /* only communicates with microscale */
      err += create_ms_job_intercomm_macro(n_jobs,job_buff_sizes,
                       ic,comm->macro);
    } else { /*=== MICROSCALE ===*/
      /* only communicates with macroscale */
      err += create_ms_job_intercomm_micro(nproc_macro,n_jobs,
                       job_buff_sizes,ic);
    }/*microscale */

  }/* end work on intercommunicator */
  return err;
}

int destroy_PGFEM_ms_job_intercomm(PGFEM_ms_job_intercomm *intercomm)
{
  int err = 0;
  if(intercomm != NULL){
    err += MPI_Comm_free(&(intercomm->comm));
    err += destroy_PGFEM_comm_info(intercomm->send_info);
    err += destroy_PGFEM_comm_info(intercomm->recv_info);
    free(intercomm->send_info);
    free(intercomm->recv_info);
  }
  return err;
}

/*==== STATIC FUNCTION DEFINITIONS ====*/
static int create_ms_job_intercomm_micro(const int nproc_macro,
                     const int n_jobs,
                     const int *job_buff_sizes,
                     PGFEM_ms_job_intercomm *ic)
{
  int err = 0;
  int nproc_inter = 0;
  int rank_inter = 0;
  err += MPI_Comm_rank(ic->comm,&rank_inter);
  err += MPI_Comm_size(ic->comm,&nproc_inter);

  /* should not have any jobs to send */
  if(n_jobs > 0){
    int rank = 0;
    PGFEM_Error_rank(&rank);
    PGFEM_printerr("[%d(w)]ERROR: have jobs on microscale!?!\n",
           rank);
    PGFEM_Abort();
  }

  /* initialize receive */
  int *n_job_recv = PGFEM_calloc(int, nproc_inter);

  /* post receives. NOTE: only commmunicate with macro processes which
     are listed first */
  MPI_Request *req_r = PGFEM_calloc(MPI_Request, nproc_macro);
  int idx = 0;
  for(int i=0; i<nproc_macro; i++){
    if(i == rank_inter) continue;
    err += MPI_Irecv(&n_job_recv[i],1,MPI_INT,i,MPI_ANY_TAG,
             ic->comm,&req_r[idx]);
    idx++;
  }
  /* finish communication */
  err += MPI_Waitall(nproc_macro,req_r,MPI_STATUS_IGNORE);
  free(req_r);

  /* build comm info based on recv data */
  err += build_PGFEM_comm_info(nproc_inter,n_job_recv,&(ic->send_info));
  err += build_PGFEM_comm_info(nproc_inter,n_job_recv,&(ic->recv_info));

  /* Communicate the buffer sizes from each process communicated with */
  req_r = PGFEM_calloc(MPI_Request, ic->recv_info->n_proc);
  for(int i=0; i<ic->recv_info->n_proc; i++){
    err += MPI_Irecv(ic->recv_info->buff_sizes + ic->recv_info->idx[i],
             ic->recv_info->n_buff[i],MPI_INT,
             ic->recv_info->proc[i],MPI_ANY_TAG,
             ic->comm,&req_r[i]);
  }

  err += MPI_Waitall(ic->recv_info->n_proc,req_r,MPI_STATUS_IGNORE);

  /* copy the buffer sizes from the receive to the send */
  for(int i=0; i<ic->recv_info->n_proc; i++){
    int *src = ic->recv_info->buff_sizes + ic->recv_info->idx[i];
    int *dest = ic->send_info->buff_sizes + ic->send_info->idx[i];
    int len = ic->recv_info->n_buff[i];
    memcpy(dest,src,len*sizeof(int));
  }

  /* Abort if work group without job */
  {
    int n_comm = 0;
    err += PGFEM_comm_info_get_n_comms(ic->recv_info,&n_comm);
    if(n_comm <= 0){
      /* print to master stderr */
      PGFEM_redirect_io_macro();
      PGFEM_printerr("[%d (inter)]ERROR: work group has no work!\n"
             "Consider using fewer groups or reducing"
             " size(comm->macro)\n\n",rank_inter);
      /* print to my (micro) stderr */
      PGFEM_redirect_io_micro();
      PGFEM_printerr("[%d (inter)]ERROR: work group has no work!\n"
             "Consider using fewer groups or reducing"
             " size(comm->macro)\n\n",rank_inter);
      PGFEM_Abort();
    }
  }

  /* cleanup */
  free(n_job_recv);
  free(req_r);

  return err;
}

static int create_ms_job_intercomm_macro(const int n_jobs,
                     const int *job_buff_sizes,
                     PGFEM_ms_job_intercomm *ic,
                     MSNET_Comm macro_comm)
{
  int err = 0;
  int nproc_macro = 0;
  int rank_macro = 0;
  err += MPI_Comm_rank(macro_comm,&rank_macro);
  err += MPI_Comm_size(macro_comm,&nproc_macro);

  int nproc_inter = 0;
  int rank_inter = 0;
  err += MPI_Comm_rank(ic->comm,&rank_inter);
  err += MPI_Comm_size(ic->comm,&nproc_inter);

  /* the number of micro workers is... */
  const int nproc_micro = nproc_inter - nproc_macro;

  /* initialize send */
  int *n_job_send = PGFEM_calloc(int, nproc_inter);

  /* get number of jobs from all procs on the macroscale */
  n_job_send[rank_inter] = n_jobs;
  MPI_Allgather(MPI_IN_PLACE,1,MPI_INT,n_job_send,1,MPI_INT,macro_comm);

  {
    /* begin with first microscale proc id */
    int idx = nproc_macro;
    /* loop through procs with lower rank */
    for(int p=0; p<rank_macro; p++){
      /* increment start_idx, resetting when gets to end */
      for(int i=0; i<n_job_send[p]; i++){
    idx++;
    if(idx >= nproc_inter) idx = nproc_macro;
      }
    }

    /*  allocate jobs from this process */
    int n_jobs_added = 0;
    while(n_jobs_added < n_jobs){
      /* wrap index space to stay in micro range */
      if(idx >= nproc_inter) idx = nproc_macro;
      n_job_send[idx]++;
      n_jobs_added++;
      idx++;
    }
  }

  /* set n_job_send[0-nproc_macro-1] = 0. Only sending to microscale! */
  memset(n_job_send,0,nproc_macro*sizeof(int));

  MPI_Request *req_s = PGFEM_calloc(MPI_Request, nproc_inter-1);
  int idx = 0;
  for(int i=0; i<nproc_inter; i++){
    if(i == rank_inter) continue;
    err += MPI_Isend(&n_job_send[i],1,MPI_INT,i,i /*tag*/,
             ic->comm,&req_s[idx]);
    idx++;
  }

  /* build comm info based on send data */
  err += build_PGFEM_comm_info(nproc_inter,n_job_send,&(ic->send_info));
  err += build_PGFEM_comm_info(nproc_inter,n_job_send,&(ic->recv_info));

  /* get the buffer sizes */
  memcpy(ic->send_info->buff_sizes,job_buff_sizes,n_jobs*sizeof(int));
  memcpy(ic->recv_info->buff_sizes,job_buff_sizes,n_jobs*sizeof(int));

  /* finish communication */
  err += MPI_Waitall(nproc_inter-1,req_s,MPI_STATUS_IGNORE);
  free(req_s);

  /* communicate buffer sizes */
  if(ic->send_info->n_proc > 0){
    req_s = PGFEM_calloc(MPI_Request, ic->send_info->n_proc);
    for(int i=0; i<ic->send_info->n_proc; i++){
      err += MPI_Isend(ic->send_info->buff_sizes + ic->send_info->idx[i],
               ic->send_info->n_buff[i],MPI_INT,
               ic->send_info->proc[i],i /*tag*/,
               ic->comm,&req_s[i]);
    }
    err += MPI_Waitall(ic->send_info->n_proc,req_s,MPI_STATUS_IGNORE);
    free(req_s);
  }

  /* reduce the number of jobs for logging */
  err += MPI_Allreduce(MPI_IN_PLACE,n_job_send+nproc_macro,nproc_micro,
               MPI_INT,MPI_SUM,macro_comm);

  if(rank_macro == 0){
    PGFEM_printf("Computing microscale contribution on %d groups\n",
         nproc_micro);
    PGFEM_printf("NJob on each group: ");
    print_array_i(PGFEM_stdout,n_job_send+nproc_macro,nproc_micro,
          1,nproc_micro);
  }

  free(n_job_send);
  return err;
}
