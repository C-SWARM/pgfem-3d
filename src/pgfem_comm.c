#include "pgfem_comm.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

void destroy_commun(COMMUN comm,long nproc)
{
  for(long i=0; i<nproc; i++){
    free(comm->SLID[i]);
    free(comm->RGID[i]);
    free(comm->SAp[i]);
    free(comm->SGRId[i]);
    free(comm->RAp[i]);
    free(comm->RGRId[i]);
  }
  free(comm->S);
  free(comm->R);
  free(comm->AS);
  free(comm->AR);
  free(comm->SLID);
  free(comm->RGID);
  free(comm->LG);
  free(comm->GL);
  free(comm->SAp);
  free(comm->SGRId);
  free(comm->RAp);
  free(comm->RGRId);
  free(comm->Nss);
  free(comm->Nrr);
  free(comm);
}

/** builds the off-process buffer for assembling the gloabal stiffness
    matrix. Posts non-blocking recieve */
int init_and_post_stiffmat_comm(double ***Lk,
				double ***receive,
				MPI_Request **req_r,
				MPI_Status **sta_r,
				const MPI_Comm mpi_comm,
				const COMMUN pgfem_comm)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);
  err += MPI_Comm_size(mpi_comm,&nproc);

  /* Allocate buffer for off-process assembly and receive */
  *Lk = PGFEM_calloc(nproc,sizeof(double));
  *receive = PGFEM_calloc(nproc,sizeof(double));
  for(int i=0; i<nproc; i++){
    size_t n_alloc_Lk = 0;
    size_t n_alloc_recv = 0;

    /* number to allocate for Lk */
    (*Lk)[i] = NULL;
    if(myrank == i || pgfem_comm->S[i] == 0){
      n_alloc_Lk = 1;
    } else {
      n_alloc_Lk = pgfem_comm->AS[i];
    }

    /* number to allocate for receive */
    (*receive)[i] = NULL;
    if(myrank == i || pgfem_comm->AR[i] == 0){
      n_alloc_recv = 1;
    } else {
      n_alloc_recv = pgfem_comm->AR[i];
    }

    /* allocate */
    (*Lk)[i] = PGFEM_calloc(n_alloc_Lk,sizeof(double));
    (*receive)[i] = PGFEM_calloc(n_alloc_recv,sizeof(double));
  }

  /* Allocate Status and Requests */
  *sta_r = NULL;
  *req_r = NULL;
  if(pgfem_comm->Nr > 0){
    *sta_r = PGFEM_calloc(pgfem_comm->Nr,sizeof(MPI_Status));
    *req_r = PGFEM_calloc(pgfem_comm->Nr,sizeof(MPI_Request));
  }

  /* post the non-blocking receive */
  for(int i=0; i<pgfem_comm->Nr; i++){
    size_t idx = pgfem_comm->Nrr[i];
    err += MPI_Irecv((*receive)[idx],pgfem_comm->AR[idx],MPI_DOUBLE,idx,
		     MPI_ANY_TAG,mpi_comm,*req_r + i);
  }

  return err;
}


/** Send the off-process data for stiffness assembly */
int send_stiffmat_comm(MPI_Status **sta_s,
		       MPI_Request **req_s,
		       double **Lk,
		       const MPI_Comm mpi_comm,
		       const COMMUN pgfem_comm)
{
  int err = 0; 
  int myrank = 0;
  int nproc = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);
  err += MPI_Comm_size(mpi_comm,&nproc);

  /* allocate status and request fields */
  *sta_s = NULL;
  *req_s = NULL;
  if(pgfem_comm->Ns > 0){
    *sta_s = PGFEM_calloc(pgfem_comm->Ns,sizeof(MPI_Status));
    *req_s = PGFEM_calloc(pgfem_comm->Ns,sizeof(MPI_Request));
  }

  /* post sends */
  for(int i=0; i<pgfem_comm->Ns; i++){
    size_t idx = pgfem_comm->Nss[i];
    err += MPI_Isend(Lk[idx],pgfem_comm->AS[idx],MPI_DOUBLE,idx,
		     myrank,mpi_comm,*req_s + i);
  }

  return err;
}

/** Finalize the stiffmat communication. Simple wrapper for the
    waitall commands. */
int finalize_stiffmat_comm(MPI_Status *sta_s,
			   MPI_Status *sta_r,
			   MPI_Request *req_s,
			   MPI_Request *req_r,
			   const COMMUN pgfem_comm)
{
  int err = 0;
  err += MPI_Waitall(pgfem_comm->Nr,req_r,sta_r);
  err += MPI_Waitall(pgfem_comm->Ns,req_s,sta_s);
  return err;
}
