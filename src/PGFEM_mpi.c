/* HEADER */
#include "PGFEM_mpi.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

/** Initialize PGFEM_mpi_comm to MPI_COMM_WORLD, or whatever is passed
    as comm_world */
int initialize_PGFEM_mpi_comm(const MPI_Comm comm_world,
			      PGFEM_mpi_comm *pgfem_mpi_comm)
{
  int err = 0;
  int myrank = 0;

  /* duplicate communicator */
  err += MPI_Comm_dup(comm_world,&pgfem_mpi_comm->world);
  pgfem_mpi_comm->macro = pgfem_mpi_comm->world;
  pgfem_mpi_comm->micro_all = pgfem_mpi_comm->world;
  pgfem_mpi_comm->micro = pgfem_mpi_comm->world;
  pgfem_mpi_comm->mm_inter = pgfem_mpi_comm->world;

  /* store rank */
  err += MPI_Comm_rank(comm_world,&myrank);
  pgfem_mpi_comm->rank_world = myrank;
  pgfem_mpi_comm->rank_macro = myrank;
  pgfem_mpi_comm->rank_micro_all = myrank;
  pgfem_mpi_comm->rank_micro = myrank;
  pgfem_mpi_comm->rank_mm_inter = myrank;

  /* set valid flag */
  pgfem_mpi_comm->valid_macro = 1;
  pgfem_mpi_comm->valid_micro_all = 1;
  pgfem_mpi_comm->valid_micro = 1;
  pgfem_mpi_comm->valid_mm_inter = 1;

  return err;
}

int destroy_PGFEM_mpi_comm(PGFEM_mpi_comm *pgfem_mpi_comm)
{
  int err = 0;
  err += MPI_Comm_free(&pgfem_mpi_comm->world);

  if(pgfem_mpi_comm->valid_macro){
    err += MPI_Comm_free(&pgfem_mpi_comm->macro);
  }

  if(pgfem_mpi_comm->valid_micro_all){
    err += MPI_Comm_free(&pgfem_mpi_comm->micro_all);
  }

  if(pgfem_mpi_comm->valid_micro){
    err += MPI_Comm_free(&pgfem_mpi_comm->micro);
  }

  if(pgfem_mpi_comm->valid_mm_inter){
    err += MPI_Comm_free(&pgfem_mpi_comm->mm_inter);
  }

  return err;
}

int PGFEM_mpi_comm_MM_split(const int macro_nproc,
			    const int micro_group_size,
			    PGFEM_mpi_comm *pgfem_mpi_comm)
{
  int err = 0;
  int nproc_world = 0;
  err += MPI_Comm_size(pgfem_mpi_comm->world,&nproc_world);

  /* error check group sizes */
  const int micro_nproc_all = nproc_world - macro_nproc;
  if(micro_nproc_all < micro_group_size
     || (micro_nproc_all%micro_group_size != 0)){
    if(pgfem_mpi_comm->rank_world == 0){
      PGFEM_printerr("ERROR: incorrect group sizes in %s!\n",__func__);
      PGFEM_printerr("WORLD = %d || MACRO = %d || GROUP = %d\n",
		     nproc_world,macro_nproc,micro_group_size);
      PGFEM_Abort();
    }
  }

  /* split into macro and micro communicators */
  {
    int rank_world = pgfem_mpi_comm->rank_world;
    const int macro_color = (rank_world < macro_nproc) ? 1 : MPI_UNDEFINED;
    err += MPI_Comm_split(pgfem_mpi_comm->world,macro_color,
			  rank_world,&(pgfem_mpi_comm->macro));

    const int micro_color = (rank_world < macro_nproc) ? MPI_UNDEFINED : 1;
    err += MPI_Comm_split(pgfem_mpi_comm->world,micro_color,
			  rank_world,&(pgfem_mpi_comm->micro_all));
  }

  /* get ranks on new communicators */
  if(pgfem_mpi_comm->macro == MPI_COMM_NULL){
    pgfem_mpi_comm->valid_macro = 0;
    pgfem_mpi_comm->rank_macro = MPI_UNDEFINED;
  } else {
    err += MPI_Comm_rank(pgfem_mpi_comm->macro,
		       &(pgfem_mpi_comm->rank_macro));
  }

  if(pgfem_mpi_comm->micro_all == MPI_COMM_NULL){
    pgfem_mpi_comm->valid_micro_all = 0;
    pgfem_mpi_comm->rank_micro_all = MPI_UNDEFINED;
  } else {
    err += MPI_Comm_rank(pgfem_mpi_comm->micro_all,
			 &(pgfem_mpi_comm->rank_micro_all));
  }

  /* split micro communicator into work groups */
  if(pgfem_mpi_comm->valid_micro_all && micro_group_size > 0){
    int rank_micro_all = pgfem_mpi_comm->rank_micro_all;
    int color = MPI_UNDEFINED;
    /* integer division for group id */
    if(rank_micro_all != MPI_UNDEFINED){
      color = rank_micro_all / micro_group_size;
    }
    err += MPI_Comm_split(pgfem_mpi_comm->micro_all,color,
			  rank_micro_all,&(pgfem_mpi_comm->micro));
    err += MPI_Comm_rank(pgfem_mpi_comm->micro,
			 &(pgfem_mpi_comm->rank_micro));
  } else {
    pgfem_mpi_comm->micro = MPI_COMM_NULL;
    pgfem_mpi_comm->valid_micro = 0;
    pgfem_mpi_comm->rank_micro = MPI_UNDEFINED;
  }

  /* create the micro-macro intercommunicator */
  {
    int rank_world = pgfem_mpi_comm->rank_world;
    int color = MPI_UNDEFINED;
    if(pgfem_mpi_comm->rank_macro != MPI_UNDEFINED
       || pgfem_mpi_comm->rank_micro == 0){
      color = 1;
    }
    err += MPI_Comm_split(pgfem_mpi_comm->world,color,
			  rank_world,&(pgfem_mpi_comm->mm_inter));
  }

  /* get rank on new communicator */
  if(pgfem_mpi_comm->mm_inter == MPI_COMM_NULL){
    pgfem_mpi_comm->valid_mm_inter = 0;
    pgfem_mpi_comm->rank_mm_inter = MPI_UNDEFINED;
  } else {
    err += MPI_Comm_rank(pgfem_mpi_comm->mm_inter,
			 &(pgfem_mpi_comm->rank_mm_inter));
  }

  /* Reduce error on mpi_comm_world */
  MPI_Allreduce(MPI_IN_PLACE,&err,1,MPI_INT,MPI_BOR,
		pgfem_mpi_comm->world);
  return err;
}

int build_PGFEM_comm_info(const int nproc,
			  const int *n_buff_proc,
			  PGFEM_comm_info **info)
{
  int err = 0;
  *info = PGFEM_calloc(1,sizeof(PGFEM_comm_info));
  PGFEM_comm_info *I = *info; /* alias */

  /* initialize values */
  I->n_proc = 0;
  I->proc = NULL;
  I->idx = NULL;
  I->n_buff = NULL;
  I->buff_sizes = NULL;

  /* compute number of procs to communicate with */
  int n_buff = 0;
  for(int i=0; i<nproc; i++){
   n_buff += n_buff_proc[i];
    if(n_buff_proc[i] > 0) I->n_proc++;
  }

  /* allocate space */
  if(I->n_proc > 0){
    I->proc = PGFEM_calloc(I->n_proc,sizeof(int));
    I->idx = PGFEM_calloc(I->n_proc + 1,sizeof(int));
    I->n_buff = PGFEM_calloc(I->n_proc,sizeof(int));
  }
  if(n_buff > 0){
    I->buff_sizes = PGFEM_calloc(n_buff,sizeof(int));
  }

  /* populate */
  int idx = 0;
  for(int i=0; i<nproc; i++){
    if(n_buff_proc[i] > 0){
      I->proc[idx] = i;
      I->n_buff[idx] = n_buff_proc[i];
      I->idx[idx+1] = I->idx[idx] + n_buff_proc[i];
      idx++;
    }
  }
  return err;
}

int destroy_PGFEM_comm_info(PGFEM_comm_info *info)
{
  int err = 0;
  if(info != NULL){
    free(info->proc);
    free(info->idx);
    free(info->n_buff);
    free(info->buff_sizes);
  }
  return err;
}

int PGFEM_comm_info_to_idx_list(const PGFEM_comm_info *info,
				int *n_comms,
				int **procs,
				int **sizes)
{
  int err = 0;
  /* compute the number of communications */
  err += PGFEM_comm_info_get_n_comms(info,n_comms);
  *procs = NULL;
  *sizes = NULL;
  if(*n_comms > 0){
    *procs = PGFEM_calloc(*n_comms,sizeof(int));
    *sizes = PGFEM_calloc(*n_comms,sizeof(int));
  }

  /* aliases */
  int *p = *procs;
  int *s = *sizes;
  int idx = 0;
  for(int i=0; i<info->n_proc; i++){
    const int proc = info->proc[i];
    const int buff_start = info->idx[i];
    const int n_buff = info->n_buff[i];
    for(int j=0; j<n_buff; j++){
      p[idx] = proc;
      s[idx] = info->buff_sizes[buff_start + j];
      idx++;
    }
  }

  /* error checking */
  if(idx != *n_comms) err++;
  return err;
}

inline int PGFEM_comm_info_get_n_comms(const PGFEM_comm_info *info,
				       int *n_comms)
{
  int err = 0;
  if(info->n_proc == 0) *n_comms = 0;
  else *n_comms = info->idx[info->n_proc];
  return err;
}

int initialize_PGFEM_server_ctx(PGFEM_server_ctx *ctx)
{
  int err = 0;
  if(ctx != NULL){
    ctx->n_comms = 0;
    ctx->procs = NULL;
    ctx->sizes = NULL;
    ctx->buffer = NULL;
    ctx->req = NULL;
    ctx->stat = NULL;
  } else err++;
  return err;
}

int build_PGFEM_server_ctx_from_PGFEM_comm_info(const PGFEM_comm_info *info,
						PGFEM_server_ctx *ctx)
{
  int err = 0;
  err += PGFEM_comm_info_to_idx_list(info,&(ctx->n_comms),
				     &(ctx->procs),&(ctx->sizes));

  /* initialize pointers */
  ctx->req = NULL;
  ctx->stat = NULL;
  ctx->buffer = NULL;

  /* allocate if necessary */
  if(ctx->n_comms > 0){
    ctx->req = PGFEM_calloc(ctx->n_comms,sizeof(MPI_Request));
    ctx->stat = PGFEM_calloc(ctx->n_comms,sizeof(MPI_Status));

    ctx->buffer = PGFEM_calloc(ctx->n_comms,sizeof(char*));
    for(int i=0; i<ctx->n_comms; i++){
      ctx->buffer[i] = PGFEM_calloc(ctx->sizes[i],sizeof(char));
      ctx->req[i] = MPI_REQUEST_NULL;
    }
  }
  return err;
}

int destroy_PGFEM_server_ctx(PGFEM_server_ctx *ctx)
{
  int err = 0;
  if(ctx != NULL){
    free(ctx->procs);
    free(ctx->sizes);
    free(ctx->req);
    free(ctx->stat);
    for(int i=0; i<ctx->n_comms; i++) free(ctx->buffer[i]);
    free(ctx->buffer);
  }
  return err;
}
