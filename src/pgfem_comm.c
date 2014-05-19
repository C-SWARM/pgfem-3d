/* HEADER */
/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 */

#include "pgfem_comm.h"
#include "allocation.h"
#include "PGFEM_io.h"

/**
 * Structure for containing an index map.
 *
 * \param n_map The number of mapped indices.
 * \param maps The index mapping stored in two consecutive
 * elements. Example: the i-th map is INDEX maps[i*2] maps to INDEX
 * maps[i*2+1].
 */
typedef struct pgfem_comm_fast_map{
  int n_maps; /**< number of maps */
  size_t *maps; /**< mapping stored in two consecutive elements */
} fast_map;

static void destroy_fast_map(fast_map *map)
{
  free(map->maps);
}

void initialize_commun(COMMUN comm)
{
  comm->S = NULL;
  comm->R = NULL;
  comm->AS = NULL;
  comm->SLID = NULL;
  comm->RGID = NULL;
  comm->LG = NULL;
  comm->SAp = NULL;
  comm->SGRId = NULL;
  comm->RAp = NULL;
  comm->RGRId = NULL;
  comm->Nss = NULL;
  comm->Nrr = NULL;
  comm->fast_LG_map = NULL;
  comm->fast_GL_map = NULL;
}

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

  /* destroy maps */
  if(comm->fast_LG_map != NULL)
    destroy_fast_map((fast_map*) comm->fast_LG_map);
  if(comm->fast_GL_map != NULL)
    destroy_fast_map((fast_map*) comm->fast_GL_map);

  free(comm->fast_LG_map);
  free(comm->fast_GL_map);

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

static int build_fast_LG_map(COMMUN pgfem_comm,
			     const long ndofd,
			     const long ngdof_owned,
			     const long start_gdof_id)
{
  int err = 0;
  pgfem_comm->fast_LG_map = PGFEM_calloc(1,sizeof(fast_map));
  fast_map *LG = pgfem_comm->fast_LG_map;
  LG->n_maps = 0;

 /* deterimine the local indices that map to global indices on the
     local domain and the corresponding number of maps */
  long *local_idx = PGFEM_calloc(ndofd,sizeof(long));
  long *global_idx = PGFEM_calloc(ndofd,sizeof(long));
  for(int i=0; i<ndofd; i++){
    global_idx[i] = pgfem_comm->LG[i] - start_gdof_id;

    /* store mapping if I own the global dof */
    if(global_idx[i] >= 0 && global_idx[i] < ngdof_owned){
      local_idx[LG->n_maps] = i;
      LG->n_maps++;
    }
  }

  /* sanity check */
  if(LG->n_maps > ngdof_owned){
    int myrank = 0;
    PGFEM_Error_rank(&myrank);
    PGFEM_printerr("[%d] ERROR: too many mappings! %s:%s:%d\n",
		   myrank,__func__,__FILE__,__LINE__);
    PGFEM_Abort();
  }

  LG->maps = PGFEM_calloc(2*LG->n_maps,sizeof(size_t));
  for(int i=0; i<LG->n_maps; i++){
    LG->maps[2*i] = local_idx[i];
    LG->maps[2*i+1] = global_idx[local_idx[i]];
  }

  free(local_idx);
  free(global_idx);

  return err;
}

static int build_fast_GL_map(COMMUN pgfem_comm,
			     const long ndofd,
			     const long ngdof_owned)
{
  int err = 0;
  pgfem_comm->fast_GL_map = PGFEM_calloc(1,sizeof(fast_map));
  fast_map *GL = pgfem_comm->fast_GL_map;
  GL->n_maps = 0;

  long *global_idx = PGFEM_calloc(ngdof_owned,sizeof(long));
  long *local_idx = PGFEM_calloc(ngdof_owned,sizeof(long));

  /* NOTE: may be able to replace. Should be identical to comm->GL? */
  for(int i=0; i< ngdof_owned; i++){
    local_idx[i] = pgfem_comm->GL[i];
    if(local_idx[i]>= 0 && local_idx[i] < ndofd){
      global_idx[GL->n_maps] = i;
      GL->n_maps++;
    }
  }

  GL->maps = PGFEM_calloc(2*GL->n_maps,sizeof(size_t));
  for(int i=0; i<GL->n_maps; i++){
    GL->maps[2*i] = global_idx[i];
    GL->maps[2*i+1] = local_idx[global_idx[i]];
  }

  free(global_idx);
  free(local_idx);

  return err;
}

int pgfem_comm_build_fast_maps(COMMUN pgfem_comm,
			       const long ndofd,
			       const long ngdof_owned,
			       const long start_gdof_id)
{
  int err = 0;
 
  err += build_fast_LG_map(pgfem_comm,ndofd,ngdof_owned,
			   start_gdof_id);
  err += build_fast_GL_map(pgfem_comm,ndofd,ngdof_owned);

  return err;
}

/* maps A to B */
static int get_mapped_values(const fast_map *map,
			     const double *A,
			     double *B)
{
  int err = 0;

  for(int i=0; i<map->n_maps; i++){
    B[map->maps[2*i+1]] = A[map->maps[2*i]];
  }

  return err;
}

int pgfem_comm_get_owned_global_dof_values(const COMMUN pgfem_comm,
					   const double *local_dofs,
					   double *global_dofs)
{
  int err = 0;

  err += get_mapped_values(pgfem_comm->fast_LG_map,
			   local_dofs,
			   global_dofs);

  return err;
}

int pgfem_comm_get_local_dof_values_from_global(const COMMUN pgfem_comm,
						const double *global_dofs,
						double *local_dofs)
{
  int err = 0;

  err += get_mapped_values(pgfem_comm->fast_GL_map,
			   global_dofs,
			   local_dofs);

  return err;
}
