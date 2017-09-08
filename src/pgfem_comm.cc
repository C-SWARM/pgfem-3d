/* HEADER */
/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 */

#include "allocation.h"
#include "pgfem_comm.h"
#include "PGFEM_io.h"

using pgfem3d::solvers::SparseSystem;

/**
 * Structure for containing an index map.
 *
 * \param n_map The number of mapped indices.
 * \param maps The index mapping stored in two consecutive
 * elements. Example: the i-th map is INDEX maps[i*2] maps to INDEX
 * maps[i*2+1].
 */
struct fast_map {
  int n_maps = 0;             /**< number of maps */
  size_t *maps = nullptr;     /**< mapping stored in two consecutive elements */
};

static void
destroy_fast_map(fast_map *map)
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
  *Lk = PGFEM_calloc(double*, nproc);
  *receive = PGFEM_calloc(double*, nproc);
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
    (*Lk)[i] = PGFEM_calloc(double, n_alloc_Lk);
    (*receive)[i] = PGFEM_calloc(double, n_alloc_recv);
  }

  /* Allocate Status and Requests */
  *sta_r = NULL;
  *req_r = NULL;
  if(pgfem_comm->Nr > 0){
    *sta_r = PGFEM_calloc(MPI_Status, pgfem_comm->Nr);
    *req_r = PGFEM_calloc(MPI_Request, pgfem_comm->Nr);
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
    *sta_s = PGFEM_calloc(MPI_Status, pgfem_comm->Ns);
    *req_s = PGFEM_calloc(MPI_Request, pgfem_comm->Ns);
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


/** Assemble non-local parts as they arrive */
int
assemble_nonlocal_stiffmat(const COMMUN pgfem_comm,
                           MPI_Status *sta_r,
                           MPI_Request *req_r,
                           SparseSystem *system,
                           double **recv)
{
  int err = 0;
  int comm_idx = 0;
  int n_received = 0;
  while (n_received < pgfem_comm->Nr){
    /* get the communication index */
    err += MPI_Waitany(pgfem_comm->Nr,req_r,&comm_idx,sta_r);

    /* convert communication index to proc_id */
    const int proc = pgfem_comm->Nrr[comm_idx];

    /* get number of rows */
    const int nrows = pgfem_comm->R[proc];

    /* allocate rows and cols to receive */
    int *row_idx = PGFEM_calloc(int, nrows);
    int *ncols = PGFEM_calloc(int, nrows);
    int *col_idx = PGFEM_calloc(int, pgfem_comm->AR[proc]);

    /* get row and column ids */
    int idx = 0;
    for(int j=0; j<pgfem_comm->R[proc]; j++){
      row_idx[j] = pgfem_comm->RGID[proc][j];
      ncols[j] = pgfem_comm->RAp[proc][j];
      for(int k=0; k<ncols[j]; k++){
        col_idx[idx] = pgfem_comm->RGRId[proc][idx];
        ++idx;
      }
    }

    /* assemble to local part of global stiffness */
    system->add(nrows, ncols, row_idx, col_idx, recv[proc]);

    /* free memory */
    free(row_idx);
    free(ncols);
    free(col_idx);

    /* increment counter */
    n_received++;
  }
  return err;
}

static int
build_fast_LG_map(COMMUN pgfem_comm,
                  const long ndofd,
                  const long ngdof_owned,
                  const long start_gdof_id)
{
  int err = 0;
  pgfem_comm->fast_LG_map = PGFEM_calloc(fast_map, 1);
  fast_map *LG = pgfem_comm->fast_LG_map;
  LG->n_maps = 0;

  /* deterimine the local indices that map to global indices on the
     local domain and the corresponding number of maps */
  long *local_idx = PGFEM_calloc(long, ndofd);
  long *global_idx = PGFEM_calloc(long, ndofd);
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

  LG->maps = PGFEM_calloc(size_t, 2*LG->n_maps);
  for(int i=0; i<LG->n_maps; i++){
    LG->maps[2*i] = local_idx[i];
    LG->maps[2*i+1] = global_idx[local_idx[i]];
  }

  free(local_idx);
  free(global_idx);

  return err;
}

static int
build_fast_GL_map(COMMUN pgfem_comm,
                  const long ndofd,
                  const long ngdof_owned)
{
  int err = 0;
  pgfem_comm->fast_GL_map = PGFEM_calloc(fast_map, 1);
  fast_map *GL = pgfem_comm->fast_GL_map;
  GL->n_maps = 0;

  long *global_idx = PGFEM_calloc(long, ngdof_owned);
  long *local_idx = PGFEM_calloc(long, ngdof_owned);

  /* NOTE: may be able to replace. Should be identical to comm->GL? */
  for(int i=0; i< ngdof_owned; i++){
    local_idx[i] = pgfem_comm->GL[i];
    if(local_idx[i]>= 0 && local_idx[i] < ndofd){
      global_idx[GL->n_maps] = i;
      GL->n_maps++;
    }
  }

  GL->maps = PGFEM_calloc(size_t, 2*GL->n_maps);
  for(int i=0; i<GL->n_maps; i++){
    GL->maps[2*i] = global_idx[i];
    GL->maps[2*i+1] = local_idx[global_idx[i]];
  }

  free(global_idx);
  free(local_idx);

  return err;
}

int
pgfem_comm_build_fast_maps(COMMUN pgfem_comm,
                           const long ndofd,
                           const long ngdof_owned,
                           const long start_gdof_id)
{
  return (build_fast_LG_map(pgfem_comm, ndofd, ngdof_owned, start_gdof_id) +
          build_fast_GL_map(pgfem_comm,ndofd,ngdof_owned));
}

/* maps A to B */
static int
get_mapped_values(const fast_map *map, const double *A, double *B)
{
  for(int i=0; i<map->n_maps; i++){
    B[map->maps[2*i+1]] = A[map->maps[2*i]];
  }
  return 0;
}

int
pgfem_comm_get_owned_global_dof_values(const COMMUN pgfem_comm,
                                       const double *local_dofs,
                                       double *global_dofs)
{
  return get_mapped_values(pgfem_comm->fast_LG_map, local_dofs, global_dofs);
}

int
pgfem_comm_get_local_dof_values_from_global(const COMMUN pgfem_comm,
                                            const double *global_dofs,
                                            double *local_dofs)
{
  return get_mapped_values(pgfem_comm->fast_GL_map, global_dofs, local_dofs);
}
