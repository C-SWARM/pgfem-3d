/* HEADER */
/*** This file contains wrapper functions for MPI calls. If the
     communication framework is changed, this is where to start */
#pragma once
#ifndef PGFEM_MPI_H
#define PGFEM_MPI_H

#include <config.h>

#include "mpi.h"

/** Future, write #define wrappers for MPI functions and
    find/replace in code */

/** Wrappers for MPI_Abort. Allows for simple change
    throughout code for different abort mechanism */
#define PGFEM_Comm_code_abort(mpi_comm,code)  MPI_Abort((mpi_comm),(code))
#define PGFEM_Comm_abort(mpi_comm) PGFEM_Comm_code_abort((mpi_comm),0)
#define PGFEM_Abort() PGFEM_Comm_code_abort(MPI_COMM_WORLD,0);

/** Get the rank on MPI_COMM_WORLD for error messages */
#define PGFEM_Error_rank(myrank) MPI_Comm_rank(MPI_COMM_WORLD,(myrank))

  /** Structure for passing around communicators etc. and managing
      macro/micro execution */
  typedef struct PGFEM_mpi_comm{
    /* communicators */
    MPI_Comm world;
    MPI_Comm macro;
    MPI_Comm micro_all;
    MPI_Comm micro;
    MPI_Comm mm_inter;
    MPI_Comm worker_inter; /* equivalent workers on different
                  microscales */

    /* stored ranks on respective communicators */
    int rank_world;
    int rank_macro;
    int rank_micro_all;
    int rank_micro;
    /* NOTE: micro procs have rank_mm_inter >= size(macro) */
    int rank_mm_inter;
    int server_id; /* rank on worker_inter */

    /* boolean valid communicator flag.
       NOTE: world is *ALWAYS* valid */
    int valid_macro;
    int valid_micro_all;
    int valid_micro;
    int valid_mm_inter;

  } PGFEM_mpi_comm;

  /** Initialize PGFEM_mpi_comm to MPI_COMM_WORLD, or whatever is passed
      as comm_world */
  int initialize_PGFEM_mpi_comm(const MPI_Comm comm_world,
                PGFEM_mpi_comm *pgfem_mpi_comm);

  /** Destroy a PGFEM_mpi_comm object */
  int destroy_PGFEM_mpi_comm(PGFEM_mpi_comm *pgfem_mpi_comm);

  /** Split up PGFEM_mpi_comm->world for macro-micro */
  int PGFEM_mpi_comm_MM_split(const int macro_nproc,
                  const int micro_group_size,
                  PGFEM_mpi_comm *pgfem_mpi_comm);

  /** Common structure for building a point-to-point communication
      graph. This is a basic building block that can be frequently
      used */
  typedef struct PGFEM_comm_info {
    int n_proc; /**< number of procs to comm with */
    int *proc; /**< proc ids to comm with [n_proc] */
    int *idx; /**< partial sum to index into buff_size [n_proc] */
    int *n_buff; /**< number of buffers to send to each proc [n_proc] */
    int *buff_sizes; /**< size of each buffer sum(n_buff) */
  } PGFEM_comm_info;

  /** build the PGFEM_comm_info object given the number of buffers
      communicated with each prcess. NOTE: the buffer sizes are
      allocated but are not computed */
  int build_PGFEM_comm_info(const int nproc,
                const int *n_buff_proc,
                PGFEM_comm_info **info);

  /** destroy a PGFEM_comm_info object */
  int destroy_PGFEM_comm_info(PGFEM_comm_info *info);

  /** expand a PGFEM_comm_info structure into a list of communications
      to carry out */
  int PGFEM_comm_info_to_idx_list(const PGFEM_comm_info *info,
                  int *n_comms,
                  int **procs,
                  int **sizes);

  /** compute the number of communications in PGFEM_comm_info */
  int PGFEM_comm_info_get_n_comms(const PGFEM_comm_info *info,
                  int *n_comms);

  /** a context for a server-type operation */
  typedef struct PGFEM_server_ctx{
    int n_comms; /**< How many communications */
    int *procs; /**< where they are going */
    int *sizes; /**< how big the messages are (in bytes) */
    int *tags; /**< tag (default MPI_ANY_TAG) */
    char **buffer; /**< buffer for communication */
    int in_process; /**< flag for if communication is in process */
    MPI_Request *req;
    MPI_Status *stat;
  } PGFEM_server_ctx;

  /** Initialize a PGFEM_server_ctx object */
  int initialize_PGFEM_server_ctx(PGFEM_server_ctx *ctx);

  /**
   * Build a server context
   */
  void build_PGFEM_server_ctx(PGFEM_server_ctx *ctx,
                  const int n_comm,
                  const int *buf_sizes);

  /** construct a PGFEM_server_ctx object based on a PGFEM_comm_info
      object */
  int build_PGFEM_server_ctx_from_PGFEM_comm_info(const PGFEM_comm_info
                          *info,
                          PGFEM_server_ctx *ctx);

  /**
   * Return the indices for communication associated with tag
   * (messages with MPI_ANY_TAG match any tag). The number of matches
   * and their indices are returned in count and indices
   * respectively. Note, indices should be an array at least as large
   * as the number of cummunications described by ctx.
   */
  int PGFEM_server_ctx_get_idx_from_tag(const PGFEM_server_ctx *ctx,
                    const int tag,
                    int *count,
                    int *indices);

  /**
   * Set the processor id for the message described at index idx.
   */
  int PGFEM_server_ctx_set_proc_at_idx(PGFEM_server_ctx *ctx,
                       const int proc,
                       const int idx);

  /**
   * Set the tag for the message described at index idx.
   */
  int PGFEM_server_ctx_set_tag_at_idx(PGFEM_server_ctx *ctx,
                      const int tag,
                      const int idx);

  /**
   * Get message info from server context at index. Note that the
   * buffer may be modified through the returned pointer.
   */
  int PGFEM_server_ctx_get_message(PGFEM_server_ctx *ctx,
                   const int idx,
                   void *buf,
                   int *n_bytes,
                   int *proc,
                   int *tag,
                   MPI_Request *req);

  /** Destroy a PGFEM_server_ctx object */
  int destroy_PGFEM_server_ctx(PGFEM_server_ctx *ctx);

#endif /* #ifndef  */
