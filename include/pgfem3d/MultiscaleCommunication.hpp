#ifndef PGFEM3D_MULTISCALE_COMMUNICATION_H
#define PGFEM3D_MULTISCALE_COMMUNICATION_H

/// @brief This file defines the multiscale communication abstractions
#include "pgfem3d/Network.hpp"
#include <cstdio>

namespace pgfem3d {
  
// The Multiscale Communicator abstraction
class MultiscaleComm {
public:
  /* communicators */
  net::PGFem3D_Comm world;
  net::PGFem3D_Comm macro;
  net::PGFem3D_Comm micro_all;
  net::PGFem3D_Comm micro;
  net::PGFem3D_Comm mm_inter;
  net::PGFem3D_Comm worker_inter; /* equivalent workers on different
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
};
  
class MultiscaleCommInfo {
public:
  MultiscaleCommInfo(int nproc, int *n_buff_proc);
  ~MultiscaleCommInfo();
  
  /** expand this objectinto a list of communications
      to carry out */
  int to_idx_list(int *n_comms,
		  int **procs,
		  int **sizes);
  
  /** compute the number of communications in this object */
  int get_n_comms(int *n_comms);
  
private:
  int n_proc;      /**< number of procs to comm with */
  int *proc;       /**< proc ids to comm with [n_proc] */
  int *idx;        /**< partial sum to index into buff_size [n_proc] */
  int *n_buff;     /**< number of buffers to send to each proc [n_proc] */
  int *buff_sizes; /**< size of each buffer sum(n_buff) */
};

class MultiscaleServerContext {
public:
  MultiscaleServerContext(const int n_comm, const int *buf_sizes);
  ~MultiscaleServerContext();
  
  /**
   * Return the indices for communication associated with tag
   * (messages with MPI_ANY_TAG match any tag). The number of matches
   * and their indices are returned in count and indices
   * respectively. Note, indices should be an array at least as large
   * as the number of cummunications described by ctx.
   */
  int get_idx_from_tag(const int tag,
		       int *count,
		       int *indices);
  
  /**
   * Set the processor id for the message described at index idx.
   */
  int set_proc_at_idx(const int proc,
		      const int idx);
  
  /**
   * Set the tag for the message described at index idx.
   */
  int set_tag_at_idx(const int tag,
		     const int idx);
  
  /**
   * Get message info from server context at index. Note that the
   * buffer may be modified through the returned pointer.
   */
  int get_message(const int idx,
		  void *buf,
		  int *n_bytes,
		  int *proc,
		  int *tag);
  
private:
  int n_comms; /**< How many communications */
  int *procs; /**< where they are going */
  int *sizes; /**< how big the messages are (in bytes) */
  int *tags; /**< tag (default MPI_ANY_TAG) */
  char **buffer; /**< buffer for communication */
  int in_process; /**< flag for if communication is in process */
};
}
  
#endif
