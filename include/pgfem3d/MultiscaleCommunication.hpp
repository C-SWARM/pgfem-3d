#ifndef PGFEM3D_MULTISCALE_COMMUNICATION_H
#define PGFEM3D_MULTISCALE_COMMUNICATION_H

/// @brief This file defines the multiscale communication abstractions
#include "pgfem3d/Communication.hpp"
#include <cstdio>

namespace pgfem3d {
  
// The Multiscale Communicator abstraction
class MultiscaleComm {
public:
  MultiscaleComm(const net::PGFem3D_Comm comm_world,
		 net::Network *n);
  ~MultiscaleComm();

  void MM_split(const int macro_nproc,
		const int micro_group_size, const int micro2_group_size, const int full_micro_np);
  
  /* communicators */
  net::PGFem3D_Comm world;
  net::PGFem3D_Comm macro;
  net::PGFem3D_Comm micro_all;
  net::PGFem3D_Comm micro;
  net::PGFem3D_Comm mm_inter;
  net::PGFem3D_Comm mm_inter_ROM;
  net::PGFem3D_Comm micro_ROM; //full micro 2 group
//  net::PGFem3D_Comm micro_ROM; //work groups within micro 2
  net::PGFem3D_Comm worker_inter; /* equivalent workers on different
				microscales */
  net::PGFem3D_Comm worker_inter_ROM;

  /* stored ranks on respective communicators */
  int rank_world;
  int rank_macro;
  int rank_micro_all;
  int rank_micro;
  int rank_micro_1;
  int rank_micro_ROM;
  /* NOTE: micro procs have rank_mm_inter >= size(macro) */
  int rank_mm_inter;
  int rank_mm_inter_ROM;
  int server_id; /* rank on worker_inter */

  /* boolean valid communicator flag.
     NOTE: world is *ALWAYS* valid */
  int valid_macro;
  int valid_micro_all;
  int valid_micro;
  int valid_mm_inter;
  int valid_mm_inter_ROM;
  int valid_micro_1;
  int valid_micro_ROM; 
private:
  net::Network *net;  // handle to the active network
};
  
class MultiscaleCommInfo {
public:
  MultiscaleCommInfo(int nproc, int *n_buff_proc);
  ~MultiscaleCommInfo();
  
  /** expand this objectinto a list of communications
      to carry out */
  void to_idx_list(int *n_comms,
		   int **procs,
		   int **sizes);
    
private:
  void get_n_comms(int *n_comms) {
    if (n_proc == 0) *n_comms = 0;
    else *n_comms = idx[n_proc];
  }
  
  int n_proc;      /**< number of procs to comm with */
  int *proc;       /**< proc ids to comm with [n_proc] */
  int *idx;        /**< partial sum to index into buff_size [n_proc] */
  int *n_buff;     /**< number of buffers to send to each proc [n_proc] */
  int *buff_sizes; /**< size of each buffer sum(n_buff) */
};

class MultiscaleServerContext {
public:
  MultiscaleServerContext(net::Network *n);
  ~MultiscaleServerContext();

  void initialize(const int n_comm, const int *buf_sizes);
  void initialize(pgfem3d::MultiscaleCommInfo *minfo);
  
  /**
   * Return the indices for communication associated with tag
   * (messages with NET_ANY_TAG match any tag). The number of matches
   * and their indices are returned in count and indices
   * respectively. Note, indices should be an array at least as large
   * as the number of cummunications described by ctx.
   */
  void get_idx_from_tag(const int tag,
			int *count,
			int *indices);
  
  
  /**
   * Set the processor id for the message described at index idx.
   */
  void set_proc_at_idx(const int proc,
		       const int idx);
  
  /**
   * Set the tag for the message described at index idx.
   */
  void set_tag_at_idx(const int tag,
		      const int idx);
  
  /**
   * Get message info from server context at index. Note that the
   * buffer may be modified through the returned pointer.
   */
  void get_message(const int idx,
		   void *buf,
		   int *n_bytes,
		   int *proc,
		   int *tag,
		   net::Request **r);

  // XXX should encapsulate all the send/recv functionality in the class
  // to protect these vars
  int n_comms; /**< How many communications */
  int *tags; /**< tag (default NET_ANY_TAG) */
  char **buffer; /**< buffer for communication */
  int *procs; /**< where they are going */
  int *sizes; /**< how big the messages are (in bytes) */
  int in_process; /**< flag for if communication is in process */

  net::Request *req;
  net::Status *stat;

private:
  net::Network *net;
};

} //end namespace pgfem3d
  
#endif
