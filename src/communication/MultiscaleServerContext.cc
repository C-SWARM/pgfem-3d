#include "pgfem3d/Communication.hpp"
#include "pgfem3d/MultiscaleCommunication.hpp"
#include "allocation.h"

using namespace pgfem3d::net;

namespace pgfem3d {

  MultiscaleServerContext::MultiscaleServerContext(net::Network *n)
  {
    n_comms = 0;
    procs = NULL;
    sizes = NULL;
    tags = NULL;
    buffer = NULL;
    req = NULL;
    stat = NULL;
    in_process = 0;
    net = n;
  }
  
  MultiscaleServerContext::~MultiscaleServerContext()
  {
    free(procs);
    free(sizes);
    free(tags);
    delete [] req;
    delete [] stat;
    for(int i=0; i<n_comms; i++)
      free(buffer[i]);
    free(buffer);
  }

  void MultiscaleServerContext::initialize(const int n_comm,
					   const int *buf_sizes)
  {
    n_comms = n_comm;
    procs = PGFEM_calloc(int, n_comm);
    sizes = PGFEM_calloc(int, n_comm);
    tags = PGFEM_calloc(int, n_comm);
    buffer = PGFEM_calloc(char*, n_comm);

    net->allocRequestArray(n_comm, &req);
    net->allocStatusArray(n_comm, &stat);
    
    /* set tags to NET_ANY_TAG, allocate individual buffers, and copy
       buffer sizes */
    for(int i=0; i<n_comm; i++){
      sizes[i] = buf_sizes[i];
      buffer[i] = static_cast<char*>(malloc(buf_sizes[i]));
      tags[i] = NET_ANY_TAG;
    }
  }

  void MultiscaleServerContext::initialize(MultiscaleCommInfo *minfo)
  {
    minfo->to_idx_list(&n_comms, &procs, &sizes);
  
    /* initialize pointers */
    req = NULL;
    stat = NULL;
    buffer = NULL;
    
    /* allocate if necessary */
    if(n_comms > 0){
      tags = PGFEM_calloc(int, n_comms);
      net->allocRequestArray(n_comms, &req);
      net->allocStatusArray(n_comms, &stat);
      
      buffer = PGFEM_calloc(char*, n_comms);
      for(int i=0; i<n_comms; i++){
	tags[i] = NET_ANY_TAG;
	buffer[i] = PGFEM_calloc(char, sizes[i]);
	req[i].setData(NULL);
      }
    }
  }
  
  void MultiscaleServerContext::get_idx_from_tag(const int tag,
						 int *count,
						 int *indices)
  {
    int c = 0;
    for(int i=0, e=n_comms; i<e; i++){
      if(tags[i] == tag || tags[i] == NET_ANY_TAG){
	indices[c] = i;
	c++;
      }
    }
    *count = c;
  }
  
  void MultiscaleServerContext::set_proc_at_idx(const int proc,
						const int idx)
  {
    if (idx >= n_comms) {
      throw;
    }
    /* can modify if comm is in process since it does not affect already
       posted comm. */
    procs[idx] = proc;
  }
  
  void MultiscaleServerContext::set_tag_at_idx(const int tag,
					       const int idx)
  {
    if (idx >= n_comms) {
      throw;
    }
    /* can modify if comm is in process since it does not affect already
       posted comm. */
    tags[idx] = tag;
  }
  

  void MultiscaleServerContext::get_message(const int idx,
					    void *buf,
					    int *n_bytes,
					    int *proc,
					    int *tag,
					    Request **r)
  {
    if (idx >= n_comms) {
      throw;
    }

    // This method is never used and this
    // assignment needs to be fixed if it ever does
    //buf = &buffer[idx];
    
    *n_bytes = sizes[idx];
    *proc = procs[idx];
    *tag = tags[idx];
    *r = req + idx;
  }
  
} // end namespace pgfem3d
  
