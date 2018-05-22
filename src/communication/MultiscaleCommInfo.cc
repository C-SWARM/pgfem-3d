#include "pgfem3d/MultiscaleCommunication.hpp"
#include "allocation.h"
#include <cassert>

using namespace pgfem3d::net;

namespace pgfem3d {

  MultiscaleCommInfo::MultiscaleCommInfo(int nproc, int *n_buff_proc)
  {
    /* initialize values */
    n_proc = 0;
    proc = NULL;
    idx = NULL;
    n_buff = NULL;
    buff_sizes = NULL;
    
    /* compute number of procs to communicate with */
    int ln_buff = 0;
    for (int i=0; i<nproc; i++){
      ln_buff += n_buff_proc[i];
      if (n_buff_proc[i] > 0) n_proc++;
    }
    
    /* allocate space */
    if (n_proc > 0){
      proc = PGFEM_calloc(int, n_proc);
      idx = PGFEM_calloc(int, n_proc + 1);
      n_buff = PGFEM_calloc(int, n_proc);
    }
    if (ln_buff > 0){
      buff_sizes = PGFEM_calloc(int, ln_buff);
    }
    
    /* populate */
    int lidx = 0;
    for (int i=0; i<nproc; i++){
      if (n_buff_proc[i] > 0){
	proc[lidx] = i;
	n_buff[lidx] = n_buff_proc[i];
	idx[lidx+1] = idx[lidx] + n_buff_proc[i];
	lidx++;
      }
    }
  }
  
  MultiscaleCommInfo::~MultiscaleCommInfo()
  {
    free(proc);
    free(idx);
    free(n_buff);
    free(buff_sizes);
  }

  void MultiscaleCommInfo::to_idx_list(int *n_comms,
				       int **procs,
				       int **sizes)
  {
    /* compute the number of communications */
    get_n_comms(n_comms);
    *procs = NULL;
    *sizes = NULL;
    if (*n_comms > 0) {
      *procs = PGFEM_calloc(int, *n_comms);
      *sizes = PGFEM_calloc(int, *n_comms);
    }
    
    /* aliases */
    int *p = *procs;
    int *s = *sizes;
    int lidx = 0;
    for (int i=0; i<n_proc; i++){
      const int lp = proc[i];
      const int buff_start = idx[i];
      const int ln_buff = n_buff[i];
      for(int j=0; j<ln_buff; j++){
	assert(p != NULL && "p can't be null");
	p[lidx] = lp;
	s[lidx] = buff_sizes[buff_start + j];
	lidx++;
      }
    }
    
    /* error checking */
    if (lidx != *n_comms) {
      throw;
    }
  }
   
  
} // end namespace pgfem3d
