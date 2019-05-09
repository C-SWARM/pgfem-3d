/* HEADER */
/**
 * AUTHORS:
 *     , Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * 2018, Ezra Kissel, Indiana Univerity, ezkissel [at] indiana.edu
 */

#include "pgfem3d/Communication.hpp"
#include "ISIR_SparseComm.hpp"
#include "PWC_SparseComm.hpp"
#include "allocation.h"
#include "utils.h"
#include "PGFEM_io.h"
#include <errno.h>
#include <cerrno>
#include <system_error>

using namespace pgfem3d::net;

namespace pgfem3d {
  
  CommHints::CommHints(const char *ipath, const char* basefname,
		       const int rank)
  {
    const char format[] = "%s/%scomm_hints_%d.in";
    
    char *tempstr = PGFEM_malloc<char>(1);  //temporary ptr allocation to prevent gcc 8.1 warning
    int nchar = snprintf(tempstr, 0, format, ipath, basefname, rank);
    free(tempstr);
    
    filestr = PGFEM_malloc<char>(nchar + 1);
    sprintf(filestr, format, ipath, basefname, rank);

    send = recv = NULL;
    nsend = nrecv = 0;
  }

  CommHints::~CommHints() {
    if (filestr)
      free(filestr);
    if (send)
      PGFEM_free(send);
    if (recv)
      PGFEM_free(recv);
  }
  
  void CommHints::read_filename(const char* fn) {
    if (!fn) {
      fn = filestr;
    }
    FILE *in = fopen(fn, "r");
    if (!in) {
      throw std::system_error(errno, std::system_category());
    } else {
      this->read(in);
      fclose(in);
    }
  }

  void CommHints::write_filename(const char* fn) {
    if (!fn) {
      fn = filestr;
    }
    FILE *out = fopen(fn, "w");
    if (!out) {
      throw std::system_error(errno, std::system_category());
    } else {
      this->write(out);
      fclose(out);
    }
  }
  
  void CommHints::read(FILE *in)
  {
    /* read hints for whom to send information */
    CHECK_SCANF(in, "%d", &nsend);
    send = PGFEM_calloc(int, nsend);
    for (int i = 0, e = nsend; i < e; i++) {
      CHECK_SCANF(in, "%d", &send[i]);
    }
    
    /* read hints from whom to receive information */
    CHECK_SCANF(in, "%d", &nrecv);
    recv = PGFEM_calloc(int, nrecv);
    for (int i = 0, e = nrecv; i < e; i++) {
      CHECK_SCANF(in, "%d", &recv[i]);
    }
  }
  
  void CommHints::write(FILE *out)
  {
    fprintf(out, "%d\t", nsend);
    for (int i = 0, e = nsend; i < e; i++) {
      fprintf(out, "%3.0d ", send[i]);
    }
    fprintf(out, "\n");
    fprintf(out, "%d\t", nrecv);
    for (int i = 0, e = nrecv; i < e; i++) {
      fprintf(out, "%3.0d ", recv[i]);
    }
    fprintf(out, "\n");
  }

  SparseComm::~SparseComm()
  {
  }
  
  SparseComm* SparseComm::Create(net::Network *n, net::PGFem3D_Comm c)
  {
    int myrank;
    n->comm_rank(c, &myrank);
    switch (n->type()) {
    case NET_ISIR:
      return new pgfem3d::ISIR_SparseComm(dynamic_cast<ISIRNetwork*>(n), c);
      break;
    case NET_PWC:
      return new pgfem3d::PWC_SparseComm(dynamic_cast<PWCNetwork*>(n), c);
      break;
    default:
      PGFEM_printerr("[%d]ERROR: Unknown network type", myrank);
      PGFEM_Abort();
    }
  }
  
  void SparseComm::get_owned_global_dof_values(const double *local_dofs,
					       double *global_dofs)
  {
    get_mapped_values(fast_LG_map, local_dofs, global_dofs);
  }
  
  void SparseComm::get_local_dof_values_from_global(const double *global_dofs,
						    double *local_dofs)
  {
    get_mapped_values(fast_GL_map, global_dofs, local_dofs);
  }
  
  int SparseComm::build_fast_maps(const long ndofd,
				  const long ngdof_owned,
				  const long start_gdof_id)
  {
    return (build_fast_LG_map(ndofd, ngdof_owned, start_gdof_id) +
	    build_fast_GL_map(ndofd,ngdof_owned));
  }

  int SparseComm::build_fast_LG_map(const long ndofd,
				    const long ngdof_owned,
				    const long start_gdof_id)
  {
    int err = 0;
    fast_LG_map = PGFEM_calloc(fast_map, 1);
    fast_LG_map->n_maps = 0;

    /* deterimine the local indices that map to global indices on the
       local domain and the corresponding number of maps */
    long *local_idx = PGFEM_calloc(long, ndofd);
    long *global_idx = PGFEM_calloc(long, ndofd);
    for(long i=0; i<ndofd; i++){
      global_idx[i] = LG[i] - start_gdof_id;

      /* store mapping if I own the global dof */
      if(global_idx[i] >= 0 && global_idx[i] < ngdof_owned){
	local_idx[fast_LG_map->n_maps] = i;
	fast_LG_map->n_maps++;
      }
    }
    
    /* sanity check */
    if(fast_LG_map->n_maps > ngdof_owned){
      PGFEM_printerr("[%d] ERROR: too many mappings! %s:%s:%d\n",
		     myrank,__func__,__FILE__,__LINE__);
      PGFEM_Abort();
    }

    fast_LG_map->maps = PGFEM_calloc(size_t, 2*fast_LG_map->n_maps);
    for(int i=0; i<fast_LG_map->n_maps; i++){
      fast_LG_map->maps[2*i] = local_idx[i];
      fast_LG_map->maps[2*i+1] = global_idx[local_idx[i]];
    }

    PGFEM_free(local_idx);
    PGFEM_free(global_idx);

    return err;
  }
  
  int SparseComm::build_fast_GL_map(const long ndofd,
				    const long ngdof_owned)
  {
    int err = 0;
    fast_GL_map = PGFEM_calloc(fast_map, 1);
    fast_GL_map->n_maps = 0;
    
    long *global_idx = PGFEM_calloc(long, ngdof_owned);
    long *local_idx = PGFEM_calloc(long, ngdof_owned);
    
    /* NOTE: may be able to replace. Should be identical to comm->GL? */
    for(long i=0; i< ngdof_owned; i++){
      local_idx[i] = GL[i];
      if(local_idx[i]>= 0 && local_idx[i] < ndofd){
	global_idx[fast_GL_map->n_maps] = i;
	fast_GL_map->n_maps++;
      }
    }
    
    fast_GL_map->maps = PGFEM_calloc(size_t, 2*fast_GL_map->n_maps);
    for(int i=0; i<fast_GL_map->n_maps; i++){
      fast_GL_map->maps[2*i] = global_idx[i];
      fast_GL_map->maps[2*i+1] = local_idx[global_idx[i]];
    }
    
    PGFEM_free(global_idx);
    PGFEM_free(local_idx);
    
    return err;
  }
  
  // Communication destructor
  CommunicationStructure::~CommunicationStructure()
  {
  }
  
  void PGFEM_Comm_code_abort(const CommunicationStructure *com, int code)
  {
    if (com && com->net->initialized()) {
      com->net->abort(com->comm, code);
    }
    
    exit(code);
  }

  // abort/error methods
  void PGFEM_Comm_Abort(const CommunicationStructure *com)
  {
    PGFEM_Comm_code_abort(com, 0);
  }
  
  void PGFEM_Abort()
  {
    PGFEM_Comm_code_abort(NULL, 0);
  }

  int PGFEM_Error_rank(int *rank)
  {
    if (net::gboot) {
      *rank = net::gboot->get_rank();
    }
    return *rank;
  }

} // namespace pgfem3d
