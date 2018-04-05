#include "PWC_SparseComm.hpp"
#include "allocation.h"
#include "incl.h"
#include "utils.h"

#include <cassert>

#ifndef UTILS_DEBUG
#define UTILS_DEBUG 0
#endif

using namespace pgfem3d;
using namespace pgfem3d::net;

namespace pgfem3d {
  
  PWC_SparseComm::PWC_SparseComm(PWCNetwork *n, PGFem3D_Comm c)
  {
    S = NULL;
    R = NULL;
    AS = NULL;
    AR = NULL;
    SLID = NULL;
    RGID = NULL;
    LG = NULL;
    GL = NULL;
    SAp = NULL;
    SGRId = NULL;
    RAp = NULL;
    RGRId = NULL;
    Ns = 0;
    Nr = 0;
    Nss = NULL;
    Nrr = NULL;

    total_ssz = 0;
    total_rsz = 0;

    total_lg_ssz = 0;
    total_lg_rsz = 0;
    
    net = n;
    comm = c;
    
    n->comm_rank(c, &myrank);
    n->comm_size(c, &nproc);
  }
  
  PWC_SparseComm::~PWC_SparseComm()
  {
    for (long i=0; i<nproc; i++) {
      free(SLID[i]);
      free(RGID[i]);
      free(SAp[i]);
      free(SGRId[i]);
      free(RAp[i]);
      free(RGRId[i]);
    }
    
    free(S);
    free(R);
    free(AS);
    free(AR);
    free(SLID);
    free(RGID);
    free(LG);
    free(GL);
    free(SAp);
    free(SGRId);
    free(RAp);
    free(RGRId);
    free(Nss);
    free(Nrr);
    
    /* destroy maps */
    if (fast_LG_map != NULL)
      free(fast_LG_map->maps);
    if (fast_GL_map != NULL)
      free(fast_GL_map->maps);
    
    free(fast_LG_map);
    free(fast_GL_map);
  }

  void PWC_SparseComm::exchange(Buffer *in, Buffer *out, long nsend,
				long nrecv, long *idxs)
  {
    CID lid = LOCAL_ID;
    
    // get remote backing buffer keys
    net->allgather(&in[myrank], sizeof(Buffer), NET_DT_BYTE,
    		   out, sizeof(Buffer), NET_DT_BYTE, comm);

    for (int i = 0; i < nsend; ++i) {
      long p = idxs[i];
      // send the receive address to our peer as a CID
      CID rid = (CID)in[p].addr;
      net->pwc(p, 0, 0, 0, lid, rid);
    }
    net->wait_n_id(nsend, lid);
    
    // wait for and update sending addresses (remote recv offsets)
    int n_received = 0;
    while (n_received < nrecv) {
      int flag;
      CID addr;
      Status stat;
      net->probe(&flag, &addr, &stat, 0);
      if (flag) {
	int p = stat.NET_SOURCE;
	out[p].addr = (uintptr_t)addr;
	n_received++;
      }
    }
  }
  
  void PWC_SparseComm::initialize()
  {
    Key skey, rkey;
    Key lg_skey, lg_rkey;
    
    send = PGFEM_calloc(double*, nproc);
    recv = PGFEM_calloc(double*, nproc);

    local_s = PGFEM_calloc(Buffer, nproc);
    local_r = PGFEM_calloc(Buffer, nproc);
    remote = PGFEM_calloc(Buffer, nproc);

    lg_S = PGFEM_calloc(double*, nproc);
    lg_R = PGFEM_calloc(double*, nproc);
    
    local_lg_S = PGFEM_calloc(Buffer, nproc);
    local_lg_R = PGFEM_calloc(Buffer, nproc);
    remote_lg_S = PGFEM_calloc(Buffer, nproc);
    remote_lg_R = PGFEM_calloc(Buffer, nproc);
    
    // how much space for send
    for (int i = 0; i < Ns; ++i) {
      long p = Nss[i];
      local_s[p].size = AS[p]*sizeof(double);
      local_lg_S[p].size = S[p]*sizeof(double);
      total_ssz += AS[p];
      total_lg_ssz += S[p];
    }

    // allocate send backing buffers
    backing_s = PGFEM_calloc_pin(double, total_ssz, net, &skey);
    backing_lg_S = PGFEM_calloc_pin(double, total_lg_ssz, net, &lg_skey);
    local_s[myrank].key = skey;
    local_lg_S[myrank].key = lg_skey;
    
    // map send arrays to backing buffers
    char *sptr = (char*)backing_s;
    char *lg_sptr = (char*)backing_lg_S;
    for (int i = 0; i < Ns; ++i) {
      long p = Nss[i];
      send[p] = (double*)sptr;
      local_s[p].addr = (uintptr_t)sptr;
      local_s[p].key = skey;

      lg_S[p] = (double*)lg_sptr;
      local_lg_S[p].addr = (uintptr_t)lg_sptr;
      local_lg_S[p].key = lg_skey;
      
      // advance backing buffer ptrs
      sptr += local_s[p].size;
      lg_sptr += local_lg_S[p].size;
    }
    
    // how much space for recv
    for (int i = 0; i < Nr; ++i) {
      long p = Nrr[i];
      local_r[p].size = AR[p]*sizeof(double);
      local_lg_R[p].size = R[p]*sizeof(double);
      total_rsz += AR[p];
      total_lg_rsz += R[p];
    }

    // allocate recv backing buffers
    backing_r = PGFEM_calloc_pin(double, total_rsz, net, &rkey);
    backing_lg_R = PGFEM_calloc_pin(double, total_lg_rsz, net, &lg_rkey);
    local_r[myrank].key = rkey;
    local_lg_R[myrank].key = lg_rkey;
    
    // map recv arrays to backing buffers
    char *rptr = (char*)backing_r;
    char *lg_rptr = (char*)backing_lg_R;
    for (int i = 0; i < Nr; ++i) {
      long p = Nrr[i];
      recv[p] = (double*)rptr;
      local_r[p].addr = (uintptr_t)rptr;
      local_r[p].key = rkey;

      lg_R[p] = (double*)lg_rptr;
      local_lg_R[p].addr = (uintptr_t)lg_rptr;
      local_lg_R[p].key = lg_rkey;
      
      // advance backing buffer ptr
      rptr += local_r[p].size;
      lg_rptr += local_lg_R[p].size;
    }

    // stiffmat metadata
    exchange(local_r, remote, Nr, Ns, Nrr);
    // LToG and GToL metadata
    exchange(local_lg_S, remote_lg_S, Ns, Nr, Nss);
    exchange(local_lg_R, remote_lg_R, Nr, Ns, Nrr);
  }
  
  void PWC_SparseComm::post_stiffmat(double ***Lk, double ***receive)
  {
    // calling code expects zeroed buffers
    memset((void*)backing_s, 0, total_ssz*sizeof(double));
    memset((void*)backing_r, 0, total_rsz*sizeof(double));

    // return buffer handles
    *Lk = send;
    *receive = recv;

    // No recv post with one-sided Photon/PWC
  }

  void PWC_SparseComm::send_stiffmat()
  {
    CID rid = myrank;
    CID lid = LOCAL_ID;
    long scount = 0;
    for (int i = 0; i < Ns; ++i) {
      long p = Nss[i];
      net->pwc(p, local_s[p].size, &local_s[p], &remote[p], lid, rid);
      scount++;
    }

    // wait for local completions
    net->wait_n_id(scount, lid);
  }

  void PWC_SparseComm::finalize_stiffmat()
  {
    net->barrier(comm);
  }
  
  void PWC_SparseComm::assemble_nonlocal_stiffmat(solvers::SparseSystem* system)
  {
    CID proc;
    int n_received = 0;
    while (n_received < Nr) {
      int flag;
      Status stat;
      net->probe(&flag, &proc, &stat, 0);
      if (flag) {
	/* get number of rows */
	const int nrows = R[proc];
	
	/* allocate rows and cols to receive */
	int *row_idx = PGFEM_calloc(int, nrows);
	int *ncols = PGFEM_calloc(int, nrows);
	int *col_idx = PGFEM_calloc(int, AR[proc]);
	
	/* get row and column ids */
	int idx = 0;
	for(int j=0; j<R[proc]; j++){
	  row_idx[j] = RGID[proc][j];
	  ncols[j] = RAp[proc][j];
	  for(int k=0; k<ncols[j]; k++){
	    col_idx[idx] = RGRId[proc][idx];
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
    }
  }

  void PWC_SparseComm::LToG(const double *f, double *Gf,
			    const long ndofd, const long *DomDof,
			    const int GDof)
  {
    long i,j,KK;
    CID lid = LOCAL_ID;
    CID rid = myrank;
    
    /* for (i=0;i<DomDof[myrank];i++) */
    /*   Gf[i] = 0.0; */
    nulld(Gf,DomDof[myrank]);
    
    if(UTILS_DEBUG){/* debug the snd-rcv */
      for(int i=0; i<nproc; i++){
	if(myrank == i){
	  for(int j=0; j<nproc; j++){
	    PGFEM_printf("[%d]: snding %ld to %d\n",myrank,S[j],j);
	    PGFEM_printf("[%d]: recieving %ld from %d\n",myrank,R[j],j);
	  }
	}
	net->barrier(comm);
      }
      net->barrier(comm);
    }
    
    /****************/
    /* Receive data */
    /****************/
    // PWC one-sided recv buffers already allocated and mapped at init
    
    /*************/
    /* Send data */
    /*************/
    for (i=0;i<Ns;i++){
      KK = Nss[i];
      for (j=0;j<S[KK];j++)
	lg_S[KK][j] = f[SLID[KK][j]];

      net->pwc(KK, local_lg_S[KK].size, &local_lg_S[KK],
	       &remote_lg_R[KK], lid, rid);
    }
    
    get_owned_global_dof_values(f, Gf);
    
    /* Wait to complete the comunicatins */
    net->wait_n_id(Ns, lid);

    CID proc;
    int n_received = 0;
    while (n_received < Nr) {
      int flag;
      Status stat;
      net->probe(&flag, &proc, &stat, 0);
      if (flag) {
	KK = proc;
	for (j=0;j<R[KK];j++)
	  Gf[RGID[KK][j] - GDof] += lg_R[KK][j];
	n_received++;
      }
    }
    net->barrier(comm);
  }
  
  void PWC_SparseComm::GToL(const double *Gr, double *r,
			    const long ndofd, const int GDof)
  {
    long i,j,KK;
    CID lid = LOCAL_ID;
    CID rid = myrank;
    
    //  for (i=0;i<ndofd;i++) r[i] = 0.0;
    nulld(r,ndofd);
    
    //  for (i=0;i<ndofd;i++) r[i] = 0.0;
    nulld(r,ndofd);
    
    /****************/
    /* Receive data */
    /****************/
    // PWC one-sided recv buffers already allocated and mapped at init
    
    /*************/
    /* Send data */
    /*************/
    for (i=0;i<Nr;i++){
      KK = Nrr[i];
      for (j=0;j<R[KK];j++)
	lg_R[KK][j] = Gr[RGID[KK][j] - GDof];

      net->pwc(KK, local_lg_R[KK].size, &local_lg_R[KK],
	       &remote_lg_S[KK], lid, rid);
    }
    
    get_local_dof_values_from_global(Gr,r);
    
    /* Wait to complete the comunicatins */
    net->wait_n_id(Nr, lid);

    CID proc;
    int n_received = 0;
    while (n_received < Ns) {
      int flag;
      Status stat;
      net->probe(&flag, &proc, &stat, 0);
      if (flag) {
	KK = proc;
	for (j=0;j<S[KK];j++)
	  r[SLID[KK][j]] = lg_S[KK][j];
	n_received++;
      }
    }
    net->barrier(comm);
  }
} //namespace pgfem3d
