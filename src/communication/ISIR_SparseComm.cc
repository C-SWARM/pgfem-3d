#include "ISIR_SparseComm.hpp"
#include "allocation.h"
#include "incl.h"
#include "utils.h"

#ifndef UTILS_DEBUG
#define UTILS_DEBUG 0
#endif

using namespace pgfem3d;
using namespace pgfem3d::net;

namespace pgfem3d {
  
  ISIR_SparseComm::ISIR_SparseComm(ISIRNetwork *n, PGFem3D_Comm c)
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
    
    net = n;
    comm = c;

    n->comm_rank(c, &myrank);
    n->comm_size(c, &nproc);
    
    sta_r = NULL;
    sta_s = NULL;
    req_r = NULL;
    req_s = NULL;
  }
  
  ISIR_SparseComm::~ISIR_SparseComm()
  {
    if (Nr > 0) {
      delete [] req_r;
      delete [] sta_r;
    }
    
    if (Ns > 0) {
      delete [] req_s;
      delete [] sta_s;
    }
    
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

    dealoc1l(send);
    dealoc1l(recv);
  }
  
  void ISIR_SparseComm::initialize()
  {
    if (Nr > 0) {
      net->allocStatusArray(Nr, &sta_r);
      net->allocRequestArray(Nr, &req_r);
    }

    if (Ns > 0) {
      net->allocStatusArray(Ns, &sta_s);
      net->allocRequestArray(Ns, &req_s);
    }
    
    send = PGFEM_calloc(double*, nproc);
    recv = PGFEM_calloc(double*, nproc);
  }
  
    /** builds the off-process buffer for assembling the gloabal stiffness
      matrix. Posts non-blocking recieve */
  void ISIR_SparseComm::post_stiffmat(double ***Lk, double ***receive)
  {
    /* Allocate buffer for off-process assembly and receive */
    for (int i = 0; i < Nr; ++i) {
      /* the current proc to setup receives for */
      long p = Nrr[i];
      size_t n_alloc_recv = 0;
      
      /* number to allocate for receive */
      recv[p] = NULL;
      if (myrank == p || AR[p] == 0) {
	n_alloc_recv = 1;
      } else {
	n_alloc_recv = AR[p];
      }
      
      /* allocate */
      recv[p] = PGFEM_calloc(double, n_alloc_recv);
    }
    
    for (int i = 0; i < Ns; ++i) {
      long p = Nss[i];
      size_t n_alloc_Lk = 0;
      
      /* number to allocate for Lk */
      send[p] = NULL;
      if (myrank == p || S[p] == 0) {
	n_alloc_Lk = 1;
      } else {
	n_alloc_Lk = AS[p];
      }
      
      send[p] = PGFEM_calloc(double, n_alloc_Lk);
    }
    
    /* point to pre-allocated buffers */
    *Lk = send;
    *receive = recv;
    
    /* post the non-blocking receives */
    for (int i = 0; i < Nr; ++i){
      long idx = Nrr[i];
      net->irecv(recv[idx], AR[idx], NET_DT_DOUBLE, idx,
		 NET_ANY_TAG, comm, &req_r[i]);
    }
  }

  void ISIR_SparseComm::send_stiffmat()
  {
    /* post sends */
    for(int i=0; i<Ns; i++){
      size_t idx = Nss[i];
      net->isend(send[idx], AS[idx], NET_DT_DOUBLE, idx,
		 myrank, comm, &req_s[i]);
    }
  }
  
  void ISIR_SparseComm::finalize_stiffmat()
  {
    net->waitall(Nr, req_r, sta_r);
    net->waitall(Ns, req_s, sta_s);
    
    /* free send/recv buffers */
    for (int i = 0; i < Ns; ++i) {
      long p = Nss[i];
      free(send[p]);
    }
    
    for (int i = 0; i < Nr; ++i) {
      long p = Nrr[i];
      free(recv[p]);
    }
  }
  
  void ISIR_SparseComm::assemble_nonlocal_stiffmat(solvers::SparseSystem* system)
  {
    int comm_idx = 0;
    int n_received = 0;
    while (n_received < Nr) {
      /* get the communication index */
      net->waitany(Nr, req_r, &comm_idx, sta_r);
      
      /* convert communication index to proc_id */
      const int proc = Nrr[comm_idx];
      
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

  void ISIR_SparseComm::LToG(const double *f, double *Gf,
			     const long ndofd, const long *DomDof,
			     const int GDof)
  {
    long i,j,KK;
    double **snd,**rcv;
    Status *sta_s,*sta_r;
    Request *req_s,*req_r;
    
    /* for (i=0;i<DomDof[myrank];i++) */
    /*   Gf[i] = 0.0; */
    nulld(Gf,DomDof[myrank]);
    
    /* Allocate receive */
    snd = PGFEM_calloc (double*, nproc);
    rcv = PGFEM_calloc (double*, nproc);
    for (i=0;i<nproc;i++) {
      if (S[i] == 0) KK = 1; else KK = S[i];
      snd[i] = PGFEM_calloc (double, KK);
      if (R[i] == 0) KK = 1; else KK = R[i];
      rcv[i] = PGFEM_calloc (double, KK);
    }
    
    /* Allocate status and request fields */
    if (Ns == 0) KK = 1; else KK = Ns;
    net->allocStatusArray(KK, &sta_s);
    net->allocRequestArray(KK, &req_s);
    
    if (Nr == 0) KK = 1; else KK = Nr;
    net->allocStatusArray(KK, &sta_r);
    net->allocRequestArray(KK, &req_r);
    
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
    for (i=0;i<Nr;i++){
      KK = Nrr[i];
      net->irecv(rcv[KK], R[KK], NET_DT_DOUBLE, KK,
		 NET_ANY_TAG, comm, &req_r[i]);
    }
    
    /*************/
    /* Send data */
    /*************/
    for (i=0;i<Ns;i++){
      KK = Nss[i];
      for (j=0;j<S[KK];j++)
	snd[KK][j] = f[SLID[KK][j]];
      
      net->isend(snd[KK], S[KK], NET_DT_DOUBLE, KK,
		 myrank, comm, &req_s[i]);
    }
    
    get_owned_global_dof_values(f, Gf);
    
    /* Wait to complete the comunicatins */
    net->waitall(Nr, req_r, sta_r);
    
    for (i=0;i<Nr;i++){
      KK = Nrr[i];
      for (j=0;j<R[KK];j++)
	Gf[RGID[KK][j] - GDof] += rcv[KK][j];
    }
    
    net->waitall(Ns, req_s, sta_s);
    /* Deallocate snd and rcv */
    for (i=0;i<nproc;i++) {
      free(snd[i]);
      free(rcv[i]);
    }
    free(snd);
    free(rcv);
    
    delete[] sta_s;
    delete[] sta_r;
    delete[] req_s;
    delete[] req_r;
  }

  void ISIR_SparseComm::GToL(const double *Gr, double *r,
			     const long ndofd, const int GDof)
  {
    long i,j,KK;
    double **snd,**rcv;
    Status *sta_s,*sta_r;
    Request *req_s,*req_r;
    
    //  for (i=0;i<ndofd;i++) r[i] = 0.0;
    nulld(r,ndofd);
    
    //  for (i=0;i<ndofd;i++) r[i] = 0.0;
    nulld(r,ndofd);
    
    /* Allocate rcv */
    snd = PGFEM_calloc (double*, nproc);
    rcv = PGFEM_calloc (double*, nproc);
    for (i=0;i<nproc;i++) {
      if (R[i] == 0) KK = 1; else KK = R[i];
      snd[i] = PGFEM_calloc (double, KK);
      if (S[i] == 0) KK = 1; else KK = S[i];
      rcv[i] = PGFEM_calloc (double, KK);
    }
    
    /* Allocate status and request fields */
    if (Nr == 0) KK = 1; else KK = Nr;
    net->allocStatusArray(KK, &sta_s);
    net->allocRequestArray(KK, &req_s);
    
    if (Ns == 0) KK = 1; else KK = Ns;
    net->allocStatusArray(KK, &sta_r);
    net->allocRequestArray(KK, &req_r);
    
    /****************/
    /* Receive data */
    /****************/
    for (i=0;i<Ns;i++){
      KK = Nss[i];
      net->irecv(rcv[KK], S[KK], NET_DT_DOUBLE, KK,
		 NET_ANY_TAG, comm, &req_r[i]);
    }
    
    /*************/
    /* Send data */
    /*************/
    for (i=0;i<Nr;i++){
      KK = Nrr[i];
      for (j=0;j<R[KK];j++)
	snd[KK][j] = Gr[RGID[KK][j] - GDof];
      
      net->isend(snd[KK], R[KK], NET_DT_DOUBLE, KK,
		 myrank, comm, &req_s[i]);
    }
    
    get_local_dof_values_from_global(Gr,r);
    
    /* Wait to complete the comunicatins */
    net->waitall(Ns, req_r, sta_r);
    
    for (i=0;i<Ns;i++){
      KK = Nss[i];
      for (j=0;j<S[KK];j++)
	r[SLID[KK][j]] = rcv[KK][j];
    }
    
    net->waitall(Nr, req_s, sta_s);
    /* Deallocate snd and rcv */
    for (i=0;i<nproc;i++) {
      free(snd[i]);
      free(rcv[i]);
    }
    
    free(snd);
    free(rcv);
    
    delete[] sta_s;
    delete[] sta_r;
    delete[] req_s;
    delete[] req_r;
  }
} //namespace pgfem3d
