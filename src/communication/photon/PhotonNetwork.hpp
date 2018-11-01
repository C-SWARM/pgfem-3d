// -------------------------------------------------------------------*- C++ -*-
// -----------------------------------------------------------------------------
#ifndef PGFEM3D_COMMUNICATION_PHOTONNETWORK_HPP
#define PGFEM3D_COMMUNICATION_PHOTONNETWORK_HPP

#include "pgfem3d/Network.hpp"
#ifdef HAVE_PGFEM_MPI
#define HAVE_MPI
#include <photon.h>
#undef HAVE_MPI
#else
#include <photon.h>
#endif

#include <stdexcept>

using pgfem3d::net::Network;
using pgfem3d::net::PGFem3D_Comm;
using pgfem3d::net::datatype_t;
using pgfem3d::net::op_t;
using pgfem3d::net::Request;
using pgfem3d::net::Status;

namespace pgfem3d {
namespace net {
namespace pwc {
      
class PhotonNetwork : public PWCNetwork {
public:  
  PhotonNetwork();
  ~PhotonNetwork();

  int type() const;
  
  int initialized();
  void finalize();
  void abort(PGFem3D_Comm comm, int code);
  
  void allocRequestArray(int count, Request *rary[]);
  void allocStatusArray(int count, Status *sary[]);

  void comm_rank(PGFem3D_Comm comm, int *rank);
  void comm_size(PGFem3D_Comm comm, int *size);
  void comm_split(PGFem3D_Comm comm, int color, int key, PGFem3D_Comm *ncomm);
  void comm_free(PGFem3D_Comm *comm);
  void comm_dup(PGFem3D_Comm comm, PGFem3D_Comm *ncomm);
  
  void barrier(PGFem3D_Comm comm);
  void bcast(void *in, int count, datatype_t dt, int root, PGFem3D_Comm comm);
  void reduce(const void *in, void *out, int count, datatype_t dt, op_t op,
	     int root, PGFem3D_Comm comm);
  void gather(const void *in, int scount, datatype_t sdt, void *out,
	     int rcount, datatype_t rdt, int root, PGFem3D_Comm comm);
  void allreduce(const void *in, void *out, int count, datatype_t dt,
		op_t op, PGFem3D_Comm comm);
  void allgather(const void *in, int scount, datatype_t sdt,
		void *out, int rcount, datatype_t rdt, PGFem3D_Comm comm);
  void allgatherv(const void *in, int scount, datatype_t sdt,
		 void *out, int *rcounts, const int *dipls,
		 datatype_t rdt, PGFem3D_Comm comm);

  void pin(const void *addr, const size_t bytes, Key *key);
  void unpin(const void *addr, const size_t bytes);
  void pwc(int dst, size_t n, Buffer *lbuf, Buffer *rbuf,
	   const CID& lid, const CID& rid);
  void gwc(int src, size_t n, Buffer *lbuf, Buffer *rbuf,
	   const CID& lid, const CID& rid);
  void probe(int *flag, CID *comp, Status *status, void (*cb)(CID));
  void progress(int *flag, CID *comp, Status *status, void (*cb)(CID));
  void wait_n(int count);
  void wait_n_id(int count, const CID& id);

  Buffer* getbuffer();
  
private:
  static void Check(int rc) {
    if (rc != PHOTON_OK) {
      throw std::exception();
    }
  }

  struct photon_config_t cfg;

  const char *fi_socket = "sockets";
  const char *fi_psm2 = "psm2";
  const char *fi_dev = NULL; 
 
  // Workspace for tracking remote buffer metadata
  Buffer *wbuf;
};

}
}
}

#endif
