#include "pgfem3d/Network.hpp"

#include <stdexcept>
#include <mpi.h>

using pgfem3d::net::Network;
using pgfem3d::net::PGFem3D_Comm;
using pgfem3d::net::datatype_t;
using pgfem3d::net::op_t;
using pgfem3d::net::Request;
using pgfem3d::net::Status;

namespace pgfem3d {
namespace net {
namespace isir {
      
class MPINetwork : public ISIRNetwork {
public:
  MPINetwork();
  ~MPINetwork();

  int type() const;
  
  int initialized();
  void finalize();
  void abort(PGFem3D_Comm comm, int code);

  int get_rank(PGFem3D_Comm comm);
  int get_nproc(PGFem3D_Comm comm);
  
  void allocRequestArray(int count, Request *rary[]);
  void allocStatusArray(int count, Status *sary[]);
  
  void barrier(PGFem3D_Comm comm);
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
  
  void isend(const void *buf, int count, datatype_t dt, int dest, int tag,
	     PGFem3D_Comm comm, Request *request);
  void irecv(void *buf, int count, datatype_t dt,
	     int source, int tag, PGFem3D_Comm comm, Request *request);
  void wait(Request *req, Status *status);
  void waitall(int count, Request *requests, Status *statuses);
  void waitany(int count, Request *requests, int *indx, Status *status);
  void waitsome(int incount, Request *requests, int *outcount,
		int array_of_indices[], Status *statuses);
  void iprobe(int source, int tag, PGFem3D_Comm comm, int *flag,
	      Status *status);
  void get_status_count(const Status *status, datatype_t dt, int *count);

private:
  static void Check(int rc) {
    if (rc != MPI_SUCCESS) {
      throw std::exception();
    }
  }
  
  static MPI_Request* ConvertRequest(Request *r) {
    if (!r) return NULL;
    MPI_Request *req;
    if (!(r->getAlloc())) {
      req = (MPI_Request*)calloc(1, sizeof(MPI_Request));
      r->setData(req);
      r->setAlloc(true);
    }
    else {
      req = (MPI_Request*)r->getData();
    }
    return req;
  }

  static MPI_Status* ConvertStatus(Status *s) {
    if (s == NET_STATUS_IGNORE)
      return MPI_STATUS_IGNORE;
    if (!s)
      return NULL;
    MPI_Status *stat;
    if (!(s->getAlloc())) {
      stat = (MPI_Status*)calloc(1, sizeof(MPI_Status));
      s->setData(stat);
      s->setAlloc(true);
    }
    else {
      stat = (MPI_Status*)s->getData();
    }
    return stat;
  }

  static void ConvertStatus(int count, MPI_Status *ms, Status *s) {
    s->setData(ms);
  }
};

}
}
}
