#ifndef PGFEM3D_NETWORK_H
#define PGFEM3D_NETWORK_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/// @brief This file defines the PGFem3D Network interface
#include "PGFem3D_options.h"
#include <cstdio>
#include <cstdint>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace pgfem3d {
namespace net {

enum net_types_t {
  NET_ISIR = 0,
  NET_PWC
};

enum net_codes_t {
  NET_SUCCESS = 0,
  NET_ERROR,
  NET_RESOURCE
};

#ifdef HAVE_MPI
// MPI compatibility
typedef MPI_Comm PGFem3D_Comm;
typedef MPI_Datatype datatype_t;
typedef MPI_Op op_t;
#else
typedef void* PGFem3D_Comm;
typedef void* datatype_t;
typedef void* op_t;
#endif
  
#ifdef HAVE_MPI
// MPI compatibility
#define NET_ANY_TAG       MPI_ANY_TAG
#define NET_WAIT_ANY      MPI_WAIT_ANY
#define NET_BOTTOM        MPI_BOTTOM
#define NET_IN_PLACE      MPI_IN_PLACE
#else
#define NET_ANY_TAG       (-1)
#define NET_WAIT_ANY      ((void*) 0)
#define NET_BOTTOM        ((void*) 0)
#define NET_IN_PLACE      ((void*) 1)
#endif

#define NET_STATUS_IGNORE ((pgfem3d::net::Status*) 0)
  
extern datatype_t datatype_byte;
extern datatype_t datatype_char;
extern datatype_t datatype_int;
extern datatype_t datatype_long;
extern datatype_t datatype_double;

#define NET_DT_BYTE ((datatype_t) (datatype_byte))
#define NET_DT_CHAR ((datatype_t) (datatype_char))
#define NET_DT_INT ((datatype_t) (datatype_int))
#define NET_DT_LONG ((datatype_t) (datatype_long))
#define NET_DT_DOUBLE ((datatype_t) (datatype_double))

extern op_t op_null;
extern op_t op_min;
extern op_t op_max;
extern op_t op_sum;
extern op_t op_bor;

#define NET_OP_NULL ((op_t) (op_null))
#define NET_OP_MAX ((op_t) (op_max))
#define NET_OP_MIN ((op_t) (op_min))
#define NET_OP_SUM ((op_t) (op_sum))
#define NET_OP_BOR ((op_t) (op_bor))

extern PGFem3D_Comm NET_COMM_WORLD;
  
// Make CIDs 8 byte values 
typedef uint64_t CID;

typedef struct key_t {
  uint64_t key0;
  uint64_t key1;
} Key;  // XXX (refactor)
  
typedef struct net_buffer_t {
  uintptr_t addr;
  size_t size;
  Key key;
} Buffer;

class NetworkObject {
public:
  NetworkObject() : data(0), _allocated(false), _idx(-1) {};
  ~NetworkObject();
  bool getAlloc() const { return _allocated; }
  void setAlloc(bool b) { _allocated = b; }
  void* getData() const { return data; }
  void setData(void *d) { data = d; }
  int getIDX() const { return _idx; }
  void setIDX(int i) { _idx = i; }
protected:
  void *data;
  bool _allocated;
  int _idx;
};
  
class Request : public NetworkObject {};
class Status : public NetworkObject
{
public:
  int NET_SOURCE;
  int NET_TAG;
  int NET_ERROR;
  size_t size;
  int count;
  int remainder;
  void *request;
};

class Network {
public:
  static Network* Create(const PGFem3D_opt& opts);
  
  virtual ~Network();

  // get the type of Network configured by PGFem3D
  virtual int type() const = 0;
  
  virtual int initialized() = 0;
  virtual void finalize() = 0;
  virtual void abort(PGFem3D_Comm comm, int code) = 0;

  virtual int get_rank(PGFem3D_Comm comm) = 0;
  virtual int get_nproc(PGFem3D_Comm comm) = 0;
  
  virtual void allocRequestArray(int count, Request *rary[]) = 0;
  virtual void allocStatusArray(int count, Status *sary[]) = 0;
  
  virtual void pin(const void *addr, const size_t bytes, Key *key) = 0;
  virtual void unpin(const void *addr, const size_t bytes) = 0;

  virtual void barrier(PGFem3D_Comm comm) = 0;
  virtual void reduce(const void *in, void *out, int count, datatype_t dt, op_t op,
		      int root, PGFem3D_Comm comm) = 0;
  virtual void gather(const void *in, int scount, datatype_t sdt, void *out,
		      int rcount, datatype_t rdt, int root, PGFem3D_Comm comm) = 0;
  virtual void allreduce(const void *in, void *out, int count, datatype_t dt,
			 op_t op, PGFem3D_Comm comm) = 0;
  virtual void allgather(const void *in, int scount, datatype_t sdt,
			 void *out, int rcount, datatype_t rdt, PGFem3D_Comm comm) = 0;
  virtual void allgatherv(const void *in, int scount, datatype_t sdt,
			  void *out, int *rcounts, const int *dipls,
			  datatype_t rdt, PGFem3D_Comm comm) = 0;
};

class ISIRNetwork : public Network {
public:
  virtual void isend(const void *buf, int count, datatype_t dt, int dest, int tag,
		     PGFem3D_Comm comm, Request *request) = 0;
  virtual void irecv(void *buf, int count, datatype_t dt,
		     int source, int tag, PGFem3D_Comm comm, Request *request) = 0;
  virtual void wait(Request *req, Status *status) = 0;
  virtual void waitall(int count, Request *requests,
		       Status *statuses) = 0;
  virtual void waitany(int count, Request *requests, int *indx,
		       Status *status) = 0;
  virtual void waitsome(int incount, Request *requests,
			int *outcount, int array_of_indices[],
			Status *statuses) = 0;
  virtual void iprobe(int source, int tag, PGFem3D_Comm comm, int *flag,
		      Status *status) = 0;
  virtual void get_status_count(const Status *status, datatype_t dt, int *count) = 0;
};

class PWCNetwork : public Network {
public:
  virtual void pwc(int dst, size_t n, Buffer *lbuf, Buffer *rbuf,
		   const CID& lid, const CID& rid) = 0;
  virtual void gwc(int src, size_t n, Buffer *lbuf, Buffer *rbuf,
		   const CID& lid, const CID& rid) = 0;
  virtual void probe(int *flag, CID *comp, Status *status, void (*cb)(CID)) = 0;
  virtual void progress(int *flag, CID *comp, Status *status, void (cb)(CID)) = 0;
  virtual void wait_n(int count) = 0;
  virtual void wait_n_id(int count, const CID& id) = 0;
};
} // namespace net
} // namespace pgfem3d

#endif
