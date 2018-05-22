#include "PGFEM_io.h"
#include "MPINetwork.hpp"
#include <cassert>

using namespace pgfem3d::net;
using namespace pgfem3d::net::isir;

MPINetwork::MPINetwork()
  : ISIRNetwork()
{
  datatype_byte = static_cast<datatype_t>(MPI_BYTE);
  datatype_char = static_cast<datatype_t>(MPI_CHAR);
  datatype_int = static_cast<datatype_t>(MPI_INT);
  datatype_long = static_cast<datatype_t>(MPI_LONG);
  datatype_double = static_cast<datatype_t>(MPI_DOUBLE);

  op_null = static_cast<op_t>(NULL);
  op_min = static_cast<op_t>(MPI_MIN);
  op_max = static_cast<op_t>(MPI_MAX);
  op_sum = static_cast<op_t>(MPI_SUM);
  op_bor = static_cast<op_t>(MPI_BOR);
  
  if (MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*)&NET_COMM_WORLD)) {
    throw PGFEM_printerr("Could not duplicate comm world communicator\n");
  }
}

MPINetwork::~MPINetwork()
{
}

int MPINetwork::type() const {
  return NET_ISIR;
}

int MPINetwork::initialized()
{
  int flag;
  MPI_Initialized(&flag);
  return (flag ? flag : 0);
}

void MPINetwork::finalize()
{
  Check(MPI_Finalize());
}

void MPINetwork::abort(PGFem3D_Comm comm, int code)
{
  Check(MPI_Abort(static_cast<MPI_Comm>(comm), code));
}

void MPINetwork::allocRequestArray(int count, Request *rary[])
{
  if (!count) {
    PGFEM_printerr("[%d]Warning: attempted to allocate 0 Requests, will allocate 1\n");
    count = 1;
  }
  assert(count>0);
  MPI_Request *req = (MPI_Request*)calloc(count, sizeof(MPI_Request));
  Request *r = new Request[count];
  for (int i=0; i<count; i++) {
    r[i].setData(&(req[i]));
    r[i].setAlloc(true);
    r[i].setIDX(i);
  }
  *rary = r;
}

void MPINetwork::allocStatusArray(int count, Status *sary[])
{
  if (!count) {
    PGFEM_printerr("[%d]Warning: attempted to allocate 0 Statuses, will allocate 1\n");
    count = 1;
  }
  assert(count>0);
  MPI_Status *stat = (MPI_Status*)calloc(count, sizeof(MPI_Status));
  Status *s = new Status[count];
  for (int i=0; i<count; i++) {
    s[i].setData(&(stat[i]));
    s[i].setAlloc(true);
    s[i].setIDX(i);
  }
  *sary = s;
}

void MPINetwork::comm_rank(PGFem3D_Comm comm, int *rank)
{
  Check(MPI_Comm_rank(static_cast<MPI_Comm>(comm), rank));
}

void MPINetwork::comm_size(PGFem3D_Comm comm, int *size)
{
  Check(MPI_Comm_size(static_cast<MPI_Comm>(comm), size));
}

void MPINetwork::comm_split(PGFem3D_Comm comm, int color, int key, PGFem3D_Comm *ncomm)
{
  Check(MPI_Comm_split(static_cast<MPI_Comm>(comm), color, key,
		       static_cast<MPI_Comm*>(ncomm)));
}

void MPINetwork::comm_free(PGFem3D_Comm *comm)
{
  Check(MPI_Comm_free(static_cast<MPI_Comm*>(comm)));
}

void MPINetwork::comm_dup(PGFem3D_Comm comm, PGFem3D_Comm *ncomm)
{
  Check(MPI_Comm_dup(static_cast<MPI_Comm>(comm),
		     static_cast<MPI_Comm*>(ncomm)));
}

void MPINetwork::barrier(PGFem3D_Comm comm)
{  
  Check(MPI_Barrier(static_cast<MPI_Comm>(comm)));
}

void MPINetwork::bcast(void *in, int count, datatype_t dt, int root, PGFem3D_Comm comm)
{
  Check(MPI_Bcast(in, count, dt, root, static_cast<MPI_Comm>(comm)));
}

void MPINetwork::reduce(const void *in, void *out, int count, datatype_t dt, op_t op,
		       int root, PGFem3D_Comm comm)
{
  Check(MPI_Reduce(in, out, count, static_cast<MPI_Datatype>(dt), static_cast<MPI_Op>(op),
		   root, static_cast<MPI_Comm>(comm)));
}

void MPINetwork::gather(const void *in, int scount, datatype_t sdt, void *out,
			int rcount, datatype_t rdt, int root, PGFem3D_Comm comm)
{
  Check(MPI_Gather(in, scount, static_cast<MPI_Datatype>(sdt), out, rcount,
		   static_cast<MPI_Datatype>(rdt), root, static_cast<MPI_Comm>(comm)));
}

void MPINetwork::allreduce(const void *in, void *out, int count, datatype_t dt,
			  op_t op, PGFem3D_Comm comm)
{
  Check(MPI_Allreduce(in, out, count, static_cast<MPI_Datatype>(dt),
		      static_cast<MPI_Op>(op), static_cast<MPI_Comm>(comm)));
}

void MPINetwork::allgather(const void *in, int scount, datatype_t sdt,
			     void *out, int rcount, datatype_t rdt, PGFem3D_Comm comm)
{
  Check(MPI_Allgather(in, scount, static_cast<MPI_Datatype>(sdt), out, rcount,
		      static_cast<MPI_Datatype>(rdt), static_cast<MPI_Comm>(comm)));
}

void MPINetwork::allgatherv(const void *in, int scount, datatype_t sdt,
			   void *out, int *rcounts, const int *dipls,
			   datatype_t rdt, PGFem3D_Comm comm)
{
  Check(MPI_Allgatherv(in, scount, static_cast<MPI_Datatype>(sdt), out, rcounts, dipls,
		       static_cast<MPI_Datatype>(rdt), static_cast<MPI_Comm>(comm)));
}

void MPINetwork::pin(const void *addr, const size_t bytes, Key *key) {}
void MPINetwork::unpin(const void *addr, const size_t bytes) {}

void MPINetwork::isend(const void *buf, int count, datatype_t dt, int dest, int tag,
		       PGFem3D_Comm comm, Request *request)
{
  Check(MPI_Isend(buf, count, static_cast<MPI_Datatype>(dt), dest, tag,
		  static_cast<MPI_Comm>(comm), ConvertRequest(request)));
}
void MPINetwork::irecv(void *buf, int count, datatype_t dt,
		      int source, int tag, PGFem3D_Comm comm, Request *request)
{
  Check(MPI_Irecv(buf, count, static_cast<MPI_Datatype>(dt), source,
		  tag, static_cast<MPI_Comm>(comm), ConvertRequest(request)));
}

void MPINetwork::send(const void *buf, int count, datatype_t dt, int dest, int tag,
		      PGFem3D_Comm comm)
{
  Check(MPI_Send(buf, count, static_cast<MPI_Datatype>(dt), dest, tag,
		 static_cast<MPI_Comm>(comm)));
}
    
void MPINetwork::recv(void *buf, int count, datatype_t dt,
		      int source, int tag, PGFem3D_Comm comm, Status *status)
{
  Check(MPI_Recv(buf, count, static_cast<MPI_Datatype>(dt), source,
		 tag, static_cast<MPI_Comm>(comm), ConvertStatus(status)));
}

void MPINetwork::wait(Request *request, Status *status)
{
  Check(MPI_Wait(ConvertRequest(request),
		 ConvertStatus(status)));
}

void MPINetwork::waitall(int count, Request *requests, Status *statuses)
{
  Check(MPI_Waitall(count, ConvertRequest(requests),
		    ConvertStatus(statuses)));
}

void MPINetwork::waitany(int count, Request *requests, int *indx,
			 Status *status)
{
  Check(MPI_Waitany(count, ConvertRequest(requests), indx, ConvertStatus(status)));
  if (!(status == NET_STATUS_IGNORE)) {
    MPI_Status *stat = (MPI_Status*)status->getData();
    status->NET_SOURCE = stat->MPI_SOURCE;
    status->NET_TAG = stat->MPI_TAG;
    status->NET_ERROR = stat->MPI_ERROR;
  }
}

void MPINetwork::waitsome(int incount, Request *requests, int *outcount,
			  int array_of_indices[], Status *statuses)
{
  MPI_Status *stat = 0;
  Check(MPI_Waitsome(incount, ConvertRequest(requests),
		     outcount, array_of_indices, stat));
  if (*outcount && stat)
    ConvertStatus(*outcount, stat, statuses);
}

void MPINetwork::iprobe(int source, int tag, PGFem3D_Comm comm, int *flag,
			Status *status)
{
  Check(MPI_Iprobe(source, tag, static_cast<MPI_Comm>(comm), flag,
		   ConvertStatus(status)));
  if (flag) {
    MPI_Status *stat = (MPI_Status*)status->getData();
    status->NET_SOURCE = stat->MPI_SOURCE;
    status->NET_TAG = stat->MPI_TAG;
    status->NET_ERROR = stat->MPI_ERROR;
  }
}

void MPINetwork::probe(int source, int tag, PGFem3D_Comm comm, Status *status)
{
  Check(MPI_Probe(source, tag, static_cast<MPI_Comm>(comm),
		  ConvertStatus(status)));
  MPI_Status *stat = (MPI_Status*)status->getData();
  status->NET_SOURCE = stat->MPI_SOURCE;
  status->NET_TAG = stat->MPI_TAG;
  status->NET_ERROR = stat->MPI_ERROR;
}

void MPINetwork::test(Request *req, int *flag, Status *status)
{
  Check(MPI_Test(ConvertRequest(req), flag, ConvertStatus(status)));
}

void MPINetwork::get_status_count(const Status *status, datatype_t dt, int *count)
{
  Check(MPI_Get_count(ConvertStatus((Status*)status), static_cast<MPI_Datatype>(dt),
		      count));
}

void MPINetwork::cancel(Request *request)
{

}
