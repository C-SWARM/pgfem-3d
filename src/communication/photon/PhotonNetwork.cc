#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "allocation.h"
#include "PGFEM_io.h"
#include "PhotonNetwork.hpp"
#include "pgfem3d/Communication.hpp"
#include <cassert>
#include <cstring>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace pgfem3d::net;
using namespace pgfem3d::net::pwc;

PhotonNetwork::PhotonNetwork()
  : PWCNetwork()
{
  memset(&cfg, 0, sizeof(cfg));
  cfg.ibv.eth_dev        = "mlnx0";
  cfg.ibv.ib_dev         = "mlx4_0+mlx5_0+qib0+hfi1_0";
  cfg.ibv.ud_gid_prefix  = "ff0e::ffff:0000:0000";
  cfg.ugni.bte_thresh    = -1;
  cfg.fi.provider        = "sockets";
  cfg.cap.max_cid_size   = -1;
  cfg.cap.small_msg_size = -1;
  cfg.cap.small_pwc_size = -1;
  cfg.cap.eager_buf_size = -1;
  cfg.cap.pwc_buf_size   = -1;
  cfg.cap.ledger_entries =  8; // This can be small as PGFem3D primarily does
                               // simple point-to-point scatter/gather RDMA right now
  cfg.cap.max_rd         = -1;
  cfg.cap.default_rd     = -1;
  cfg.cap.num_cq         = -1;
  cfg.cap.use_rcq        =  1;
  cfg.attr.comp_order    = PHOTON_ORDER_DEFAULT;
  cfg.meta_exch          = PHOTON_EXCH_MPI;
  cfg.backend            = PHOTON_BACKEND_DEFAULT;
  cfg.coll               = PHOTON_COLL_IFACE_PWC;

  int rc = photon_init(&cfg);
  if (rc != PHOTON_OK) {
    int myrank;
    PGFEM_Error_rank(&myrank);
    PGFEM_printerr("[%d]ERROR: Photon network failed to initialize\n", myrank);
    PGFEM_Abort();
  }

  datatype_byte = static_cast<datatype_t>(photon_byte);
  datatype_char = static_cast<datatype_t>(photon_char);
  datatype_int = static_cast<datatype_t>(photon_int);
  datatype_long = static_cast<datatype_t>(photon_long);
  datatype_double = static_cast<datatype_t>(photon_double);
  
  op_null = static_cast<op_t>(NULL);
  op_min = static_cast<op_t>(photon_op_min);
  op_max = static_cast<op_t>(photon_op_max);
  op_sum = static_cast<op_t>(photon_op_sum);
  op_bor = static_cast<op_t>(photon_op_bor);

  NET_COMM_WORLD = static_cast<PGFem3D_Comm>(PHOTON_COMM_WORLD);
}

PhotonNetwork::~PhotonNetwork()
{
}

int PhotonNetwork::type() const {
  return NET_PWC;
}

int PhotonNetwork::initialized()
{
  return photon_initialized();
}

void PhotonNetwork::finalize()
{
  photon_finalize();
}

void PhotonNetwork::abort(PGFem3D_Comm comm, int code)
{
}

void PhotonNetwork::allocRequestArray(int count, Request *rary[])
{
}

void PhotonNetwork::allocStatusArray(int count, Status *sary[])
{
}

void PhotonNetwork::comm_rank(PGFem3D_Comm comm, int *rank)
{
  *rank = cfg.address;
}

void PhotonNetwork::comm_size(PGFem3D_Comm comm, int *size)
{
  *size = cfg.nproc;
}

void PhotonNetwork::comm_split(PGFem3D_Comm comm, int color, int key, PGFem3D_Comm *ncomm)
{
}

void PhotonNetwork::comm_free(PGFem3D_Comm *comm)
{
}

void PhotonNetwork::comm_dup(PGFem3D_Comm comm, PGFem3D_Comm *ncomm)
{
}

void PhotonNetwork::barrier(PGFem3D_Comm comm)
{  
  photon_rid req;
  Check(photon_collective_init(static_cast<photonComm>(comm), PHOTON_COLL_BARRIER,
			       (photon_cid){0}, &req, PHOTON_REQ_PWC_NO_LCE));
  Check(photon_collective_join(req, NULL, NULL, 0, 0, photon_datatype_null,
			       photon_datatype_null, 0, PHOTON_OP_NULL));
}

void PhotonNetwork::bcast(void *in, int count, datatype_t dt, int root, PGFem3D_Comm comm)
{
  photon_rid req;
  Check(photon_collective_init(static_cast<photonComm>(comm), PHOTON_COLL_BCAST,
			       (photon_cid){0}, &req, PHOTON_REQ_PWC_NO_LCE));
  Check(photon_collective_join(req, (void*)in, NULL, count, count, dt, dt, root,
			       PHOTON_OP_NULL));
}

void PhotonNetwork::reduce(const void *in, void *out, int count, datatype_t dt, op_t op,
		       int root, PGFem3D_Comm comm)
{
  photon_rid req;
  Check(photon_collective_init(static_cast<photonComm>(comm), PHOTON_COLL_REDUCE,
			       (photon_cid){0}, &req, PHOTON_REQ_PWC_NO_LCE));
  Check(photon_collective_join(req, (void*)in, out, count, count, dt, dt, root, op));
}

void PhotonNetwork::gather(const void *in, int scount, datatype_t sdt, void *out,
			int rcount, datatype_t rdt, int root, PGFem3D_Comm comm)
{
  photon_rid req;
  Check(photon_collective_init(static_cast<photonComm>(comm), PHOTON_COLL_GATHER,
			       (photon_cid){0}, &req, PHOTON_REQ_PWC_NO_LCE));
  Check(photon_collective_join(req, (void*)in, out, scount, rcount, sdt, rdt, root,
			       PHOTON_OP_NULL));
}

void PhotonNetwork::allreduce(const void *in, void *out, int count, datatype_t dt,
			  op_t op, PGFem3D_Comm comm)
{
  photon_rid req;
  Check(photon_collective_init(static_cast<photonComm>(comm), PHOTON_COLL_ALLREDUCE,
			       (photon_cid){0}, &req, PHOTON_REQ_PWC_NO_LCE));
  Check(photon_collective_join(req, (void*)in, out, count, count, dt, dt, 0, op));
}

void PhotonNetwork::allgather(const void *in, int scount, datatype_t sdt,
			     void *out, int rcount, datatype_t rdt, PGFem3D_Comm comm)
{
  photon_rid req;
  Check(photon_collective_init(static_cast<photonComm>(comm), PHOTON_COLL_ALLGATHER,
			       (photon_cid){0}, &req, PHOTON_REQ_PWC_NO_LCE));
  Check(photon_collective_join(req, (void*)in, out, scount, rcount, sdt, rdt, 0,
			       PHOTON_OP_NULL));
}

// current Photon collective interface does not support vectored versions
void PhotonNetwork::allgatherv(const void *in, int scount, datatype_t sdt,
			       void *out, int *rcounts, const int *dipls,
			       datatype_t rdt, PGFem3D_Comm comm)
{
  MPI_Allgatherv(in, scount, static_cast<MPI_Datatype>(sdt), out, rcounts, dipls,
		 static_cast<MPI_Datatype>(rdt), static_cast<MPI_Comm>(comm));
}

void PhotonNetwork::pin(const void *addr, const size_t bytes, Key *key)
{
  Check(photon_register_buffer((void*)addr, bytes));
  if (key) {
    const struct photon_buffer_priv_t *priv;
    Check(photon_get_buffer_private((void*)addr, bytes, &priv));
    key->key0 = priv->key0;
    key->key1 = priv->key1;
  }
}

void PhotonNetwork::unpin(const void *addr, const size_t bytes)
{
  Check(photon_unregister_buffer((void*)addr, bytes));
}

void PhotonNetwork::pwc(int dst, size_t n, Buffer *lbuf, Buffer *rbuf,
			const CID& lid, const CID& rid)
{
  struct photon_buffer_t lb = {0};
  if (lbuf) {
    lb.addr = lbuf->addr;
    lb.size = n;
    lb.priv.key0 = lbuf->key.key0;
    lb.priv.key1 = lbuf->key.key1;
  }
  
  struct photon_buffer_t rb = {0};
  if (rbuf) {
    rb.addr = rbuf->addr;
    rb.size = n;
    rb.priv.key0 = rbuf->key.key0;
    rb.priv.key1 = rbuf->key.key1;
  }
  
  photon_cid li;
  li.u64 = lid;
  li.size = 0;

  photon_cid ri;
  ri.u64 = rid;
  ri.size = 0;

  Check(photon_put_with_completion(dst, n, &lb, &rb, li, ri, 0));
}

void PhotonNetwork::gwc(int src, size_t n, Buffer *lbuf, Buffer *rbuf,
			const CID& lid, const CID& rid)
{
  struct photon_buffer_t lb = {0};
  if (lbuf) {
    lb.addr = lbuf->addr;
    lb.size = n;
    lb.priv.key0 = lbuf->key.key0;
    lb.priv.key1 = lbuf->key.key1;
  }
    
  struct photon_buffer_t rb = {0};
  if (rbuf) {
    rb.addr = rbuf->addr;
    rb.size = n;
    rb.priv.key0 = rbuf->key.key0;
    rb.priv.key1 = rbuf->key.key1;
  }
  
  photon_cid li;
  li.u64 = lid;
  li.size = 0;

  photon_cid ri;
  ri.u64 = rid;
  ri.size = 0;
  
  Check(photon_get_with_completion(src, n, &lb, &rb, li, ri, 0));
}

void PhotonNetwork::probe(int *flag, CID *comp, Status *status, void (*cb)(CID))
{
  photon_cid rid;
  int src;
  Check(photon_probe_completion(PHOTON_ANY_SOURCE, flag, NULL, &rid, &src, NULL,
				PHOTON_PROBE_LEDGER));
  if (comp) {
    *comp = rid.u64;
    if (cb) {
      cb(*comp);
    }
  }
  
  if (status) {
    status->NET_SOURCE = src;
  }
}

void PhotonNetwork::progress(int *flag, CID *comp, Status *status, void (*cb)(CID))
{
  photon_cid rid;
  int src;
  Check(photon_probe_completion(PHOTON_ANY_SOURCE, flag, NULL, &rid, &src, NULL,
				PHOTON_PROBE_EVQ));
  if (comp) {
    *comp = rid.u64;
    if (cb) {
      cb(*comp);
    }
  }
  
  if (status) {
    status->NET_SOURCE = src;
  }
}

void PhotonNetwork::wait_n(int count)
{
  do {
    photon_cid rid;
    int src, flag;
    Check(photon_probe_completion(PHOTON_ANY_SOURCE, &flag, NULL, &rid, &src, NULL,
				  PHOTON_PROBE_EVQ));
    if (flag)
      --count;
  } while (count);
}

void PhotonNetwork::wait_n_id(int count, const CID& id)
{
  do {
    photon_cid rid;
    int src, flag;
    Check(photon_probe_completion(PHOTON_ANY_SOURCE, &flag, NULL, &rid, &src, NULL,
				  PHOTON_PROBE_EVQ));
    if (flag && rid.u64 == id)
      --count;
  } while (count);
}
