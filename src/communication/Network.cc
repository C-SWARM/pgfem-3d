#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFEM_io.h"
#include "pgfem3d/Communication.hpp"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdlib>

#ifdef HAVE_PHOTON
#include "photon/PhotonNetwork.hpp"
#endif

#ifdef HAVE_PGFEM_MPI
#include "mpi/MPINetwork.hpp"
#endif

using pgfem3d::net::NetworkObject;
using pgfem3d::net::Network;

namespace pgfem3d::net {
  datatype_t datatype_byte;
  datatype_t datatype_char;
  datatype_t datatype_int;
  datatype_t datatype_long;
  datatype_t datatype_double;

  op_t op_null;
  op_t op_min;
  op_t op_max;
  op_t op_sum;
  op_t op_bor;

  PGFem3D_Comm NET_COMM_WORLD;
}

NetworkObject::~NetworkObject()
{
  if (data && !_idx)
    free(data);
}

Network::~Network()
{
}

Network*
Network::Create(const PGFem3D_opt& opts)
{
  int myrank;
  PGFEM_Error_rank(&myrank);
  
  switch (opts.network) {
  case NETWORK_ISIR:
#ifdef HAVE_PGFEM_MPI
    return new isir::MPINetwork();
#else
    PGFEM_printerr("[%d]ERROR: ISIR Network not available in this build\n", myrank);
    PGFEM_Abort();
#endif
    break;
  case NETWORK_PWC:
#ifdef HAVE_PHOTON
    return new pwc::PhotonNetwork();
#else
    PGFEM_printerr("[%d]ERROR: PWC Network not available in this build\n", myrank);
    PGFEM_Abort();
#endif
    break;
  case NETWORK_ENV:
    char* network_type;
    network_type = std::getenv("PGFEM3D_NET");
    if (network_type != nullptr && strcmp(network_type, "pwc") == 0) {
      PGFEM_printf("[%d] Trying to use PWC Network from env\n", myrank);
#ifdef HAVE_PHOTON
      return new pwc::PhotonNetwork();
#else
      PGFEM_printerr("[%d]ERROR: PWC Network not available in this build\n", myrank);
      PGFEM_Abort();
#endif
    } else {
      PGFEM_printf("[%d] Trying to use ISIR Network by default\n", myrank);
#ifdef HAVE_MPI
      return new isir::MPINetwork();
#else
      PGFEM_printerr("[%d]ERROR: ISIR Network not available in this build\n", myrank);
      PGFEM_Abort();
#endif
    }
    break;
  default:
    PGFEM_printerr("[%d]ERROR: Unknown network type\n", myrank);
    PGFEM_Abort();
  }
}
