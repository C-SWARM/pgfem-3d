#ifndef PGFEM3D_BOOT_H
#define PGFEM3D_BOOT_H

/// @brief This file defines the PGFem3D Boot interface
#include <cstdio>

// XXX unecessary once all explicit MPI calls removed
// from pgfem3d
#include <mpi.h>

namespace pgfem3d {
namespace net {

  enum net_boot_t {
    NET_BOOT_DEFAULT = 0,
    NET_BOOT_MPI,
    NET_BOOT_PMI
  };
  
  class Boot {
  public:
    Boot() : Boot(NET_BOOT_DEFAULT) {}
    Boot(net_boot_t type);
    ~Boot();

    net_boot_t type() const { return _type; }
    int get_rank() const { return rank; }
    int get_nproc() const { return nproc; }
    void get_processor_name(char *name, int *len);

  protected:
    net_boot_t _type;
    int rank;
    int nproc;
    char *processor_name;
    int pname_len;
  };

  extern Boot *gboot;
  extern int NET_MAX_PROCESSOR_NAME;
}
}

#endif
