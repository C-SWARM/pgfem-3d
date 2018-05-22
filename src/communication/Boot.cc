#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "PGFEM_io.h"
#include "pgfem3d/Boot.hpp"
#include "pgfem3d/Communication.hpp"
#include <string.h>

using namespace pgfem3d;
using namespace pgfem3d::net;

int pgfem3d::net::NET_MAX_PROCESSOR_NAME;

// maintain a global handle to the initial boot class
Boot* pgfem3d::net::gboot = 0;

Boot::Boot(net_boot_t type) {
  if (type == NET_BOOT_DEFAULT){
    _type = NET_BOOT_MPI;
  }
  else {
    _type = type;
  }
  
  switch (_type) {
  case NET_BOOT_MPI:
    {
#ifdef HAVE_MPI
      int init;
      if (MPI_Initialized(&init)) {
	throw PGFEM_printerr("MPI initialization failed\n");
      }
      
      if (!init) {
	int level, thread_level = MPI_THREAD_SINGLE;
	if (MPI_Init_thread(NULL, NULL, thread_level, &level)) {
	  throw PGFEM_printerr("MPI initialization failed\n");
	}
	
	if (level != thread_level) {
	  throw PGFEM_printerr("MPI thread level failed: requested %d, received %d\n",
			       thread_level, level);
	}
      }
      
      if (MPI_Comm_rank(MPI_COMM_WORLD, &rank)) {
	throw PGFEM_printerr("Could not get rank\n");
      }
      
      if (MPI_Comm_size(MPI_COMM_WORLD, &nproc)) {
	throw PGFEM_printerr("Could not get the number of ranks\n");
      }

      NET_MAX_PROCESSOR_NAME = MPI_MAX_PROCESSOR_NAME;
      processor_name = new char[NET_MAX_PROCESSOR_NAME];
      if (!processor_name) {
	throw PGFEM_printerr("Could not allocate space for processor name\n");
      }
      
      if (MPI_Get_processor_name(processor_name, &pname_len)) {
	throw PGFEM_printerr("Could not get the processor name\n");
      }
#else
      PGFEM_printerr("MPI bootstrap not supported in current configuration\n");
      PGFEM_Abort();
#endif
    }
    break;
  case NET_BOOT_PMI:
    PGFEM_printerr("PMI bootstrap currently not supported\n");
    break;
  default:
    break;
  }
  
  // save a handle to the initial instantiated boot class
  if (!gboot)
    gboot = this;
}

Boot::~Boot() {
  switch (_type) {
  case NET_BOOT_MPI: {
#ifdef HAVE_MPI
    int fin;
    if (MPI_Finalized(&fin)) {
      PGFEM_printerr("MPI finalize failed\n");
    }
    if (!fin) {
      if (MPI_Finalize()) {
	PGFEM_printerr("MPI finalize failed\n");
      }
    }
    if (processor_name) {
      delete [] processor_name;
    }
#endif
    break;
  }
  case NET_BOOT_PMI:
    PGFEM_printerr("PMI bootstrap currently not supported\n");
    break;
  default:
    break;
  }
}

void Boot::get_processor_name(char *name, int *len)
{
  memcpy(name, processor_name, pname_len);
  *len = pname_len;
}
