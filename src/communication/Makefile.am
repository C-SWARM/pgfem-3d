SUBDIRS = 
noinst_LTLIBRARIES     = libcommunication.la

libcommunication_la_SOURCES  = Communication.cc           \
			       ISIR_SparseComm.cc         \
			       PWC_SparseComm.cc          \
			       MultiscaleComm.cc          \
			       MultiscaleCommInfo.cc      \
			       MultiscaleServerContext.cc \
			       Network.cc                 \
			       Boot.cc

libcommunication_la_CPPFLAGS = -Wall -Werror \
                               -I$(top_srcdir)/include   \
                               -I$(top_builddir)/include \
                               $(LIBGCM_CPPFLAGS)        \
                               $(LAPACK_INCLUDE)

if HAVE_PGFEM_MPI
libcommunication_la_SOURCES += mpi/MPINetwork.cc
endif

if HAVE_PHOTON
libcommunication_la_SOURCES += photon/PhotonNetwork.cc
endif
