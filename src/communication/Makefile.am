SUBDIRS = 
noinst_LTLIBRARIES     = libcommunication.la

libcommunication_la_SOURCES  = Communication.cc           \
			       ISIR_SparseComm.cc         \
			       PWC_SparseComm.cc

libcommunication_la_CPPFLAGS = -Wall -Werror \
                               -I$(top_srcdir)/include   \
                               -I$(top_builddir)/include \
                               $(LIBGCM_CPPFLAGS)        \
                               $(LAPACK_INCLUDE)         \
			       $(MSNET_INCLUDE)
