SUBDIRS = hypre
noinst_LTLIBRARIES     = libsolvers.la
libsolvers_la_SOURCES  = Solver.cc SparseSystem.cc
libsolvers_la_CPPFLAGS = -Wall -Werror \
                         -I$(top_srcdir)/include \
                         -I$(top_builddir)/include \
                         $(LIBGCM_CPPFLAGS) \
                         $(LAPACK_INCLUDE) \
                         $(HYPRE_INCLUDE)  \
                         $(MSNET_INCLUDE)

libsolvers_la_LIBADD   = hypre/libhypre.la

if HAVE_TRILINOS
SUBDIRS += trilinos
libsolvers_la_LIBADD += trilinos/libtrilinos.la
libsolvers_la_CPPFLAGS += -Wno-strict-aliasing \
			-Wno-unused-function \
			-Wno-unused-local-typedefs \
			$(TRILINOS_INCLUDE)
endif
