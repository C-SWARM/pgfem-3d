SUBDIRS = hypre trilinos
noinst_LTLIBRARIES     = libsolvers.la
libsolvers_la_SOURCES  = Solver.cc SparseSystem.cc
libsolvers_la_CPPFLAGS = -Wall -Werror \
                         -I$(top_srcdir)/include \
                         -I$(top_builddir)/include \
                         $(LIBGCM_CPPFLAGS) \
                         $(LAPACK_INCLUDE) \
                         $(HYPRE_INCLUDE)
libsolvers_la_LIBADD   = hypre/libhypre.la trilinos/libtrilinos.la
