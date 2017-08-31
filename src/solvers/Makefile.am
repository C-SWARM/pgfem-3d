SUBDIRS = hypre
noinst_LTLIBRARIES     = libsolvers.la
libsolvers_la_SOURCES  = Solver.cc SparseSystem.cc
libsolvers_la_CPPFLAGS = -Wall -Werror \
                         -I$(top_srcdir)/include \
                         -I$(top_builddir)/include \
                         $(hypre_include)
libsolvers_la_LIBADD   = hypre/libhypre.la 
