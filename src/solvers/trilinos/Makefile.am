noinst_LTLIBRARIES = libtrilinos.la

libtrilinos_la_SOURCES = Trilinos.cc

libtrilinos_la_CPPFLAGS = -Wall -Werror \
                       -I$(top_srcdir)/include \
                       -I$(top_builddir)/include \
                       $(LIBGCM_CPPFLAGS) \
                       $(LAPACK_INCLUDE) \
                       $(TRILINOS_INCLUDE)
