noinst_LTLIBRARIES = libhypre.la

libhypre_la_SOURCES = Hypre.cc \
                      Preconditioner.cc Euclid.cc ParaSails.cc Pilut.cc \
                      Boomer.cc ScaleDiag.cc Jacobi.cc \
                      Solver.cc GMRes.cc BCGStab.cc AMG.cc Flex.cc Hybrid.cc \
                      boomerAMGInterface.cc

libhypre_la_CPPFLAGS = -Wall -Werror \
                       -I$(top_srcdir)/include \
                       -I$(top_builddir)/include \
                       $(LIBGCM_CPPFLAGS) \
                       $(LAPACK_INCLUDE) \
                       $(HYPRE_INCLUDE)
