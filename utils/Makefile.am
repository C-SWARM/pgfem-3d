bin_PROGRAMS = gen_restart_from_NP2NP gen_co_from_NP2NP

gen_restart_from_NP2NP_SOURCES = gen_restart_from_NP2NP.cc
gen_co_from_NP2NP_SOURCES = gen_co_from_NP2NP.cc

# Universal flags for all binaries
AM_CPPFLAGS = -Wall -Werror \
              -I$(top_srcdir)/include \
              -I$(top_srcdir)/FEMLIB/include \
              -I$(top_srcdir)/lib/VTK_IO/include \
              -I$(top_srcdir)/lib/fastmathparser/include \
              -I$(top_srcdir)/energy_equation/include \
              -I$(top_builddir)/include \
              $(HYPRE_INCLUDE) \
              $(SUITESPARSE_INCLUDE) \
              $(TTL_INCLUDE) \
              $(LIBGCM_CPPFLAGS) \
              $(LAPACK_INCLUDE) \
              $(MSNET_INCLUDE)

AM_LDFLAGS  = $(VTK_LDFLAGS)

LDADD       = $(top_builddir)/src/libpgfem3d.la


