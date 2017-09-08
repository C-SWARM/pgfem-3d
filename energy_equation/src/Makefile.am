noinst_LTLIBRARIES = libenergy_equation.la

libenergy_equation_la_CPPFLAGS = -Wall -Werror \
                                 -I$(srcdir)/../include \
                                 -I$(top_srcdir)/FEMLIB/include \
                                 -I$(top_srcdir)/include \
                                 -I$(top_builddir)/include \
                                 $(LIBGCM_CPPFLAGS) \
                                 $(TTL_INCLUDE) \
                                 $(LAPACK_INCLUDE)

libenergy_equation_la_SOURCES  = energy_equation.cc
