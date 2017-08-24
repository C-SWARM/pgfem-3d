noinst_LTLIBRARIES = libenergy_equation.la

libenergy_equation_la_CPPFLAGS = -I$(srcdir)/../include \
                                 -I$(top_srcdir)/include \
                                 -I$(top_builddir)/include \
                                 $(hypre_include) \
                                 $(cnstvm_include) \
                                 $(ttl_include) \
                                 $(mkl_include) 

libenergy_equation_la_SOURCES  = energy_equation.cc
