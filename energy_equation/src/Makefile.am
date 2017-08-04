noinst_LIBRARIES = libPGFem3D_Energy_equation.a
libPGFem3D_Energy_equation_a_SOURCES = energy_equation.cc

AM_CPPFLAGS = -I$(top_srcdir)/include \
              $(hypre_include) \
              $(cnstvm_include) \
              $(mkl_include) \
              -I$(srcdir)/../include

AM_CXXFLAGS = -fpermissive
