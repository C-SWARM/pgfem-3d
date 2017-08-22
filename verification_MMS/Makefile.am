noinst_PROGRAMS = PGFem3D verification_tf_MMS

PGFem3D_SOURCES = ../src/super_main.cc ../src/dynamics.cc MMS.cc
verification_tf_MMS_SOURCES = ../src/dynamics.cc verification_tf_MMS.cc MMS.cc

# Universal flags for all binaries
# AM_CXXFLAGS = 

AM_CPPFLAGS = -Wall -Werror \
              -I$(top_srcdir)/include \
              -I$(top_builddir)/include \
              -I$(top_srcdir)/lib/VTK_IO/include \
              $(hypre_include)

LDADD       = $(top_builddir)/src/libpgfem3d.a \
              $(top_builddir)/lib/VTK_IO/src/libPGFem3D_to_VTK.a \
              $(top_builddir)/energy_equation/src/libenergy_equation.a \
              $(suitesparse_lib) \
              $(hypre_lib) \
              $(VTK_LIBS) \
              $(cnstvm_lib) \
              $(mkl_lib)

AM_LDFLAGS  = -Wl,--allow-multiple-definition \
              $(VTK_LDFLAGS)
