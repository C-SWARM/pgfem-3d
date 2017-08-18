noinst_PROGRAMS = shear tension_disp tension_pressure UQ_study \
                  crystal_plasticity crystal_plasticity_verification_test

shear_SOURCES = shear.cc
tension_disp_SOURCES = tension_disp.cc
tension_pressure_SOURCES = tension_pressure.cc
UQ_study_SOURCES = UQ_study.cc
crystal_plasticity_SOURCES = crystal_plasticity.cc
crystal_plasticity_verification_test_SOURCES = crystal_plasticity_verification_test.cc

# Universal flags for all binaries
# AM_CXXFLAGS = 

AM_CPPFLAGS = -Wall -Werror \
              -I$(top_srcdir)/include \
              -I$(top_builddir)/include \
              -I$(top_srcdir)/lib/VTK_IO/include

LDADD       = $(top_builddir)/src/libpgfem3d.a \
              $(top_builddir)/lib/VTK_IO/src/libPGFem3D_to_VTK.a \
              $(top_builddir)/energy_equation/src/libenergy_equation.a \
              $(suitesparse_lib) \
              $(hypre_lib) \
              $(cnstvm_lib) \
              $(mkl_lib)

AM_LDFLAGS  = -Wl,--allow-multiple-definition
