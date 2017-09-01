bin_PROGRAMS = shear \
               tension_disp \
               tension_pressure \
               UQ_study \
               crystal_plasticity \
               crystal_plasticity_verification_test

shear_SOURCES                                = shear.cc
tension_disp_SOURCES                         = tension_disp.cc
tension_pressure_SOURCES                     = tension_pressure.cc
UQ_study_SOURCES                             = UQ_study.cc
crystal_plasticity_SOURCES                   = crystal_plasticity.cc
crystal_plasticity_verification_test_SOURCES = crystal_plasticity_verification_test.cc

AM_CPPFLAGS = -Wall -Werror \
              -I$(top_srcdir)/include \
              -I$(top_srcdir)/lib/VTK_IO/include \
              -I$(top_builddir)/include \
               $(HYPRE_INCLUDE)

AM_LDFLAGS  = $(VTK_LDFLAGS)
LDADD       = $(top_builddir)/src/libpgfem3d.la
