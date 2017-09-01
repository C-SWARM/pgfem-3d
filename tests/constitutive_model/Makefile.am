bin_PROGRAMS = BPA_driver \
               unit_test_43_driver \
               IVD_driver \
               J2D_driver \
               CH_layer_driver \
               comm_hints_driver

# Universal flags for all binaries
AM_CPPFLAGS = -Wall -Werror \
              -I$(top_srcdir)/include \
              -I$(top_srcdir)/lib/VTK_IO/include \
              -I$(top_srcdir)/energy_equation/include \
              -I$(top_builddir)/include \
              $(VTK_IO_include) \
              $(energy_equation_include) \
              $(HYPRE_INCLUDE) \
              $(SUITESPARSE_INCLUDE) \
              $(TTL_INCLUDE) \
              $(GCM_INCLUDE) \
              $(LAPACK_LIBS)

# https://gitlab-cswarm.crc.nd.edu/Constitutive_Model/Generalizsed_constitutive_model/issues/9
AM_LDFLAGS  = $(LAPACK_LDFLAGS) \
              $(HYPRE_LDFLAGS) \
              $(SUITESPARSE_LDFLAGS) \
              $(VTK_LDFLAGS) \
              -Wl,--allow-multiple-definition

LDADD       = $(top_builddir)/src/libpgfem3d.la

comm_hints_driver_SOURCES   = comm_hints_driver.cc
BPA_driver_SOURCES          = plasticity_model_BPA_driver.cc
IVD_driver_SOURCES          = cm_iso_viscous_damage_driver.cc
J2D_driver_SOURCES          = cm_j2_plasticity_driver.cc
CH_layer_driver_SOURCES     = CH_layer_driver.cc
unit_test_43_driver_SOURCES = unit_test_43_driver.cc

