noinst_LIBRARIES = libPGFem3D_Energy_equation.a
libPGFem3D_Energy_equation_a_SOURCES = energy_equation.cc

AM_INCLUDES = -I../../include \
              -I../../FEMLIB/include \
              $(hypre_include) \
              $(cnstvm_include) \
              $(mkl_include)

AM_CXXFLAGS = -I../include $(AM_INCLUDES)
