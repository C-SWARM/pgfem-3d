noinst_LIBRARIES = libPGFem3D_FEMLIB.a
libPGFem3D_FEMLIB_a_SOURCES = femlib.cc

AM_INCLUDES = -I../../include \
              $(hypre_include) \
              $(cnstvm_include) \
              $(mkl_include)

AM_CXXFLAGS = -I../include $(AM_INCLUDES)
