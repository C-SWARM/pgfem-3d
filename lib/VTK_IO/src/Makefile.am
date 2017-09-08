noinst_LTLIBRARIES = libPGFem3D_to_VTK.la

libPGFem3D_to_VTK_la_SOURCES = PGFem3D_to_VTK.cc

libPGFem3D_to_VTK_la_CPPFLAGS = -I$(top_srcdir)/include \
                                -I$(srcdir)/../include \
                                -I$(top_builddir)/include \
                                $(GCM_INCLUDE) \
                                $(LAPACK_INCLUDE) \
                                $(VTK_CPPFLAGS)

libPGFem3D_to_VTK_la_CXXFLAGS = $(VTK_CXXFLAGS)
