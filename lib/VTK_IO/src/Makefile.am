noinst_LTLIBRARIES = libPGFem3D_to_VTK.la

libPGFem3D_to_VTK_la_SOURCES = PGFem3D_to_VTK.cc

libPGFem3D_to_VTK_la_CPPFLAGS = -I$(top_srcdir)/include \
                                -I$(srcdir)/../include \
                                -I$(top_builddir)/include \
                                $(LIBGCM_CPPFLAGS) \
                                $(LAPACK_INCLUDE) \
				$(MSNET_INCLUDE) \
                                $(VTK_CPPFLAGS)

libPGFem3D_to_VTK_la_CXXFLAGS = $(VTK_CXXFLAGS)
