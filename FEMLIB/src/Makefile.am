noinst_LTLIBRARIES = libPGFem3D_FEMLIB.la

libPGFem3D_FEMLIB_la_CPPFLAGS = -Wall -Werror \
                                -I$(srcdir)/../include \
                                -I$(top_srcdir)/include \
                                -I$(top_builddir)/include \
                                $(LIBGCM_CPPFLAGS) \
                                $(LAPACK_INCLUDE)  \
				$(MSNET_INCLUDE)

libPGFem3D_FEMLIB_la_SOURCES  = femlib.cc

