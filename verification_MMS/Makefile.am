noinst_PROGRAMS = PGFem3D \
                  verification_tf_MMS

PGFem3D_SOURCES             = super_main.cc dynamics.cc MMS.cc
verification_tf_MMS_SOURCES = verification_tf_MMS.cc dynamics.cc MMS.cc

AM_CPPFLAGS = -Wall -Werror \
              -DVERIFICATION_USING_MMS=1 \
              -I$(top_srcdir)/include \
              -I$(top_srcdir)/FEMLIB/include \
              -I$(top_srcdir)/lib/VTK_IO/include \
              -I$(top_builddir)/include \
              $(LIBGCM_CPPFLAGS) \
              $(LAPACK_INCLUDE) \
              $(HYPRE_INCLUDE)  \
	      $(MSNET_INCLUDE)

LDADD       = $(top_builddir)/src/libpgfem3d.la
