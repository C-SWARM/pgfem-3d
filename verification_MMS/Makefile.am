noinst_PROGRAMS = PGFem3D \
                  verification_tf_MMS

PGFem3D_SOURCES             = super_main.cc dynamics.cc MMS.cc
verification_tf_MMS_SOURCES = verification_tf_MMS.cc dynamics.cc MMS.cc

AM_CPPFLAGS = -Wall -Werror \
              -DVERIFICATION_USING_MMS=1 \
              -I$(top_srcdir)/include \
              -I$(top_srcdir)/lib/VTK_IO/include \
              -I$(top_builddir)/include \
              $(HYPRE_INCLUDE)

LDADD       = $(top_builddir)/src/libpgfem3d.la

# This is a hack to get around the fact that there are conflicting symbols in
# the Generalizsed_constitutive_model.
#
# https://gitlab-cswarm.crc.nd.edu/Constitutive_Model/Generalizsed_constitutive_model/issues/9
#
AM_LDFLAGS  = $(HYPRE_LDFLAGS) -Wl,--allow-multiple-definition 
