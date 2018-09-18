noinst_LTLIBRARIES = libtrilinos.la

libtrilinos_la_SOURCES = Trilinos.cc

libtrilinos_la_CPPFLAGS = -Wall -Werror \
			-Wno-strict-aliasing \
			-Wno-unused-function \
			-Wno-unused-local-typedefs \
			-Wno-deprecated-declarations \
			-Wno-unused-variable \
			-Wno-sign-compare \
			-I$(top_srcdir)/include \
			-I$(top_builddir)/include \
			$(LIBGCM_CPPFLAGS) \
			$(LAPACK_INCLUDE) \
			$(TRILINOS_INCLUDE)

