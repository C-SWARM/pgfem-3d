noinst_LIBRARIES = libenergy_equation.a

# libenergy_equation_a_CXXFLAGS = 
libenergy_equation_a_CPPFLAGS = -I$(srcdir)/../include \
                                -I$(top_srcdir)/include \
                                -I$(top_builddir)/include \
                                $(hypre_include) \
                                $(cnstvm_include) \
                                $(mkl_include) 

libenergy_equation_a_SOURCES  = energy_equation.cc
