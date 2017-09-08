EXC = main_test
OBJ_TEST=$(shell ls ../../src/*.o | grep -v 'super_main.o' | grep -v '[a-z]*'_driver.o)

UQ_stdy_inc=../include

PGFem3D_dir=../..
pgfem_incl_dir = $(PGFem3D_dir)/include

# SUITESPARSE
SUITESPARSE_lib	= -L${PGFEM3D_EXT_LIB}/UMFPACK/UMFPACK/Lib -lumfpack \
		  -L${PGFEM3D_EXT_LIB}/UMFPACK/AMD/Lib -lamd

# HYPRE
HYPRE_inc       = -I$(PGFEM3D_EXT_LIB)/hypre-2.4.0b/src/hypre/include
HYPRE_lib	= -L$(PGFEM3D_EXT_LIB)/hypre-2.4.0b/src/lib -lHYPRE

# CNSTVM
CNSTVM_dir = $(shell grep 'cnstvm_lib' ../../Makefile* | grep -v '../Makefile.in:cnstvm_lib = ' | grep -o '[^-L]*' | grep 'Generalizsed_constitutive_model' | sed 's/\/lib //')
CNSTVM_lib = -L$(CNSTVM_dir)/lib -lConstitutiveModel
CNSTVM_incl = -I$(CNSTVM_dir)/utils/include \
              -I$(CNSTVM_dir)/material/include \
              -I$(CNSTVM_dir)/elasticity/include \
              -I$(CNSTVM_dir)/crystal_plasticity/include

# Compressed VTK I/O
VTK_IO_dir	= $(PGFem3D_dir)/lib/VTK_IO
VTK_IO_incl	= -I$(VTK_IO_dir)/include
VTK_IO_lib	= -L$(VTK_IO_dir)/src -lPGFem3D_to_VTK

energy_equation_include = -I$(TEST_DIR)/energy_equation/include
energy_equation_lib     = -L$(TEST_DIR)/energy_equation/src -lPGFem3D_Energy_equation

# VTK
VTK_dir		= /opt/crc/vtk/5.10.1/gcc
VTK_incl	= -I$(VTK_dir)/include/vtk-6.0
VTK_lib		= -Wl,-rpath=$(VTK_dir)/lib \
		  -L$(VTK_dir)/lib -lvtkIOXML-6.0 -lvtkIOXMLParser-6.0 \
		  -lvtkIOCore-6.0 -lvtkzlib-6.0 \
		  -lvtkCommonExecutionModel-6.0 -lvtkCommonDataModel-6.0 \
		  -lvtkCommonMisc-6.0 -lvtkCommonSystem-6.0 \
		  -lvtkCommonTransforms-6.0 -lvtkCommonMath-6.0 \
		  -lvtkIOGeometry-6.0 -lvtkCommonCore-6.0 \
	 	  -lvtksys-6.0 -lvtkexpat-6.0

# Intel MKL: MKLROOT must be defined or define it here
MKL_incl	= -I$(MKLROOT)/include
MKL_lib		= $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
	     	  -Wl,--start-group \
	     	  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	     	  $(MKLROOT)/lib/intel64/libmkl_sequential.a \
	     	  $(MKLROOT)/lib/intel64/libmkl_core.a \
	     	  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	     	  -Wl,--end-group -lpthread


CXX = mpicxx
LINK = mpicxx

RM = rm -rf

CXXFLAGS = -Wall -std=c++14 -fpermissive -Ofast -g 

INCLUDES = -I$(pgfem_incl_dir) \
            $(HYPRE_inc) \
	    $(CNSTVM_incl) \
	    $(VTK_IO_incl) \
	    $(MKL_incl) \
	    -I$(UQ_stdy_inc) \
            $(energy_equation_include)

LIBS	= $(VTK_IO_lib) \
	  $(SUITESPARSE_lib) \
	  $(HYPRE_lib) \
          $(CNSTVM_lib) \
	  $(VTK_lib) \
	  $(MKL_lib) \
          $(energy_equation_lib)

OBJ = $(SRC:.cc=.o)

all:
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) get_max_disp.cc 		
	$(LINK) -o get_max_disp $(CXXFLAGS) get_max_disp.o  $(OBJ_TEST) $(LIBS)

clean:  
	rm -f *.out *.o *.mod core* get_max_disp
