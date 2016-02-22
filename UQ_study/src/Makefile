EXC = main_test
OBJ_TEST=$(shell ls ../../src/*.o | grep -v 'super_main.o' | grep -v plasticity_model_BPA_driver.o | grep -v cm_iso_viscous_damage_driver.o | grep -v unit_test_43_driver.o | grep -v cm_j2_plasticity_driver.o)

UQ_stdy_inc=../include

PGFem3D_dir=../..
pfem_incl_dir = $(PGFem3D_dir)/include

# SUITESPARSE
SUITESPARSE_lib	= -L${PGFEM3D_EXT_LIB}/UMFPACK/UMFPACK/Lib -lumfpack \
		  -L${PGFEM3D_EXT_LIB}/UMFPACK/AMD/Lib -lamd

# HYPRE
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


CC = mpicc
LINK = mpicxx

RM = rm -rf

CFLAGS = -Wall -std=c99 -O3 -g

INCLUDES = -I$(pfem_incl_dir) \
                  $(CNSTVM_incl) \
		  $(VTK_IO_incl) \
		  $(MKL_incl) \
		  -I$(UQ_stdy_inc)

LIBS	= $(VTK_IO_lib) \
	  $(SUITESPARSE_lib) \
	  $(HYPRE_lib) \
          $(CNSTVM_lib) \
	  $(VTK_lib) \
	  $(MKL_lib)

OBJ = $(SRC:.c=.o)

all:
	$(CC) -c $(INCLUDES) $(CFLAGS) get_max_disp.c 		
	$(LINK) -o get_max_disp $(CFLAGS) get_max_disp.o  $(OBJ_TEST) $(LIBS)

clean:  
	rm -f *.out *.o *.mod core* get_max_disp