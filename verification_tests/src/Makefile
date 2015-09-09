EXC = main_test
OBJ_TEST=$(shell ls ../../src/*.o | grep -v 'super_main.o' | grep -v 'plasticity_model_BPA_driver.o')

POST_PROCESSING=../include

pfem_incl_dir = $(TEST_DIR)/include

# SUITESPARSE
SUITESPARSE_lib	= -L${PGFEM3D_EXT_LIB}/UMFPACK/UMFPACK/Lib -lumfpack \
		  -L${PGFEM3D_EXT_LIB}/UMFPACK/AMD/Lib -lamd

# HYPRE
HYPRE_lib	= -L$(PGFEM3D_EXT_LIB)/hypre-2.4.0b/src/lib -lHYPRE

# CNSTVM
CNSTVM_dir = $(HOME)/temp/Constitutive_Model_implicit
CNSTVM_lib = -L$(CNSTVM_dir)/src -lCM
CNSTVM_incl = -I$(CNSTVM_dir)/include

# Compressed VTK I/O
VTK_IO_dir	= $(TEST_DIR)/lib/VTK_IO
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
		  -I$(POST_PROCESSING)
	  
LIBS	= $(VTK_IO_lib) \
	  $(SUITESPARSE_lib) \
	  $(HYPRE_lib) \
          $(CNSTVM_lib) \
	  $(VTK_lib) \
	  $(MKL_lib)
	  	  		  		  		
OBJ = $(SRC:.c=.o)

all:
#	$(CC) -c $(INCLUDES) $(CFLAGS) post_processing.c 	
	$(CC) -c $(INCLUDES) $(CFLAGS) shear.c 		
	$(LINK) -o shear $(CFLAGS) shear.o  $(OBJ_TEST) $(LIBS)
	$(CC) -c $(INCLUDES) $(CFLAGS) tension_disp.c 		
	$(LINK) -o tension_disp $(CFLAGS) tension_disp.o  $(OBJ_TEST) $(LIBS)
	$(CC) -c $(INCLUDES) $(CFLAGS) tension_pressure.c 		
	$(LINK) -o tension_pressure $(CFLAGS) tension_pressure.o $(OBJ_TEST) $(LIBS)		
	$(CC) -c $(INCLUDES) $(CFLAGS) UQ_study.c 		
	$(LINK) -o UQ_study $(CFLAGS) UQ_study.o $(OBJ_TEST) $(LIBS)	
	$(CC) -c $(INCLUDES) $(CFLAGS) crystal_plasticity.c
	$(LINK) -o crystal_plasticity $(CFLAGS) crystal_plasticity.o $(OBJ_TEST) $(LIBS)
 							
clean:  
	rm -f *.out *.o *.mod core* shear tension_disp tension_pressure
