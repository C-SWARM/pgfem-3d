#!/bin/bash

EXC = main_test
OBJ_TEST=$(shell ls ../../src/*.o | grep -v 'super_main.o' | grep -v '[a-z]*'_driver.o)

POST_PROCESSING=../include

pgfem_incl_dir = $(pgfem_root_dir)/include

# SUITESPARSE
suitesparse_dir := $(shell grep -o 'with-suitesparse-dir=[^ ]*' ../../config.log | sed -e 's/^.*=//' )
SUITESPARSE_lib	= -L${suitesparse_dir}/UMFPACK/Lib -lumfpack \
            		  -L${suitesparse_dir}/AMD/Lib -lamd

# HYPRE
hypre_dir := $(shell grep -o 'with-hypre-dir=[^ ]*' ../../config.log | sed -e 's/^.*=//' )
HYPRE_inc	= -I$(hypre_dir)/include
HYPRE_lib	= -L$(hypre_dir)/lib -lHYPRE

# CNSTVM
CNSTVM_dir := $(shell grep -o 'with-cnstvm-dir=[^ ]*' ../../config.log | sed -e 's/^.*=//' )
CNSTVM_lib = -L$(CNSTVM_dir)/lib -lConstitutiveModel
CNSTVM_incl = -I$(CNSTVM_dir)/utils/include \
              -I$(CNSTVM_dir)/material/include \
              -I$(CNSTVM_dir)/elasticity/include \
              -I$(CNSTVM_dir)/crystal_plasticity/include

# Compressed VTK I/O
VTK_IO_dir	= $(pgfem_root_dir)/lib/VTK_IO
VTK_IO_incl	= -I$(VTK_IO_dir)/include
VTK_IO_lib	= -L$(VTK_IO_dir)/src -lPGFem3D_to_VTK

energy_equation_include = -I$(pgfem_root_dir)/energy_equation/include
energy_equation_lib     = -L$(pgfem_root_dir)/energy_equation/src -lPGFem3D_Energy_equation

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
	     	  -Wl,--end-group -lpthread -lm -ldl


CXX = mpicxx
LINK = mpicxx

CXXFLAGS = -Wall -std=c++14 -fpermissive -Ofast -g

INCLUDES = -I$(pgfem_incl_dir) \
	    $(HYPRE_inc) \
	    $(CNSTVM_incl) \
	    $(VTK_IO_incl) \
	    $(MKL_incl) \
	   -I$(POST_PROCESSING) \
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
#	$(CXX) -c $(INCLUDES) $(CXXFLAGS) post_processing.cc 	
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) shear.cc 		
	$(LINK) -o shear $(CXXFLAGS) shear.o  $(OBJ_TEST) $(LIBS)
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) tension_disp.cc 		
	$(LINK) -o tension_disp $(CXXFLAGS) tension_disp.o  $(OBJ_TEST) $(LIBS)
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) tension_pressure.cc 		
	$(LINK) -o tension_pressure $(CXXFLAGS) tension_pressure.o $(OBJ_TEST) $(LIBS)		
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) UQ_study.cc 		
	$(LINK) -o UQ_study $(CXXFLAGS) UQ_study.o $(OBJ_TEST) $(LIBS)	
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) crystal_plasticity.cc
	$(LINK) -o crystal_plasticity $(CXXFLAGS) crystal_plasticity.o $(OBJ_TEST) $(LIBS)
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) crystal_plasticity_verification_test.cc
	$(LINK) -o crystal_plasticity_verification_test $(CXXFLAGS) crystal_plasticity_verification_test.o $(OBJ_TEST) $(LIBS) 							
clean:  
	rm -f *.out *.o *.mod core* shear tension_disp tension_pressure UQ_study crystal_plasticity  crystal_plasticity_verification_test
