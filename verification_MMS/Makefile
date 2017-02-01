OBJ_TEST = $(shell ls ../src/*.o | grep -v 'super_main.o' | grep -v 'dynamics.o' | grep -v '[a-z]*'_driver.o)

pfem_incl_dir = $(TEST_DIR)/include

# SUITESPARSE
SUITESPARSE_lib	= -L${PGFEM3D_EXT_LIB}/UMFPACK/UMFPACK/Lib -lumfpack \
		  -L${PGFEM3D_EXT_LIB}/UMFPACK/AMD/Lib -lamd

# HYPRE
HYPRE_lib	= -L$(PGFEM3D_EXT_LIB)/hypre-2.4.0b/src/lib -lHYPRE

# CNSTVM
CNSTVM_dir = $(shell grep 'cnstvm_lib' ../Makefile* | grep -v '../Makefile.in:cnstvm_lib = ' | grep -o '[^-L]*' | grep 'Generalizsed_constitutive_model' | sed 's/\/lib //')
CNSTVM_lib = -L$(CNSTVM_dir)/lib -lConstitutiveModel
CNSTVM_incl = -I$(CNSTVM_dir)/utils/include \
              -I$(CNSTVM_dir)/material/include \
              -I$(CNSTVM_dir)/elasticity/include \
              -I$(CNSTVM_dir)/crystal_plasticity/include

# Compressed VTK I/O
VTK_IO_dir	= $(TEST_DIR)/lib/VTK_IO
VTK_IO_incl	= -I$(VTK_IO_dir)/include
VTK_IO_lib	= -L$(VTK_IO_dir)/src -lPGFem3D_to_VTK
#VTK_IO_lib	= $(VTK_IO_dir)/src/libPGFem3D_to_VTK.a

# VTK
VTK_dir	 = /afs/crc.nd.edu/group/cswarm/VTK/7.0.0/intel/15.0/mvapich2/2.1/static-lib
VTK_incl = -I$(VTK_dir)/include/vtk-7.0
VTK_lib= -Wl,-rpath=$vtk_dir/lib -L$vtk_dir/lib \
         -lvtkIOXML-7.0 -lvtkIOXMLParser-7.0 -lvtkIOCore-7.0 -lvtkzlib-7.0 \
         -lvtkCommonExecutionModel-7.0 -lvtkCommonDataModel-7.0 \
         -lvtkCommonMisc-7.0 -lvtkCommonSystem-7.0 -lvtkCommonTransforms-7.0 \
         -lvtkCommonMath-7.0 -lvtkIOGeometry-7.0 -lvtkCommonCore-7.0 \
         -lvtksys-7.0 -lvtkexpat-7.0


VTK_lib	 = -Wl,-rpath=$(VTK_dir)/lib \
		  -L$(VTK_dir)/lib -lvtkIOXML-7.0 -lvtkIOXMLParser-7.0 \
		  -lvtkIOCore-7.0 -lvtkzlib-7.0 \
		  -lvtkCommonExecutionModel-7.0 -lvtkCommonDataModel-7.0 \
		  -lvtkCommonMisc-7.0 -lvtkCommonSystem-7.0 \
		  -lvtkCommonTransforms-7.0 -lvtkCommonMath-7.0 \
		  -lvtkIOGeometry-7.0 -lvtkCommonCore-7.0 \
	 	  -lvtksys-7.0 -lvtkexpat-7.0

# Intel MKL: MKLROOT must be defined or define it here
MKL_incl  = -I$(MKLROOT)/include
MKL_lib   = $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
	          -Wl,--start-group \
	          $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	          $(MKLROOT)/lib/intel64/libmkl_sequential.a \
	          $(MKLROOT)/lib/intel64/libmkl_core.a \
	          $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	          -Wl,--end-group -lpthread

CXX = mpicxx
LINK = mpicxx

RM = rm -rf

DFLAGS = -DVERIFICATION_USING_MMS=1

CXXFLAGS = -Wall -std=c++14 -fpermissive -Ofast -g 
CXXFLAGS += $(DFLAGS)

INCLUDES = -I$(pfem_incl_dir) \
		  $(CNSTVM_incl) \
		  $(VTK_IO_incl) \
		  $(MKL_incl)
	  
LIBS	= $(VTK_IO_lib) \
	  $(SUITESPARSE_lib) \
	  $(HYPRE_lib) \
	  $(CNSTVM_lib) \
	  $(VTK_lib) \
	  $(MKL_lib) 

#SRC = verification_tf_MMS.cc
OBJ = $(SRC:.cc=.o)

all:	
#	$(CXX) -c $(INCLUDES) $(CXXFLAGS) MMS.cc 		
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) ../src/dynamics.cc
	$(LINK) -o PGFem3D $(CXXFLAGS) $(OBJS) $(LIBS) ../src/super_main.o MMS.o dynamics.o $(OBJ_TEST) $(LIBS)
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) verification_tf_MMS.cc 		
	$(LINK) -o verification_tf_MMS $(CXXFLAGS) verification_tf_MMS.o MMS.o dynamics.o $(OBJ_TEST) $(LIBS)


 							
clean:  
	rm -f *.out *.o *.mod core* PGFem3D verification_tf_MMS
