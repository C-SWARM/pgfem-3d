SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)
OLIB = libPMTL_ACG.a
PROG = test_mtl_acg

CC = mpicc
CXX = mpicxx
LINK = mpicxx

PMTL_lib = -L/opt/crc/scilib/mtl/boost_1_58_0/build_boost_mvapich2/lib -lboost_mpi -lboost_serialization \
	         -L/opt/crc/scilib/mtl/src/PMTL-all-4.0.9555-Linux/usr/lib -lpmtl4_security_mvapich2.2
PMTL_inc = -I/opt/crc/scilib/mtl/boost_1_58_0/build_boost_mvapich2/include

INCDIRS = -I../include $(PMTL_inc)
LIBs = $(PMTL_ACG_lib) $(PMTL_lib)

CFLAGS = -Wall -O3 -g -DMTL_HAS_MPI

all: $(OBJs)
	$(LINK) $(CFLAGS) $(OBJs) $(LIBs) -o $(PROG)  

%.o:%.cc
	$(CXX) $(CFLAGS) $(INCDIRS) -c $<

clean :
	rm -rf *.o
	