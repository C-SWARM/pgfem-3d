SRC = $(shell ls *.c)
OBJs = $(SRC:.c=.o)
PROG = test_mtl_acg

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_pmtl_acg = $(DIR_TOP)/PMTL_ACG_IO

PMTL_ACG_lib = -L$(DIR_pmtl_acg)/lib -lPMTL_ACG
PMTL_lib = -L/opt/crc/scilib/mtl/boost_1_58_0/build_boost_mvapich2/lib -lboost_mpi -lboost_serialization \
	         -L/opt/crc/scilib/mtl/src/PMTL-all-4.0.9555-Linux/usr/lib -lpmtl4_security_mvapich2.2
PMTL_inc = -I/opt/crc/scilib/mtl/boost_1_58_0/build_boost_mvapich2/include

CC = mpicxx
CXX = mpicxx
LINK = mpicxx

INCs = -I$(DIR_pmtl_acg)/include $(PMTL_inc)

LIBs = $(PMTL_ACG_lib) $(PMTL_lib)
#CFLAGS = -Wall -std=c99 -O3 -g -DMTL_HAS_MPI
CFLAGS = -Wall -O3 -g -DMTL_HAS_MPI

all: $(OBJs)
	$(LINK) $(CFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.c
	$(CC) $(CFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o
	
