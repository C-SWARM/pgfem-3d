SRC = $(shell ls *.c)
OBJs = $(SRC:.c=.o)
PROG = test_mtl_acg

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_mtl_acg = $(DIR_TOP)/MTL_ACG_IO

MTL_ACG_lib = -L$(DIR_mtl_acg)/lib -lMTL_ACG

CC = gcc
CXX = g++
LINK = g++

AR       = ar
ARFLAGS = -rcs

INCs = -I$(DIR_mtl_acg)/include

LIBs = $(MTL_ACG_lib)
CFLAGS = -Wall -std=c99 -O3 -g

all: $(OBJs)
	$(LINK) $(CFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.c
	$(CC) $(CFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o
	
