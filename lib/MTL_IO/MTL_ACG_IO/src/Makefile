SRC = $(shell ls *.cc)

OBJ = $(SRC:.cc=.o)
OLIB = libMTL_ACG.a

CC = gcc
CXX = g++

INCDIRS = -I../include
CFLAGS = 

all: $(OBJ) $(OLIB)

%.o:%.cc
	$(CXX) $(CFLAGS) $(INCDIRS) -c $<

$(OLIB): $(OBJ)
	$(AR) $(ARFLAGS) $(OLIB) $(OBJ)	
	rm -rf $(OBJ)
	mv $(OLIB) ../lib

clean :
	rm -rf *.o
	
