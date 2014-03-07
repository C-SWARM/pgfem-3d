include ../../../Makefile.path

CC = mpicxx
CFLAGS = -Wall -O3

AR		= ar
ARFLAGS		= -rcs

INCLUDES = $(VTK_IO_incl) $(VTK_incl) -I$(pfem_incl_dir)

SRC = PGFem3D_to_VTK.cxx
OBJ = $(SRC:.cxx=.o)
LIB = $(VTK_IO_dir)/lib/libPGFem3D_to_VTK.a

all: $(OBJ) $(LIB)

%.o: %.cxx
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

$(LIB): $(OBJ)
	$(AR) $(ARFLAGS) $(LIB) $(OBJ)

clean:
	rm -rf *.o *~ 
realclean: clean
	rm -rf $(LIB)
tags:
	etags $(SRC)
depend:
	makedepend -- $(INCLUDES) $(CFLAGS) -- $(SRC)


# DO NOT DELETE
