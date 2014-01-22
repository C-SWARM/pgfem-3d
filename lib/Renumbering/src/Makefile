# get paths
include ../../../Makefile.path

libdir		= $(RENUMBERING_dir)/lib

CC		= mpicc
CFLAGS		= -m64 -O3 -ansi

AR		= ar
ARFLAGS		= -rcs

INCLUDES	= $(RENUMBERING_incl) $(PARMETIS_incl)

OBJ		= check_num_procs.o remove_diagonal.o renumber.o\
		  switch_domains.o renumber_utils.o domain_container.o\
		  comm_profiler.o

LIB		= $(libdir)/librenumbering.a

#########################################
#########################################

vpath %.c $(srcdir)

all: lib $(OBJ)

lib: $(OBJ)
	$(AR) $(ARFLAGS) $(LIB) $(OBJ)

%.o : %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $<


clean:
	rm -rf $(LIB) *.o *~
