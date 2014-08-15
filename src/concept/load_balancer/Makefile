# This is a makefile to compile the load balancer concept code

CC = gcc
CFLAGS = -Wall -g -O0
CSTD =  -std=c99
PROGRAM = load_bal

SRC = load.c partition.c stats.c load_balancer_utils.c load_balancer.c lb_transfer.c main.c
LIB = -lm

# should not need to modify below
OBJ = $(SRC:.c=.o)
depfile = Makefile.depend

all: $(OBJ) $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(CC) -o $@ $^ $(LIB)

%.o : %.c
	$(CC) $(CFLAGS) $(CSTD) -c $<

main.o : main.c # cannot compile with -std=c99 becasue of getopt
	$(CC) $(CFLAGS) -std=gnu99 -c $<
clean:
	rm -rf *.o *~ $(depfile) TAGS

realclean: clean
	rm -rf $(PROGRAM)

depend: $(depfile)
	makedepend -f $(depfile) $(SRC)
	rm $(depfile).bak

$(depfile):
	touch $(depfile)

tags:
	etags $(SRC)

# DO NOT DELETE
-include $(depfile)
