# This is a makefile to compile the load balancer concept code

CC = gcc
CFLAGS = -Wall -std=c99 -g -O0
PROGRAM = load_bal

SRC = load.c partition.c stats.c load_balancer_utils.c load_balancer.c

#main.c load_balancer.c load_list.c load_balancer_utils.c load.c partition.c
LIB = -lm

# should not need to modify below
OBJ = $(SRC:.c=.o)
depfile = Makefile.depend

all: $(OBJ) #$(PROGRAM)

$(PROGRAM): $(OBJ)
	$(CC) -o $@ $^ $(LIB)

%.o : %.c
	$(CC) $(CFLAGS) -c $<

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
