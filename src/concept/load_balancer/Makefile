# This is a makefile to compile the load balancer concept code

CC = gcc
CFLAGS = -Wall -Werror -std=c99 -g -O0
PROGRAM = load_bal

SRC = main.c load_balancer.c load_list.c load_balancer_utils.c
LIB = -lm

# should not need to modify below
OBJ = $(SRC:.c=.o)
depfile = Makefile.depend

all: $(PROGRAM)

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
