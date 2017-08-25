/* HEADER */
/**
 * Simple program to test the Comm_hints interface.
 *
 * AUTHORS:
 *  Matthew Mosby, University of Notre Dame <mmosby1@nd.edu>
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "comm_hints.h"
#include <stdio.h>
#include <stdlib.h>

static void print_test_comm_hints_file(const char *fn)
{
  FILE *out = fopen(fn, "w");
  if (!out) abort();

  /* print dummy hints */
  fprintf(out, "3\t1  2  3\n");
  fprintf(out, "5\t4  5  6  7  8\n");

  fclose(out);
}

int main(int argc, char **argv)
{
  int err = 0;
  char *filename = Comm_hints_filename("./", "test_", 0);
  print_test_comm_hints_file(filename);

  Comm_hints *hints = Comm_hints_construct();
  err += Comm_hints_read_filename(hints, filename);
  Comm_hints_nsend(hints);
  Comm_hints_nrecv(hints);
  Comm_hints_send_list(hints);
  Comm_hints_recv_list(hints);
  err += Comm_hints_write(hints, stdout);
  err += Comm_hints_destroy(hints);
  err -= remove(filename); /* returns -1 on error */
  free(filename);
  return err;
}
