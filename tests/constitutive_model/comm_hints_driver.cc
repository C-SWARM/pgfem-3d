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

#include "pgfem3d/Communication.hpp"
#include <stdio.h>
#include <stdlib.h>

using namespace pgfem3d;

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
  CommHints *hints;
  hints = new CommHints("./", "test", 0);
  char *filename = hints->get_filename();
  print_test_comm_hints_file(filename);

  delete hints;
  
  hints = new CommHints();
  try {
    hints->read_filename(filename);
    hints->get_nsend();
    hints->get_nrecv();
    hints->get_send_list();
    hints->get_recv_list();
    hints->write(stdout);
  } catch (int e) {
    return e;
  }

  delete hints;
  
  return 0;
}
