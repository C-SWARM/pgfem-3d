/* HEADER */
/**
 * Define Comm_hints structure.
 *
 * AUTHORS:
 *  Matthew Mosby, University of Notre Dame <mmosby1@nd.edu>
 */

#include "comm_hints.h"
#include "allocation.h"
#include "utils.h"

struct COMM_HINTS {
  int *send;
  int *recv;
  int nsend;
  int nrecv;
};

Comm_hints* Comm_hints_construct()
{
  Comm_hints *hints = PGFEM_malloc<Comm_hints>();
  hints->send = hints->recv = NULL;
  hints->nsend = hints->nrecv = 0;
  return hints;
}

int Comm_hints_destroy(Comm_hints *hints)
{
  if (hints) {
    free(hints->send);
    free(hints->recv);
    free(hints);
  }
  return 0;
}

int Comm_hints_read(Comm_hints *hints,
                    FILE *in)
{
  int err = 0;

  /* read hints for whom to send information */
  CHECK_SCANF(in, "%d", &hints->nsend);
  hints->send = PGFEM_calloc(int, hints->nsend);
  for (int i = 0, e = hints->nsend; i < e; i++) {
    CHECK_SCANF(in, "%d", hints->send + i);
  }

  /* read hints from whom to receive information */
  CHECK_SCANF(in, "%d", &hints->nrecv);
  hints->recv = PGFEM_calloc(int, hints->nrecv);
  for (int i = 0, e = hints->nrecv; i < e; i++) {
    CHECK_SCANF(in, "%d", hints->recv + i);
  }

  return err;
}

int Comm_hints_read_filename(Comm_hints *hints,
                             const char *fn)
{
  int err = 0;
  FILE *in = fopen(fn, "r");
  if (!in) {
    err++;
  } else {
    err += Comm_hints_read(hints, in);
    fclose(in);
  }
  return err;
}

int Comm_hints_write(const Comm_hints *hints,
                     FILE *out)
{
  int err = 0;
  fprintf(out, "%d\t", hints->nsend);
  for (int i = 0, e = hints->nsend; i < e; i++) {
    fprintf(out, "%3.0d ", hints->send[i]);
  }
  fprintf(out, "\n");
  fprintf(out, "%d\t", hints->nrecv);
  for (int i = 0, e = hints->nrecv; i < e; i++) {
    fprintf(out, "%3.0d ", hints->recv[i]);
  }
  fprintf(out, "\n");
  return err;
}

int Comm_hints_write_filename(const Comm_hints *hints,
                              const char *fn)
{
  int err = 0;
  FILE *out = fopen(fn, "w");
  if (!out) {
    err++;
  } else {
    err += Comm_hints_write(hints, out);
    fclose(out);
  }
  return err;
}


int Comm_hints_nsend(const Comm_hints *hints) { return hints->nsend; }
int Comm_hints_nrecv(const Comm_hints *hints) { return hints->nrecv; }
const int* Comm_hints_send_list(const Comm_hints *hints) { return hints->send; }
const int* Comm_hints_recv_list(const Comm_hints *hints) { return hints->recv; }

char* Comm_hints_filename(const char *ipath,
                          const char *base_filename,
                          const int rank)
{
  char *str = NULL;
  const char format[] = "%s/%scomm_hints_%d.in";
  int nchar = snprintf(str, 0, format, ipath, base_filename, rank);
  str = PGFEM_malloc<char>(nchar + 1);
  sprintf(str, format, ipath, base_filename, rank);
  return str;
}
