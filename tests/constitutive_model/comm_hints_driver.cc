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
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace pgfem3d;

static std::vector<int> send = { 1, 2, 3 };
static std::vector<int> recv = { 4, 5, 6, 7, 8 };

static void print_test_comm_hints_file(std::string fn)
{
  const char* cfn = fn.c_str();
  if (FILE *out = fopen(cfn, "w")) {
    int nsend = send.size();
    fprintf(out, "%d\t", nsend);
    for (auto s : send) {
      fprintf(out, "%d\t", s);
    }
    fprintf(out, "\n");

    int nrecv = recv.size();
    fprintf(out, "%d\t", nrecv);
    for (auto r : recv) {
      fprintf(out, "%d\t", r);
    }
    fprintf(out, "\n");
    fclose(out);
  }
  else {
    std::cerr << "Failed to open file " << fn << "\n";
    abort();
  }
}

int main(int argc, char **argv)
{
  auto fn = CommHints::BuildFilename("./", "test", 0);
  print_test_comm_hints_file(fn);


  CommHints hints("./", "test", 0);

  if (hints.nSend() != send.size())  {
    abort();
  }

  if (hints.nRecv() != recv.size()) {
    abort();
  }

  auto& sends = hints.sends();
  for (size_t i = 0, e = hints.nSend(); i < e; ++i) {
    if (sends[i] != send[i]) { abort(); }
  }

  auto& recvs = hints.recvs();
  for (size_t i = 0, e = hints.nRecv(); i < e; ++i) {
    if (recvs[i] != recv[i]) { abort(); }
  }

  std::cout << hints.nSend() << " ";
  for (auto s : hints.sends()) {
    std::cout << s << " ";
  }
  std::cout << "\n";

  std::cout << hints.nRecv() << " ";
  for (auto r : hints.recvs()) {
    std::cout << r << " ";
  }
  std::cout << "\n";

  try {
    CommHints hints("", "", 0);
    return 1;
  }
  catch (std::system_error&) {
    return 0;
  }
}
