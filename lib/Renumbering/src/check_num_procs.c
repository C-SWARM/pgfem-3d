#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "check_num_procs.h"

void check_num_procs(MPI_Comm comm)
/* Function checs to see that nproc is a power of two and returns an
   error message and exits if it is not. */
{
  double result;
  int nproc,result2;

  MPI_Comm_size(comm,&nproc);

  result = log((double) nproc) / log(2);
  result2 = (int) result;

  if(result2 != result)
    {
      printf("ERROR::Number of processors must be a power of 2!\nExiting program.\n");
      MPI_Abort(comm,2);
      exit(2);
    }
}
