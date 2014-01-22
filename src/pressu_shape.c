#include "pressu_shape.h"
#include <stdlib.h>

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

void pressu_shape (long npres,double ksi,double eta,double zet,double *Psi)
{
  if (npres == 1){
    Psi[0] = 1.0;
  }
  else{
    PGFEM_printf ("Incorrect number of pressure nodes\n"); abort();
  }
}
