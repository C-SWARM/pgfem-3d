#include "skyline.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif


long skyline(int ncol, int *Ap, int *Ai, int start)
{
  int i;
  long sky;

  sky = 0;
  for(i=0;i<ncol;i++)
    sky += i+start - Ai[Ap[i]] +1;

/*   PGFEM_printf("skyline = %d\n",sky); */
  return sky;
}
