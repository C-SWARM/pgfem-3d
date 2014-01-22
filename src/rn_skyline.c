#include "rn_skyline.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

long rn_skyline(int ncol, int *Ap, int *Ai, int *order, int start)
{
  long sky=0;
  int i,j;

  int *ap, **AA;

  ap = (int*) PGFEM_calloc (ncol, sizeof(int));
  AA = (int**) PGFEM_calloc (ncol, sizeof(int*));
  for(i=0;i<ncol;i++)
    {
      ap[i] = Ap[i+1] - Ap[i];
      AA[i] = (int*) PGFEM_calloc (ap[i], sizeof(int));
    }

  /* Copy re-numbered values of Ai into AA */
  for(i=0;i<ncol;i++){
    for(j=0;j<ap[i];j++)
      AA[i][j] = order[ Ai[ j+Ap[i] ] ];
  }

  /* Sort the column indicies */
  for(i=0;i<ncol;i++)
    qsort(AA[i],ap[i],sizeof(int),compare_int);

  /* Determine the skyline */
  for(i=0;i<ncol;i++)
    sky += order[i+start] - AA[i][0] + 1;

  return sky;
}
