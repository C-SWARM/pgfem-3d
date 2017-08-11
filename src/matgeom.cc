#include "matgeom.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif


MATGEOM build_matgeom (long nc,long np)
     /*
       nc - number of volume fractions
       np - number of psi angle
     */
{
  MATGEOM pom;
  long i;

  pom = PGFEM_calloc (MATGEOM_1, 1);

  pom->cf = PGFEM_calloc (double, nc);
  pom->cd = PGFEM_calloc (double, nc);
  pom->cm = PGFEM_calloc (double, nc);
  pom->ee = PGFEM_calloc (double*, np);

  for (i=0;i<np;i++) {
    pom->ee[i] = PGFEM_calloc (double, 9);
  }

  return (pom);
}

void destroy_matgeom(MATGEOM mg, long np)
{
  for(long i=0; i<np; i++){
    free(mg->ee[i]);
  }
  free(mg->ee);
  free(mg->cf);
  free(mg->cd);
  free(mg->cm);
  free(mg);
}
