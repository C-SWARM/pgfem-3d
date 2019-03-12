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
    PGFEM_free(mg->ee[i]);
  }
  PGFEM_free(mg->ee);
  PGFEM_free(mg->cf);
  PGFEM_free(mg->cd);
  PGFEM_free(mg->cm);
  PGFEM_free(mg);
}
