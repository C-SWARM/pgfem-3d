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

  pom = (MATGEOM) PGFEM_calloc (1, sizeof(MATGEOM_1));
  
  pom->cf = (double *) PGFEM_calloc (nc,sizeof(double));
  pom->cd = (double *) PGFEM_calloc (nc,sizeof(double));
  pom->cm = (double *) PGFEM_calloc (nc,sizeof(double));
  pom->ee = (double **) PGFEM_calloc (np,sizeof(double*));
  
  for (i=0;i<np;i++) pom->ee[i] = (double *) PGFEM_calloc (9,sizeof(double));

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
