#include "hommat.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

HOMMAT* build_hommat (long i)
     /*
       i - pocet homogennich materialu
     */
{
  HOMMAT *pom;
  long ii;
  
  pom = (HOMMAT*) PGFEM_calloc (i, sizeof(HOMMAT));
  
  for (ii=0;ii<i;ii++){	 
    pom[ii].M = (double*) PGFEM_calloc (9,sizeof(double));
    pom[ii].L = (double*) PGFEM_calloc (9,sizeof(double));
  }
  if (pom == NULL){
    PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    abort ();
  }

  return (pom);
}

void destroy_hommat(HOMMAT* hm, long nm)
{
  for(long i=0; i<nm; i++){
    free(hm[i].M);
    free(hm[i].L);
  }
  free(hm);
}
