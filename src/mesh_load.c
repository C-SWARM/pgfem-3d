#include "mesh_load.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

ZATNODE* build_zatnode (long ndofn,long nln)
     /*
       Fx, Fy, Fz
     */
{
  ZATNODE *pom;
  long i;
  
  if (nln == 0) pom = (ZATNODE*) PGFEM_calloc (1, sizeof(ZATNODE));
  else          pom = (ZATNODE*) PGFEM_calloc (nln, sizeof(ZATNODE)); 
  
  for (i=0;i<nln;i++){  
    pom[i].load = (double *) PGFEM_calloc (ndofn,sizeof(double));
  }
  
  if (pom==NULL){
    PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    abort ();
  }
  
  return (pom);
}

void destroy_zatnode(ZATNODE* zn, long nln)
{
  for(long i=0; i<nln; i++){
    free(zn[i].load);
  }
  free(zn);
}

ZATELEM* build_zatelem (long ndofn,long nle)
     /*
       Fx, Fy, Fz
     */
{
  ZATELEM *pom;
  long i;
  
  if (nle == 0) pom = (ZATELEM*) PGFEM_calloc (1, sizeof(ZATELEM));
  else          pom = (ZATELEM*) PGFEM_calloc (nle, sizeof(ZATELEM)); 
  
  for (i=0;i<nle;i++){  
    pom[i].load = (double *) PGFEM_calloc (ndofn,sizeof(double));
  }
  
  if (pom==NULL){
    PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    abort ();
  }

  return (pom);
}

void destroy_zatelem(ZATELEM* ze, long nle)
{
  for(long i=0; i<nle; i++){
    free(ze[i].load);
  }
  free(ze);
}
