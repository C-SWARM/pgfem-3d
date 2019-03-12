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

  if (nln == 0) pom = PGFEM_calloc (ZATNODE, 1);
  else          pom = PGFEM_calloc (ZATNODE, nln);

  for (i=0;i<nln;i++){
    pom[i].load = PGFEM_calloc (double, ndofn);
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
    PGFEM_free(zn[i].load);
  }
  PGFEM_free(zn);
}

ZATELEM* build_zatelem (long ndofn,long nle)
     /*
       Fx, Fy, Fz
     */
{
  ZATELEM *pom;
  long i;

  if (nle == 0) pom = PGFEM_calloc (ZATELEM, 1);
  else          pom = PGFEM_calloc (ZATELEM, nle);

  for (i=0;i<nle;i++){
    pom[i].load = PGFEM_calloc (double, ndofn);
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
    PGFEM_free(ze[i].load);
  }
  PGFEM_free(ze);
}
