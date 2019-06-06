#include "mesh_load.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#include "pgfem3d/Communication.hpp"

using pgfem3d::PGFEM_Abort;

ZATNODE* build_zatnode(long ndofn, long nln) {
  ZATNODE *pom;

  if (nln == 0) {
    pom = PGFEM_calloc (ZATNODE, 1);
  }
  else {
    pom = PGFEM_calloc (ZATNODE, nln);
  }

  if (pom == nullptr) {
    PGFEM_printf("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    PGFEM_Abort();
  }

  for (long i = 0; i < nln; ++i) {
    pom[i].load = PGFEM_calloc (double, ndofn);
  }

  return pom;
}

void destroy_zatnode(ZATNODE* zn, long nln) {
  for(long i = 0; i < nln; ++i){
    PGFEM_free(zn[i].load);
  }
  PGFEM_free(zn);
}

ZATELEM* build_zatelem (long ndofn,long nle) {
  ZATELEM *pom;

  if (nle == 0) {
    pom = PGFEM_calloc (ZATELEM, 1);
  }
  else {
    pom = PGFEM_calloc (ZATELEM, nle);
  }

  if (pom == nullptr) {
    PGFEM_printf("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    PGFEM_Abort();
  }

  for (long i = 0; i < nle; ++i) {
    pom[i].load = PGFEM_calloc (double, ndofn);
  }

  return pom;
}

void destroy_zatelem(ZATELEM* ze, long nle) {
  for(long i = 0; i < nle; ++i) {
    PGFEM_free(ze[i].load);
  }
  PGFEM_free(ze);
}
