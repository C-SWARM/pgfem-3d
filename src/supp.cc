#include "supp.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

void destroy_supp(SUPP sup)
{
  PGFEM_free(sup->free);
  PGFEM_free(sup->supp);
  PGFEM_free(sup->lnpd);
  PGFEM_free(sup->lepd);
  PGFEM_free(sup->lbepd);
  PGFEM_free(sup->defl);
  PGFEM_free(sup->defl_d);
  PGFEM_free(sup->F0);
  PGFEM_free(sup->N0);
  PGFEM_free(sup);
}
