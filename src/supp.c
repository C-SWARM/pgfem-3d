#include "supp.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

void destroy_supp(SUPP sup)
{
  free(sup->free);
  free(sup->supp);
  free(sup->lnpd);
  free(sup->lepd);
  free(sup->lbepd);
  free(sup->defl);
  free(sup->defl_d);
  free(sup->F0);
  free(sup->N0);
  free(sup);
}
