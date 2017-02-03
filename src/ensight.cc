#include "ensight.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

void destroy_ensight(ENSIGHT ensight)
{
  if(ensight != NULL){
    free(ensight->Sm);
    free(ensight->Sp);
    free(ensight->No);
    free(ensight->Vp);
    free(ensight->Cp);
    free(ensight);
  }
}
