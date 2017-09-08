#include "initialize_damage.h"
#include <stdlib.h>

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

void initialize_damage(const int ne,
               const ELEMENT *elem,
               const HOMMAT *hommat,
               EPS *eps,
               const int analysis)
{
  int n_params = sizeof(damage_params)/sizeof(double);
  double *dam_params = PGFEM_calloc (double, std::max(n_params, 4));

  for(int i=0; i<ne; i++){/* for each element */
    const int mat = elem[i].mat[2];

    /* Currently get damage parameters from the extra parameters in
       HOMMAT. Eventually, we need new input routines to read damage
       parameters directly from file. */

    dam_params[0] = hommat[mat].e1;
    dam_params[1] = hommat[mat].e2;
    dam_params[2] = hommat[mat].e3;
    dam_params[3] = hommat[mat].e4;

    /* Currently only have 1 damage function/evolution equation, so
       hard set the flags */
    const int eq_flag = 0;

    long n_ip;
    int_point(elem[i].toe,&n_ip);
    for(int j=0; j<n_ip; j++){
      init_damagef(&eps[i].dam[j],eq_flag,dam_params,n_params);
    }

    if(analysis == STABILIZED){
      int_point(10,&n_ip);
      for(int j=0; j<n_ip; j++){
    init_damagef(&eps[i].st[j].dam,eq_flag,dam_params,n_params);
      }
    }

  }
  free(dam_params);
}
