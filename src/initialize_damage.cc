#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "initialize_damage.h"
#include "allocation.h"
#include "elem3d.h"
#include "enumerations.h"
#include <cstdlib>

void initialize_damage(const int ne,
                       const Element *elem,
                       const HOMMAT *hommat,
                       EPS *eps,
                       const int analysis)
{
  int n_params = sizeof(DamageParams)/sizeof(double);
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
  PGFEM_free(dam_params);
}
