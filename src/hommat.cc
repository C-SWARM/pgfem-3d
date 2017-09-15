#include "hommat.h"
#include "PGFEM_io.h"
#include "allocation.h"
#include "constitutive_model.h"

HOMMAT* build_hommat (long i)
     /*
       i - number of homogeneous materials
     */
{
  HOMMAT *pom;
  long ii;

  pom = PGFEM_calloc (HOMMAT, i);

  for (ii=0;ii<i;ii++){
    /* internal allocation */
    pom[ii].M = PGFEM_calloc (double, 9);
    pom[ii].L = PGFEM_calloc (double, 9);

    /* initialize variables */
    pom[ii].m10 = 0.0;
    pom[ii].m01 = 0.0;
    pom[ii].E = 0.0;
    pom[ii].G = 0.0;
    pom[ii].nu = 0.0;
    pom[ii].density = 0.0;
    pom[ii].e1 = 0.0;
    pom[ii].e2 = 0.0;
    pom[ii].e3 = 0.0;
    pom[ii].e4 = 0.0;
    pom[ii].devPotFlag = -1; /* poisoned */
    pom[ii].volPotFlag = -1; /* poisoned */
    pom[ii].mat_id = -1;
    pom[ii].param = NULL;
  }

  return (pom);
}

void destroy_hommat(HOMMAT* hm, long nm)
{
  for(long i=0; i<nm; i++){
    free(hm[i].M);
    free(hm[i].L);

    if(hm[i].param==NULL)
      continue;
      
    for(int ia = 0; ia < nm; ia++)
    {
      delete hm[i].param;
      hm[i].param = NULL;
    }    
  }
  free(hm);
}

double hommat_get_kappa(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}
