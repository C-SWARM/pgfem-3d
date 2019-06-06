#include "matgeom.h"
#include <cassert>

MATGEOM build_matgeom (long nc,long np)
     /*
       nc - number of volume fractions
       np - number of psi angle
     */
{
  MATGEOM pom = new MATGEOM_1{};
  pom->cf = new double[nc]{};
  pom->cd = new double[nc]{};
  pom->cm = new double[nc]{};
  pom->ee = new double*[np]{};

  for (long i = 0; i < np; ++i) {
    pom->ee[i] = new double[9]{};
  }

  return pom;
}

void destroy_matgeom(MATGEOM mg, long np)
{
  assert(mg);
  for(long i = 0; i < np; ++i) {
    delete [] mg->ee[i];
  }
  delete [] mg->ee;
  delete [] mg->cf;
  delete [] mg->cd;
  delete [] mg->cm;
  delete [] mg;
}
