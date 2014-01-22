#include "macroscopic_load_AL.h"

void macroscopic_load_AL (long TYPE,double lam,EPS *eps)
{
  /* Plane strain deviatoric tension */
  if (TYPE == 1){
    eps[0].F[0][0] = 1./(1. - (lam)*eps[0].load);
    eps[0].F[1][1] = 1. - (lam)*eps[0].load;
    eps[0].F[2][2] = 1.;
  }
  /* Plane strain simple shear */
  if (TYPE == 2){
    eps[0].F[0][0] = 1.;
    eps[0].F[1][1] = 1.;
    eps[0].F[2][2] = 1.;
    eps[0].F[0][1] = (lam)*eps[0].load;
  }
  /* Deviatoric tension */
  if (TYPE == 3){
    eps[0].F[0][0] = 1./((1. - (lam)*eps[0].load)*(1. - (lam)*eps[0].load));
    eps[0].F[1][1] = 1. - (lam)*eps[0].load1;
    eps[0].F[2][2] = 1. - (lam)*eps[0].load1;
  }
  /* tension in x3 + shear x1x2 */
  if (TYPE == 4){
    eps[0].F[0][0] = 1.;
    eps[0].F[1][1] = 1.;
    eps[0].F[2][2] = 1. + (lam)*eps[0].load;
    eps[0].F[0][1] = (lam)*eps[0].load1;
  }
}
