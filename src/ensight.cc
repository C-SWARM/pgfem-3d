#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "ensight.h"

Ensight::~Ensight()
{
  delete [] Sm;
  delete [] Sp;
  delete [] No;
  delete [] Vp;
  delete [] Cp;
}
