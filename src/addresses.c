#include "addresses.h"

/** Computes addresses of the diagonal matrix elements. */
void addresses (long *adr,long ndof)
{
  long i,j,k;
  
  j = adr[0]; adr[0] = 0;
  for (i=1;i<=ndof;i++){
    k = adr[i];
    adr[i] = adr[i-1]+j;
    j = k;
  }
}

/** Computes addresses of the diagonal matrix elements. */
void addresses_nonsym (long *adr,long ndof)
{
  long i,j;
  
  adr[0] = 0;
  for (i=1;i<ndof;i++){
    j = adr[i];
    adr[i] = adr[i-1]+j;
  }
}
