#include "remove_diagonal.h"

void remove_diagonal(int ncol, int start, int *Ap, int *Ai, int *Apd, int *Aid)
{
  int i,j,k;

  k = 0;
  Apd[0] = 0;
  for(i=0;i<ncol;i++)
    {
      for(j=Ap[i];j<Ap[i+1];j++)
	{
	  if(Ai[j] == i+start) continue;
	  Aid[k] = Ai[j];
	  k++;
	}
      Apd[i+1] = k;
    }
}
