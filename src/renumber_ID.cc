#include "renumber_ID.h"

void renumber_ID (int ndofn,int nn,NODE *node,int *g_order,MPI_Comm comm)
{
  long i,k;
  int myrank;
  MPI_Comm_rank(comm,&myrank);

  for (i=0;i<nn;i++){
    if(node[i].Dom != myrank) continue;
    for (k=0;k<ndofn;k++){
      if (node[i].Gid[k] <= 0) continue; /* Fixed or prescribed deflection */
      else{
	node[i].Gid[k] = g_order[ node[i].Gid[k] -1 ] +1;
      }
    }
  }
}
