#include "renumber_ID.h"

void renumber_ID (int ndofn,int nn,NODE *node,int *g_order,MPI_Comm comm, const int mp_id)
{
  long i,k;
  int myrank;
  MPI_Comm_rank(comm,&myrank);

  for (i=0;i<nn;i++){
    if(node[i].Dom != myrank) continue;
    for (k=0;k<ndofn;k++){
      if (node[i].id_map[mp_id].Gid[k] <= 0) continue; /* Fixed or prescribed deflection */
      else{
	node[i].id_map[mp_id].Gid[k] = g_order[ node[i].id_map[mp_id].Gid[k] -1 ] +1;
      }
    }
  }
}
