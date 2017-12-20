#include "build_distribution.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

/* Create nodal distribution array */
void build_distribution(long *DomDof, int *Dist, int nproc)
{
  int i;

  Dist[0] = 0;
  for(i=1;i<=nproc;i++) Dist[i] = DomDof[i-1] + Dist[i-1];

  /*****************************************/
  /*          FUNCTION TESTING             */
  /*****************************************/

/*   MPI_Comm_rank(comm,&myrank); */
/*   if(myrank == 0) */
/*     { */
/*       for(i=0; i<=nproc; i++) */
/* 	PGFEM_printf("%d\n",Dist[i]); */
/*     } */
}
