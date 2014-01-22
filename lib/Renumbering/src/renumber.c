#include <stdio.h>
#include <stdlib.h>

#include "parmetis.h"
#include "remove_diagonal.h"
#include "renumber.h"

void renumber(int *dist, int *Ap, int *Ai, int *l_order,
	      int *sizes, int *order,MPI_Comm comm)
{
  int *Apd, *Aid, *nsend, l_ncol, start,
      numflag, myrank, nproc, options[3], i;

  /* Initialize parameters */
  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&myrank);
  l_ncol = dist[myrank +1] - dist[myrank];
  start = dist[myrank];
  numflag = 0;

  /* Options::Defined at compile time by macros. */
  if(_MOPT1>0){options[0] = 1; options[1] = _MOPT2; options[2] = _MOPT3;}
  else {options[0] = 0; options[1] = 0; options[2] = 0;}

  /* Initialize and create diagonal free graph */
  Apd = (int*) calloc (l_ncol+1,sizeof(int));
  Aid = (int*) calloc (Ap[l_ncol] - l_ncol,sizeof(int));

  /* Remove Diagonal */
#ifdef RENUM_DEBUG
  if(myrank == 0) printf("Removing diagonal..."); 
#endif /* RENUM_DEBUG */

  remove_diagonal(l_ncol,start,Ap,Ai,Apd,Aid);

#ifdef RENUM_DEBUG
  if(myrank == 0) printf("Diagonal removed.\n"); 
#endif /* RENUM_DEBUG */


  /* ParMETIS Routine */
#ifdef RENUM_DEBUG
  if(myrank == 0) printf("Entering ParMETIS routine..."); 
#endif /* RENUM_DEBUG */

  ParMETIS_V3_NodeND(dist,Apd,Aid,&numflag,options,l_order,sizes,&comm);

#ifdef RENUM_DEBUG
  if(myrank == 0) printf("ParMETIS finished.\n"); 
#endif /* RENUM_DEBUG */


  /* Allgatherv l_order into order */
  nsend = (int*) calloc (nproc,sizeof(int));
  for(i=0;i<nproc;i++)
    nsend[i] = dist[i+1] - dist[i];
  MPI_Allgatherv(l_order,l_ncol,MPI_INT,order,nsend,dist,MPI_INT,comm);

  free(Apd); free(Aid); free(nsend);
}
