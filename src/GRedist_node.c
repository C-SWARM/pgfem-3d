#include "GRedist_node.h"
#include <stdlib.h>

#include "PGFEM_io.h"
#include "allocation.h"
#include "incl.h"
#include "matice.h"

#ifndef DOF_DEBUG
#define DOF_DEBUG 0
#endif

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

long GRedist_node (const int nproc,
		   const int myrank,
		   const long nn,
		   const long ndofn,
		   NODE *node,
		   const MPI_Comm Comm)
{
  long NBN = 0; /* return value */

  long GN = 0;
  long need = 0;
  for (int i=0;i<nn;i++){
    /* Number of global nodes OWNED by the domain */
    if (node[i].Gnn >= 0 && node[i].Dom == myrank && node[i].Pr == -1) GN++;

    /* Number of global nodes NEEDED by the domain */
    if (node[i].Gnn >= 0 && node[i].Dom != myrank && node[i].Pr == -1) need++;
  }

  /* allocate space for OWNED node/dof ids */
  long *hu1 = NULL;
  long *ID = NULL;
  if(GN > 0){
    hu1 = aloc1l (GN);
    ID = aloc1l (GN*ndofn);
  }

  /* allocate space for NEEDED node/dof ids */
  long *hu3 = NULL;
  long *LoNo = NULL;
  if(need > 0){
    hu3 = aloc1l (need);
    LoNo = aloc1l (need);
  }

  /* Get information to send/recieve */
  {
    int l = 0;
    int k = 0;
    for (int i=0;i<nn;i++){
      if(node[i].Gnn >= 0){
	/* global node owned by THIS domain */
	if (node[i].Dom == myrank && node[i].Pr == -1){
	  /* Array of global node nnumbers on domain */
	  hu1[k] = node[i].Gnn;
	  for (int j=0;j<ndofn;j++){
	    /* ID contains the Gid numbers for the global nodes on the
	       domain*/
	    ID[k*ndofn+j] = node[i].Gid[j];
	  }
	  k++;
	}
	/* global node owned by OTHER domain */
	else if (node[i].Dom != myrank && node[i].Pr == -1){
	  /* list of global node numbers needed from other domains */
	  hu3[l] = node[i].Gnn;
	  /* list of associated local node numbers */
	  LoNo[l] = i;
	  l++;
	}
      }
      /* do nothing for purely local nodes (Gnn < 0) */
    }/* end i < nn */
  }

  /* Gather number of boundary nodes owned by each domain */
  long *BN = aloc1l (nproc);
  int *BNint = aloc1i(nproc);
  int *displ = aloc1i (nproc);
  MPI_Allgather (&GN,1,MPI_LONG,BN,1,MPI_LONG,Comm);

  displ[0] = 0;
  for (int i=0;i<nproc;i++){
    BNint[i] = BN[i];
    NBN += BN[i];
    if(i > 0){
      displ[i] = displ[i-1] + BN[i-1];
    }
  }
  
  /* Gather nodes on boundaries between domains */
  long *GNn = NULL;
  long *GID = NULL;
  if(NBN > 0){
    GNn = aloc1l (NBN);
    GID = aloc1l (NBN*ndofn);
  }
  MPI_Allgatherv (hu1,GN,MPI_LONG,GNn,BNint,displ,MPI_LONG,Comm);
  free (hu1);
  free(BNint);

  /* hu2 contains the number of dofs on each domain from the boundary
     nodes */
  int *hu2 = aloc1i (nproc);
  for (int i=0;i<nproc;i++) {
    hu2[i] = BN[i]*ndofn;
    displ[i] *= ndofn;
  }

  /* Gather global ID numbers for nodes on boundaries between domains */
  MPI_Allgatherv (ID,GN*ndofn,MPI_LONG,GID,hu2,displ,MPI_LONG,Comm);
  free(ID);
  free(hu2);
  free(displ);
  free (BN);

  // Need to sort GNn with key.  Copy into val_key structure and sort,
  // then copy back.
  val_key *new_GNn = NULL;
  long *new_GID = NULL;
  if(NBN > 0){
    new_GNn = PGFEM_calloc(NBN,sizeof(*new_GNn));
    new_GID = aloc1l (NBN*ndofn);
  }

  // copy GNn into sort container
  for (int j=0; j<NBN; j++) {
    new_GNn[j].val = GNn[j];
    new_GNn[j].key = j;

    // Make copy of original GID
    for (int k=0; k<ndofn; k++)
      new_GID[j*ndofn+k] = GID[j*ndofn+k];
  }

  // sort GNn container w/ key
  qsort(new_GNn,NBN,sizeof(val_key),compare_val_w_key);
  
  for (int i=0; i<NBN; i++){
    for(int k=0; k<ndofn; k++){
      GID[i*ndofn + k] = new_GID[new_GNn[i].key*ndofn+k];
    }
    GNn[i] = new_GNn[i].val;
  }
  free(new_GNn);
  free(new_GID);

  if (PFEM_DEBUG){
    /* Check global node numbers */
    int k = 0;
    for (int i=0;i<NBN;i++){
      if (GNn[i] != k){
	PGFEM_printf ("Error in global node numbers\n");
	PGFEM_Comm_abort (Comm);
      } else k++;
    }
  }
  free(GNn);

  /* RENUMBER GLOBAL ID ON DOMAINS */
  for (int i=0;i<need;i++){
    for (int j=0;j<ndofn;j++){
      if(node[LoNo[i]].id[j] <= 0){
	/* BC overrides periodicity */
      	node[LoNo[i]].Gid[j] = node[LoNo[i]].id[j];
      } else {
	node[LoNo[i]].Gid[j] = GID[hu3[i]*ndofn+j];
      }
    }
  }
  free (GID);
  free (hu3);
  free (LoNo);

  for (int i=0;i<nn;i++){
    if (node[i].Pr == -1) continue;
    for (int j=0;j<ndofn;j++){
      if(node[i].id[j] <= 0){
	/* BC overrides periodicity */
	node[i].Gid[j] = node[i].id[j];
      } else {
	node[i].Gid[j] = node[node[i].Pr].Gid[j];
      }
    }
  }
  
  return (NBN);
}
