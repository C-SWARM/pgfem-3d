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
  int *displ = NULL;
  int *hu2 = NULL;
  int m = 0;
  int *BNint = NULL;
  long *BN = NULL;
  long *hu1 = NULL;
  long *ID = NULL;
  long *GNn = NULL;
  long i = 0;
  long j = 0;
  long k = 0;
  long l = 0;
  long GN = 0;
  long NBN = 0;
  long *GID = NULL;
  long need=0;
  long *hu3 = NULL;
  long *LoNo = NULL;
  long *Ipom = NULL;

  for (i=0;i<nn;i++){
    if (node[i].Gnn >= 0 && node[i].Dom == myrank && node[i].Pr == -1) GN++;
    if (node[i].Gnn >= 0 && node[i].Dom != myrank && node[i].Pr == -1) need++;
  }

  if(GN > 0){
    hu1 = aloc1l (GN);
    ID = aloc1l (GN*ndofn);
  }

  if(need > 0){
    hu3 = aloc1l (need);
    LoNo = aloc1l (need);
  }

  for (i=0;i<nn;i++){
    if (node[i].Gnn >= 0 
	&& node[i].Dom == myrank 
	&& node[i].Pr == -1){
      hu1[k] = node[i].Gnn; /* Array of global node nnumbers on domain */
      for (j=0;j<ndofn;j++){
	ID[k*ndofn+j] = node[i].Gid[j];
	/* ID contains the Gid numbers for the global nodes on the domain*/
      }
      k++;
    }
    if (node[i].Gnn >= 0 
	&& node[i].Dom != myrank 
	&& node[i].Pr == -1){
      hu3[l] = node[i].Gnn;
      LoNo[l] = i;
      /* hu3 contains the global node numbers on the other
	 domains. LoNo contains the local node numbers of the global
	 nodes on other domains. */
      l++;
    }
  }/* end i < nn */

  /* Gather number of domain boundary nodes owned by each domain */
  BN = aloc1l (nproc);
  BNint = aloc1i(nproc);
  displ = aloc1i (nproc);
  MPI_Allgather (&GN,1,MPI_LONG,BN,1,MPI_LONG,Comm);

  for (i=0;i<nproc;i++){
    BNint[i] = BN[i];
    NBN += BN[i];
    if(i > 0){
      displ[i] = displ[i-1] + BN[i-1];
    }
  }

  
  /* Gather nodes on boundaries between domains */
  if(NBN > 0){
    GNn = aloc1l (NBN);
    GID = aloc1l (NBN*ndofn);
  }
  MPI_Allgatherv (hu1,GN,MPI_LONG,GNn,BNint,displ,MPI_LONG,Comm);

  hu2 = aloc1i (nproc);
  for (i=0;i<nproc;i++) {
    hu2[i] = BN[i]*ndofn;
    displ[i] *= ndofn;
  }
  /* hu2 contains the number of dofs on each domain from the boundary nodes */

  /* Gather global ID numbers for nodes on boundaries between domains */
  MPI_Allgatherv (ID,GN*ndofn,MPI_LONG,GID,hu2,displ,MPI_LONG,Comm);
  Ipom = aloc1l (ndofn);


  // Need to sort GNn with key.  Copy into val_key structure and sort,
  // then copy back.
  val_key *new_GNn = NULL;
  long *new_GID = NULL;
  if(NBN > 0){
    new_GNn = (val_key*) PGFEM_calloc(NBN,sizeof(val_key));
    new_GID = aloc1l (NBN*ndofn);
  }

  // copy GNn into sort container
  for (j=0; j<NBN; j++) {
    new_GNn[j].val = GNn[j];
    new_GNn[j].key = j;

    // Make copy of original GID
    for (k=0; k<ndofn; k++)
      new_GID[j*ndofn+k] = GID[j*ndofn+k];
  }

  // sort GNn container w/ key
  qsort(new_GNn,NBN,sizeof(val_key),compare_val_w_key);
  
  for (i=0; i<NBN; i++){
    for(k=0; k<ndofn; k++){
      GID[i*ndofn + k] = new_GID[new_GNn[i].key*ndofn+k];
    }
    GNn[i] = new_GNn[i].val;
  }

  if (PFEM_DEBUG){
    /* Check global node numbers */
    k = 0;
    for (i=0;i<NBN;i++){
      if (GNn[i] != k){
	PGFEM_printf ("Error in global node numbers\n");
	PGFEM_Comm_code_abort (Comm,m);
      } else k++;
    }
  }

  /* RENUMBER GLOBAL ID ON DOMAINS */
  for (i=0;i<need;i++){
    for (j=0;j<ndofn;j++){
      if(node[LoNo[i]].id[j] <= 0){ /* ADDED 1/29/2013 MM */
	/* BC overrides periodicity */
      	node[LoNo[i]].Gid[j] = node[LoNo[i]].id[j];
      } else {
	node[LoNo[i]].Gid[j] = GID[hu3[i]*ndofn+j];
      }
      /* if(myrank==0)
	 PGFEM_printf("node[%d].Gid[%d]::%d\n",LoNo[i],j,GID[hu3[i]*ndofn+j]-1); */
    }
  }
  
  //  if (periodic == 1){
  for (i=0;i<nn;i++){
    if (node[i].Pr == -1) continue;
    for (j=0;j<ndofn;j++){
      if(node[i].id[j] <= 0){ /* ADDED 1/29/2013 MM*/
	/* BC overrides periodicity */
	node[i].Gid[j] = node[i].id[j];
      } else {
	node[i].Gid[j] = node[node[i].Pr].Gid[j];
      }
    }
  }
  //  }
  
  if (PFEM_DEBUG){  /* PRINT */
    
    PGFEM_printf ("Give || Process [%d]\n",myrank);
    for (i=0;i<GN;i++){
      PGFEM_printf ("GNn = %ld || ID =",hu1[i]);
      for (j=0;j<ndofn;j++)  PGFEM_printf ("  %ld",ID[i*ndofn+j]);
      PGFEM_printf ("\n");
    }
    
    PGFEM_printf ("\n");
    PGFEM_printf ("Need || Process [%d]\n",myrank);
    for (i=0;i<need;i++) PGFEM_printf ("GNn %ld - LNn %ld\n",hu3[i],LoNo[i]);
    
    if (myrank == 0){
      PGFEM_printf ("\n");
      PGFEM_printf ("BN :  ");
      for (i=0;i<nproc;i++)
	PGFEM_printf ("%ld  ",BN[i]); PGFEM_printf ("\n");
      PGFEM_printf ("Displ :  ");
      for (i=0;i<nproc;i++)
	PGFEM_printf ("%d  ",displ[i]); PGFEM_printf ("\n");
      PGFEM_printf ("Nodes :  ");
      for (i=0;i<NBN;i++)
	PGFEM_printf ("%ld  ",GNn[i]); PGFEM_printf ("\n");
      
      PGFEM_printf ("\n");
      for (i=0;i<NBN;i++){
	PGFEM_printf ("[%ld] ::  ",GNn[i]);
	for (j=0;j<ndofn;j++) PGFEM_printf ("  %ld",GID[i*ndofn+j]);
	PGFEM_printf ("\n");
      }
    }

    if(myrank == 0){
      for(i=0;i<nn;i++){
	for(j=0;j<ndofn;j++){
	  if(node[i].Gid[j] > 0){
	    PGFEM_printf("node[%ld].Gid[%ld] :: %ld\n",i,j,node[i].Gid[j]);
	  }
	}
      }
    }
  } /* PFEM_DEBUG */
  
  if(DOF_DEBUG){
    char jmeno[50];
    FILE *out;
    sprintf (jmeno,"%s%d.out","./",myrank);
    if ((out = fopen(jmeno,"a")) == NULL ){ /* append (after code_num_node */
      PGFEM_printf("Output file is not possible to open on processor [%d]\n",myrank);
      PGFEM_printf("Check the output file and run program again\n");
      return (0);
    } 
    
    PGFEM_fprintf(out,"\nnode : Gnn || dof[0] dof[1] dof[2] :: Gdof[0] Gdof[1] Gdof[2]\n");
    for (i=0;i<nn;i++){
      PGFEM_fprintf (out,"%ld : %ld || %ld %ld %ld :: %ld %ld %ld\n",
	       i,node[i].Gnn,node[i].id[0],node[i].id[1],
	       node[i].id[2],node[i].Gid[0],node[i].Gid[1],
	       node[i].Gid[2]);
    }
    
    fclose (out);
  }/* end DOF_DEBUG */

  dealoc1l (BN);
  dealoc1i (BNint);
  dealoc1i (displ);
  dealoc1i (hu2);
  dealoc1l (hu1);
  dealoc1l (ID);
  dealoc1l (GNn);
  dealoc1l (GID);
  dealoc1l (hu3);
  dealoc1l (LoNo);
  dealoc1l (Ipom);
  free(new_GNn);
  free(new_GID);
  
  return (NBN);
}
