/*
  SetGlobalNodeNumbers.c || Given a distributed set of nodes and (a)
  the number of nodes on each domain, (b) the number of global
  (intreficial) nodes on each domain, (c) the total number of nodes,
  assign a global number to each node in the set.
*/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "SetGlobalNodeNumbers.h"
#include "PGFEM_io.h"
#include "allocation.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

void SetGlobalNodeNumbers(int nNodesDom, Node *node, CommunicationStructure *com)
{
  int myrank = com->rank;
  int nproc = com->nproc;

  /* iterators */
  int i,j,k;

  /* Determine number of nodes owned on the domain */
  int nOwn, nGlobalOwn, nGlobalNeed, nPeriodic;
  nOwn = nNodesDom;
  nGlobalOwn = nGlobalNeed = 0;
  for(i=0; i<nNodesDom; i++){
    if(node[i].Dom != myrank || node[i].Pr != -1)
      nOwn--;
    if(node[i].Dom != myrank)
      nGlobalNeed++;
    if(node[i].Gnn >= 0 && node[i].Dom == myrank && node[i].Pr == -1)
      nGlobalOwn++;
  }

  nGlobalOwn *= 2;
  nPeriodic = nNodesDom - nOwn + nGlobalNeed;

  int *distNode = PGFEM_calloc (int, nproc+1);
  int *nodeMappingPtr = PGFEM_calloc (int, nproc+1);
  int *nodeMappingCnt = PGFEM_calloc (int, nproc);

  com->net->allgather(&nOwn,1,NET_DT_INT,&distNode[1],1,NET_DT_INT,com->comm);
  com->net->allgather(&nGlobalOwn,1,NET_DT_INT,&nodeMappingPtr[1],1,NET_DT_INT,com->comm);
  
  for(i=1; i<=nproc; i++){
    distNode[i] += distNode[i-1];
    nodeMappingCnt[i-1] = nodeMappingPtr[i];
    nodeMappingPtr[i] += nodeMappingPtr[i-1];
  }

  /*
    Allocate boundaryNodeMapping.  Will store old and new node numbers
    old = boundaryNodeMapping[nodeMappingPtr[rank]+2*i]
    new = boundaryNodeMapping[nodeMappingPtr[rank]+2*i+1]
  */
  int *boundaryNodeMapping = PGFEM_calloc(int, nodeMappingPtr[nproc]);

  /* Allocate quick lookup node container */
  int *boundaryLookup = PGFEM_calloc(int, (nGlobalNeed > 0) ? nGlobalNeed : 1);
  int *periodicLookup = PGFEM_calloc(int, (nPeriodic > 0) ? nPeriodic : 1);
  int boundaryIdx = 0;
  int periodicIdx = 0;

  j = k = 0;
  for(i=0; i<nNodesDom; i++){

    /* Skip periodic and unowned nodes after adding to lookup tables */
    if(node[i].Dom != myrank || node[i].Pr != -1){
      if(node[i].Dom != myrank){
        boundaryLookup[boundaryIdx] = i;
        boundaryIdx++;
      }
      if(node[i].Pr != -1){
        periodicLookup[periodicIdx] = i;
        periodicIdx++;
      }
      continue;
    }

    /* Record old Gnn and new Gnn */
    if(node[i].Dom == myrank && node[i].Gnn >= 0){
      boundaryNodeMapping[nodeMappingPtr[myrank]+2*k] = node[i].Gnn;
      boundaryNodeMapping[nodeMappingPtr[myrank]+2*k+1] = j + distNode[myrank];
      k++;
    }

    /* Set global number  */
    node[i].Gnn = j + distNode[myrank];
    j++;
  }

  com->net->allgatherv(&boundaryNodeMapping[nodeMappingPtr[myrank]],
		       nodeMappingPtr[myrank+1]-nodeMappingPtr[myrank],
		       NET_DT_INT,boundaryNodeMapping,nodeMappingCnt,
		       nodeMappingPtr,NET_DT_INT,com->comm);
  
  /* Loop through the boundary nodes and update global numbers */
  int domain, oldGnn, newGnn{};
  for(i=0; i<nGlobalNeed; i++){
    /* determine domain which owns the node and what the old number was */
    domain = node[boundaryLookup[i]].Dom;
    oldGnn =  node[boundaryLookup[i]].Gnn;

    /* Search for new Gnn */
    for(j=nodeMappingPtr[domain]; j<nodeMappingPtr[domain+1]; j += 2){

      if(boundaryNodeMapping[j] == oldGnn){
        newGnn = boundaryNodeMapping[j+1];
        break;
      }

      if(j == nodeMappingPtr[domain+1] -1){
        PGFEM_printf("ERROR: could not find match for node %d on domain %d owned by domain %d!\n%s",
                     boundaryLookup[i],myrank,domain,"Exiting.\n");
        abort();
      }
    } /* search for newGnn */

    node[boundaryLookup[i]].Gnn = newGnn;
  } /* update Gnn on boundary nodes */

  /* Loop through the periodic nodes and update global numbers */
  for(i=0; i<nPeriodic; i++){
    node[periodicLookup[i]].Gnn = node[ node[periodicLookup[i]].Pr ].Gnn;
  }

  /*
    if(myrank == 0){
    for(i=0; i<nproc; i++){
    PGFEM_printf("Mappings for process %d:\n",i);
    for(j=nodeMappingPtr[i]; j<nodeMappingPtr[i+1]; j+=2){
    PGFEM_printf("%d -> %d\n",boundaryNodeMapping[j],boundaryNodeMapping[j+1]);
    }
    }
    }
  */

  if(myrank == 0){
    PGFEM_printf("Summary of %s:\n",__func__);
    PGFEM_printf("Total number of nodes: %d\n",distNode[nproc]);
    PGFEM_printf("Total number of boundary nodes: %d\n",nodeMappingPtr[nproc]/2);
    PGFEM_printf("These numbers should match those printed above.\n");
  }

  abort();
}
