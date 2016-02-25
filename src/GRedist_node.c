#include "GRedist_node.h"
#include <stdlib.h>
#include <string.h>

#include "PGFEM_io.h"
#include "allocation.h"
#include "matice.h"

#ifndef NDEBUG
#define PFEM_DEBUG 1
#else
#define PFEM_DEBUG 0
#endif

/**
 * Original implementation of GRedist_node
 */
static long fallback_GRedist_node(const int nproc,
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

  /* allocate buffers which contian information in the following order:
     own_buf[i] = {Gnn,Gid[1],...,Gid[n]}[i].
     need_buf[i] = {Gnn,node_idx}[i].

     'own_buf' is condensed across all buffers and sorted by
     Gnn. 'need_buf' is local and is sorted by Gnn as well. After
     sorting, dofs are assigned to the local 'need' nodes.
  */
  long *own_buf = NULL;
  long *need_buf = NULL;
  size_t own_buf_elem_size = (ndofn+1)*sizeof(long);
  size_t need_buf_elem_size = 2*sizeof(long);
  if(GN > 0) own_buf = malloc(GN*own_buf_elem_size);
  if(need > 0) need_buf = malloc(need*need_buf_elem_size);

  {
    size_t own_idx = 0;
    size_t need_idx = 0;
    for (size_t i=0;i<nn;i++){
      if(node[i].Gnn >= 0){
	/* global node owned by THIS domain */
	if (node[i].Dom == myrank && node[i].Pr == -1){
	  own_buf[own_idx++] = node[i].Gnn;
	  memcpy(own_buf + own_idx,node[i].Gid,ndofn*sizeof(long));
	  own_idx += ndofn;
	}
	/* global node owned by OTHER domain */
	else if (node[i].Dom != myrank && node[i].Pr == -1){
	  need_buf[need_idx++] = node[i].Gnn;
	  need_buf[need_idx++] = i;
	}
      }
      /* do nothing for purely local nodes (Gnn < 0) */
    }/* end i < nn */
  }

  /* Gather number of boundary nodes owned by each domain */
  long *BN = aloc1l (nproc);
  MPI_Allgather (&GN,1,MPI_LONG,BN,1,MPI_LONG,Comm);

  /* Compute recvcount and displ arrays */
  int *recvcount = PGFEM_calloc(nproc,sizeof(int));
  int *displ = PGFEM_calloc(nproc,sizeof(int));
  displ[0] = 0;
  for (size_t i=0;i<nproc;i++){
    recvcount[i] = BN[i]*(ndofn + 1);
    NBN += BN[i];
    if(i > 0){
      displ[i] = displ[i-1] + recvcount[i-1];
    }
  }

  /* gather the global nodes and their associated global dof ids */
  long *Gnn_Gid = NULL;
  if(NBN > 0) Gnn_Gid = malloc(NBN*own_buf_elem_size);
  MPI_Allgatherv(own_buf,recvcount[myrank],MPI_LONG,
		 Gnn_Gid,recvcount,displ,MPI_LONG,Comm);

  /* sort the list of Gnn and their associated Gid by Gnn */
  qsort(Gnn_Gid,NBN,own_buf_elem_size,compare_long);

  /* sort need_buf by Gnn */
  qsort(need_buf,need,need_buf_elem_size,compare_long);

  if (PFEM_DEBUG){
    /* Check global node numbers, should be ordered and contiguous */
    size_t k = 0;
    for (size_t i=0; i<NBN*(ndofn+1); i+=ndofn+1){
      if (Gnn_Gid[i] != k){
	PGFEM_printf ("Error in global node numbers (%ld)\n",k);
	//PGFEM_Comm_abort (Comm);
      } else k++;
    }
  }

  /* RENUMBER GLOBAL ID ON DOMAINS */
  for (size_t i=0;i<need;i++){
    const size_t Gnn = need_buf[i*2];
    const size_t nod = need_buf[i*2+1];
    const size_t Gnn_Gid_idx = Gnn*(ndofn+1) + 1;
    for (size_t j=0;j<ndofn;j++){
      if(node[nod].id[j] <= 0){
	/* BC overrides periodicity */
      	node[nod].Gid[j] = node[nod].id[j];
      } else {
	node[nod].Gid[j] = Gnn_Gid[Gnn_Gid_idx + j];
      }
    }
  }

  for (size_t i=0;i<nn;i++){
    if (node[i].Pr == -1) continue;
    for (size_t j=0;j<ndofn;j++){
      if(node[i].id[j] <= 0){
	/* BC overrides periodicity */
	node[i].Gid[j] = node[i].id[j];
      } else {
	node[i].Gid[j] = node[node[i].Pr].Gid[j];
      }
    }
  }

  free(BN);
  free(Gnn_Gid);
  free(need_buf);
  free(own_buf);
  free(recvcount);
  free(displ);

  return (NBN);
}

/**
 * New implementation of GRedist_node that avoids global communication
 *
 * \return total number of global nodes
 */
static long comm_hints_GRedist_node(const int nproc,
                                    const int myrank,
                                    const long nnode,
                                    const long ndofn,
                                    NODE *nodes,
                                    const Comm_hints *hints,
                                    const MPI_Comm comm)
{
  long owned_gnn = 0;
  long total_gnn = -1;

  return total_gnn;
}

long GRedist_node (const int nproc,
		   const int myrank,
		   const long nn,
		   const long ndofn,
		   NODE *node,
                   const Comm_hints *hints,
		   const MPI_Comm Comm)
{
  if (!hints) return fallback_GRedist_node(nproc, myrank, nn, ndofn, node, Comm);
  else return comm_hints_GRedist_node(nproc, myrank, nn, ndofn, node, hints, Comm);
}
