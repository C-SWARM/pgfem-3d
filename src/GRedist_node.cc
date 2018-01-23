#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "GRedist_node.h"
#include "PGFEM_io.h"
#include "allocation.h"
#include "matice.h"
#include <cstdlib>
#include <cstring>

#ifndef NDEBUG
#define PFEM_DEBUG 1
#else
#define PFEM_DEBUG 0
#endif

#define GREDIST_TAG 999

/**
 * Original implementation of GRedist_node
 */
static long fallback_GRedist_node(const int nproc,
                                  const int myrank,
                                  const long nn,
                                  const long ndofn,
                                  Node *node,
                                  const MPI_Comm Comm,
                                  const long mp_id)
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
  if(GN > 0) own_buf = static_cast<long*>(malloc(GN*own_buf_elem_size));
  if(need > 0) need_buf = static_cast<long*>(malloc(need*need_buf_elem_size));

  {
    size_t own_idx = 0;
    size_t need_idx = 0;
    for (auto i=0;i<nn;i++){
      if(node[i].Gnn >= 0){
        /* global node owned by THIS domain */
        if (node[i].Dom == myrank && node[i].Pr == -1){
          own_buf[own_idx++] = node[i].Gnn;
          memcpy(own_buf + own_idx,node[i].id_map[mp_id].Gid,ndofn*sizeof(long));
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
  int *recvcount = PGFEM_calloc(int,nproc);
  int *displ = PGFEM_calloc(int,nproc);
  displ[0] = 0;
  for (auto i=0;i<nproc;i++){
    recvcount[i] = BN[i]*(ndofn + 1);
    NBN += BN[i];
    if(i > 0){
      displ[i] = displ[i-1] + recvcount[i-1];
    }
  }

  /* gather the global nodes and their associated global dof ids */
  long *Gnn_Gid = NULL;
  if(NBN > 0) Gnn_Gid = static_cast<long*>(malloc(NBN*own_buf_elem_size));
  MPI_Allgatherv(own_buf,recvcount[myrank],MPI_LONG,
                 Gnn_Gid,recvcount,displ,MPI_LONG,Comm);

  /* sort the list of Gnn and their associated Gid by Gnn */
  qsort(Gnn_Gid,NBN,own_buf_elem_size,compare_long);

  /* sort need_buf by Gnn */
  if (need_buf) qsort(need_buf,need,need_buf_elem_size,compare_long);

  if (PFEM_DEBUG){
    /* Check global node numbers, should be ordered and contiguous */
    long k = 0;
    for (auto i=0; i<NBN*(ndofn+1); i+=ndofn+1){
      if (Gnn_Gid[i] != k){
        PGFEM_printf ("Error in global node numbers (%ld)\n",k);
        //PGFEM_Comm_abort (Comm);
      } else k++;
    }
  }

  /* RENUMBER GLOBAL ID ON DOMAINS */
  for (int i=0;i<need;i++){
    const size_t Gnn = need_buf[i*2];
    const size_t nod = need_buf[i*2+1];
    const size_t Gnn_Gid_idx = Gnn*(ndofn+1) + 1;
    for (long j=0;j<ndofn;j++){
      if(node[nod].id_map[mp_id].id[j] <= 0){
        /* BC overrides periodicity */
        node[nod].id_map[mp_id].Gid[j] = node[nod].id_map[mp_id].id[j];
      } else {
        node[nod].id_map[mp_id].Gid[j] = Gnn_Gid[Gnn_Gid_idx + j];
      }
    }
  }

  for (auto i=0;i<nn;i++){
    if (node[i].Pr == -1) continue;
    for (auto j=0;j<ndofn;j++){
      if(node[i].id_map[mp_id].id[j] <= 0){
        /* BC overrides periodicity */
        node[i].id_map[mp_id].Gid[j] = node[i].id_map[mp_id].id[j];
      } else {
        node[i].id_map[mp_id].Gid[j] = node[node[i].Pr].id_map[mp_id].Gid[j];
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
                                    Node *nodes,
                                    const Comm_hints *hints,
                                    const MPI_Comm comm,
                                    const long mp_id)
{
  long owned_gnn = 0;
  long total_gnn = 0;
  int owned_range[2] = {0};

  /* get the reduced list of shared global nodes */
  int n_shared = 0;
  const Node *shared = NULL;
  nodes_filter_shared_nodes(nnode, nodes, &n_shared, &shared);

  /* get the owned Gnn and Gid and pack into a buffer to
     communicate. The buffer stores information in the follwoing
     order:

     owned_Gnn_Gid[i] = {Gnn, Gid[1], ..., Gid[ndofn]}[i]

     This buffer is constructed in Gnn order and does not require
     sorting.
  */
  long *owned_Gnn_Gid = NULL;
  int len_owned_Gnn_Gid = 0;
  if (!nodes_get_shared_idx_range(n_shared, shared, myrank, owned_range)) {
    owned_gnn = owned_range[1] - owned_range[0];
    len_owned_Gnn_Gid = owned_gnn * (ndofn + 1);
    owned_Gnn_Gid = static_cast<long*>(malloc(len_owned_Gnn_Gid * sizeof(*owned_Gnn_Gid)));

    /* loop through list and decrement owned_gnn if duplicates are
       found (e.g., periodic nodes). Additionally construct the
       owned_Gnn_Gid buffer */
    int idx = 0;
    int i, e;
    for (i = owned_range[0], e = owned_range[1] - 1; i < e; i++) {
      owned_Gnn_Gid[idx++] = shared[i].Gnn;
      memcpy(owned_Gnn_Gid + idx, shared[i].id_map[mp_id].Gid, ndofn * sizeof(*owned_Gnn_Gid));
      idx += ndofn;
      if (shared[i].Gnn == shared[i+1].Gnn) owned_gnn--;
    }
    owned_Gnn_Gid[idx++] = shared[i].Gnn;
    memcpy(owned_Gnn_Gid + idx, shared[i].id_map[mp_id].Gid, ndofn * sizeof(*owned_Gnn_Gid));
  }

  /* initialize communication of the owned node information based on
     the communication hints -- NOTE: This is a SCATTER
     operation. ***Therefore, we use the receive hints to post
     sends.*** */
  const int nsend = Comm_hints_nrecv(hints);
  const int *send = Comm_hints_recv_list(hints);
  MPI_Request *req = static_cast<MPI_Request*>(malloc(nsend * sizeof(*req)));
  for (int i = 0; i < nsend; i++) {
    MPI_Isend(owned_Gnn_Gid, len_owned_Gnn_Gid, MPI_LONG,
                     send[i], GREDIST_TAG, comm, req + i);
  }

  /* reduce the total number of boundary nodes */
  MPI_Allreduce(&owned_gnn, &total_gnn, 1, MPI_LONG, MPI_SUM, comm);

  /* busy loop: Probe for message, post matching receive, do work --
     NOTE: This is a SCATTER operation. ***Therefore, we use the send
     hints to post receives.*** */
  const int nrecv = Comm_hints_nsend(hints);
  const int *recv = Comm_hints_send_list(hints);
  int *finished = static_cast<int*>(calloc(nrecv, sizeof(*finished)));
  int remaining = nrecv;
  while (remaining > 0) {
    int idx = 0;
    int msg_waiting = 0;
    MPI_Status stat;
    /* Probe for waiting message */
    for (idx = 0; idx < nrecv; idx++) {
      if (finished[idx]) continue;
      MPI_Iprobe(recv[idx], GREDIST_TAG, comm, &msg_waiting, &stat);
      if (msg_waiting) break;
    }

    /* if we didn't find a message waiting, then restart the busy
       loop */
    if (!msg_waiting) continue;

    /* determine the size of the incomming message and allocate an
       appropriately sized buffer */
    int msg_count = 0;
    MPI_Get_count(&stat, MPI_LONG, &msg_count);
    long *recv_Gnn_Gid = static_cast<long*>(malloc(msg_count * sizeof(*recv_Gnn_Gid)));

    /* We want to ensure that we are posting the *identically
       matching* receive, so we use the tag instead of MPI_ANY_TAG */
    MPI_Request breq;
    MPI_Irecv(recv_Gnn_Gid, msg_count, MPI_LONG,
                     recv[idx], stat.MPI_TAG, comm, &breq);

    /* While waiting for communication to complete, find the range of
       nodes we need to work on that are owned by the current
       domain */
    int range[2] = {0};
    if (nodes_get_shared_idx_range(n_shared, shared, recv[idx], range)) {
      PGFEM_printerr("ERROR: got bad hints on proc [%d] for proc [%d]! No matching nodes\n",
                     myrank, recv[idx]);
      PGFEM_Abort();
    }

    /* complete the receive communication */
    MPI_Wait(&breq, MPI_STATUS_IGNORE);

    /*
      Now we loop through the nodes that need information from the
      other domain and search for matches. We take advantage of the
      fact that both the key space (range) and search space
      (recv_Gnn_Gid) are sorted by Gnn. Therefore, we use binary
      search over a decreasing domain to find matches.
    */
    int size = ndofn + 1;
    const long *match = recv_Gnn_Gid;
    long key[3] = {0};
    for (int i = range[0], e = range[1]; i < e; i++) {
      int nmemb = (msg_count - (match - recv_Gnn_Gid)) / size;
      key[0] = shared[i].Gnn;
      match = static_cast<const long*>(bsearch(key, match, nmemb, size * sizeof(*match), compare_long));
      if (!match) {
        PGFEM_printerr("[%d] ERROR: did not find match!\n", myrank);
        PGFEM_Abort();
      }

      /* set Gid */
      for (int j = 0; j < ndofn; j++) {
        if (shared[i].id_map[mp_id].id[j] <= 0) {
          /* BC overrides periodicity */
          shared[i].id_map[mp_id].Gid[j] = shared[i].id_map[mp_id].id[j];
        } else {
          shared[i].id_map[mp_id].Gid[j] = match[j + 1];
        }
      }
    }

    /* mark job as done */
    remaining--;
    finished[idx] = 1;

    /* loop cleanup */
    free(recv_Gnn_Gid);
  }

  /* re-sort the nodes by their local number */
  nodes_sort_loc_id(nnode, nodes);

  /* complete non-blocking sends */
  MPI_Waitall(nsend, req, MPI_STATUS_IGNORE);

  /* cleanup */
  free(owned_Gnn_Gid);
  free(finished);
  free(req);
  return total_gnn;
}

long GRedist_node (const int nproc,
                   const int myrank,
                   const long nn,
                   const long ndofn,
                   Node *node,
                   const Comm_hints *hints,
                   const MPI_Comm Comm,
                   const int mp_id)
{
  if (!hints) return fallback_GRedist_node(nproc, myrank, nn, ndofn, node, Comm, mp_id);
  else return comm_hints_GRedist_node(nproc, myrank, nn, ndofn, node, hints, Comm, mp_id);
}
