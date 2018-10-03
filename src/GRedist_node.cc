#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "GRedist_node.h"
#include "PGFEM_io.h"
#include "allocation.h"
#include "matice.h"
#include <cstdlib>
#include <cstring>
#include <cassert>

#ifndef NDEBUG
#define PFEM_DEBUG 1
#else
#define PFEM_DEBUG 0
#endif

#define GREDIST_TAG 999

using namespace pgfem3d;
using namespace pgfem3d::net;

static long fallback_GRedist_node_ISIR(const int nproc,
				       const int myrank,
				       const long nn,
				       const long ndofn,
				       Node *node,
				       const CommunicationStructure *com,
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
  com->net->allgather(&GN,1,NET_DT_LONG,BN,1,NET_DT_LONG,com->comm);

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
  //assert(NBN > 0 && "NBN must be > 0");
  if(NBN > 0) Gnn_Gid = static_cast<long*>(malloc(NBN*own_buf_elem_size));
  com->net->allgatherv(own_buf,recvcount[myrank],NET_DT_LONG,
		       Gnn_Gid,recvcount,displ,NET_DT_LONG,com->comm);

  /* sort the list of Gnn and their associated Gid by Gnn */
  if (Gnn_Gid != NULL)
    qsort(Gnn_Gid,NBN,own_buf_elem_size,compare_long);

  /* sort need_buf by Gnn */
  if (need_buf) qsort(need_buf,need,need_buf_elem_size,compare_long);

  if (PFEM_DEBUG){
    /* Check global node numbers, should be ordered and contiguous */
    long k = 0;
    for (auto i=0; i<NBN*(ndofn+1); i+=ndofn+1){
      assert(Gnn_Gid != NULL && "Gnn_Gid can't be NULL");
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
	assert(Gnn_Gid != NULL && "Gnn_Gid can't be NULL");
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
  if(NBN > 0) free(Gnn_Gid);
  free(need_buf);
  free(own_buf);
  free(recvcount);
  free(displ);

  return (NBN);
}

static long fallback_GRedist_node_PWC(const int nproc,
				      const int myrank,
				      const long nn,
				      const long ndofn,
				      Node *node,
				      const CommunicationStructure *com,
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
  com->net->allgather(&GN,1,NET_DT_LONG,BN,1,NET_DT_LONG,com->comm);

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
  com->net->allgatherv(own_buf,recvcount[myrank],NET_DT_LONG,
		       Gnn_Gid,recvcount,displ,NET_DT_LONG,com->comm);

  /* sort the list of Gnn and their associated Gid by Gnn */
  if (Gnn_Gid != NULL)
    qsort(Gnn_Gid,NBN,own_buf_elem_size,compare_long);

  /* sort need_buf by Gnn */
  if (need_buf) qsort(need_buf,need,need_buf_elem_size,compare_long);

  if (PFEM_DEBUG){
    /* Check global node numbers, should be ordered and contiguous */
    long k = 0;
    for (auto i=0; i<NBN*(ndofn+1); i+=ndofn+1){
      assert(Gnn_Gid != NULL && "Gnn_Gid can't be NULL");
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
	assert(Gnn_Gid != NULL && "Gnn_Gid can't be NULL");
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
  if(NBN > 0) free(Gnn_Gid);
  free(need_buf);
  free(own_buf);
  free(recvcount);
  free(displ);

  return (NBN);
}

/**
 * Original implementation of GRedist_node
 */
static long fallback_GRedist_node(const int nproc,
                                  const int myrank,
                                  const long nn,
                                  const long ndofn,
                                  Node *node,
                                  const CommunicationStructure *com,
                                  const long mp_id)
{
  switch (com->net->type()) {
  case NET_ISIR:
    return fallback_GRedist_node_ISIR(nproc, myrank, nn, ndofn, node, com, mp_id);
    break;
  case NET_PWC:
    return fallback_GRedist_node_PWC(nproc, myrank, nn, ndofn, node, com, mp_id);
    break;
  default:
    PGFEM_Abort();    
  }
  return 1;
}

/**
 * New implementation of GRedist_node that avoids global communication
 *
 * \return total number of global nodes
 */
static long comm_hints_GRedist_node_ISIR(const int nproc,
					 const int myrank,
					 const long nnode,
					 const long ndofn,
					 Node *nodes,
					 const CommunicationStructure *com,
					 const long mp_id)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(com->net);
  
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
  const int nsend = com->hints->get_nrecv();
  const int *send = com->hints->get_recv_list();
  Request *req;
  net->allocRequestArray(nsend, &req);
  for (int i = 0; i < nsend; i++) {
    net->isend(owned_Gnn_Gid, len_owned_Gnn_Gid, NET_DT_LONG,
		    send[i], GREDIST_TAG, com->comm, &req[i]);
  }

  /* reduce the total number of boundary nodes */
  net->allreduce(&owned_gnn, &total_gnn, 1, NET_DT_LONG, NET_OP_SUM, com->comm);

  /* busy loop: Probe for message, post matching receive, do work --
     NOTE: This is a SCATTER operation. ***Therefore, we use the send
     hints to post receives.*** */
  const int nrecv = com->hints->get_nsend();
  const int *recv = com->hints->get_send_list();
  int *finished = static_cast<int*>(calloc(nrecv, sizeof(*finished)));
  int remaining = nrecv;
  while (remaining > 0) {
    int idx = 0;
    int msg_waiting = 0;
    Status stat;
    /* Probe for waiting message */
    for (idx = 0; idx < nrecv; idx++) {
      if (finished[idx]) continue;
      net->iprobe(recv[idx], GREDIST_TAG, com->comm, &msg_waiting, &stat);
      if (msg_waiting) break;
    }

    /* if we didn't find a message waiting, then restart the busy
       loop */
    if (!msg_waiting) continue;

    /* determine the size of the incomming message and allocate an
       appropriately sized buffer */
    int msg_count = 0;
    net->get_status_count(&stat, NET_DT_LONG, &msg_count);
    long *recv_Gnn_Gid = static_cast<long*>(malloc(msg_count * sizeof(*recv_Gnn_Gid)));

    /* We want to ensure that we are posting the *identically
       matching* receive, so we use an explicit tag */
    Request breq;
    net->irecv(recv_Gnn_Gid, msg_count, NET_DT_LONG,
	       recv[idx], stat.NET_TAG, com->comm, &breq);

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
    net->wait(&breq, NET_STATUS_IGNORE);
    
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
  net->waitall(nsend, req, NET_STATUS_IGNORE);

  /* cleanup */
  free(owned_Gnn_Gid);
  free(finished);
  delete [] req;
  return total_gnn;
}

static long comm_hints_GRedist_node_PWC(const int nproc,
					const int myrank,
					const long nnode,
					const long ndofn,
					Node *nodes,
					const CommunicationStructure *com,
					const long mp_id)
{
  PWCNetwork *net = static_cast<PWCNetwork*>(com->net);
  CID lid = 0xcafebeee;

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
  Buffer sbuffer;
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
  
    /* pin send buffer */
    sbuffer.addr = reinterpret_cast<uintptr_t> (owned_Gnn_Gid);
    sbuffer.size = len_owned_Gnn_Gid * sizeof(*owned_Gnn_Gid);
    net->pin(reinterpret_cast<void *>(sbuffer.addr), sbuffer.size, &sbuffer.key);
  }  

  /* Send allocation information */
  const int nsend = com->hints->get_nrecv();
  const int *send = com->hints->get_recv_list();
  for (int i = 0; i < nsend; i++) {
    CID rid = (CID)len_owned_Gnn_Gid;
    net->pwc(send[i], 0, 0, 0, lid, rid);
  }

  /* wait for local PWC completions from above 
     receive the size of recv buffer*/
  const int nrecv = com->hints->get_nsend();
  long *len_recv_Gnn_Gid = PGFEM_calloc (long, nproc);
  for (int i = 0; i < nproc; i++) {
    len_recv_Gnn_Gid[i] = -1;
  }

  net->wait_n_id(nsend, lid);

  int t_count = 0;
  while (t_count < nrecv) {
    int flag;
    CID val;
    Status stat;
    net->probe(&flag, &val, &stat, 0);
    if (flag) {
      int p = stat.NET_SOURCE;
      /* the long values exchanged are encoded in the PWC RIDs */
      len_recv_Gnn_Gid[p] = (long)val;
      t_count++;
    }
  }

  /* allocate and pin recv buffers */
  long **RECI = PGFEM_calloc (long*, nproc);
  Buffer *rbuffers = PGFEM_calloc (Buffer, nproc);
  for (int i = 0; i < nproc; i++) {
    long JJ = 0;
    if (len_recv_Gnn_Gid[i] <= 0) JJ = 1; else JJ = len_recv_Gnn_Gid[i];
    RECI[i] = PGFEM_calloc_pin (long, JJ,
				net, &rbuffers[i].key);
    rbuffers[i].addr = reinterpret_cast<uintptr_t> (RECI[i]);
    rbuffers[i].size = sizeof(long)*JJ;
  }

  /* exchange receive buffers */
  for (int i = 0; i < nproc; i++) {
    net->gather(&rbuffers[i], sizeof(Buffer), NET_DT_BYTE,
	net->getbuffer(), sizeof(Buffer), NET_DT_BYTE, i, com->comm);
  }

  /* initialize communication of the owned node information based on
     the communication hints -- NOTE: This is a SCATTER
     operation. ***Therefore, we use the receive hints to post
     sends.*** */
  for (int i = 0; i < nsend; i++) {
    CID rid = (CID)len_owned_Gnn_Gid;
    net->pwc(send[i], sbuffer.size, &sbuffer, &net->getbuffer()[send[i]], lid, rid);
  }

  /* reduce the total number of boundary nodes */
  net->allreduce(&owned_gnn, &total_gnn, 1, NET_DT_LONG, NET_OP_SUM, com->comm);

  /* Wait to complete the communications */
  net->wait_n_id(nsend, lid);

  /* busy loop: Probe for message, post matching receive, do work --
     NOTE: This is a SCATTER operation. ***Therefore, we use the send
     hints to post receives.*** */
  int remaining = nrecv;
  while (remaining > 0) {
    int flag;
    CID val;
    Status stat;
    /* Probe for finished message */
    net->probe(&flag, &val, &stat, 0);
    /* if we didn't find a message waiting, then restart the busy
       loop */
    if (flag) {
      int source = stat.NET_SOURCE;
      long *recv_Gnn_Gid = RECI[source];
      long msg_count = len_recv_Gnn_Gid[source];

      /* While waiting for communication to complete, find the range of
         nodes we need to work on that are owned by the current
         domain */
      int range[2] = {0};
      if (nodes_get_shared_idx_range(n_shared, shared, source, range)) {
        PGFEM_printerr("ERROR: got bad hints on proc [%d] for proc [%d]! No matching nodes\n",
                       myrank, source);
        PGFEM_Abort();
      }

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
    }
  }

  /* re-sort the nodes by their local number */
  nodes_sort_loc_id(nnode, nodes);

  /* cleanup */
  net->unpin(reinterpret_cast<void *> (sbuffer.addr), sbuffer.size);
  for (int i = 0; i < nproc; i++) {
    net->unpin(reinterpret_cast<void *> (rbuffers[i].addr), rbuffers[i].size);
  }

  dealoc1l(rbuffers);
  dealoc2l(RECI,nproc);
  free(owned_Gnn_Gid);
  return total_gnn;
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
				    const CommunicationStructure *com,
                                    const long mp_id)
{
  switch (com->net->type()) {
  case NET_ISIR:
    return comm_hints_GRedist_node_ISIR(nproc, myrank, nnode, ndofn, nodes, com, mp_id);
    break;
  case NET_PWC:
    return comm_hints_GRedist_node_PWC(nproc, myrank, nnode, ndofn, nodes, com, mp_id);
    break;
  default:
    PGFEM_Abort();    
  }
  return 1;
}

long GRedist_node (const int nproc,
                   const int myrank,
                   const long nn,
                   const long ndofn,
                   Node *node,
		   const CommunicationStructure *com,
                   const int mp_id)
{
  if (!(com->hints)) return fallback_GRedist_node(nproc, myrank, nn, ndofn, node, com, mp_id);
  else return comm_hints_GRedist_node(nproc, myrank, nn, ndofn, node, com, mp_id);
}
