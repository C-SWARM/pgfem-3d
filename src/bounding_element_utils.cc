#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "bounding_element_utils.h"
#include "cohesive_element_utils.h"
#include "enumerations.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "incl.h"
#include "index_macros.h"
#include "matice.h"
#include "transform_coordinates.h"
#include "utils.h"
#include <cmath>

using namespace pgfem3d;
using namespace pgfem3d::net;

#ifndef DEBUG_BE_UTILS
#define DEBUG_BE_UTILS 0
#endif

int bounding_element_set_local_ids(const int n_be,
                                   BoundingElement *b_elems,
                                   const Element *elem)
{
  /* loop through bounding elements and get local node id's from
     volume element connectivity */
  int err = 0;
  for(int i=0; i<n_be; i++){
    const BoundingElement *ptr_be = &b_elems[i];
    const long *be_nod = ptr_be->nodes;
    const int be_nne = ptr_be->nnodes;

    const Element *ptr_ve = &elem[ptr_be->vol_elem_id];
    const long *ve_nod = ptr_ve->nod;
    const int ve_nne = ptr_ve->toe;

    /* this is what we are after */
    long *loc_nod = ptr_be->loc_nodes;

    for(int j=0; j<be_nne; j++){
      for(int k=0; k<ve_nne; k++){
        if(be_nod[j] == ve_nod[k]){
          loc_nod[j] = k;
          break;
        }
      }
    }

  }/* for each bounding element */
  return err;
}/* bounding_element_compute_local_ids */

int bounding_element_reverse_mapping(const int n_be,
                                     const BoundingElement *b_elems,
                                     Element *elem)
{
  int err = 0;
  if(n_be == 0){ /* no bounding elements on domain, exit */
    return err;
  }

  /* loop through b_elems and store id,vol_elem_id pairs. sort pairs
     by vol_elem_id. create list of b_elems assocaiated with each
     unique vol_elem_id. loop through elems and set n_be, allocate and
     populate lists. */

  val_key *map = PGFEM_calloc(val_key, n_be);
  for(int i=0; i<n_be; i++){
    map[i].val = b_elems[i].vol_elem_id;
    map[i].key = i;
  }

  qsort(map,n_be,sizeof(val_key),compare_val_w_key);

  int idx = 1;
  int n_vol_elems = 1;
  while(idx<n_be){
    if(map[idx].val > map[idx-1].val){
      n_vol_elems ++;
    }
    idx++;
  }

  int *vol_elem_ids = PGFEM_calloc(int, n_vol_elems);
  int *n_be_on_elem = PGFEM_calloc(int, n_vol_elems);

  idx = 0;
  vol_elem_ids[idx] = map[0].val;
  n_be_on_elem[idx] = 1;
  for(int i=1; i<n_be; i++){
    if(map[i].val > map[i-1].val){
      idx ++;
      vol_elem_ids[idx] = map[i].val;
    }
    n_be_on_elem[idx]++;
  }

  for(int i=0; i<n_vol_elems; i++){
    /* compute where to index from in the mapping */
    int map_start = 0;
    for(int j=0; j<i; j++){
      map_start += n_be_on_elem[j];
    }

    /* set values on the element */
    Element *ptr_elem = &elem[vol_elem_ids[i]];
    ptr_elem->n_be = n_be_on_elem[i];
    ptr_elem->be_ids = PGFEM_calloc(long, ptr_elem->n_be);
    for(int j=0; j<ptr_elem->n_be; j++){
      ptr_elem->be_ids[j] = map[map_start+j].key;
    }
  }

  free(map);
  free(vol_elem_ids);
  free(n_be_on_elem);

  return err;
}

static int _bounding_element_communicate_damage_ISIR(const int n_be,
						     BoundingElement *b_elems,
						     const int ne,
						     const EPS *eps,
						     const CommunicationStructure *com)
{
  int err = 0;
  int myrank = com->rank;
  int nproc = com->nproc;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(com->net);
  
  /* no bounding elements on domain, return without doing anything. */
  if(n_be == 0){
    return err;
  }

  /* loop through the boundary elements and determine how many to send
     to each process. The number sent to each process is also the
     number to recieve from each process since the master-slave
     relationship is one-to-one. */
  int *n_SR = PGFEM_calloc(int, nproc);
  /* src_dest tells where to send stuff to and where it came
     from: dest_dom dest_id src_id (src_id is always local) */
  const int n_sd = 3;
  int *src_dest = PGFEM_calloc(int, n_be*n_sd);
  for(int src_id=0; src_id<n_be; src_id++){
    const BoundingElement *pbe = &b_elems[src_id];
    if(pbe->master){
      src_dest[src_id*n_sd + 0] = pbe->slave_dom;
      src_dest[src_id*n_sd + 1] = pbe->slave_be_id;
      src_dest[src_id*n_sd + 2] = src_id;
      n_SR[pbe->slave_dom]++;
    } else {
      src_dest[src_id*n_sd + 0] = pbe->master_dom;
      src_dest[src_id*n_sd + 1] = pbe->master_be_id;
      src_dest[src_id*n_sd + 2] = src_id;
      n_SR[pbe->master_dom]++;
    }
  }

  /* sort src_dest by dest_dom */
  qsort(src_dest,n_be,n_sd*sizeof(int),compare_int);

  int **receive_id = PGFEM_calloc(int*, nproc);
  int **send_id = PGFEM_calloc(int*, nproc);
  double **receive_val = PGFEM_calloc(double*, nproc);
  double **send_val = PGFEM_calloc(double*, nproc);
  int *p_src_dest = &src_dest[0];
  int n_proc_comm = 0;
  for(int i=0; i<nproc; i++){
    const int n_val = n_SR[i];
    if(n_val > 0){ /* I will send & recieve from this process */
      /* allocate space to recieve */
      receive_id[i] = PGFEM_calloc(int, n_val);
      receive_val[i] = PGFEM_calloc(double, n_val);

      /* allocate and populate stuff to send */
      send_id[i] = PGFEM_calloc(int, n_val); /* this is the destination id */
      send_val[i] = PGFEM_calloc(double, n_val);

      int *p_send_id = send_id[i];
      double *p_send_val = send_val[i];

      for(int j=0; j<n_val; j++){
        p_send_id[j] = p_src_dest[1];
        const int ve_id = b_elems[p_src_dest[2]].vol_elem_id;
        p_send_val[j] = eps[ve_id].dam[0].w;
        p_src_dest += n_sd;
      }

      n_proc_comm++;
    } else {
      receive_id[i] = NULL;
      receive_val[i] = NULL;
      send_id[i] = NULL;
      send_val[i] = NULL;
    }
  }

  /* remove myrank from number of comm procs */
  if(n_SR[myrank] > 0){
    n_proc_comm--;
  }

  /*==================================*/
  /*  x     post send and recieve      */
  /*==================================*/
  Status *sta_s;
  Status *sta_r;
  Request *req_s;
  Request *req_r;

  net->allocStatusArray(2*n_proc_comm, &sta_s);
  net->allocStatusArray(2*n_proc_comm, &sta_r);
  net->allocRequestArray(2*n_proc_comm, &req_s);
  net->allocRequestArray(2*n_proc_comm, &req_r);
  
  int count = 0;
  for(int i=0; i<nproc; i++){
    if(n_SR[i] > 0 && i != myrank){
      /* post receives */
      net->irecv(receive_id[i], n_SR[i], NET_DT_INT, i,
		      NET_ANY_TAG, com->comm, &req_r[2*count+0]);
      net->irecv(receive_val[i], n_SR[i], NET_DT_DOUBLE, i,
		      NET_ANY_TAG, com->comm, &req_r[2*count+1]);

      /* post sends */
      net->isend(send_id[i], n_SR[i], NET_DT_INT, i,
		      NET_ANY_TAG, com->comm, &req_s[2*count+0]);
      net->isend(send_val[i], n_SR[i], NET_DT_DOUBLE, i,
		      NET_ANY_TAG, com->comm, &req_s[2*count+1]);
      count++;
    }
  }

  /*==================================*/
  /*       local computations         */
  /*==================================*/
  if(n_SR[myrank] > 0){ /* have local elements */
    p_src_dest = src_dest;
    int pos = 0;
    while(*p_src_dest != myrank
          && pos < n_be*n_sd){
      pos+= n_sd;
      p_src_dest += n_sd;
    }

    for(int i=0; i<n_SR[myrank]; i++){
      BoundingElement *p_dest_be = &b_elems[p_src_dest[1]];
      const BoundingElement *p_src_be = &b_elems[p_src_dest[2]];
      p_src_dest += n_sd;
      p_dest_be->other_val = eps[p_src_be->vol_elem_id].dam[0].w;
    }
  }
  /*==================================*/
  /*       Finish communication       */
  /*==================================*/
  net->waitall(2*n_proc_comm, req_s, sta_s);
  net->waitall(2*n_proc_comm, req_r, sta_r);

  for(int i=0; i<nproc; i++){
    if(n_SR[i] > 0 && i != myrank){
      for(int j=0; j<n_SR[i]; j++){
        BoundingElement *p_dest_be = &b_elems[receive_id[i][j]];
        p_dest_be->other_val = receive_val[i][j];
      }
    }
  }

  /* debug, print communication */
  static int repeat = 0;
  if(DEBUG_BE_UTILS && !repeat){
    FILE *out;
    char fname[100];
    sprintf(fname,"%s_%05d.log",__func__,myrank);
    out = fopen(fname,"w");
    PGFEM_fprintf(out,"=== src_dest structure ===\n");
    for(int i=0; i<n_be; i++){
      PGFEM_fprintf(out,"dest dom: %3d || dest id: %5d || src id: %5d\n",
                    src_dest[i*n_sd],src_dest[i*n_sd+1],src_dest[i*n_sd+2]);
    }
    fclose(out);
  }
  repeat++;

  free(n_SR);
  free(src_dest);

  for(int i=0; i<nproc; i++){
    free(receive_id[i]);
    free(receive_val[i]);
    free(send_id[i]);
    free(send_val[i]);
  }
  free(receive_id);
  free(receive_val);
  free(send_id);
  free(send_val);
  
  delete [] sta_s;
  delete [] sta_r;
  delete [] req_s;
  delete [] req_r;

  return err;
}

static int _bounding_element_communicate_damage_PWC(const int n_be,
						    BoundingElement *b_elems,
						    const int ne,
						    const EPS *eps,
						    const CommunicationStructure *com)
{
  int err = 0;
  int myrank = com->rank;
  int nproc = com->nproc;
  PWCNetwork *net = static_cast<PWCNetwork*>(com->net);
  CID lid = 0xcafebabe;  

  /* no bounding elements on domain, return without doing anything. */
  if(n_be == 0){
    return err;
  }

  /* loop through the boundary elements and determine how many to send
     to each process. The number sent to each process is also the
     number to recieve from each process since the master-slave
     relationship is one-to-one. */
  int *n_SR = PGFEM_calloc(int, nproc);
  /* src_dest tells where to send stuff to and where it came
     from: dest_dom dest_id src_id (src_id is always local) */
  const int n_sd = 3;
  int *src_dest = PGFEM_calloc(int, n_be*n_sd);
  for(int src_id=0; src_id<n_be; src_id++){
    const BoundingElement *pbe = &b_elems[src_id];
    if(pbe->master){
      src_dest[src_id*n_sd + 0] = pbe->slave_dom;
      src_dest[src_id*n_sd + 1] = pbe->slave_be_id;
      src_dest[src_id*n_sd + 2] = src_id;
      n_SR[pbe->slave_dom]++;
    } else {
      src_dest[src_id*n_sd + 0] = pbe->master_dom;
      src_dest[src_id*n_sd + 1] = pbe->master_be_id;
      src_dest[src_id*n_sd + 2] = src_id;
      n_SR[pbe->master_dom]++;
    }
  }

  /* sort src_dest by dest_dom */
  qsort(src_dest,n_be,n_sd*sizeof(int),compare_int);

  int **receive_id = PGFEM_calloc(int*, nproc);
  int **send_id = PGFEM_calloc(int*, nproc);
  double **receive_val = PGFEM_calloc(double*, nproc);
  double **send_val = PGFEM_calloc(double*, nproc);
  Buffer sbuf_id, sbuf_val, rbuf_id, rbuf_val;
  
  Buffer *rbuffers = PGFEM_calloc (Buffer, nproc*2);
  Buffer *sbuffers = PGFEM_calloc (Buffer, nproc*2);

  /* Determing total backing buffer size for rbuf and sbuf */
  size_t n_size = 0;
  for (int i=0; i<nproc; i++) {
    const int n_val = n_SR[i];
    if (n_val > 0) {
      n_size += n_val;
    }
  }
  sbuf_id.addr = reinterpret_cast<uintptr_t> (PGFEM_calloc_pin(int, n_size, net,
							       &sbuf_id.key));
  sbuf_id.size = sizeof(int)*n_size;
  sbuf_val.addr = reinterpret_cast<uintptr_t> (PGFEM_calloc_pin(double, n_size, net,
								&sbuf_val.key));
  sbuf_val.size = sizeof(double)*n_size;

  rbuf_id.addr = reinterpret_cast<uintptr_t> (PGFEM_calloc_pin(int, n_size, net,
							       &rbuf_id.key));
  rbuf_id.size = sizeof(int)*n_size;
  rbuf_val.addr = reinterpret_cast<uintptr_t> (PGFEM_calloc_pin(double, n_size, net,
								&rbuf_val.key));
  rbuf_val.size = sizeof(double)*n_size;
  
  int *p_src_dest = &src_dest[0];
  int n_proc_comm = 0;
  int n_id_offset = 0;
  int n_val_offset = 0;
  for (int i=0; i<nproc; i++) {
    const int n_val = n_SR[i];
    if(n_val > 0){ /* I will send & recieve from this process */
      /* allocate space to recieve */
      receive_id[i] = reinterpret_cast<int*> (rbuf_id.addr + n_id_offset);
      receive_val[i] = reinterpret_cast<double*> (rbuf_val.addr + n_val_offset);

      rbuffers[2*i].addr = reinterpret_cast<uintptr_t> (receive_id[i]);
      rbuffers[2*i].size = sizeof(int) * n_val;
      rbuffers[2*i].key = rbuf_id.key;
      
      rbuffers[2*i+1].addr = reinterpret_cast<uintptr_t> (receive_val[i]);
      rbuffers[2*i+1].size = sizeof(double) * n_val;
      rbuffers[2*i+1].key = rbuf_val.key;
      
      /* allocate and populate stuff to send */
      send_id[i] = reinterpret_cast<int*> (sbuf_id.addr + n_id_offset);
      send_val[i] = reinterpret_cast<double*> (sbuf_val.addr + n_val_offset);
      
      sbuffers[2*i].addr = reinterpret_cast<uintptr_t> (send_id[i]);
      sbuffers[2*i].size = sizeof(int) * n_val;
      sbuffers[2*i].key = sbuf_id.key;

      sbuffers[2*i+1].addr = reinterpret_cast<uintptr_t> (send_val[i]);
      sbuffers[2*i+1].size = sizeof(double) * n_val;
      sbuffers[2*i+1].key = sbuf_val.key;

      int *p_send_id = send_id[i];
      double *p_send_val = send_val[i];

      for(int j=0; j<n_val; j++){
        p_send_id[j] = p_src_dest[1];
        const int ve_id = b_elems[p_src_dest[2]].vol_elem_id;
        p_send_val[j] = eps[ve_id].dam[0].w;
        p_src_dest += n_sd;
      }

      n_proc_comm++;
      n_id_offset = sizeof(int)*n_val;
      n_val_offset = sizeof(double)*n_val;
    } else {
      receive_id[i] = NULL;
      receive_val[i] = NULL;

      rbuffers[2*i].addr = 0;
      rbuffers[2*i].size = 0;
      rbuffers[2*i+1].addr = 0;
      rbuffers[2*i+1].size = 0;

      send_id[i] = NULL;
      send_val[i] = NULL;

      sbuffers[2*i].addr = 0;
      sbuffers[2*i].size = 0;
      sbuffers[2*i+1].addr = 0;
      sbuffers[2*i+1].size = 0;
    }
  }

  /* Allocate local write buffers */  
  Buffer *wbuffers = new Buffer[2 * nproc];

  /* Exchange receive buffers */
  for (int i = 0; i < nproc; i++) {
    net->gather(&rbuffers[2*i], sizeof(Buffer)*2, NET_DT_BYTE,
		wbuffers, sizeof(Buffer)*2, NET_DT_BYTE, i, com->comm);
  }

  /* remove myrank from number of comm procs */
  if(n_SR[myrank] > 0){
    n_proc_comm--;
  }

  /*==================================*/
  /*  x     post send and recieve      */
  /*==================================*/
  int num_send = 0;
  for(int i=0; i<nproc; i++){
    if(n_SR[i] > 0 && i != myrank){
      /* send with pwc */
      CID rid = (CID)n_SR[i];
      net->pwc(i, sbuffers[2*i].size, &sbuffers[2*i], &wbuffers[2*i], lid, rid);
      net->pwc(i, sbuffers[2*i+1].size, &sbuffers[2*i+1], &wbuffers[2*i+1], lid, rid);
      num_send += 2;
    }
  }

  /*==================================*/
  /*       local computations         */
  /*==================================*/
  if(n_SR[myrank] > 0){ /* have local elements */
    p_src_dest = src_dest;
    int pos = 0;
    while(*p_src_dest != myrank
          && pos < n_be*n_sd){
      pos+= n_sd;
      p_src_dest += n_sd;
    }

    for(int i=0; i<n_SR[myrank]; i++){
      BoundingElement *p_dest_be = &b_elems[p_src_dest[1]];
      const BoundingElement *p_src_be = &b_elems[p_src_dest[2]];
      p_src_dest += n_sd;
      p_dest_be->other_val = eps[p_src_be->vol_elem_id].dam[0].w;
    }
  }

  /*==================================*/
  /*       Finish communication       */
  /*==================================*/
  net->wait_n_id(num_send, lid);

  int t_count = 0;
  while (t_count < num_send) {
    int flag;
    CID val;
    Status stat;
    net->probe(&flag, &val, &stat, 0);
    if (flag) {
      int p = stat.NET_SOURCE;
      for(int j=0; j<n_SR[p]; j++){
        BoundingElement *p_dest_be = &b_elems[receive_id[p][j]];
        p_dest_be->other_val = receive_val[p][j];
      }
      t_count++;
    }
  }

  net->barrier(com->comm);

  /* debug, print communication */
  static int repeat = 0;
  if(DEBUG_BE_UTILS && !repeat){
    FILE *out;
    char fname[100];
    sprintf(fname,"%s_%05d.log",__func__,myrank);
    out = fopen(fname,"w");
    PGFEM_fprintf(out,"=== src_dest structure ===\n");
    for(int i=0; i<n_be; i++){
      PGFEM_fprintf(out,"dest dom: %3d || dest id: %5d || src id: %5d\n",
                    src_dest[i*n_sd],src_dest[i*n_sd+1],src_dest[i*n_sd+2]);
    }
    fclose(out);
  }
  repeat++;

  net->unpin(reinterpret_cast<void *> (sbuf_id.addr), sbuf_id.size);
  net->unpin(reinterpret_cast<void *> (sbuf_val.addr), sbuf_val.size);
  net->unpin(reinterpret_cast<void *> (rbuf_id.addr), rbuf_id.size);
  net->unpin(reinterpret_cast<void *> (rbuf_val.addr), rbuf_val.size);

  dealoc1l((void*)sbuf_id.addr);
  dealoc1l((void*)sbuf_val.addr);
  dealoc1l((void*)rbuf_id.addr);
  dealoc1l((void*)rbuf_val.addr);
  
  dealoc1l(rbuffers);
  dealoc1l(sbuffers);

  delete [] wbuffers;

  free(n_SR);
  free(src_dest);

  free(receive_id);
  free(receive_val);
  free(send_id);
  free(send_val);

  return err;
}

int bounding_element_communicate_damage(const int n_be,
                                        BoundingElement *b_elems,
                                        const int ne,
                                        const EPS *eps,
                                        const CommunicationStructure *com)
{
  switch(com->net->type()) {
  case NET_ISIR:
    return _bounding_element_communicate_damage_ISIR(n_be, b_elems, ne, eps, com);
    break;
  case NET_PWC:
    return _bounding_element_communicate_damage_PWC(n_be, b_elems, ne, eps, com);
    break;
  default:
    PGFEM_Abort();
  }
  return 1;
}

int bounding_element_compute_resulting_traction(const int n_be,
                                                const BoundingElement *b_elems,
                                                const Element *elems,
                                                const Node *nodes,
                                                const EPS *eps,
                                                const SIG *sig,
                                                const int ndofd,
                                                const long *DomDof,
                                                const int GDof,
						const CommunicationStructure *com,
                                                const int analysis,
                                                double *res_trac)
{
  int err = 0;
  int myrank = com->rank;

  /* Check that bounding elements exist and we should proceed,
     otehrwise exit early */
  int g_be = n_be;
  com->net->allreduce(NET_IN_PLACE,&g_be,1,NET_DT_INT,NET_OP_SUM,com->comm);
  if(g_be <= 0){
    return err;
  }

  double *L_trac = aloc1(ndofd);
  double *L_trac_sgn = aloc1(ndofd);
  double *G_trac = aloc1(DomDof[myrank]);
  double *G_trac_sgn = aloc1(DomDof[myrank]);
  nulld(res_trac,3);

  /*** THIS FUNCTION ASSUMES LINEAR ELEMENTS ***/

  /*
   * F -> eps[ve_id].il[ip].F
   * S -> eps[ve_id].il[ip].o (sym formmat)
   */

  for(int be_id=0; be_id<n_be; be_id++){
    const int ip = 0;
    const BoundingElement *p_be = &b_elems[be_id];
    const long *L_be_dofs = p_be->L_dof_ids;
    const int nn_be = p_be->nnodes;
    const int ve_id = p_be->vol_elem_id;
    const long *b_nodes = p_be->nodes;
    const double *F = eps[ve_id].il[ip].F;
    const double *Ssym =  sig[ve_id].il[ip].o;
    double *S = aloc1(9);
    S[0] = Ssym[0]; S[4] = Ssym[1]; S[8] = Ssym[2];
    S[idx_2(1,2)] = S[idx_2(2,1)] = Ssym[3];
    S[idx_2(0,2)] = S[idx_2(2,0)] = Ssym[4];
    S[idx_2(0,1)] = S[idx_2(1,0)] = Ssym[5];

    /* linear element, so constant values, integrate at centroid */
    const double ksi = 1./3.;
    const double eta = 1./3.;
    const double wt = 1.0;
    double jj = 0.0;
    double *normal = aloc1(3);
    double sgn = 1.0;

    /* Compute Jacobian of the transformation and the normal to the
       element surface */
    {
      double *x = aloc1(nn_be);
      double *y = aloc1(nn_be);
      double *z = aloc1(nn_be);
      double *N2d = aloc1(nn_be);
      double *e1 = aloc1(3);
      double *e2h = aloc1(3);
      double *e2 = aloc1(3);
      double *xl = aloc1(nn_be);
      double *yl = aloc1(nn_be);
      double *zl = aloc1(nn_be);
      double *Nx = aloc1(nn_be);
      double *Ny = aloc1(nn_be);

      switch(analysis){
       case DISP:
        nodecoord_total(nn_be,b_nodes,nodes,x,y,z);
        break;
       default:
        nodecoord_updated(nn_be,b_nodes,nodes,x,y,z);
        break;
      }
      base_vec(nn_be,ksi,eta,x,y,z,e1,e2,e2h,normal,myrank);
      transform_coordinates(nn_be,x,y,z,e1,e2,normal,0,xl,yl,zl);
      jj = dN3_xy(ksi,eta,nn_be,xl,yl,zl,Nx,Ny);

      free(x);
      free(y);
      free(z);
      free(N2d);
      free(e1);
      free(e2);
      free(e2h);
      free(xl);
      free(yl);
      free(zl);
      free(Nx);
      free(Ny);
    }

    /* compute sgn (+-1) according to master/slave */
    sgn = normal[0] + normal[1] + normal[2];
    if(sgn < 0.0){
      sgn = -1.0;
    } else {
      sgn = 1.0;
    }

    /* integrate traction on bounding element */
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        for(int k=0; k<3; k++){
          const int ik = idx_2(i,k);
          const int kj = idx_2(k,j);
          L_trac[L_be_dofs[i]-1] += F[ik]*S[kj]*normal[j]*jj*wt;
          L_trac_sgn[L_be_dofs[i]-1] += sgn*F[ik]*S[kj]*normal[j]*jj*wt;
          res_trac[i] += F[ik]*S[kj]*normal[j]*jj*wt;
        }
      }
    }

    free(S);
  }/* for each be on local domain */

  /* sum up tractions on all domains */
  com->net->allreduce(NET_IN_PLACE,res_trac,3,NET_DT_DOUBLE,NET_OP_SUM,com->comm);

  /* get GLOBAL trac/trac_sgn */
  LToG(L_trac,G_trac,ndofd,com);
  LToG(L_trac_sgn,G_trac_sgn,ndofd,com);

  /* compute norms */
  double *norms = aloc1(3);
  for(int i=0; i<DomDof[myrank]; i++){
    norms[0] += G_trac[i]*G_trac[i];
    norms[1] += G_trac_sgn[i]*G_trac_sgn[i];
  }
  com->net->allreduce(NET_IN_PLACE,norms,2,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  norms[0] = sqrt(norms[0]);
  norms[1] = sqrt(norms[1]);
  norms[2] = 2.*norms[0]/norms[1];

  if(myrank == 0){
    PGFEM_printf("[1]||t(+) + t(-)|| = %12.5e :: [2]||t(+) - t(-)|| = %12.5e\n"
                 "2*[1]/[2] = %12.5e\n",norms[0],norms[1],norms[2]);
    PGFEM_printf("Resulting tractions on bounding elements: ");
    print_array_d(PGFEM_stdout,res_trac,3,1,3);
  }

  free(L_trac);
  free(L_trac_sgn);
  free(G_trac);
  free(G_trac_sgn);
  free(norms);

  return err;
}

/** Compute the permutation to transform the global stiffness into a
    block-ordered matrix for preconditioning. NOTE THIS WILL PROBABLY
    BOTTLENECK FOR LARGE SYSTEMS */
int compute_block_permutation(const int n_be,
                              const int ndofn,
                              const BoundingElement *b_elems,
                              const Element *elems,
                              const Node *nodes,
                              const long *DomDof,
			      const CommunicationStructure *com,
                              long *perm,
                              const int mp_id)
{
  int err = 0;
  int myrank = com->rank;
  int nproc = com->nproc;
  
  /* compute bounds of owned rows. These are stored as global
     variables, but I am trying to phase them out... */
  int l_bound = 0;
  int u_bound = 0;
  for(int i=0; i<myrank; i++){
    l_bound += DomDof[i];
  }
  u_bound = l_bound + DomDof[myrank] - 1;

  int n_dof = 0;
  int *recvcount = aloc1i(nproc);
  int *displs = aloc1i(nproc);
  for(int i=0; i<nproc; i++){
    n_dof += DomDof[i];
    recvcount[i] = 2*DomDof[i];
    if(i>0){
      displs[i] = displs[i-1] + recvcount[i-1];
    }
  }

  /* mark the global dofs associated with constraints */
  int n_bnd_u_dofs = 0;
  int n_lm_dofs = 0;
  int *is_lm_dof = aloc1i(n_dof);
  int *dof_ids = aloc1i(2*n_dof); /* contains new_id, orig_id */
  for(int i=0; i<n_be; i++){
    const BoundingElement *p_be = &b_elems[i];
    const Element *p_ve = &elems[p_be->vol_elem_id];
    const int ndofe = get_ndof_on_bnd_elem(nodes,p_be,elems,ndofn);
    const int ndof_ve = get_ndof_on_elem_nodes(p_ve->toe,p_ve->nod,nodes,ndofn);
    long *dof = aloc1l(ndofe);
    get_dof_ids_on_bnd_elem(1,ndofn,nodes,p_be,elems,dof,mp_id);
    for(int j=0; j<ndofe; j++){
      const int id = dof[j] - 1;
      if(id >= 0){

        /* count number of dofs on domain */
        if(l_bound <= id && id <= u_bound
           && !is_lm_dof[id] /* don't double-count */){
          if(j<ndof_ve) n_bnd_u_dofs++;
          else n_lm_dofs++;
        }

        /* set marker */
        if(j < ndof_ve){ /* disp dofs */
          is_lm_dof[id] = 0;
        } else { /* lm dofs */
          is_lm_dof[id] = 3;
        }
      }
    }
    free(dof);
  }

  for(int i=0; i<n_dof; i++){
    dof_ids[2*i] = i + is_lm_dof[i]*n_dof; /* new_id */
    dof_ids[2*i+1] = i; /* old id */
  }

  /* sort global array. Implement parallel sorting algorithm in the
     future, but for now we will just gather and sort locally */
  int *g_dof_ids = aloc1i(2*n_dof);
  com->net->allreduce(dof_ids,g_dof_ids,2*n_dof,NET_DT_INT,NET_OP_MAX,com->comm);
  if(myrank == 0){
    com->net->reduce(NET_IN_PLACE,&n_bnd_u_dofs,1,NET_DT_INT,NET_OP_SUM,0,com->comm);
    com->net->reduce(NET_IN_PLACE,&n_lm_dofs,1,NET_DT_INT,NET_OP_SUM,0,com->comm);
  } else {
    int tmp = 0;
    com->net->reduce(&n_bnd_u_dofs,&tmp,1,NET_DT_INT,NET_OP_SUM,0,com->comm);
    com->net->reduce(&n_lm_dofs,&tmp,1,NET_DT_INT,NET_OP_SUM,0,com->comm);
  }
  qsort(g_dof_ids,n_dof,2*sizeof(int),compare_int);

  /* store permutation in 'perm' */
  for(int i=0; i<DomDof[myrank]; i++){
    perm[i] = g_dof_ids[2*(i+l_bound)+1];
  }

  free(recvcount);
  free(displs);
  free(is_lm_dof);
  free(dof_ids);
  free(g_dof_ids);

  if(1 /* DEBUG */){
    if(myrank == 0){
      PGFEM_printf("n_bnd_u_dofs: %d n_lm_dofs: %d\n",n_bnd_u_dofs,n_lm_dofs);
    }
    char fname[50];
    sprintf(fname,"perm_%d.out",myrank);
    FILE *out = fopen(fname,"w");
    PGFEM_fprintf(out,"%d\t%d\n",l_bound,u_bound);
    for(int i=0; i<DomDof[myrank]; i++){
      PGFEM_fprintf(out,"%d\t%ld\n",i+l_bound,perm[i]);
    }
    fclose(out);
  }

  return err;
}/* compute_block_permutation() */
