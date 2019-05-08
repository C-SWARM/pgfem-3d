#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "generate_dof_ids.h"
#include "GRedist_node.h"
#include "PGFEM_io.h"
#include "allocation.h"
#include "matice.h"
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cassert>

using namespace pgfem3d;
using namespace pgfem3d::net;

static constexpr int ndn = 3;

/*=== STATIC FUNCTIONS ===*/
static long generate_local_dof_ids_on_elem(const long nnode,
                                           const long ndofn,
                                           const long start_id,
                                           const long elem_id,
                                           long *visited_node_dof,
                                           Node *nodes,
                                           Element *elems,
                                           BoundingElement *b_elems,
                                           const CommunicationStructure *com,
                                           const int mp_id);

static long generate_local_dof_ids_on_coel(const long nnode,
                                           const long ndofn,
                                           const long start_id,
                                           const long coel_id,
                                           long *visited_node_dof,
                                           Node *nodes,
                                           COEL *coel,
                                           const CommunicationStructure *com,
                                           const int mp_id);

static long generate_global_dof_ids_on_elem(const long ndofn,
                                            const long start_id,
                                            const long elem_id,
                                            long *visited_node_dof,
                                            Node *nodes,
                                            Element *elems,
                                            BoundingElement *b_elems,
                                            const CommunicationStructure *com,
                                            const int mp_id);

static long generate_global_dof_ids_on_coel(const long ndofn,
                                            const long start_id,
                                            const long coel_id,
                                            long *visited_node_dof,
                                            Node *nodes,
                                            COEL *coel,
                                            const CommunicationStructure *com,
                                            const int mp_id);

static void distribute_global_dof_ids_on_bounding_elements
(const int n_belem,
 const int ndof_be,
 BoundingElement *b_elems,
 const CommunicationStructure *com);

/*=== API FUNCTIONS ===*/
long generate_local_dof_ids(const long nelem,
                            const long ncoel,
                            const long nnode,
                            const long ndofn,
                            Node *nodes,
                            Element *elems,
                            COEL *coel,
                            BoundingElement *b_elems,
                            const CommunicationStructure *com,
                            const int mp_id)
{
  int myrank = com->rank;

  /* set periodicity */
  for(long i=0; i<nnode; i++){
    Node *ptr_node = &nodes[ i ];
    ptr_node->Pr = -1;
    for(long j=0; j<i; j++){
      const Node *ptr_pnode = &nodes[ j ];
      if(i<=j || ptr_node->Gnn == -1 || ptr_pnode->Gnn == -1
         || ptr_node->Gnn != ptr_pnode->Gnn){
        continue;
      } else {
        ptr_node->Pr = j;
        break;
      }
    }
  }

  /* allocate boolean array for whether or not we have assigned the
     node dof id (index = node_id * ndofn + dof_idx */
  long *visited_node_dof = PGFEM_calloc(long, ndofn*nnode);

  long ndof = 1;
  for(long i=0; i<nelem; i++){
    ndof += generate_local_dof_ids_on_elem(nnode,ndofn,ndof,
                                           i,visited_node_dof,
                                           nodes,elems,b_elems,
                                           com,mp_id);
  }

  /* Nodes on cohesive elements may not be counted if the cohesive
     element coincides with the communication boundary. Note that we
     still need to provide numbers for ndofn on each node if we own
     the element. See similar note in Psparse_ApAi*/
  for(long i=0; i<ncoel; i++){
    ndof += generate_local_dof_ids_on_coel(nnode,ndofn,ndof,
                                           i,visited_node_dof,
                                           nodes,coel,com,mp_id);
  }

  /* some periodic nodes may not actually be a part of an element,
     make sure to number them correctly */
  int err = 0;
  for(long i=0; i<nnode; i++){
    for(long j=0; j<ndofn; j++){
      const long idx = i*ndofn+j;
      if(!visited_node_dof[idx]){
        /* periodic nodes have different pressure ids! */
        if(j<ndn){
          if(nodes[i].Pr == -1){
            PGFEM_printerr("[%d] ERROR: Hanging node (%d)!\n",myrank,i);
            err = 1;
            break;
          } else { /* node is periodic */
            const long pidx = nodes[i].Pr*ndofn+j;
            if(!visited_node_dof[pidx]){
              PGFEM_printerr("[%d] ERROR: master periodic node (%ld)"
                             " has not been numbered!\n",myrank,nodes[i].Pr);
              err = 1;
              break;
            } else {
              nodes[i].id_map[mp_id].id[j] = nodes[nodes[i].Pr].id_map[mp_id].id[j];
              visited_node_dof[idx] = 1;
            }
          }
        } else { /* pressure node */
          PGFEM_printerr("[%d] TESTING: Periodic pressure node (%d)\n",myrank,i);
          nodes[i].id_map[mp_id].id[j] = ndof;
          ndof ++;
        }
      }
    }
  }

#ifndef NDEBUG
  /* sanity check */
  for(long i=0; i<ndofn*nnode; i++){
    if(!visited_node_dof[i]){
      PGFEM_printerr("[%d]ERROR: did not visit all dofs!\n",myrank);
      PGFEM_Abort();
    }
  }
#endif

  PGFEM_free(visited_node_dof);
  if(err){
    PGFEM_Abort();
  }

  ndof --;
  return ndof;
}/* generate_local_dof_ids() */

long generate_global_dof_ids(const long nelem,
                             const long ncoel,
                             const long nnode,
                             const long ndofn,
                             Node *nodes,
                             Element *elems,
                             COEL *coel,
                             BoundingElement *b_elems,
                             const CommunicationStructure *com,
                             const int mp_id)
{
  long *visited_node_dof = PGFEM_calloc(long, nnode*ndofn);
  long ndof = 1;

  for(long i=0; i<nelem; i++){
    ndof += generate_global_dof_ids_on_elem(ndofn,ndof,i,visited_node_dof,
                                            nodes,elems,b_elems,
                                            com,mp_id);
  }

  /* Nodes on cohesive elements may not be counted if the cohesive
     element coincides with the communication boundary. Note that we
     still need to provide numbers for ndofn on each node if we own
     the element. See similar note in Psparse_ApAi*/
  for(long i=0; i<ncoel; i++){
    ndof += generate_global_dof_ids_on_coel(ndofn,ndof,i,visited_node_dof,
                                            nodes,coel,com,mp_id);
  }

  /* sanity check */
#ifndef NDEBUG
  int myrank = com->rank;
  for(long i=0; i<ndofn*nnode; i++){
    if(!visited_node_dof[i]){
      PGFEM_printerr("[%d]ERROR: did not visit all dofs!\n",myrank);
      PGFEM_Abort();
    }
  }
#endif

  PGFEM_free(visited_node_dof);
  ndof --;

  /* Abort with message if no global dofs on the domain */
  assert(ndof > 0);
  return ndof;
}/* generate_global_dof_ids() */

void renumber_global_dof_ids(const long nelem,
                             const long ncoel,       /* UNUSED */
                             const int n_belem,
                             const long nnode,
                             const long ndofn,
                             const long *n_G_dof_on_dom,
                             Node *nodes,
                             Element *elems,
                             COEL *coel,
                             BoundingElement *b_elems,
                             const CommunicationStructure *com,
                             const int mp_id)
{
  int myrank = com->rank;

  /* compute how much to add to each Gdof id */
  long dof_id_adder = 0;
  for(int i=0; i<myrank; i++){
    dof_id_adder += n_G_dof_on_dom[i];
  }

  /* loop through each entity and increment Gdofs if I own them */
  /*=== ELEMENTS ===*/
  for(long i=0; i<nelem; i++){
    for(int j=0; j<elems[i].n_dofs; j++){
      elems[i].G_dof_ids[j] += dof_id_adder;
    }
  }

  /*=== COHESIVE ELEMENTS ===*/
  /* currently no local dofs on cohesive elements */

  /*=== BOUNDING ELEMENTS ===*/
  for(int i=0; i<n_belem; i++){
    if(!b_elems[i].periodic || b_elems[i].master_dom == myrank){
      for(int j=0; j<b_elems[i].n_dofs; j++){
        b_elems[i].G_dof_ids[j] += dof_id_adder;
      }
    }
  }

  /*=== NODES ===*/
  for(long i=0; i<nnode; i++){
    if(nodes[i].Dom == myrank){
      for(long j=0; j<ndofn; j++){
        if(nodes[i].id_map[mp_id].Gid[j] > 0){
          nodes[i].id_map[mp_id].Gid[j] += dof_id_adder;
        }
      }
    }
  }

} /* renumber_global_dof_ids */

long distribute_global_dof_ids(const long nelem,    /* UNUSED */
                               const long ncoel,    /* UNUSED */
                               const int n_belem,
                               const long nnode,
                               const long ndofn,
                               const int ndof_be,
                               Node *nodes,
                               Element *elems,
                               COEL *coel,
                               BoundingElement *b_elems,
                               const CommunicationStructure *com,
                               const int mp_id)
{
  int myrank = com->rank;
  int nproc = com->nproc;

  /* Distrubute the global dofs on the nodes and return the number of
     communication boundary nodes */
  long n_bnd_nodes = GRedist_node(nproc, myrank, nnode,
                                  ndofn, nodes, com, mp_id);

  /* Under current formulation, element dofs are always owned by the
     current domain and no communication is required. */

  distribute_global_dof_ids_on_bounding_elements(n_belem,ndof_be,
                                                 b_elems,com);

  return n_bnd_nodes;
} /* distribute_global_dof_ids */

/*=== STATIC FUNCTIONS ===*/
static long generate_local_dof_ids_on_elem(const long nnode,
                                           const long ndofn,
                                           const long start_id,
                                           const long elem_id,
                                           long *visited_node_dof,
                                           Node *nodes,
                                           Element *elems,
                                           BoundingElement *b_elems,
                                           const CommunicationStructure *com,
                                           const int mp_id)
{
  int myrank = com->rank;

  Element *ptr_elem = &elems[elem_id];
  const long *elem_nod = ptr_elem->nod;
  const long nne = ptr_elem->toe;
  long id_adder = 0;
  for(long i=0; i<nne; i++){
    const long node_id = elem_nod[i];
    Node *ptr_node = &nodes[node_id];

    for(long j=0; j<ndofn; j++){
      const long dof_idx = node_id*ndofn + j;

      /* set dof id only if the dof has not been previously visited */
      if(!visited_node_dof[dof_idx]){
        if(ptr_node->id_map[mp_id].id[j] < 0){ /* prescribed dof */
          visited_node_dof[dof_idx] = 1;
        } else if(ptr_node->id_map[mp_id].id[j] == 1){ /* fixed dof */
          ptr_node->id_map[mp_id].id[j] = 0;
          visited_node_dof[dof_idx] = 1;
        } else if(ptr_node->id_map[mp_id].id[j] >= 0){
          /* periodic node/dofs. Do not set pressure periodic! */
          if(ptr_node->Pr != -1 && j < ndn){
            const long pdof_idx = ptr_node->Pr*ndofn + j;
            /* if the dof id has not been set yet on the master node,
               set it now */
            if(!visited_node_dof[pdof_idx]){
              nodes[ptr_node->Pr].id_map[mp_id].id[j] = start_id + id_adder;
              visited_node_dof[pdof_idx] = 1;
              id_adder ++;
            }
            ptr_node->id_map[mp_id].id[j] = nodes[ptr_node->Pr].id_map[mp_id].id[j];
            visited_node_dof[dof_idx] = 1;
          } else {
            ptr_node->id_map[mp_id].id[j] = start_id + id_adder;
            id_adder ++;
            visited_node_dof[dof_idx] = 1;
          }
        }
      } else { /* have already set the current dof id */
        continue;
      }
    } /* for each dof on the node */
  } /* for each node on the element */

  /* Assign element dof ids (if there are any) */
  for(int i=0; i<ptr_elem->n_dofs; i++){
    ptr_elem->L_dof_ids[i] = start_id + id_adder;
    id_adder ++;
  }

  /* Assign dof ids on boundary elements. Currently unassigned if == 0 */
  for(int i=0; i<ptr_elem->n_be; i++){
    const long be_id = ptr_elem->be_ids[i];
    BoundingElement *ptr_be = &b_elems[be_id];
    if(!ptr_be->periodic){
      for(int j=0; j<ptr_be->n_dofs; j++){
        if(ptr_be->L_dof_ids[j] == 0){
          ptr_be->L_dof_ids[j] = start_id + id_adder;
          id_adder ++;
        }
      }
    } else { /* periodic element */
      if( ptr_be->master_dom != myrank
          || (ptr_be->master_dom == myrank
              && ptr_be->master_be_id == be_id) ){
        for(int j=0; j<ptr_be->n_dofs; j++){
          if(ptr_be->L_dof_ids[j] == 0){
            ptr_be->L_dof_ids[j] = start_id + id_adder;
            id_adder ++;
          }
        }
      } else {
        /* master is on this domain but is not this element. */
        BoundingElement *ptr_mbe = &b_elems[ptr_be->master_be_id];
        if(ptr_mbe->n_dofs != ptr_be->n_dofs){
          PGFEM_printerr("ERROR: different number of dofs on periodic"
                         " bounding elements!\n");
          PGFEM_Abort();
        }
        /* Number the master dofs first */
        for(int j=0; j<ptr_mbe->n_dofs; j++){
          if(ptr_mbe->L_dof_ids[j] == 0){
            ptr_mbe->L_dof_ids[j] = start_id + id_adder;
            id_adder ++;
          }
        }
        memcpy(ptr_be->L_dof_ids,ptr_mbe->L_dof_ids,
               ptr_mbe->n_dofs*sizeof(long));

      }
    }/* periodic element */
  }/* each bounding element on vol elem */

  return id_adder;
}/* generate_local_dof_ids_on_elem() */

static long generate_local_dof_ids_on_coel(const long nnode,
                                           const long ndofn,
                                           const long start_id,
                                           const long coel_id,
                                           long *visited_node_dof,
                                           Node *nodes,
                                           COEL *coel,
                                           const CommunicationStructure *com,
                                           const int mp_id)
{
  COEL *ptr_coel = &coel[coel_id];
  const long *elem_nod = ptr_coel->nod;
  const long nne = ptr_coel->toe;

  long id_adder = 0;
  for(long i=0; i<nne; i++){
    const long node_id = elem_nod[i];
    Node *ptr_node = &nodes[node_id];
    for(long j=0; j<ndofn; j++){
      const long dof_idx = node_id*ndofn + j;

      /* set dof id only if the dof has not been previously visited */
      if(!visited_node_dof[dof_idx]){
        if(ptr_node->id_map[mp_id].id[j] < 0){ /* prescribed dof */
          visited_node_dof[dof_idx] = 1;
        } else if(ptr_node->id_map[mp_id].id[j] == 1){ /* fixed dof */
          ptr_node->id_map[mp_id].id[j] = 0;
          visited_node_dof[dof_idx] = 1;
        } else if(ptr_node->id_map[mp_id].id[j] >= 0){
          /* periodic node/dofs. Do not set pressure periodic! */
          if(ptr_node->Pr != -1 && j < ndn){
            const long pdof_idx = ptr_node->Pr*ndofn + j;
            /* if the dof id has not been set yet on the master node,
               set it now */
            if(!visited_node_dof[pdof_idx]){
              nodes[ptr_node->Pr].id_map[mp_id].id[j] = start_id + id_adder;
              visited_node_dof[pdof_idx] = 1;
              id_adder ++;
            }
            ptr_node->id_map[mp_id].id[j] = nodes[ptr_node->Pr].id_map[mp_id].id[j];
            visited_node_dof[dof_idx] = 1;
          } else {
            ptr_node->id_map[mp_id].id[j] = start_id + id_adder;
            id_adder ++;
            visited_node_dof[dof_idx] = 1;
          }
        }
      } else { /* have already set the current dof id */
        continue;
      }
    } /* for each dof on the node */
  } /* for each node on the element */

  return id_adder;
} /* generate_local_dof_ids_on_coel */

static long generate_global_dof_ids_on_elem(const long ndofn,
                                            const long start_id,
                                            const long elem_id,
                                            long *visited_node_dof,
                                            Node *nodes,
                                            Element *elems,
                                            BoundingElement *b_elems,
                                            const CommunicationStructure *com,
                                            const int mp_id)
{
  int myrank = com->rank;

  long id_adder = 0;
  Element *ptr_elem = &elems[elem_id];
  const long *elem_nod = ptr_elem->nod;
  const long nne = ptr_elem->toe;

  for(long i=0; i<nne; i++){
    const long node_id = elem_nod[i];
    Node *ptr_node = &nodes[node_id];
    if(ptr_node->Dom != myrank){
      for(long j=0; j<ndofn; j++){
        const long dof_idx = node_id*ndofn + j;
        visited_node_dof[dof_idx] = 1;
      }
      continue;
    } else {
      for(long j=0; j<ndofn; j++){
        const long dof_idx = node_id*ndofn + j;
        if(!visited_node_dof[dof_idx]){
          if(ptr_node->id_map[mp_id].id[j] <= 0){ /* prescribed or supported dof */
            ptr_node->id_map[mp_id].Gid[j] = ptr_node->id_map[mp_id].id[j];
            visited_node_dof[dof_idx] = 1;
          }
          /* periodic node/dofs. Do not set pressure periodic! */
          else if(ptr_node->Pr != -1 && j < ndn){
            Node *ptr_pnode = &nodes[ptr_node->Pr];
            const long pdof_idx = ptr_node->Pr*ndofn + j;
            if(!visited_node_dof[pdof_idx]){ /* have not numbered master node yet */
              ptr_pnode->id_map[mp_id].Gid[j] = start_id + id_adder;
              visited_node_dof[pdof_idx] = 1;
              id_adder ++;
            }
            ptr_node->id_map[mp_id].Gid[j] = ptr_pnode->id_map[mp_id].Gid[j];
            visited_node_dof[dof_idx] = 1;
          } else { /* all other nodes */
            ptr_node->id_map[mp_id].Gid[j] = start_id + id_adder;
            visited_node_dof[dof_idx] = 1;
            id_adder ++;
          }
        }
      } /* for each dof on node */
    }/* node is owned by current domain */
  }/* for each node on element */

  /* number element-wise dof(s) */
  for(int i=0; i<ptr_elem->n_dofs; i++){
    ptr_elem->G_dof_ids[i] = start_id + id_adder;
    id_adder ++;
  }

  /* number bounding element dofs */
  for(int i=0; i<ptr_elem->n_be; i++){
    const long be_id = ptr_elem->be_ids[i];
    BoundingElement *ptr_be = &b_elems[ be_id ];
    if(!ptr_be->periodic){
      for(int j=0; j<ptr_be->n_dofs; j++){
        if(ptr_be->G_dof_ids[j] == 0){
          ptr_be->G_dof_ids[j] = start_id + id_adder;
          id_adder ++;
        }
      }
    } else { /* periodic bounding element */
      if( ptr_be->master_dom != myrank ){
        /* master on different domain, do not number G_id */
        continue;
      } else if (ptr_be->master_dom == myrank && ptr_be->master_be_id == be_id){
        /* master is current bounding element */
        for(int j=0; j<ptr_be->n_dofs; j++){
          if(ptr_be->G_dof_ids[j] == 0){
            ptr_be->G_dof_ids[j] = start_id + id_adder;
            id_adder ++;
          }
        }
      } else {
        /* master is current domain, but different element. Number
           master element and copy G_dof_ids to current element */
        BoundingElement *ptr_mbe = &b_elems[ptr_be->master_be_id];
        if(ptr_mbe->n_dofs != ptr_be->n_dofs){
          PGFEM_printerr("ERROR: non-matching n_dofs on bnd elems in %s\n",__func__);
          PGFEM_Abort();
        }
        for(int j=0; j<ptr_mbe->n_dofs; j++){
          if(ptr_mbe->G_dof_ids[j] == 0){
            ptr_mbe->G_dof_ids[j] = start_id + id_adder;
            id_adder ++;
          }
        }
        memcpy(ptr_be->G_dof_ids,ptr_mbe->G_dof_ids,
               ptr_be->n_dofs*sizeof(long));
      }
    } /* periodic b_elem */
  } /* for each b_elem on vol_elem */

  return id_adder;
}/* generate_global_dof_ids_on_elem() */

static long generate_global_dof_ids_on_coel(const long ndofn,
                                            const long start_id,
                                            const long coel_id,
                                            long *visited_node_dof,
                                            Node *nodes,
                                            COEL *coel,
                                            const CommunicationStructure *com,
                                            const int mp_id)
{
  int myrank = com->rank;

  long id_adder = 0;
  COEL *ptr_coel = &coel[coel_id];
  const long *elem_nod = ptr_coel->nod;
  const long nne = ptr_coel->toe;

  for(long i=0; i<nne; i++){
    const long node_id = elem_nod[i];
    Node *ptr_node = &nodes[node_id];
    if(ptr_node->Dom != myrank){
      for(long j=0; j<ndofn; j++){
        const long dof_idx = node_id*ndofn + j;
        visited_node_dof[dof_idx] = 1;
      }
      continue;
    } else {
      for(long j=0; j<ndofn; j++){
        const long dof_idx = node_id*ndofn + j;
        if(!visited_node_dof[dof_idx]){
          if(ptr_node->id_map[mp_id].id[j] <= 0){ /* prescribed or supported dof */
            ptr_node->id_map[mp_id].Gid[j] = ptr_node->id_map[mp_id].id[j];
            visited_node_dof[dof_idx] = 1;
          }
          /* periodic node/dofs. Do not set pressure periodic! */
          else if(ptr_node->Pr != -1 && j < ndn){
            Node *ptr_pnode = &nodes[ptr_node->Pr];
            const long pdof_idx = ptr_node->Pr*ndofn + j;
            if(!visited_node_dof[pdof_idx]){ /* have not numbered master node yet */
              ptr_pnode->id_map[mp_id].Gid[j] = start_id + id_adder;
              visited_node_dof[pdof_idx] = 1;
              id_adder ++;
            }
            ptr_node->id_map[mp_id].Gid[j] = ptr_pnode->id_map[mp_id].Gid[j];
            visited_node_dof[dof_idx] = 1;
          } else { /* all other nodes */
            ptr_node->id_map[mp_id].Gid[j] = start_id + id_adder;
            visited_node_dof[dof_idx] = 1;
            id_adder ++;
          }
        }
      } /* for each dof on node */
    }/* node is owned by current domain */
  }/* for each node on element */

  return id_adder;
}/* generate_global_dof_ids_on_coel() */


static void distribute_global_dof_ids_on_bounding_elements(const int n_belem,
                                                           const int ndof_be,
                                                           BoundingElement *b_elems,
                                                           const CommunicationStructure *com)
{
  int myrank = com->rank;
  int nproc = com->nproc;

  /* get number of global elems/dofs */
  int n_Gbelem_on_dom = 0;
  for(int i=0; i<n_belem; i++){
    const BoundingElement *p_be = &b_elems[i];
    if(p_be->periodic /* non periodic are local to domain */
       && p_be->master_dom == myrank
       && p_be->master_be_id == i){
      n_Gbelem_on_dom += 1;
    }
  }

  /* allocate space */
  int *NG_bel = PGFEM_calloc(int, nproc);
  int n_Gdof_on_dom = n_Gbelem_on_dom*(ndof_be+1);
  long *Gdof_on_dom = NULL;
  if(n_Gdof_on_dom > 0){
    Gdof_on_dom = PGFEM_calloc(long, n_Gdof_on_dom);
  }

  /*
   * Populate Gdof_on_dom with be_id Gdof_id_0 ... Gdof_id_n.  NOTE:
   * Gdof_on_dom is sorted by be_id by construction and is unique with
   * regards to be_id
   */
  {
    int idx = 0;
    for(int i=0; i<n_belem; i++){
      const BoundingElement *p_be = &b_elems[i];
      if(p_be->periodic /* non periodic are local to domain */
         && p_be->master_dom == myrank
         && p_be->master_be_id == i){
        long *ptr = &Gdof_on_dom[idx*(ndof_be+1)];
	assert(ptr != NULL && "ptr can't be null");
        *ptr = i; /* store elem id */
        ptr += 1; /* increment pointer */
        memcpy(ptr,b_elems[i].G_dof_ids,ndof_be*sizeof(long));
        idx++;
      }
    }
  }

  /* compute number of global dofs from all domains */
  com->net->allgather(&n_Gbelem_on_dom,1,NET_DT_INT,NG_bel,1,NET_DT_INT,com->comm);
  int n_Gbelem_on_all_dom = 0;
  int *allgatherv_disp = PGFEM_calloc(int, nproc);
  int *allgatherv_count = PGFEM_calloc(int, nproc);
  for(int i=0; i<nproc; i++){
    allgatherv_disp[i] = n_Gbelem_on_all_dom*(ndof_be + 1);
    allgatherv_count[i] = NG_bel[i]*(ndof_be + 1);
    n_Gbelem_on_all_dom += NG_bel[i];
  }
  int n_Gdof_on_all_dom = n_Gbelem_on_all_dom*(ndof_be+1);
  long *Gdof_on_all_dom = NULL;
  if(n_Gdof_on_all_dom > 0){
    PGFEM_calloc(long, n_Gdof_on_all_dom);
    com->net->allgatherv(Gdof_on_dom,n_Gdof_on_dom,NET_DT_LONG,Gdof_on_all_dom,
			 allgatherv_count,allgatherv_disp,NET_DT_LONG,com->comm);
  }

  /* Now I have the global dofs from all domains. NOTE:
     Gdof_on_all_dom is sorted (on each domain) by be_id by
     construction */
  for(int i=0; i<n_belem; i++){
    BoundingElement *p_be = &b_elems[i];
    if(p_be->periodic && p_be->master_dom != myrank){
      const int master_dom = p_be->master_dom;
      const int master_be_id = p_be->master_be_id;
      const int n_elem = NG_bel[master_dom];
      const size_t size = (ndof_be+1)*sizeof(long);
      const long *ptr_base = Gdof_on_all_dom + allgatherv_disp[master_dom];
      const long *ptr_match = static_cast<const long*>(bsearch(&master_be_id,
                                                               ptr_base,
                                                               n_elem,
                                                               size,
                                                               compare_long));

      if(ptr_match == NULL){
        PGFEM_printerr("[%d]ERROR: did not find matching"
                       " master element in %s!\n",myrank,__func__);
        PGFEM_Abort();
      }

      /* ptr_match is incremented to ommit the master_be_id from the copy */
      memcpy(p_be->G_dof_ids,ptr_match+1,ndof_be*sizeof(long));
    }
  }

  PGFEM_free(NG_bel);
  PGFEM_free(Gdof_on_dom);
  PGFEM_free(Gdof_on_all_dom);
  PGFEM_free(allgatherv_disp);
  PGFEM_free(allgatherv_count);
}
