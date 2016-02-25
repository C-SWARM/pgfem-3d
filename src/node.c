/** This file defines some routines related to the node structure such
    as allocation, deallocation, reading and writing */
#include "node.h"
#include "allocation.h"
#include <assert.h>

NODE* build_node (const long nn,
		  const long ndofn)
     /*
       ne - Number of nodes
       ndofn - number of DOF of one node
     */
{
  NODE *pom = (NODE*) PGFEM_calloc (nn, sizeof(NODE));
  
  for (int ii=0;ii<nn;ii++){	 
    pom[ii].id = (long*) PGFEM_calloc (ndofn,sizeof(long));
    pom[ii].Gid = (long*) PGFEM_calloc (ndofn,sizeof(long));

    pom[ii].ndofn = ndofn;
  }
  
  return (pom);
}

void destroy_node(const long nn,
		  NODE* nod)
{
  for(long i=0; i<nn; i++){
    free(nod[i].id);
    free(nod[i].Gid);
  }
  free(nod);
}


long read_nodes (FILE *in,
		 const long nn,
		 NODE *node,
		 const int legacy,
		 MPI_Comm comm)
     /*
       in   - Input file
       nn   - Number of nodes
       node - Structure type of NODE
       
       returns: total number of nodes over all domains, counting nodes
       on boundries only once
       %%%%%%%%%%%%%%%% TESTED 6.12.99 %%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%% MODIFIED 7.20.05 %%%%%%%%%%%%%%%
     */
{
  int myrank = 0;
  MPI_Comm_rank(comm,&myrank);

  long Gtnn = 0;
  long tnn = nn;
  NODE *p_node = NULL;

  for (long i=0; i<nn; i++){
    {
      long id = 0;
      long Gnn = 0;
      long Dom = 0;
      fscanf (in,"%ld %ld %ld",&Gnn,&Dom,&id);
      p_node = &node[id];
      p_node->loc_id = id;
      p_node->Gnn = Gnn;
      p_node->Dom = Dom;
    }
    
    /* Input file error checking */
    if (p_node->Gnn < 0 && p_node->Dom != myrank){
      PGFEM_printerr("[%d] ERROR: incorrect node domain info (node %ld)!"
	      " %s:%s:%d\n",myrank,i,__func__,__FILE__,__LINE__);
      PGFEM_Abort();
    } 

    /* If we get a global node that doesn't live on this domain,
       subtract it from tnn */
    if ( p_node->Gnn != -1 && p_node->Dom != myrank ){
      tnn--;
    }

    if (legacy) {
      fscanf (in,"%lf %lf %lf %ld",&p_node->x1,&p_node->x2,&p_node->x3,
	      &p_node->pr);
    } else {
      fscanf (in,"%lf %lf %lf %d %d %ld",
	      &p_node->x1,&p_node->x2,&p_node->x3,
	      &p_node->model_type,&p_node->model_id,&p_node->pr);
    }

    p_node->x1_fd = p_node->x1;
    p_node->x2_fd = p_node->x2;
    p_node->x3_fd = p_node->x3;

    /* error check read */
    if(ferror(in)){
      PGFEM_printerr("[%d]ERROR:fscanf returned error"
	      " reading node %ld!\n",myrank,i);
      PGFEM_Abort();
    } else if(feof(in)){
      PGFEM_printerr("[%d]ERROR:prematurely reached end of input file!\n",
	      myrank);
      PGFEM_Abort();
    }
  }

  /* Gather tnn from all domains */
  MPI_Allreduce(&tnn,&Gtnn,1,MPI_LONG,MPI_SUM,comm);  

  return Gtnn;
}

void write_node_fname(const char *filename,
		      const int nnodes,
		      const NODE *nodes)
{
  FILE *ofile = fopen(filename,"w");
  if(ofile == NULL){
    PGFEM_printerr("ERROR: cannot open file %s in %s\n",filename,__func__);
    PGFEM_Abort();
  }

  write_node(ofile,nnodes,nodes);

  fclose(ofile);
}

void write_node(FILE *ofile,
		const int nnodes,
		const NODE *nodes)
{
  /* write header describing format */
  PGFEM_fprintf(ofile,"  Gnn DOM   Lnn       X            Y            Z"
	  "              X_fd         Y_fd         Z_fd         Lid::Gid ...\n");
  PGFEM_fprintf(ofile,"===================================================="
	  "================================================================\n");
  for(int i=0; i<nnodes; i++){
    const NODE *p_node = &nodes[i];
    PGFEM_fprintf(ofile,"%5ld %3ld %5d ",p_node->Gnn,p_node->Dom,i);
    PGFEM_fprintf(ofile,"%12.5e %12.5e %12.5e    %12.5e %12.5e %12.5e    ",
	    p_node->x1,p_node->x2,p_node->x3,
	    p_node->x1_fd,p_node->x2_fd,p_node->x3_fd);
    for(int j=0; j<p_node->ndofn; j++){
      PGFEM_fprintf(ofile,"%5ld::%-5ld ",p_node->id[j],p_node->Gid[j]);
    }
    PGFEM_fprintf(ofile,"\n");
  }
}

static int node_comp_loc_id(const void *a,
                            const void *b)
{
  return (((NODE*) a)->loc_id - ((NODE*) b)->loc_id);
}

static int node_comp_own(const void *a,
                         const void *b)
{
  return (((NODE*) a)->Dom - ((NODE*) b)->Dom);
}

static int node_comp_Gnn(const void *a,
                         const void *b)
{
  return (((NODE*) a)->Gnn - ((NODE*) b)->Gnn);
}

static int node_comp_own_Gnn(const void *a,
                             const void *b)
{
  int own = node_comp_own(a,b);
  if (own) return own;
  else return node_comp_Gnn(a,b);
}

static int node_comp_Gnn_loc(const void *a,
                             const void *b)
{
  int Gnn = node_comp_Gnn(a,b);
  if (Gnn) return Gnn;
  else return node_comp_loc_id(a,b);
}

static int node_comp_own_Gnn_loc(const void *a,
                                 const void *b)
{
  int own_gnn = node_comp_own_Gnn(a,b);
  if (own_gnn) return own_gnn;
  else return node_comp_loc_id(a,b);
}

void nodes_sort_loc_id(const int nnode,
                       NODE *nodes)
{
  qsort(nodes, nnode, sizeof(*nodes), node_comp_loc_id);
}

void nodes_sort_own_Gnn_loc(const int nnode,
                            NODE *nodes)
{
  qsort(nodes, nnode, sizeof(*nodes), node_comp_own_Gnn_loc);
}

/**
 * Compute the index range of Global Nodes owned by the specified domain.
 *
 * For valid results, requires the nodes to be sorted by
 * `nodes_sort_own_Gnn`. Index range may include duplicate global node
 * numbersin the case of periodicity.
 *
 * \return non-zero if no shared/global nodes are owned by the
 * specified domain. On success, `range` specifies matches in
 * [range[0], range[1]).
 */
int nodes_get_Gnn_idx_range(const int nnode,
                            const NODE *nodes,
                            const int dom,
                            int range[2])
{
  int err = 0;
  /* create a node for comparison */
  NODE comp_node = {0};
  comp_node.Dom = dom;

  /* search for a node with matching ownership */
  NODE *ptr_lb = bsearch(&comp_node, nodes, nnode,
                         sizeof(*nodes), node_comp_Gnn);

  /* exit early if no match found */
  if (!ptr_lb) return 1;

  /* linearly search for bounds */
  NODE *ptr_ub = ptr_lb;
  while (ptr_ub->Dom == dom) ++ptr_ub;
  if (ptr_lb->Gnn < 0) {
    /* matched node is purely local -> search forward for first owned
       boundary node */
    while (ptr_lb->Gnn < 0 && ptr_lb->Dom == dom) ++ptr_lb;
  } else {
    /* matched node is on the boundary -> search backward for first
       owned boundary node */
    while (ptr_lb->Gnn >= 0 && ptr_lb->Dom == dom) --ptr_lb;

    /* Lower-bound is inclusive -> increment pointer */
    ++ptr_lb;
  }

  assert(ptr_lb != ptr_ub);
  if (ptr_lb == ptr_ub) err++;

  /* perform pointer arithmetic w.r.t nodes to get the index range */
  range[0] = ptr_lb - nodes;
  range[1] = ptr_ub - nodes;

  return err;
}

int nodes_filter_shared_nodes(const int nnode,
                              NODE *nodes,
                              int *n_shared,
                              NODE **shared)
{
  int err = 0;

  /* sort by Gnn */
  qsort(nodes, nnode, sizeof(*nodes), node_comp_Gnn_loc);

  /* perform linear search from the end to find the beginning of the
     shared nodes list. We start from the end as typically there are
     fewer boundary nodes than local nodes. */
  NODE *ptr = &nodes[nnode-1];
  while (ptr->Gnn >= 0) --ptr;
  ++ptr; /* lower bound is inclusive, increment pointer */

  /* compute number of shared nodes and starting index */
  *n_shared = nnode - (ptr - nodes);
  *shared = ptr;
  assert(*n_shared >= 0); /* check for implementation error */

  /* re-sort shared nodes by own->Gnn->loc_id */
  qsort(ptr, *n_shared, sizeof(*ptr), node_comp_own_Gnn_loc);

  return err;
}
