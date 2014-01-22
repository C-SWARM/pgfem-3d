/** This file defines some routines related to the node structure such
    as allocation, deallocation, reading and writing */
#include "node.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

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
