// BuildConnectivity.c Builds nodal connectivity for the distributed mesh.
// Returns connectivity in CSR format.  Currently only linear tetras
// are suported (every node in element touches every other node)

void BuildConnectivity(int nNodes, int nElements, ELEMENT *elem, NODE *node, int **Ap, int **Ai, MPI_Comm comm)
{
  int nproc, myrank;
  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&myrank);

  // Get nodes in each element and sort by local ID number
  int **allNodes;
  allNodes = (int**) calloc (nElements,sizeof(int*));
  for(i=0; i<nElements; i++)
    allNodes[i] = (int*) calloc (4,sizeof(int));

  GetAllElementNodes(elem,nElements,allNodes);

}

void GetAllElementNodes(ELEMENT *elem, int nElements, int **allNodes)
{
  int i;
  for(i=0; i<nElements; i++){
    if(elem[i].toe != 4){
      PGFEM_printf("ERROR: only tetra elements are supported! Exiting.\n");
      abort();
    }
    elemnodes(i,4,allNodes[i],elem);
  }
}
