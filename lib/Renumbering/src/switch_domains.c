#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "renumber_utils.h"
#include "sort_container.h"
#include "domain_container.h"
#include "comm_profiler.h"
#include "switch_domains.h"

void switch_domains(int ncol, int *l_order, int *dist, MPI_Comm Old, MPI_Comm *New)
{
  int i, j, nproc, myrank, new_rank;

  sort_container *dom_it, *Domains;

  /* Determine number of procs and old rank */
  MPI_Comm_size(Old,&nproc);
  MPI_Comm_rank(Old,&myrank);

  /* Sort l_order */
  qsort(l_order,ncol,sizeof(int),compare_integer);

  /* Determine the number of entries of l_order which land within each domain. */
  dom_it = (sort_container*) calloc (nproc,sizeof(sort_container));

  j = 0;
  for(i=0;i<nproc;i++)
    {
      dom_it[i].value = 0;
      dom_it[i].index = i;
      for(j;j<ncol;j++)
	{
	  if(dist[i] <= l_order[j] && l_order[j] < dist[i+1])
	    dom_it[i].value++;
	  else
	    break;
	}
    }

  /* Determine the ranking of these domains */
  qsort(dom_it,nproc,sizeof(sort_container),compare_sort_container_reverse);

  /* Create a new MPI_Datatype */
  MPI_Datatype mpi_domain_info;
  MPI_Type_contiguous(2,MPI_INT,&mpi_domain_info);
  MPI_Type_commit(&mpi_domain_info);

  MPI_Request s_req, *r_req;
  MPI_Status s_status, *r_status;
  r_req = (MPI_Request*) calloc (nproc-1,sizeof(MPI_Request));
  r_status = (MPI_Status*) calloc (nproc-1,sizeof(MPI_Status));

  /* Send domain information to Proc 0 */
  if(myrank != 0) MPI_Isend(dom_it,nproc,mpi_domain_info,0,myrank,Old,&s_req);

  /* Proc 0 recieves information */
  if(myrank == 0){
    Domains = (sort_container*) calloc (nproc*nproc,sizeof(sort_container));
    memcpy(Domains,dom_it,nproc*sizeof(sort_container));
  }
  
  if(myrank == 0)
    for(i=1;i<nproc;i++)
      MPI_Irecv(&Domains[i*nproc],nproc,mpi_domain_info,i,MPI_ANY_TAG,Old,&r_req[i-1]);

  if(myrank != 0) MPI_Wait(&s_req,&s_status);
  else MPI_Waitall(nproc-1,r_req,r_status);

  /* Free up requests and statuses */
  free(r_req); free(r_status);

  /* Debugging print statment */
#ifdef RENUM_DEBUG
  if(myrank == 0)
    {
      printf("Collected preference vector (split by processor):\n");
      for(i=0;i<nproc;i++)
	{
	  for(j=0;j<nproc;j++)
	    printf("%d(%d)  ",Domains[i*nproc+j].index,Domains[i*nproc+j].value);
	  printf("\n");
	}
    }
#endif /* RENUM_DEBUG */

  /* Now we have all of the domain preferences and values on processor
     0. We will determine conflicts and uncovered domains, and then
     based on the preference order and number of elements on the
     domain, we will determine the new processor ranks. */

  int *new_ranks;
  if(myrank == 0)
    {
      Domain_container *Dom_cont;
      Dom_cont = (Domain_container*) calloc (nproc,sizeof(Domain_container));
      initialize_domain_container(nproc,Domains,Dom_cont);
      new_ranks = (int*) calloc (nproc,sizeof(int));

      for(i=0;i<nproc;i++)
	new_ranks[Dom_cont[i].owner] = i;

      /* Print renumbering profile */
#ifdef RENUM_DEBUG
      comm_profiler(nproc,dist,new_ranks,Domains);
#endif /* RENUM_DEBUG */

      free(Domains); free(Dom_cont);
    }

  MPI_Scatter(new_ranks,1,MPI_INT,&new_rank,1,MPI_INT,0,Old);
  MPI_Comm_split(Old,0,new_rank,New);

}/* End of switch_domains */
