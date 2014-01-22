#include <stdio.h>
#include <stdlib.h>

#include "domain_container.h"

void initialize_domain_container(int n, sort_container *sort, Domain_container *Domain)
{
  int i,j,conflict,loop_cnt;
  int *procs;

  procs = (int*) calloc (n,sizeof(int));

  conflict = loop_cnt = 0;
  while(!conflict){

#ifdef RENUM_DEBUG
    printf("Conflict loop %d\n",loop_cnt);
#endif /* RENUM_DEBUG */

    for(i=0; i<n; i++) /* Processor */
      {
	for(j=0;j<n;j++) /* Preference */
	  {
	    if(!Domain[sort[i*n+j].index].covered && procs[i] == 0)
	      {
		procs[i] += 1; /* mark proc i as used */
		Domain[sort[i*n+j].index].owner = i;
		Domain[sort[i*n+j].index].preference = j;
		Domain[sort[i*n+j].index].value = sort[i*n+j].value;
		Domain[sort[i*n+j].index].covered = 1; /* true */
		break;
	      }
	    else /* Conflict resolution */
	      {
		if(procs[i] == 0 && sort[i*n+j].value > Domain[sort[i*n+j].index].value)
		  {
		    procs[Domain[sort[i*n+j].index].owner] -= 1; /* mark original owner as unused */
		    procs[i] += 1; /* mark new owner as used */
		    Domain[sort[i*n+j].index].owner = i;
		    Domain[sort[i*n+j].index].preference = j;
		    Domain[sort[i*n+j].index].value = sort[i*n+j].value;
		  } else continue;
	      } /* if(!Domain[sort[i*n+j].index].covered) */
	  } /* Preference */
      } /* Processor */
    
    conflict = 1;
    for(i=0;i<n;i++)
      if(procs[i] != 1){
	conflict = 0;	
	break;
      }

    loop_cnt++;
  } /* while(!conflict) */

  /* Error checking */
  for(i=0;i<n;i++)
    {
      if(!Domain[i].covered) printf("Domain [%d] is uncovered!\n",i);
      if(procs[i] != 1) printf("Processor [%d] is either unused or overused. (%d)\n",i,procs[i]);
    }

  /* Testing */
#ifdef RENUM_DEBUG
  for(i=0;i<n;i++)
    printf("[%d]:: Owner: %d Preference: %d Value: %d Covered: %d\n",
	  i,Domain[i].owner,Domain[i].preference,Domain[i].value,Domain[i].covered);
#endif /* RENUM_DEBUG */
}
