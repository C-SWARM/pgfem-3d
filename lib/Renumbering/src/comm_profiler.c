#include <stdio.h>
#include <stdlib.h>

#include "renumber_utils.h"
#include "comm_profiler.h"

void comm_profiler(int ndom, int *dist, int *new_ranks, sort_container *info)
{
  int i,j,num_send,num_recv,tot_send,tot_recv;
  char filename[20];
  FILE *out;

  sprintf(filename,"rn_comm.report");
  out = fopen(filename,"w");
    if(out == NULL){
      printf("Could not open %s. Report will not be created.\n",filename);
      return;
    }

  /* re-sort info to the original configuration (i.e. by index) */
    for(i=0;i<ndom;i++)
      qsort(&info[i*ndom],ndom,sizeof(sort_container),compare_sort_container_revert);

  fprintf(out,"Number of domains:\t%d\n",ndom);
  fprintf(out,"Number of columns:\t%d\n",dist[ndom]);

  fprintf(out,"\n############################################################\n\n");

  tot_send = tot_recv = 0;

  for(i=0;i<ndom;i++){
    fprintf(out,"Domain %d:\n",i);
    fprintf(out,"Number of columns: %d\n\n",dist[i+1]-dist[i]);

    fprintf(out,"Number of columns on each domain:\n");
    for(j=0;j<ndom;j++)
      fprintf(out,"\tDomain %d: %d\n",j,info[i*ndom+j].value);

    tot_send += num_send = dist[i+1]-dist[i] - info[i*ndom+new_ranks[i]].value;
    tot_recv += num_recv = dist[new_ranks[i]+1]-dist[new_ranks[i]]
      - info[i*ndom+new_ranks[i]].value;

    fprintf(out,"\nPreferred domain (optimized for system): %d\n",new_ranks[i]);
    fprintf(out,"Total number of rows to send:\t %d\n",num_send);
    fprintf(out,"Total number of rows to recieve: %d\n",num_recv);

    fprintf(out,"\n############################################################\n\n");
  } /* end for(i=0;i<ndom;i++) */

  fprintf(out,"Total communication statistics (All domains):\n");
  fprintf(out,"Rows sent:\t%d\n",tot_send);
  fprintf(out,"Rows received:\t%d\n",tot_recv);

  fclose(out);
}
