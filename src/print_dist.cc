#include "print_dist.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

void print_dist(char *name,int length, int *dist)
{
  int i;
  FILE *out;
  out = PGFEM_fopen(name,"w");

  if(out == NULL)
    PGFEM_printf("ERROR: Cannot open file %s. File will not be created.\n",name);
  else
    {
      PGFEM_fprintf(out,"%d\n",length);
      
      for(i=0;i<length;i++)
	PGFEM_fprintf(out,"%d ",dist[i]);
      PGFEM_fclose(out);
    }
}
