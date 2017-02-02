#include "metis_printing.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif


void print_loc_reorder(char *name, int *order, int lndof, int offset)
{
  FILE *out;
  int i;
  out = PGFEM_fopen(name,"w");
  for(i=0;i<lndof;i++)
    PGFEM_fprintf(out,"%d::%d::%d\n",i,i+offset,order[i]);
  PGFEM_fclose(out);
}


void print_reorder(char *name, int *order, int gndof)
{
  FILE *out;
  int i;
  out=PGFEM_fopen(name,"w");
  for(i=0;i<gndof;i++)
    PGFEM_fprintf(out,"%d::%d\n",i,order[i]);
  PGFEM_fclose(out);
}


void print_sizes(char *name, int *sizes, int nproc)
{
  FILE *out;
  int i;
  out=PGFEM_fopen(name,"w");
  for(i = 0; i < 2*nproc-1; i++)
    PGFEM_fprintf(out,"%d\n",sizes[i]);
  PGFEM_fclose(out);
}
