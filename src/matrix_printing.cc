#include "matrix_printing.h"
#include <stdlib.h>

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

/*****************************************/
/*   MATRIX PRINTING (COLUMN & MATLAB)   */
/*****************************************/

void print_columns(char *name, int ncol, int *Ap, int *Ai, int start)
{
  int i,j;
  FILE *out;
  out = PGFEM_fopen(name,"w");
  {
    for(i=0;i<ncol;i++)
      {
	PGFEM_fprintf(out,"%d :: ",i+start);
	for(j=Ap[i];j<Ap[i+1];j++)
	  PGFEM_fprintf(out,"%d ",Ai[j]);
	PGFEM_fprintf(out,"\n");
      }
    PGFEM_fclose(out);
  }
}


void print_matrix(char *name, int ncol, int *Ap, int *Ai, int start)
{
  int i,j;
  FILE *out;
  out = PGFEM_fopen(name,"w");
  if(out == NULL)
    PGFEM_printf("ERROR: Unable to open file %s. File will not be created.\n",name);
  else
    {
      for(i=0;i<ncol;i++)
	{
	  for(j=Ap[i];j<Ap[i+1];j++)
	    PGFEM_fprintf(out,"%d %d %d\n",i+start+1,Ai[j]+1,1);
	}
      PGFEM_fclose(out);
    }
}


void print_rn_matrix(char *name, int ncol, int *Ap, int *Ai, int *order, int start)
{
  int i,j;
  FILE *out;
  out = PGFEM_fopen(name,"w");
  if(out == NULL)
    PGFEM_printf("ERROR: Unable to open file %s. File will not be created.\n",name);
  else
    {
      for(i=0;i<ncol;i++)
	{
	  for(j=Ap[i];j<Ap[i+1];j++)
	    PGFEM_fprintf(out,"%d %d %d\n",order[i+start]+1,order[Ai[j]]+1,1);
	}
      PGFEM_fclose(out);
    }
}


void print_CSR(char *name, int ncol, int *Ap, int *Ai)
{
  int i;
  FILE *out;
  out = PGFEM_fopen(name,"w");
  if (out == NULL)
    PGFEM_printf("ERROR: Unable to open file %s. File will not be created.\n",name);
  else
    {
      /* Print header (allocation info) */
      PGFEM_fprintf(out,"%d %d\n",ncol,Ap[ncol]);
      
      /* Print Ap */
      for(i=0;i<=ncol;i++)
	PGFEM_fprintf(out,"%d ",Ap[i]);
      PGFEM_fprintf(out,"\n");
      
      /* Print Ai */
      for(i=0;i<Ap[ncol];i++)
	PGFEM_fprintf(out,"%d ",Ai[i]);
      
      PGFEM_fclose(out);
    }
}

void printVector(double *vector, int length,char *filename)
{
  FILE *out;
  int i;
  out = PGFEM_fopen(filename,"w");

  for(i=0; i<length; i++) {
    PGFEM_fprintf(out,"%8.8f\n",vector[i]);
  }
  PGFEM_fclose(out);
}
