#include "read_cryst_plast.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

void read_cryst_plast (FILE *in1,
		       long nmat,
		       CRPL *crpl,
		       const int plc)
{
  long i,j,k,nss;
  
  for (i=0;i<nmat;i++) {
    if (plc == 0)
      fscanf (in1,"%ld %lf %lf %lf %lf %lf %lf %lf",&crpl[i].nss,&crpl[i].a,&crpl[i].m,&crpl[i].TH,&crpl[i].To,&crpl[i].Tso,&crpl[i].Go,&crpl[i].mm);
    /* PLC */
    if (plc == 1)
      fscanf (in1,"%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&crpl[i].nss,&crpl[i].a,&crpl[i].m,&crpl[i].To,
	      &crpl[i].b,&crpl[i].fo,&crpl[i].nu,&crpl[i].l1,&crpl[i].l2,&crpl[i].c1,&crpl[i].c2,&crpl[i].c3,&crpl[i].to);
    
    nss = crpl[i].nss;
    
    crpl[i].P = (double**) PGFEM_calloc (nss,sizeof(double*));
    for (j=0;j<nss;j++) crpl[i].P[j] = (double*) PGFEM_calloc (6,sizeof(double));
    
    for (j=0;j<nss;j++) for (k=0;k<6;k++) fscanf (in1,"%lf",&crpl[i].P[j][k]);
  }
}
