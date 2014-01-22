#include "sig.h"

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

SIG* build_sig_el (const long ne)
     /*
       ELASTIC      
     */
{
  SIG *pom;
  long i;
  
  pom = (SIG*) PGFEM_calloc (ne, sizeof(SIG));
  
  for (i=0;i<ne;i++){
    pom[i].el.o  = (double *) PGFEM_calloc (6,sizeof(double ));
    pom[i].el.f  = (double *) PGFEM_calloc (6,sizeof(double ));
    pom[i].el.d  = (double *) PGFEM_calloc (6,sizeof(double ));
    pom[i].el.m  = (double *) PGFEM_calloc (6,sizeof(double ));
  }
  
  if (pom == NULL){
    PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    abort ();
  }

  return (pom);
}

void destroy_sig_el(SIG* sig,
		    const long ne)
{
  for(long i=0; i<ne; i++){
    free(sig[i].el.o);
    free(sig[i].el.f);
    free(sig[i].el.d);
    free(sig[i].el.m);
    free(sig[i].p);
    free(sig[i].d_p);
  }
  free(sig);
}

SIG* build_sig_il (const long ne,
		   const int analysis,
		   ELEMENT *elem)
/*
  INELASTIC
*/
{
  SIG *pom;
  long i,j,II,nne;
  
  pom = (SIG*) PGFEM_calloc (ne,sizeof(SIG));
  
  for (i=0;i<ne;i++){
    
    nne = elem[i].toe;
    
    /* Integration */
    int_point (nne,&II);
    
    switch(analysis){ /* this can be cleaned up a bit more */
    case FS_CRPL:
    case FINITE_STRAIN:
    case STABILIZED:
    case MINI:
    case MINI_3F:
    case DISP:    
      pom[i].el.o  = (double *) PGFEM_calloc (6,sizeof(double ));
      pom[i].il   = (IL0_sig *) PGFEM_calloc (II,sizeof(IL0_sig));
      for (j=0;j<II;j++)
	pom[i].il[j].o    = (double *) PGFEM_calloc (6,sizeof(double));
      break;

    default:
      pom[i].il   = (IL0_sig *) PGFEM_calloc (II,sizeof(IL0_sig));
      pom[i].el.o  = (double *) PGFEM_calloc (6,sizeof(double ));
      
      if (analysis == TP_ELASTO_PLASTIC){
	pom[i].d_il = (IL1_sig *) PGFEM_calloc (II,sizeof(IL1_sig));
	
	pom[i].el.f  = (double *) PGFEM_calloc (6,sizeof(double ));
	pom[i].el.d  = (double *) PGFEM_calloc (6,sizeof(double ));
	pom[i].el.m  = (double *) PGFEM_calloc (6,sizeof(double ));
      }
      
      for (j=0;j<II;j++){
	pom[i].il[j].o    = (double *) PGFEM_calloc (6,sizeof(double));
	if (analysis == TP_ELASTO_PLASTIC){
	  pom[i].d_il[j].o  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].f    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].f  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].m    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].m  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].d    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].d  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].a    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].a  = (double *) PGFEM_calloc (6,sizeof(double));
	}
      }
      break;
    } /* switch(analysis) */
  }/* i < ne */
  
  if (pom == NULL){
    PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    abort ();
  }
 
  return pom;
}

void destroy_sig_il(SIG* sig,
		    const ELEMENT *elem,
		    const long ne,
		    const int analysis)
{
  if(elem == NULL){
    PGFEM_printf("Must destroy sig_il before element\n");
    abort();
  }

  long nip;

  for(long i=0; i<ne; i++){
    int_point (elem[i].toe,&nip);

    switch(analysis){
    default:
      for(long j=0; j<nip; j++){
	free(sig[i].il[j].o);
      }

      free(sig[i].el.o);
      free(sig[i].il);
      break;

    case ELASTIC:
    case TP_ELASTO_PLASTIC:
      for(long j=0; j<nip; j++){
	free(sig[i].il[j].o);
	if(analysis == TP_ELASTO_PLASTIC){ 
	  free(sig[i].d_il[j].o);

	  free(sig[i].il[j].f);
	  free(sig[i].d_il[j].f);

	  free(sig[i].il[j].m);
	  free(sig[i].d_il[j].m);

	  free(sig[i].il[j].d);
	  free(sig[i].d_il[j].d);

	  free(sig[i].il[j].a);
	  free(sig[i].d_il[j].a);
	}
      }

      if(analysis == TP_ELASTO_PLASTIC){
	free(sig[i].d_il);
	free(sig[i].el.f);
	free(sig[i].el.d);
	free(sig[i].el.m);
      }

      free(sig[i].el.o);
      free(sig[i].il);
      break;
    }/* switch(analysis) */

    free(sig[i].p);
    free(sig[i].d_p);
    free(sig[i].pn_1);
  } /* for each elem */
  free(sig);
}
