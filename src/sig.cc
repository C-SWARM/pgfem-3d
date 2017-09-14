#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "sig.h"
#include "allocation.h"
#include "elem3d.h"
#include "enumerations.h"

SIG* build_sig_el (const long ne)
/*
  ELASTIC
*/
{
  SIG *pom = PGFEM_calloc (SIG, ne);

  for (int i=0;i<ne;i++){
    pom[i].el.o  = PGFEM_calloc (double, 6);
    pom[i].el.f  = PGFEM_calloc (double, 6);
    pom[i].el.d  = PGFEM_calloc (double, 6);
    pom[i].el.m  = PGFEM_calloc (double, 6);
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
                   Element *elem)
/*
  INELASTIC
*/
{
  SIG *pom;
  long i,j,II,nne;

  pom = PGFEM_calloc (SIG, ne);

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
      pom[i].el.o = PGFEM_calloc (double, 6);
      pom[i].il = PGFEM_calloc (IL0_sig, II);
      for (j=0;j<II;j++)
        pom[i].il[j].o = PGFEM_calloc (double, 6);
      break;

     default:
      pom[i].il = PGFEM_calloc (IL0_sig, II);
      pom[i].el.o = PGFEM_calloc (double, 6);

      if (analysis == TP_ELASTO_PLASTIC){
        pom[i].d_il = PGFEM_calloc (IL1_sig, II);
        pom[i].el.f = PGFEM_calloc (double, 6);
        pom[i].el.d = PGFEM_calloc (double, 6);
        pom[i].el.m = PGFEM_calloc (double, 6);
      }

      for (j=0;j<II;j++){
        pom[i].il[j].o = PGFEM_calloc (double, 6);
        if (analysis == TP_ELASTO_PLASTIC){
          pom[i].d_il[j].o = PGFEM_calloc (double, 6);
          pom[i].il[j].f   = PGFEM_calloc (double, 6);
          pom[i].d_il[j].f = PGFEM_calloc (double, 6);
          pom[i].il[j].m   = PGFEM_calloc (double, 6);
          pom[i].d_il[j].m = PGFEM_calloc (double, 6);
          pom[i].il[j].d   = PGFEM_calloc (double, 6);
          pom[i].d_il[j].d = PGFEM_calloc (double, 6);
          pom[i].il[j].a   = PGFEM_calloc (double, 6);
          pom[i].d_il[j].a = PGFEM_calloc (double, 6);
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
                    const Element *elem,
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
