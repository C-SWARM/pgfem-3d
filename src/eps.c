#include "eps.h"
#include "enumerations.h"
#include "allocation.h"
#include "elem3d.h"

#include <string.h>

static const int periodic = 0;

EPS* build_eps_il (const long ne,
		   const ELEMENT *elem,
		   const int analysis)
/*
  INELASTIC
*/
{
  EPS *pom;
  long i,j,k,M,N,P,II,JJ,nne;
  
  pom = (EPS*) PGFEM_calloc (ne, sizeof(EPS));
  
  for (i=0;i<ne;i++){
    /* initialize ALL variables */
    EPS *p_pom = &pom[i];

    /* el */
    p_pom->el.o = NULL;
    p_pom->el.f = NULL;
    p_pom->el.m = NULL;
    p_pom->el.i = NULL;
    p_pom->el.d = NULL;
    p_pom->el.eq = 0;
    p_pom->el.eq_m = 0;
    p_pom->el.eq_i = 0;

    /* IL0 */
    p_pom->il = NULL;

    /* IL1 */
    p_pom->d_il = NULL;

    /* IL2 */
    p_pom->st = NULL;

    p_pom->dam = NULL;
    p_pom->T = NULL;
    p_pom->d_T = NULL;
    p_pom->GD = 0;
    p_pom->load1 = 0;
    p_pom->load = 0;
    p_pom->eff = 0;
    p_pom->type = 0;
    p_pom->F = NULL;
    p_pom->Fn = NULL;
    p_pom->P = NULL;
    p_pom->S = NULL;
    p_pom->Fe = NULL;
    p_pom->Fp = NULL;
    p_pom->FB = NULL;
    
    nne = elem[i].toe;
    /* Integration */
    int_point (nne,&II); int_point (10,&JJ);

    if(analysis == MINI
       || analysis == MINI_3F){ /* linear plus bubble */
      int_point (5,&JJ);
    }
  
    switch(analysis){ /* can be cleaned up a bit */
    default:
      {
	pom[i].el.o = (double *) PGFEM_calloc (6,sizeof(double ));
	pom[i].pl.o = (double *) PGFEM_calloc (6,sizeof(double));
	pom[i].il = (IL0_eps *) PGFEM_calloc (II,sizeof(IL0_eps)); 
	pom[i].st = (IL2_eps *) PGFEM_calloc (JJ,sizeof(IL2_eps));

	/* volumetric damage structure */
	pom[i].dam = (damage*) PGFEM_calloc (II,sizeof(damage));
      
	/* Pressure integration part */
	if (analysis == STABILIZED
	    || analysis == MINI){
	  for (j=0;j<JJ;j++)
	    pom[i].st[j].Fpp = (double *) PGFEM_calloc (9,sizeof(double));
	}

	/* the following seems to be a bug so I have bracketed it in
	   an 'if' statement. Looking through all revisions, it has
	   always been this way. Only the 0th element is referenced
	   throughtout the code so I think the intention was only to
	   allocate the first element anyhow... MM 2/20/2013 */
	if(i == 0){
	  pom[0].Dp = (double **) PGFEM_calloc (3,sizeof(double*));
	  for (j=0;j<3;j++)
	    pom[0].Dp[j] = (double *) PGFEM_calloc (3,sizeof(double));
	}

	for (j=0;j<II;j++){
	  pom[i].il[j].o = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].il[j].F = (double *) PGFEM_calloc (9,sizeof(double));
	
	  if (periodic == 1) {
	    pom[i].il[j].Fe = (double *) PGFEM_calloc (9,sizeof(double));
	    pom[i].il[j].Fe1 = (double *) PGFEM_calloc (9,sizeof(double));
	  }
	
	  if (analysis == FS_CRPL){
	    pom[i].il[j].Fp = (double *) PGFEM_calloc (9,sizeof(double));
	    pom[i].il[j].UU = (double *) PGFEM_calloc (9,sizeof(double));
	  
	    pom[i].il[j].dUU_Tr   = (double **) PGFEM_calloc (3,sizeof(double*));
	    pom[i].il[j].dUU_Tr_n = (double **) PGFEM_calloc (3,sizeof(double*));
	    pom[i].il[j].dUU_Fr   = (double ****) PGFEM_calloc (3,sizeof(double***));
	    pom[i].il[j].dUU_Fr_n = (double ****) PGFEM_calloc (3,sizeof(double***));
	  
	    for (M=0;M<3;M++){
	    
	      pom[i].il[j].dUU_Tr[M]   = (double *) PGFEM_calloc (3,sizeof(double));
	      pom[i].il[j].dUU_Tr_n[M] = (double *) PGFEM_calloc (3,sizeof(double));
	      pom[i].il[j].dUU_Fr[M]   = (double ***) PGFEM_calloc (3,sizeof(double**));
	      pom[i].il[j].dUU_Fr_n[M] = (double ***) PGFEM_calloc (3,sizeof(double**));
	    
	      for (N=0;N<3;N++){
		pom[i].il[j].dUU_Fr[M][N]=
		  (double **) PGFEM_calloc (3,sizeof(double*));
		pom[i].il[j].dUU_Fr_n[M][N] = 
		  (double **) PGFEM_calloc (3,sizeof(double*));
	      
		for (P=0;P<3;P++){
		  pom[i].il[j].dUU_Fr[M][N][P] =
		    (double *) PGFEM_calloc (3,sizeof(double));
		  pom[i].il[j].dUU_Fr_n[M][N][P] =
		    (double *) PGFEM_calloc (3,sizeof(double));
		}
	      }
	    }
	  }/* end analysis == FS_CRPL */
	}/* j < II */
      
	if (periodic == 1 && i == 0){
	
	  pom[0].F  = (double **) PGFEM_calloc (3,sizeof(double*));
	  pom[0].Fn = (double **) PGFEM_calloc (3,sizeof(double*));
	  pom[0].FB = (double **) PGFEM_calloc (3,sizeof(double*));
	  pom[0].P  = (double **) PGFEM_calloc (3,sizeof(double*));
	  pom[0].S  = (double **) PGFEM_calloc (3,sizeof(double*));
	
	  pom[0].Fe  = (double **) PGFEM_calloc (3,sizeof(double*));
	  pom[0].Fp  = (double **) PGFEM_calloc (3,sizeof(double*));
	
	  for (k=0;k<3;k++){
	    pom[0].F[k]  = (double *) PGFEM_calloc (3,sizeof(double));
	    pom[0].Fn[k] = (double *) PGFEM_calloc (3,sizeof(double));
	    pom[0].FB[k] = (double *) PGFEM_calloc (3,sizeof(double));
	    pom[0].P[k]  = (double *) PGFEM_calloc (3,sizeof(double));
	    pom[0].S[k]  = (double *) PGFEM_calloc (3,sizeof(double));
	  
	    pom[0].Fe[k]  = (double *) PGFEM_calloc (3,sizeof(double));
	    pom[0].Fp[k]  = (double *) PGFEM_calloc (3,sizeof(double));
	  }
	}/*end periodic == 1 */
	break;
      }/* end default case */
    case ELASTIC:
    case TP_ELASTO_PLASTIC:
      pom[i].il    = (IL0_eps *) PGFEM_calloc (II,sizeof(IL0_eps));
      pom[i].el.o  = (double *) PGFEM_calloc (6,sizeof(double ));
      if (analysis == TP_ELASTO_PLASTIC){
	pom[i].d_il  = (IL1_eps *) PGFEM_calloc (II,sizeof(IL1_eps));

	pom[i].el.f  = (double *) PGFEM_calloc (6,sizeof(double ));
	pom[i].el.d  = (double *) PGFEM_calloc (6,sizeof(double ));
	pom[i].el.m  = (double *) PGFEM_calloc (6,sizeof(double ));
	pom[i].el.i  = (double *) PGFEM_calloc (6,sizeof(double ));
      }
      
      for (j=0;j<II;j++){
	pom[i].il[j].o    = (double *) PGFEM_calloc (6,sizeof(double));
	if (analysis == TP_ELASTO_PLASTIC){
	  pom[i].d_il[j].o  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].f    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].f  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].m    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].m  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].i    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].i  = (double *) PGFEM_calloc (6,sizeof(double));
	  
	  pom[i].il[j].d    = (double *) PGFEM_calloc (6,sizeof(double));
	  pom[i].d_il[j].d  = (double *) PGFEM_calloc (6,sizeof(double));
	}
      }
      break;
    } /* switch(analysis) */
  }/* end i < ne */

  return pom;
}

static void copy_IL0_eps(IL0_eps *restrict dest,
			 const IL0_eps *restrict src)
{

}

static void copy_IL1_eps(IL1_eps *restrict dest,
			 const IL1_eps *restrict src)
{

}

static void copy_IL2_eps(IL2_eps *restrict dest,
			 const IL2_eps *restrict src)
{

}

void copy_eps(EPS *restrict dest,
	      const EPS *restrict src,
	      const long ne,
	      const ELEMENT *elem,
	      const int analysis)
{

}

static void destroy_IL0_eps(IL0_eps *p_il0)
{
  free(p_il0->o);
  free(p_il0->f);
  free(p_il0->m);
  free(p_il0->d);
  free(p_il0->i);
  free(p_il0->F);
  free(p_il0->Fp);
  free(p_il0->UU);
  free(p_il0->Fe);
  free(p_il0->Fe1);
  free(p_il0->GA);
  free(p_il0->GA1);
  free(p_il0->PLC_B);

  /* multi-dim pointers */
  if(p_il0->dUU_Fr!=NULL) dealoc4(p_il0->dUU_Fr,3,3,3);
  if(p_il0->dUU_Fr_n!=NULL) dealoc4(p_il0->dUU_Fr_n,3,3,3);
  if(p_il0->dUU_Tr!=NULL) dealoc2(p_il0->dUU_Tr,3);
  if(p_il0->dUU_Tr_n!=NULL) dealoc2(p_il0->dUU_Tr_n,3);
}

static void destroy_IL1_eps(IL1_eps *p_il1)
{
  free(p_il1->o);
  free(p_il1->f);
  free(p_il1->m);
  free(p_il1->d);
  free(p_il1->i);
}

static void destroy_IL2_eps(IL2_eps *p_il2)
{
  free(p_il2->Fpp);
}

void destroy_eps_il(EPS* eps,
		    const ELEMENT *elem,
		    const long ne,
		    const int analysis)
{
  for(long i=0; i<ne; i++){
    const int nne = elem[i].toe;
    long n_ip = 0;
    long n_ipq = 0;

    /* number of integration points */
    int_point (nne,&n_ip);
    switch(analysis){
    case MINI: case MINI_3F: int_point (5,&n_ipq); break;
    default: int_point (10,&n_ipq); break;
    }

    EPS *p_eps = &eps[i];
    /* elastic */
    free(p_eps->el.o);
    free(p_eps->el.f);
    free(p_eps->el.m);
    free(p_eps->el.i);
    free(p_eps->el.d);

    /* plastic ? */
    free(p_eps->pl.o);

    /* inelastic */
    if(p_eps->il != NULL){
      for(int j=0; j<n_ip; j++){
	destroy_IL0_eps(&p_eps->il[j]);
      }
    }

    /* inelastic 2 */
    if(p_eps->d_il != NULL){
      for(int j=0; j<n_ip; j++){
	destroy_IL1_eps(&p_eps->d_il[j]);
      }
    }

    /* stabilized */
    if(p_eps->st != NULL){
      for(int j=0; j<n_ipq; j++){
	destroy_IL2_eps(&p_eps->st[j]);
      }
    }

    /* remaining */
    free(p_eps->il);
    free(p_eps->d_il);
    free(p_eps->st);
    free(p_eps->dam);
    free(p_eps->T);
    free(p_eps->d_T);
    if(p_eps->F != NULL) dealoc2(p_eps->F,3);
    if(p_eps->Fn != NULL) dealoc2(p_eps->Fn,3);
    if(p_eps->P != NULL) dealoc2(p_eps->P,3);
    if(p_eps->S != NULL) dealoc2(p_eps->S,3);
    if(p_eps->Fe != NULL) dealoc2(p_eps->Fe,3);
    if(p_eps->Fp != NULL) dealoc2(p_eps->Fp,3);
    if(p_eps->FB != NULL) dealoc2(p_eps->FB,3);
    if(p_eps->Dp != NULL) dealoc2(p_eps->Dp,3);

  }/* for each element */

  free(eps);
}
