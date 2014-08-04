#include "eps.h"
#include "enumerations.h"
#include "allocation.h"
#include "utils.h"
#include "elem3d.h"

#include <string.h>

#define cast_V4 (void****)
#define cast_const_V4 (const void****)
#define cast_V2 (void**)
#define cast_const_V2 (const void**)

static const int periodic = 0;
static const size_t SYM_TENSOR = 6;
static const size_t NDN = 3;
static const size_t TENSOR_2 = 9; /* NDN*NDN */

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
	pom[i].el.o = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double ));
	pom[i].pl.o = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	pom[i].il = (IL0_eps *) PGFEM_calloc (II,sizeof(IL0_eps)); 
	pom[i].st = (IL2_eps *) PGFEM_calloc (JJ,sizeof(IL2_eps));

	/* volumetric damage structure */
	pom[i].dam = (damage*) PGFEM_calloc (II,sizeof(damage));
      
	/* Pressure integration part */
	if (analysis == STABILIZED
	    || analysis == MINI){
	  for (j=0;j<JJ;j++)
	    pom[i].st[j].Fpp = (double *) PGFEM_calloc (TENSOR_2,sizeof(double));
	}

	/* the following seems to be a bug so I have bracketed it in
	   an 'if' statement. Looking through all revisions, it has
	   always been this way. Only the 0th element is referenced
	   throughtout the code so I think the intention was only to
	   allocate the first element anyhow... MM 2/20/2013 */
	if(i == 0){
	  pom[0].Dp = (double **) PGFEM_calloc (NDN,sizeof(double*));
	  for (j=0;j<NDN;j++)
	    pom[0].Dp[j] = (double *) PGFEM_calloc (NDN,sizeof(double));
	}

	for (j=0;j<II;j++){
	  pom[i].il[j].o = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  pom[i].il[j].F = (double *) PGFEM_calloc (TENSOR_2,sizeof(double));
	
	  if (periodic == 1) {
	    pom[i].il[j].Fe = (double *) PGFEM_calloc (TENSOR_2,sizeof(double));
	    pom[i].il[j].Fe1 = (double *) PGFEM_calloc (TENSOR_2,sizeof(double));
	  }
	
	  if (analysis == FS_CRPL){
	    pom[i].il[j].Fp = (double *) PGFEM_calloc (TENSOR_2,sizeof(double));
	    pom[i].il[j].UU = (double *) PGFEM_calloc (TENSOR_2,sizeof(double));
	  
	    pom[i].il[j].dUU_Tr   = (double **) PGFEM_calloc (NDN,sizeof(double*));
	    pom[i].il[j].dUU_Tr_n = (double **) PGFEM_calloc (NDN,sizeof(double*));
	    pom[i].il[j].dUU_Fr   = (double ****) PGFEM_calloc (NDN,sizeof(double***));
	    pom[i].il[j].dUU_Fr_n = (double ****) PGFEM_calloc (NDN,sizeof(double***));
	  
	    for (M=0;M<NDN;M++){
	    
	      pom[i].il[j].dUU_Tr[M]   = (double *) PGFEM_calloc (NDN,sizeof(double));
	      pom[i].il[j].dUU_Tr_n[M] = (double *) PGFEM_calloc (NDN,sizeof(double));
	      pom[i].il[j].dUU_Fr[M]   = (double ***) PGFEM_calloc (NDN,sizeof(double**));
	      pom[i].il[j].dUU_Fr_n[M] = (double ***) PGFEM_calloc (NDN,sizeof(double**));
	    
	      for (N=0;N<NDN;N++){
		pom[i].il[j].dUU_Fr[M][N]=
		  (double **) PGFEM_calloc (NDN,sizeof(double*));
		pom[i].il[j].dUU_Fr_n[M][N] = 
		  (double **) PGFEM_calloc (NDN,sizeof(double*));
	      
		for (P=0;P<NDN;P++){
		  pom[i].il[j].dUU_Fr[M][N][P] =
		    (double *) PGFEM_calloc (NDN,sizeof(double));
		  pom[i].il[j].dUU_Fr_n[M][N][P] =
		    (double *) PGFEM_calloc (NDN,sizeof(double));
		}
	      }
	    }
	  }/* end analysis == FS_CRPL */
	}/* j < II */
      
	if (periodic == 1 && i == 0){
	
	  pom[0].F  = (double **) PGFEM_calloc (NDN,sizeof(double*));
	  pom[0].Fn = (double **) PGFEM_calloc (NDN,sizeof(double*));
	  pom[0].FB = (double **) PGFEM_calloc (NDN,sizeof(double*));
	  pom[0].P  = (double **) PGFEM_calloc (NDN,sizeof(double*));
	  pom[0].S  = (double **) PGFEM_calloc (NDN,sizeof(double*));
	
	  pom[0].Fe  = (double **) PGFEM_calloc (NDN,sizeof(double*));
	  pom[0].Fp  = (double **) PGFEM_calloc (NDN,sizeof(double*));
	
	  for (k=0;k<NDN;k++){
	    pom[0].F[k]  = (double *) PGFEM_calloc (NDN,sizeof(double));
	    pom[0].Fn[k] = (double *) PGFEM_calloc (NDN,sizeof(double));
	    pom[0].FB[k] = (double *) PGFEM_calloc (NDN,sizeof(double));
	    pom[0].P[k]  = (double *) PGFEM_calloc (NDN,sizeof(double));
	    pom[0].S[k]  = (double *) PGFEM_calloc (NDN,sizeof(double));
	  
	    pom[0].Fe[k]  = (double *) PGFEM_calloc (NDN,sizeof(double));
	    pom[0].Fp[k]  = (double *) PGFEM_calloc (NDN,sizeof(double));
	  }
	}/*end periodic == 1 */
	break;
      }/* end default case */
    case ELASTIC:
    case TP_ELASTO_PLASTIC:
      pom[i].il    = (IL0_eps *) PGFEM_calloc (II,sizeof(IL0_eps));
      pom[i].el.o  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double ));
      if (analysis == TP_ELASTO_PLASTIC){
	pom[i].d_il  = (IL1_eps *) PGFEM_calloc (II,sizeof(IL1_eps));

	pom[i].el.f  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double ));
	pom[i].el.d  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double ));
	pom[i].el.m  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double ));
	pom[i].el.i  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double ));
      }
      
      for (j=0;j<II;j++){
	pom[i].il[j].o    = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	if (analysis == TP_ELASTO_PLASTIC){
	  pom[i].d_il[j].o  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  
	  pom[i].il[j].f    = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  pom[i].d_il[j].f  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  
	  pom[i].il[j].m    = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  pom[i].d_il[j].m  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  
	  pom[i].il[j].i    = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  pom[i].d_il[j].i  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  
	  pom[i].il[j].d    = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	  pom[i].d_il[j].d  = (double *) PGFEM_calloc (SYM_TENSOR,sizeof(double));
	}
      }
      break;
    } /* switch(analysis) */
  }/* end i < ne */

  return pom;
}

static void copy_IL0_eps(IL0_eps *dest,
			 const IL0_eps *src,
			 const int analysis)
{
  static const size_t d = sizeof(double);
  memcpy(dest->o,src->o,SYM_TENSOR*d);
  memcpy(dest->F,src->F,TENSOR_2*d);

  dest->Un_1 = src->Un_1;
  dest->Un = src->Un;
  dest->Jn_1 = src->Jn_1;
  dest->Y = src->Y;

  if(periodic == 1){
    memcpy(dest->Fe,src->Fe,TENSOR_2*d);
    memcpy(dest->Fe1,src->Fe1,TENSOR_2*d);
  }

  if(analysis == TP_ELASTO_PLASTIC){
    memcpy(dest->f,src->f,SYM_TENSOR*d);
    memcpy(dest->m,src->m,SYM_TENSOR*d);
    memcpy(dest->d,src->d,SYM_TENSOR*d);
    memcpy(dest->i,src->i,SYM_TENSOR*d);
  }

  if (analysis == FS_CRPL){
    dest->lam = src->lam;
    dest->eff = src->eff;
    dest->GAMA = src->GAMA;

    memcpy(dest->Fp,src->Fp,TENSOR_2*d);
    memcpy(dest->UU,src->UU,TENSOR_2*d);

    copy_4mat(cast_V4 dest->dUU_Fr,
	      cast_const_V4 src->dUU_Fr,
	      NDN,NDN,NDN,NDN,d);
    copy_2mat(cast_V2 dest->dUU_Tr,
	      cast_const_V2 src->dUU_Tr,
	      NDN,NDN,d);
    copy_4mat(cast_V4 dest->dUU_Fr_n,
	      cast_const_V4 src->dUU_Fr_n,
	      NDN,NDN,NDN,NDN,d);
    copy_2mat(cast_V2 dest->dUU_Tr_n,
	      cast_const_V2 src->dUU_Tr_n,
	      NDN,NDN,d);

    /** UNUSED    
	dest->GA = NULL;
	dest->GA1 = NULL;
	dest->PLC_B = NULL;
    */
  }
}

static void copy_IL1_eps(IL1_eps *dest,
			 const IL1_eps *src,
			 const int analysis)
{
  static const size_t d = sizeof(double);
  if(analysis == TP_ELASTO_PLASTIC){
    memcpy(dest->o,src->o,SYM_TENSOR*d);
    memcpy(dest->f,src->f,SYM_TENSOR*d);
    memcpy(dest->m,src->m,SYM_TENSOR*d);
    memcpy(dest->d,src->d,SYM_TENSOR*d);
    memcpy(dest->i,src->i,SYM_TENSOR*d);
  }
}

static void copy_IL2_eps(IL2_eps *dest,
			 const IL2_eps *src,
			 const int analysis)
{
  static const size_t d = sizeof(double);

  copy_damage(&(dest->dam),&(src->dam));
  dest->Un = src->Un;
  dest->Un_1 = src->Un_1;
  dest->Jn_1 = src->Jn_1;

  if(analysis == STABILIZED
     || analysis == MINI){
    memcpy(dest->Fpp,src->Fpp,TENSOR_2*d);
  }
}

static void copy_eps_local(EPS *dest,
			   const EPS *src,
			   const int analysis)
{
  static const size_t d = sizeof(double);

  switch(analysis){
  default:
    memcpy(dest->el.o,src->el.o,SYM_TENSOR*d);
    memcpy(dest->pl.o,src->pl.o,SYM_TENSOR*d);
    memcpy(dest->pl.eq,src->pl.eq,2*d);
    /* skipping periodic stuff since deprecated (periodic == 0) */
    break;
  case ELASTIC:
    memcpy(dest->el.o,src->el.o,SYM_TENSOR*d);
    break;
  case TP_ELASTO_PLASTIC:
    memcpy(dest->el.o,src->el.o,SYM_TENSOR*d);
    memcpy(dest->el.f,src->el.f,SYM_TENSOR*d);
    memcpy(dest->el.d,src->el.d,SYM_TENSOR*d);
    memcpy(dest->el.m,src->el.m,SYM_TENSOR*d);
    memcpy(dest->el.i,src->el.i,SYM_TENSOR*d);
    break;
  }
}

/**
 * Copy one EPS object into another. if src == dest, then no copy is
 * performed. In this function, all pointers are to corresponding
 * entries in the list.
 */
static void copy_eps(EPS *dest,
		     const EPS *src,
		     const ELEMENT *elem,
		     const int analysis)
{
  static const size_t d = sizeof(double);
  if(dest == src) return;
  long pt_I = 0;
  long pt_J = 0;

  /* Get number of integration points */
  int_point(elem->toe,&pt_I);
  switch(analysis){
  default: int_point(10,&pt_J); break;
  case MINI: case MINI_3F: int_point(5,&pt_J); break;
  }

  /* begin copying. See build_eps_il */
  copy_eps_local(dest,src,analysis);

  for(int i=0; i<pt_I; i++){
    copy_damage(dest->dam+i,src->dam+i);
    copy_IL0_eps((dest->il + i),(src->il + i),analysis);
    copy_IL1_eps((dest->d_il + i),(src->d_il + i),analysis);
  }

  for(int j=0; j<pt_J; j++){
    copy_IL2_eps(dest->st+j,src->st+j,analysis);
  }
}

void copy_eps_list(EPS *dest,
		   const EPS *src,
		   const int ne,
		   const ELEMENT *elem,
		   const int analysis)
{
  static const size_t d = sizeof(double);

  switch(analysis){
  default:
    copy_2mat(cast_V2 dest[0].Dp,
	      cast_const_V2 src[0].Dp,
	      NDN,NDN,d);
    break;
  case ELASTIC: break;
  case TP_ELASTO_PLASTIC: break;
  }

  /* for each element */
  for(int n=0; n<ne; n++){
    copy_eps(dest+n,src+n,elem+n,analysis);
  }
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
  if(p_il0->dUU_Fr!=NULL) dealoc4(p_il0->dUU_Fr,NDN,NDN,NDN);
  if(p_il0->dUU_Fr_n!=NULL) dealoc4(p_il0->dUU_Fr_n,NDN,NDN,NDN);
  if(p_il0->dUU_Tr!=NULL) dealoc2(p_il0->dUU_Tr,NDN);
  if(p_il0->dUU_Tr_n!=NULL) dealoc2(p_il0->dUU_Tr_n,NDN);
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
    if(p_eps->F != NULL) dealoc2(p_eps->F,NDN);
    if(p_eps->Fn != NULL) dealoc2(p_eps->Fn,NDN);
    if(p_eps->P != NULL) dealoc2(p_eps->P,NDN);
    if(p_eps->S != NULL) dealoc2(p_eps->S,NDN);
    if(p_eps->Fe != NULL) dealoc2(p_eps->Fe,NDN);
    if(p_eps->Fp != NULL) dealoc2(p_eps->Fp,NDN);
    if(p_eps->FB != NULL) dealoc2(p_eps->FB,NDN);
    if(p_eps->Dp != NULL) dealoc2(p_eps->Dp,NDN);

  }/* for each element */

  free(eps);
}
