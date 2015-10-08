#include "eps.h"
#include "constitutive_model.h"
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
  int err = 0;
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

        /* Generalized constitutive modeling interface */
  if(analysis==CM)      
    pom[i].model = PGFEM_calloc(II,sizeof(*(pom[i].model)));
      
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
	  if(analysis==CM) 
          err += constitutive_model_construct(pom[i].model + j);
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

static size_t sizeof_IL0_eps(const int analysis)
{
  static IL0_eps src;
  size_t s = 0;
  s = (SYM_TENSOR*sizeof(*(src.o)) + TENSOR_2*sizeof(*(src.F))
       + sizeof(src.Un_1) + sizeof(src.Un) + sizeof(src.Jn_1)
       + sizeof(src.Y));

  if(analysis == TP_ELASTO_PLASTIC){
    s += SYM_TENSOR*(sizeof(*(src.f)) + sizeof(*(src.m))
		     + sizeof(*(src.d)) + sizeof(*(src.i)));
  }

  if (analysis == FS_CRPL){
    s += (sizeof(src.lam) + sizeof(src.eff) + sizeof(src.GAMA)
	  + TENSOR_2*(sizeof(*(src.Fp)) + sizeof(*(src.UU)))
	  + NDN*NDN*(sizeof(**(src.dUU_Tr)) + sizeof(**(src.dUU_Tr_n))
		     + NDN*NDN*(sizeof(****(src.dUU_Fr))
				+ sizeof(****(src.dUU_Fr_n))))
	  );
  }

  return s;
}

static size_t sizeof_IL1_eps(const int analysis)
{
  static IL1_eps src;
  size_t s = 0;
  if(analysis == TP_ELASTO_PLASTIC){
    s = SYM_TENSOR*(sizeof(*(src.o)) + sizeof(*(src.f))
		    + sizeof(*(src.m)) + sizeof(*(src.d))
		    + sizeof(*(src.i)));
  }
  return s;
}

static size_t sizeof_IL2_eps(const int analysis)
{
  static IL2_eps src;
  size_t s = 0;
  s = (sizeof_damage(&(src.dam)) + sizeof(src.Un_1)
       + sizeof(src.Un) + sizeof(src.Jn_1));

  if(analysis == STABILIZED
     || analysis == MINI){
    s += TENSOR_2*sizeof(*(src.Fpp));
  }
  return s;
}

static size_t sizeof_eps_local(const int analysis)
{
  static EPS src;
  size_t s = 0;
  switch(analysis){
  default:
    s = (SYM_TENSOR*(sizeof(*(src.el.o)) + sizeof(*(src.pl.o)))
	 + 2*sizeof(*(src.pl.eq)));
    /* skipping periodic stuff since deprecated (periodic == 0) */
    break;
  case ELASTIC:
    s = SYM_TENSOR*sizeof(*(src.el.o));
    break;
  case TP_ELASTO_PLASTIC:
    s = SYM_TENSOR*(sizeof(*(src.el.o)) + sizeof(*(src.el.f))
		    + sizeof(*(src.el.d)) + sizeof(*(src.el.m))
		    + sizeof(*(src.el.i)));
    break;
  }
  s += sizeof((src.el.eq)) + sizeof((src.el.eq_m)) + sizeof((src.el.eq_i));
  return s;
}

/** return the size of a particular EPS object */
static size_t sizeof_eps(const EPS *eps,
			 const ELEMENT *elem,
			 const int analysis)
{
  size_t s = 0;
  long pt_I = 0;
  long pt_J = 0;

  /* Get number of integration points */
  int_point(elem->toe,&pt_I);
  switch(analysis){
  default: int_point(10,&pt_J); break;
  case MINI: case MINI_3F: int_point(5,&pt_J); break;
  }

  s += sizeof_eps_local(analysis);
  if(pt_I > 0){
    s += pt_I * (sizeof_IL0_eps(analysis)
		 + sizeof_IL1_eps(analysis)
		 + sizeof_damage((eps->dam)));
  }

  if(pt_J > 0){
    s += pt_J * sizeof_IL2_eps(analysis);
  }

  return s;
}

size_t sizeof_eps_list(const EPS *src,
		       const int ne,
		       const ELEMENT *elem,
		       const int analysis)
{
  size_t s = 0;
  switch(analysis){
  default:
    s += NDN*NDN*sizeof(double);
    break;
  case ELASTIC: break;
  case TP_ELASTO_PLASTIC: break;
  }

  for(int i=0; i<ne; i++){
    s += sizeof_eps(src+i,elem+i,analysis);
  }
  return s;
}


static void pack_IL0_eps(const IL0_eps *src,
			 const int analysis,
			 char *buffer,
			 size_t *pos)
{
  pack_data(src->o,buffer,pos,SYM_TENSOR,sizeof(*(src->o)));
  pack_data(src->F,buffer,pos,TENSOR_2,sizeof(*(src->F)));
  pack_data(&(src->Un_1),buffer,pos,1,sizeof(src->Un_1));
  pack_data(&(src->Un),buffer,pos,1,sizeof(src->Un));
  pack_data(&(src->Jn_1),buffer,pos,1,sizeof(src->Jn_1));
  pack_data(&(src->Y),buffer,pos,1,sizeof(src->Y));

  if(analysis == TP_ELASTO_PLASTIC){
    pack_data(src->f,buffer,pos,SYM_TENSOR,sizeof(*(src->f)));
    pack_data(src->m,buffer,pos,SYM_TENSOR,sizeof(*(src->m)));
    pack_data(src->d,buffer,pos,SYM_TENSOR,sizeof(*(src->d)));
    pack_data(src->i,buffer,pos,SYM_TENSOR,sizeof(*(src->i)));
  }

  if (analysis == FS_CRPL){
    pack_data(&(src->lam),buffer,pos,1,sizeof(src->lam));
    pack_data(&(src->eff),buffer,pos,1,sizeof(src->eff));
    pack_data(&(src->GAMA),buffer,pos,1,sizeof(src->GAMA));
    pack_data(src->Fp,buffer,pos,TENSOR_2,sizeof(*(src->Fp)));
    pack_data(src->UU,buffer,pos,TENSOR_2,sizeof(*(src->UU)));
    pack_2mat(cast_const_V2 src->dUU_Tr,NDN,NDN,sizeof(**(src->dUU_Tr)),buffer,pos);
    pack_2mat(cast_const_V2 src->dUU_Tr_n,NDN,NDN,sizeof(**(src->dUU_Tr_n)),buffer,pos);
    pack_4mat(cast_const_V4 src->dUU_Fr,NDN,NDN,NDN,NDN,sizeof(****(src->dUU_Fr)),buffer,pos);
    pack_4mat(cast_const_V4 src->dUU_Fr_n,NDN,NDN,NDN,NDN,sizeof(****(src->dUU_Fr_n)),buffer,pos);
  }
}

static void unpack_IL0_eps(IL0_eps *dest,
			   const int analysis,
			   const char *buffer,
			   size_t *pos)
{
  unpack_data(buffer,dest->o,pos,SYM_TENSOR,sizeof(*(dest->o)));
  unpack_data(buffer,dest->F,pos,TENSOR_2,sizeof(*(dest->F)));
  unpack_data(buffer,&(dest->Un_1),pos,1,sizeof(dest->Un_1));
  unpack_data(buffer,&(dest->Un),pos,1,sizeof(dest->Un));
  unpack_data(buffer,&(dest->Jn_1),pos,1,sizeof(dest->Jn_1));
  unpack_data(buffer,&(dest->Y),pos,1,sizeof(dest->Y));

  if(analysis == TP_ELASTO_PLASTIC){
    unpack_data(buffer,dest->f,pos,SYM_TENSOR,sizeof(*(dest->f)));
    unpack_data(buffer,dest->m,pos,SYM_TENSOR,sizeof(*(dest->m)));
    unpack_data(buffer,dest->d,pos,SYM_TENSOR,sizeof(*(dest->d)));
    unpack_data(buffer,dest->i,pos,SYM_TENSOR,sizeof(*(dest->i)));
  }

  if (analysis == FS_CRPL){
    unpack_data(buffer,&(dest->lam),pos,1,sizeof(dest->lam));
    unpack_data(buffer,&(dest->eff),pos,1,sizeof(dest->eff));
    unpack_data(buffer,&(dest->GAMA),pos,1,sizeof(dest->GAMA));
    unpack_data(buffer,dest->Fp,pos,TENSOR_2,sizeof(*(dest->Fp)));
    unpack_data(buffer,dest->UU,pos,TENSOR_2,sizeof(*(dest->UU)));
    unpack_2mat(cast_V2 dest->dUU_Tr,NDN,NDN,sizeof(**(dest->dUU_Tr)),buffer,pos);
    unpack_2mat(cast_V2 dest->dUU_Tr_n,NDN,NDN,sizeof(**(dest->dUU_Tr_n)),buffer,pos);
    unpack_4mat(cast_V4 dest->dUU_Fr,NDN,NDN,NDN,NDN,sizeof(****(dest->dUU_Fr)),buffer,pos);
    unpack_4mat(cast_V4 dest->dUU_Fr_n,NDN,NDN,NDN,NDN,sizeof(****(dest->dUU_Fr_n)),buffer,pos);
  }
}

static void pack_IL1_eps(const IL1_eps *src,
			 const int analysis,
			 char *buffer,
			 size_t *pos)
{
  if(analysis == TP_ELASTO_PLASTIC){
    pack_data(src->o,buffer,pos,SYM_TENSOR,sizeof(*(src->o)));
    pack_data(src->f,buffer,pos,SYM_TENSOR,sizeof(*(src->f)));
    pack_data(src->m,buffer,pos,SYM_TENSOR,sizeof(*(src->m)));
    pack_data(src->d,buffer,pos,SYM_TENSOR,sizeof(*(src->d)));
    pack_data(src->i,buffer,pos,SYM_TENSOR,sizeof(*(src->i)));
  }
}

static void unpack_IL1_eps(IL1_eps *dest,
			   const int analysis,
			   const char *buffer,
			   size_t *pos)
{
  if(analysis == TP_ELASTO_PLASTIC){
    unpack_data(buffer,dest->o,pos,SYM_TENSOR,sizeof(*(dest->o)));
    unpack_data(buffer,dest->f,pos,SYM_TENSOR,sizeof(*(dest->f)));
    unpack_data(buffer,dest->m,pos,SYM_TENSOR,sizeof(*(dest->m)));
    unpack_data(buffer,dest->d,pos,SYM_TENSOR,sizeof(*(dest->d)));
    unpack_data(buffer,dest->i,pos,SYM_TENSOR,sizeof(*(dest->i)));
  }
}

static void pack_IL2_eps(const IL2_eps *src,
			 const int analysis,
			 char *buffer,
			 size_t *pos)
{
  pack_damage(&(src->dam),buffer,pos);
  pack_data(&(src->Un_1),buffer,pos,1,sizeof(src->Un_1));
  pack_data(&(src->Un),buffer,pos,1,sizeof(src->Un));
  pack_data(&(src->Jn_1),buffer,pos,1,sizeof(src->Jn_1));

  if(analysis == STABILIZED
     || analysis == MINI){
    pack_data(src->Fpp,buffer,pos,TENSOR_2,sizeof(*(src->Fpp)));
  }
}

static void unpack_IL2_eps(IL2_eps *dest,
			   const int analysis,
			   const char *buffer,
			   size_t *pos)
{
  unpack_damage(&(dest->dam),buffer,pos);
  unpack_data(buffer,&(dest->Un_1),pos,1,sizeof(dest->Un_1));
  unpack_data(buffer,&(dest->Un),pos,1,sizeof(dest->Un));
  unpack_data(buffer,&(dest->Jn_1),pos,1,sizeof(dest->Jn_1));

  if(analysis == STABILIZED
     || analysis == MINI){
    unpack_data(buffer,dest->Fpp,pos,TENSOR_2,sizeof(*(dest->Fpp)));
  }
}

static void pack_eps_local(const EPS *src,
			   const int analysis,
			   char *buffer,
			   size_t *pos)
{
  pack_data(&(src->el.eq),buffer,pos,1,sizeof((src->el.eq)));
  pack_data(&(src->el.eq_m),buffer,pos,1,sizeof((src->el.eq_m)));
  pack_data(&(src->el.eq_i),buffer,pos,1,sizeof((src->el.eq_i)));
  switch(analysis){
  default:
    pack_data(src->el.o,buffer,pos,SYM_TENSOR,sizeof(*(src->el.o)));
    pack_data(src->pl.o,buffer,pos,SYM_TENSOR,sizeof(*(src->pl.o)));
    pack_data(src->pl.eq,buffer,pos,2,sizeof(*(src->pl.eq)));
    break;
  case ELASTIC:
    pack_data(src->el.o,buffer,pos,SYM_TENSOR,sizeof(*(src->el.o)));
    break;
  case TP_ELASTO_PLASTIC:
    pack_data(src->el.o,buffer,pos,SYM_TENSOR,sizeof(*(src->el.o)));
    pack_data(src->el.f,buffer,pos,SYM_TENSOR,sizeof(*(src->el.f)));
    pack_data(src->el.d,buffer,pos,SYM_TENSOR,sizeof(*(src->el.d)));
    pack_data(src->el.m,buffer,pos,SYM_TENSOR,sizeof(*(src->el.m)));
    pack_data(src->el.i,buffer,pos,SYM_TENSOR,sizeof(*(src->el.i)));
    break;
  } 
}

static void unpack_eps_local(EPS *dest,
			     const int analysis,
			     const char *buffer,
			     size_t *pos)
{
  unpack_data(buffer,&(dest->el.eq),pos,1,sizeof((dest->el.eq)));
  unpack_data(buffer,&(dest->el.eq_m),pos,1,sizeof((dest->el.eq_m)));
  unpack_data(buffer,&(dest->el.eq_i),pos,1,sizeof((dest->el.eq_i)));
  switch(analysis){
  default:
    unpack_data(buffer,dest->el.o,pos,SYM_TENSOR,sizeof(*(dest->el.o)));
    unpack_data(buffer,dest->pl.o,pos,SYM_TENSOR,sizeof(*(dest->pl.o)));
    unpack_data(buffer,dest->pl.eq,pos,2,sizeof(*(dest->pl.eq)));
    break;
  case ELASTIC:
    unpack_data(buffer,dest->el.o,pos,SYM_TENSOR,sizeof(*(dest->el.o)));
    break;
  case TP_ELASTO_PLASTIC:
    unpack_data(buffer,dest->el.o,pos,SYM_TENSOR,sizeof(*(dest->el.o)));
    unpack_data(buffer,dest->el.f,pos,SYM_TENSOR,sizeof(*(dest->el.f)));
    unpack_data(buffer,dest->el.d,pos,SYM_TENSOR,sizeof(*(dest->el.d)));
    unpack_data(buffer,dest->el.m,pos,SYM_TENSOR,sizeof(*(dest->el.m)));
    unpack_data(buffer,dest->el.i,pos,SYM_TENSOR,sizeof(*(dest->el.i)));
    break;
  } 
}

/** Pack a single EPS object */
static void pack_eps(const EPS *src,
		     const ELEMENT *elem,
		     const int analysis,
		     char *buffer,
		     size_t *pos)
{
  long pt_I = 0;
  long pt_J = 0;

  /* Get number of integration points */
  int_point(elem->toe,&pt_I);
  switch(analysis){
  default: int_point(10,&pt_J); break;
  case MINI: case MINI_3F: int_point(5,&pt_J); break;
  }

  pack_eps_local(src,analysis,buffer,pos);

  for(int i=0; i<pt_I; i++){
    pack_damage((src->dam)+i,buffer,pos);
    pack_IL0_eps((src->il)+i,analysis,buffer,pos);
    pack_IL1_eps((src->d_il)+i,analysis,buffer,pos);
  }

  for(int j=0; j<pt_J; j++){
    pack_IL2_eps((src->st)+j,analysis,buffer,pos);
  }
}

/** Unpack a single EPS object */
static void unpack_eps(EPS *dest,
		       const ELEMENT *elem,
		       const int analysis,
		       const char *buffer,
		       size_t *pos)
{
  long pt_I = 0;
  long pt_J = 0;

  /* Get number of integration points */
  int_point(elem->toe,&pt_I);
  switch(analysis){
  default: int_point(10,&pt_J); break;
  case MINI: case MINI_3F: int_point(5,&pt_J); break;
  }

  unpack_eps_local(dest,analysis,buffer,pos);

  for(int i=0; i<pt_I; i++){
    unpack_damage((dest->dam)+i,buffer,pos);
    unpack_IL0_eps((dest->il)+i,analysis,buffer,pos);
    unpack_IL1_eps((dest->d_il)+i,analysis,buffer,pos);
  }

  for(int j=0; j<pt_J; j++){
    unpack_IL2_eps((dest->st)+j,analysis,buffer,pos);
  }
}

void pack_eps_list(const EPS *src,
		   const int ne,
		   const ELEMENT *elem,
		   const int analysis,
		   char *buffer,
		   size_t *pos)
{
  switch(analysis){
  default:
    pack_2mat(cast_const_V2 src[0].Dp,NDN,NDN,
	      sizeof(**(src[0].Dp)),buffer,pos);
    break;
  case ELASTIC: break;
  case TP_ELASTO_PLASTIC: break;
  }

  for(int n=0; n<ne; n++){
    pack_eps(src+n,elem+n,analysis,buffer,pos);
  }
}

void unpack_eps_list(EPS *dest,
		     const int ne,
		     const ELEMENT *elem,
		     const int analysis,
		     const char *buffer,
		     size_t *pos)
{
  switch(analysis){
  default:
    unpack_2mat(cast_V2 dest[0].Dp,NDN,NDN,
	      sizeof(**(dest[0].Dp)),buffer,pos);
    break;
  case ELASTIC: break;
  case TP_ELASTO_PLASTIC: break;
  }

  for(int n=0; n<ne; n++){
    unpack_eps(dest+n,elem+n,analysis,buffer,pos);
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
  int err = 0;
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

    /* constitutive_model */
    if(analysis==CM)
    {  
      if (p_eps->model != NULL) {
        for (int j = 0; j < n_ip; j++) {
          err += constitutive_model_destroy((p_eps->model) + j);
        }
      }
    }

    /* remaining */
    free(p_eps->il);
    free(p_eps->d_il);
    free(p_eps->st);
    free(p_eps->dam);
    free(p_eps->T);
    free(p_eps->d_T);
    free(p_eps->model);
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
