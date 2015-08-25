/*******************************************
 * HOMOGENIZATION of Composite media    3D *
 * Karel Matous                            *
 * December 2000                       (c) *
 *******************************************/

#include "homogen.h"
#include <math.h>

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef INCL_H
#include "incl.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

void Mat_3D_orthotropic (const long nmat,
			 MATERIAL *mater,
			 const int analysis)
/*
  
 */
{
  long i;
  double **L,**M;
  
  L = aloc2 (6,6);
 
  for (i=0;i<nmat;i++){
    M = aloc2 (6,6);
    
    switch(analysis){
    case STABILIZED:
    case MINI:
    case MINI_3F:
    case DISP:
    case TF:
    case CM:
      M[0][0] = 1./mater[i].Ex;            M[0][1] = -mater[i].nyz/mater[i].Ex;  M[0][2] = -mater[i].nyz/mater[i].Ex;
      M[1][0] = -mater[i].nyz/mater[i].Ex; M[1][1] = 1./mater[i].Ex;             M[1][2] = -mater[i].nyz/mater[i].Ex;
      M[2][0] = -mater[i].nyz/mater[i].Ex; M[2][1] = -mater[i].nyz/mater[i].Ex;  M[2][2] = 1./mater[i].Ex;
      
      M[3][3] = 1./mater[i].Gyz;  M[4][4] = 1./mater[i].Gyz;  M[5][5] = 1./mater[i].Gyz;
      break;
    default:
      M[0][0] = 1./mater[i].Ex;            M[0][1] = -mater[i].nxy/mater[i].Ey;  M[0][2] = -mater[i].nxz/mater[i].Ez;
      M[1][0] = -mater[i].nxy/mater[i].Ey; M[1][1] = 1./mater[i].Ey;             M[1][2] = -mater[i].nyz/mater[i].Ez;
      M[2][0] = -mater[i].nxz/mater[i].Ez; M[2][1] = -mater[i].nyz/mater[i].Ez;  M[2][2] = 1./mater[i].Ez;
      
      M[3][3] = 1./mater[i].Gyz;  M[4][4] = 1./mater[i].Gxz;  M[5][5] = 1./mater[i].Gxy;
      break;
    }
    
    mater[i].M[0] = M[0][0];
    mater[i].M[1] = M[0][1];
    mater[i].M[2] = M[0][2];
    mater[i].M[3] = M[1][1];
    mater[i].M[4] = M[1][2];
    mater[i].M[5] = M[2][2];
    mater[i].M[6] = M[3][3];
    mater[i].M[7] = M[4][4];
    mater[i].M[8] = M[5][5];
    
    inv_I (M,L,6);
    
    mater[i].L[0] = L[0][0];
    mater[i].L[1] = L[0][1];
    mater[i].L[2] = L[0][2];
    mater[i].L[3] = L[1][1];
    mater[i].L[4] = L[1][2];
    mater[i].L[5] = L[2][2];
    mater[i].L[6] = L[3][3];
    mater[i].L[7] = L[4][4];
    mater[i].L[8] = L[5][5];
    
    dealoc2 (M,6);
  }
  dealoc2 (L,6);
}

void Mat_trans_isotropic (long nmat,MATERIAL *mater)
/*
       
 */
{
  long ii;
  double  k,l,m;
  double **L,**M;
  
  L = aloc2 (6,6);
  M = aloc2 (6,6);
  
  for (ii=0;ii<nmat;ii++){
    
    k=-1/(1/(mater[ii].Gyz)-4/(mater[ii].Ey)+4*mater[ii].nxy*mater[ii].nxy/(mater[ii].Ex));
    l=2*k*mater[ii].nxy;
    m=mater[ii].Ex+l*l/k;
    
    /* Matice tuhosti */                   
    
    L[0][0] = m;
    L[0][1] = l;
    L[1][0] = l;
    L[0][2] = l;
    L[2][0] = l;
    L[1][1] = k+mater[ii].Gyz;
    L[2][2] = k+mater[ii].Gyz;
    L[1][2] = k-mater[ii].Gyz;
    L[2][1] = k-mater[ii].Gyz;
    L[3][3] = mater[ii].Gyz;
    L[4][4] = mater[ii].Gxy;
    L[5][5] = mater[ii].Gxy;
    
    /* Matice poddajnosti */
    
    inv_I(L,M,6);
    
    mater[ii].L[0] = L[0][0];
    mater[ii].L[1] = L[0][1];
    mater[ii].L[2] = L[1][1];
    mater[ii].L[3] = L[1][2];
    mater[ii].L[4] = L[3][3];
    mater[ii].L[5] = L[4][4];
    
    mater[ii].M[0] = M[0][0];
    mater[ii].M[1] = M[0][1];
    mater[ii].M[2] = M[1][1];
    mater[ii].M[3] = M[1][2];
    mater[ii].M[4] = M[3][3];
    mater[ii].M[5] = M[4][4];
  }
  
  dealoc2 (M,6);  dealoc2 (L,6);
}

void hom_matrices (long ***a,
		   const long ne,
		   const long nmat,
		   const long nc,
		   ELEMENT *elem,
		   MATERIAL *mater,
		   MATGEOM matgeom,
		   HOMMAT *hommat,
		   const long SHAPE,
		   const int analysis)
/*
  creates material matrices of the homogeneous medium
*/
{
  long i,j,k,nn;
  
  nn=0;
  for (i=0;i<nmat;i++){
    for (j=0;j<nmat;j++){
      for (k=0;k<nc;k++){
	if (a[i][j][k]==1){
	  Overall_Mat (i,j,k,nn,mater,matgeom,hommat,SHAPE);
	  hommat[nn].e1 =  hommat[nn].e2
	    = hommat[nn].e3 = hommat[nn].e4 =  0.0;
	  switch(analysis){
	  case DISP: /* only DISP gets damage parameters */
	    /* use extra entries for extra variables 
	       (vol damage: Yin, p1, p2, mu) */
	    hommat[nn].e1 = mater[i].Gxz;
	    hommat[nn].e2 = mater[i].Gxy;
	    hommat[nn].e3 = mater[i].nxz;
	    hommat[nn].e4 = mater[i].nxy;
	    /* intentional drop through for isotropic properties */

	    /* ISOTROPIC MAT PROPS */
	  case STABILIZED:
	  case MINI:
	  case MINI_3F:
	  case TF:
	  case CM:
	    
	    /* Mooney - Rivlin */
	    /* This is a material input from file LOOK OUT */
	    hommat[nn].m10 = mater[i].Ey;
	    hommat[nn].m01 = mater[i].Ez; 
	    hommat[nn].E = mater[i].Ex;
	    hommat[nn].G = mater[i].Gyz;
	    hommat[nn].nu = mater[i].nyz;

	    /* Potential functions for isotropic materials */
	    hommat[nn].devPotFlag = mater[i].devPotFlag;
	    hommat[nn].volPotFlag = mater[i].volPotFlag;
	    break;

	  default:
	    break;
	  }

	  a[i][j][k] = nn;
	  nn++;
	}
      }
    }
  }

  /* REAPLCED 12/10/2012 MM */  
  assign_homogenized_material(ne,elem,a,analysis);
}

void assign_homogenized_material(const long ne,
				 ELEMENT *elem,
				 long ***a,
				 const int analysis)
{
  /* Copied from hom_matrices for more general use 12/10/2012 MM */
  long i;
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);


  for (i=0;i<ne;i++){
    switch(analysis){
    case ELASTIC:
    case TP_ELASTO_PLASTIC:
      break;
    default:
      if (elem[i].mat[0] != elem[i].mat[1]){ 
	if (err_rank == 0){
	  PGFEM_printf("Incorrect material input\n");
	}
	PGFEM_Abort();
      }
      break;
    }
    elem[i].mat[2] = a[elem[i].mat[0]][elem[i].mat[1]][elem[i].hom[0]];
  }
}

void funkce_Wf (long jj,MATERIAL *mater,double **Lf,double **Mf,double **Lm,double **Mm,double **Wf)
/*
       
 */
{
  long i,j;
  double  k,l,m,**P,**A,**C,**I;
  
  k=-1/(1/(mater[jj].Gyz)-4/(mater[jj].Ey)+4*mater[jj].nxy*mater[jj].nxy/(mater[jj].Ex));
  l=2*k*mater[jj].nxy;
  m=mater[jj].Ex+l*l/k;
  
  P = aloc2 (6,6);
  A = aloc2 (6,6);
  C = aloc2 (6,6);
  I = aloc2 (6,6);
  
  P[1][1] = (k+4*mater[jj].Gyz)/(8*mater[jj].Gyz*(k+mater[jj].Gyz));
  P[2][2] = (k+4*mater[jj].Gyz)/(8*mater[jj].Gyz*(k+mater[jj].Gyz));
  P[3][3] = (k+2*mater[jj].Gyz)/(2*mater[jj].Gyz*(k+mater[jj].Gyz));
  P[4][4] = 1/(2*mater[jj].Gxy);
  P[5][5] = 1/(2*mater[jj].Gxy);
  P[1][2] = -k/(8*mater[jj].Gyz*(k+mater[jj].Gyz));
  P[2][1] = -k/(8*mater[jj].Gyz*(k+mater[jj].Gyz));
  
  for(i=0;i<6;i++){
    for(j=0;j<6;j++){
      A[i][j] = Lf[i][j]-Lm[i][j];
    }
  }
  
  nas_AB (P,A,C,6,6,6);
  
  I[0][0] = 1.;  I[1][1] = 1.;  I[2][2] = 1.;  I[3][3] = 1.;  I[4][4] = 1.;  I[5][5] = 1.;
  
  for(i=0;i<6;i++){
    for(j=0;j<6;j++){
      A[i][j] = I[i][j]+C[i][j];
    }
  }
  
  inv_I(A,C,6);
  nas_AB (Lf,C,A,6,6,6);
  nas_AB (A,Mm,Wf,6,6,6);
  
  dealoc2 (P,6);  dealoc2 (A,6);  dealoc2 (C,6);  dealoc2 (I,6);
}

void funkce_Bf (long kk,MATGEOM matgeom,double **Wf,double **Bf)
/*
       
 */
{
  long i,j;
  double **A,**I,**B;
  
  A = aloc2 (6,6);  I = aloc2 (6,6);  B = aloc2 (6,6);
  
  I[0][0] = 1-matgeom->cf[kk];  I[1][1] = 1-matgeom->cf[kk];  I[2][2] = 1-matgeom->cf[kk];
  I[3][3] = 1-matgeom->cf[kk];  I[4][4] = 1-matgeom->cf[kk];  I[5][5] = 1-matgeom->cf[kk];
  
  for(i=0;i<6;i++){
    for(j=0;j<6;j++){
      A[i][j] = matgeom->cf[kk]*Wf[i][j];
    }
  }
  
  for(i=0;i<6;i++){
    for(j=0;j<6;j++){
      B[i][j] = I[i][j]+A[i][j];
    }
  }
  
  inv_I  (B,A,6);
  nas_AB (Wf,A,Bf,6,6,6);
  
  dealoc2 (A,6);  dealoc2 (I,6);  dealoc2 (B,6);
}
void funkce_Bm (long kk,MATGEOM matgeom,double **Wf,double **Bm)
/*
       
 */
{
  long i,j;
  double **A,**I,**B;
  
  A = aloc2 (6,6);  I = aloc2 (6,6);  B = aloc2 (6,6);
  
  I[0][0] = 1-matgeom->cf[kk];  I[1][1] = 1-matgeom->cf[kk];  I[2][2] = 1-matgeom->cf[kk];
  I[3][3] = 1-matgeom->cf[kk];  I[4][4] = 1-matgeom->cf[kk];  I[5][5] = 1-matgeom->cf[kk];
  
  for(i=0;i<6;i++){
    for(j=0;j<6;j++){
      A[i][j] = matgeom->cf[kk]*Wf[i][j];
    }
  }
  
  for(i=0;i<6;i++){
    for(j=0;j<6;j++){
      B[i][j] = I[i][j]+A[i][j];
    }
  }
  
  inv_I(B,Bm,6);
  
  dealoc2 (A,6);  dealoc2 (I,6);  dealoc2 (B,6);
}

void mat_vrs (long nn,long kk,HOMMAT *hommat,MATERIAL *mater,MATGEOM matgeom,double **Bf,double **Bm,double **Mf,double **Mm)
/*
       
 */	  
{
  double **A,**B;
  
  A = aloc2 (6,6);  B = aloc2 (6,6);
  
  nas_AB (Mf,Bf,A,6,6,6);
  nas_AB (Mm,Bm,B,6,6,6);
  
  hommat[nn].M[0] = matgeom->cf[kk]*A[0][0] + (1-matgeom->cf[kk])*B[0][0];
  hommat[nn].M[1] = matgeom->cf[kk]*A[0][1] + (1-matgeom->cf[kk])*B[0][1];
  hommat[nn].M[2] = matgeom->cf[kk]*A[1][1] + (1-matgeom->cf[kk])*B[1][1];
  hommat[nn].M[3] = matgeom->cf[kk]*A[1][2] + (1-matgeom->cf[kk])*B[1][2];
  hommat[nn].M[4] = matgeom->cf[kk]*A[3][3] + (1-matgeom->cf[kk])*B[3][3];
  hommat[nn].M[5] = matgeom->cf[kk]*A[4][4] + (1-matgeom->cf[kk])*B[4][4];
  
  dealoc2 (A,6);  dealoc2 (B,6);
}

void mori_tanaka (long ii,long jj,long kk,long nn,MATERIAL *mater,MATGEOM matgeom,HOMMAT *hommat)
/*
       
 */
{
  double **Lf,**Mf,**Lm,**Mm,**Wf,**Bf,**Bm;
  
  Lf = aloc2 (6,6);  Mf = aloc2 (6,6);  Lm = aloc2 (6,6);  Mm = aloc2 (6,6);
  Bf = aloc2 (6,6);  Bm = aloc2 (6,6);  Wf = aloc2 (6,6);
  
  /* Matice tuhosti */                                 /* Matice poddajnosti */
  
  Lf[0][0] = mater[ii].L[0]; Lm[0][0]=mater[jj].L[0];    Mf[0][0] = mater[ii].M[0]; Mm[0][0] = mater[jj].M[0];
  Lf[0][1] = mater[ii].L[1]; Lm[0][1]=mater[jj].L[1];    Mf[0][1] = mater[ii].M[1]; Mm[0][1] = mater[jj].M[1];
  Lf[0][2] = mater[ii].L[1]; Lm[0][2]=mater[jj].L[1];    Mf[0][2] = mater[ii].M[1]; Mm[0][2] = mater[jj].M[1];
  
  Lf[1][0] = mater[ii].L[1]; Lm[1][0]=mater[jj].L[1];    Mf[1][0] = mater[ii].M[1]; Mm[1][0] = mater[jj].M[1];
  Lf[1][1] = mater[ii].L[2]; Lm[1][1]=mater[jj].L[2];    Mf[1][1] = mater[ii].M[2]; Mm[1][1] = mater[jj].M[2];
  Lf[1][2] = mater[ii].L[3]; Lm[1][2]=mater[jj].L[3];    Mf[1][2] = mater[ii].M[3]; Mm[1][2] = mater[jj].M[3];
  
  Lf[2][0] = mater[ii].L[1]; Lm[2][0]=mater[jj].L[1];    Mf[2][0] = mater[ii].M[1]; Mm[2][0] = mater[jj].M[1];
  Lf[2][1] = mater[ii].L[3]; Lm[2][1]=mater[jj].L[3];    Mf[2][1] = mater[ii].M[3]; Mm[2][1] = mater[jj].M[3];
  Lf[2][2] = mater[ii].L[2]; Lm[2][2]=mater[jj].L[2];    Mf[2][2] = mater[ii].M[2]; Mm[2][2] = mater[jj].M[2];
  
  Lf[3][3] = mater[ii].L[4]; Lm[3][3]=mater[jj].L[4];    Mf[3][3] = mater[ii].M[4]; Mm[3][3] = mater[jj].M[4];
  Lf[4][4] = mater[ii].L[5]; Lm[4][4]=mater[jj].L[5];    Mf[4][4] = mater[ii].M[5]; Mm[4][4] = mater[jj].M[5];
  Lf[5][5] = mater[ii].L[5]; Lm[5][5]=mater[jj].L[5];    Mf[5][5] = mater[ii].M[5]; Mm[5][5] = mater[jj].M[5];
  
  funkce_Wf (jj,mater,Lf,Mf,Lm,Mm,Wf);
  funkce_Bf (kk,matgeom,Wf,Bf);
  funkce_Bm (kk,matgeom,Wf,Bm);
  mat_vrs (nn,kk,hommat,mater,matgeom,Bf,Bm,Mf,Mm);
  
  /* FREE */
  dealoc2(Lf,6);  dealoc2(Mf,6);  dealoc2(Lm,6);  dealoc2(Mm,6);  dealoc2(Wf,6);  dealoc2(Bf,6);  
  dealoc2(Bm,6);
}

void Stiffness_Matrix_3D (long ii,long ip,ELEMENT *elem,HOMMAT *hommat,double **D,long TYPE)
/*
  elem[ii].mat[0] -> Matrix
  elem[ii].mat[1] -> Fiber
  elem[ii].mat[2] -> Homogeneous
*/
{
  long i,j,k;
  
  nulld2 (D,6,6);
  
  if (TYPE == 0){
    D[0][0] = hommat[elem[ii].mat[2]].L[0];
    D[0][1] = hommat[elem[ii].mat[2]].L[1];
    D[0][2] = hommat[elem[ii].mat[2]].L[2];
    
    D[1][0] = hommat[elem[ii].mat[2]].L[1];
    D[1][1] = hommat[elem[ii].mat[2]].L[3];
    D[1][2] = hommat[elem[ii].mat[2]].L[4];
    
    D[2][0] = hommat[elem[ii].mat[2]].L[2];
    D[2][1] = hommat[elem[ii].mat[2]].L[4];
    D[2][2] = hommat[elem[ii].mat[2]].L[5];
    
    D[3][3] = hommat[elem[ii].mat[2]].L[6];
    D[4][4] = hommat[elem[ii].mat[2]].L[7];
    D[5][5] = hommat[elem[ii].mat[2]].L[8];
  }
  if (TYPE == 1){
    k = 0;
    for (i=0;i<6;i++) for (j=i;j<6;j++) { D[i][j] = D[j][i] = elem[ii].L[ip][k];k++;}
  }
}

/****************************************************************************************************************/
/**********************  TRANSFORMATION FILED ANALYSIS FOR THREE PHASE MEDIUM SYSTEM  ***************************/
/****************************************************************************************************************/

TFA* build_tfa (long i)
/*
  i - pocet homogennich materialu
*/
{
  TFA *pom;
  long ii;
  
  pom = (TFA*) PGFEM_calloc (i, sizeof(TFA));
  
  for (ii=0;ii<i;ii++){	 
    pom[ii].Am = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].A2 = (double*) PGFEM_calloc (36,sizeof(double));
    
    pom[ii].Dmm = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Dmd = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Dmb = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Dbm = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Dbb = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Dbd = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Ddm = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Ddd = (double*) PGFEM_calloc (36,sizeof(double));
    pom[ii].Ddb = (double*) PGFEM_calloc (36,sizeof(double));
  }
  if (pom == NULL){
    PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__);
    abort ();
  }
  return (pom);
}

void Overall_Mat (long ii,long jj,long kk,long nn,MATERIAL *mater,MATGEOM matgeom,HOMMAT *hommat,long TYPE)
/*
       
 */
{
  long i,j;
  double **S,**Lo,**Mo,**I,**P,E,G,K,nu,si,th,**SI,**Lst,**Mm,**Lm,**M2,**L2,**Am,**A2;
  double **A,**B,**C;
  
  Lo = aloc2 (6,6); Mo = aloc2 (6,6); I = aloc2 (6,6); P = aloc2 (6,6); Lst = aloc2 (6,6); SI = aloc2 (6,6); A = aloc2 (6,6); B = aloc2 (6,6);
  M2 = aloc2 (6,6); L2 = aloc2 (6,6); Lm = aloc2 (6,6); Mm = aloc2 (6,6); A2 = aloc2 (6,6); Am = aloc2 (6,6); C = aloc2 (6,6); S = aloc2 (6,6);
  
  for (i=0;i<6;i++) for (j=0;j<6;j++){if (i == j) I[i][j] = 1.0; else I[i][j] = 0.0;}
  
  Lm[0][0] = mater[ii].L[0];  L2[0][0] = mater[jj].L[0];  Mm[0][0] = mater[ii].M[0];  M2[0][0] = mater[jj].M[0];
  Lm[0][1] = mater[ii].L[1];  L2[0][1] = mater[jj].L[1];  Mm[0][1] = mater[ii].M[1];  M2[0][1] = mater[jj].M[1];
  Lm[1][0] = mater[ii].L[1];  L2[1][0] = mater[jj].L[1];  Mm[1][0] = mater[ii].M[1];  M2[1][0] = mater[jj].M[1];
  Lm[0][2] = mater[ii].L[2];  L2[0][2] = mater[jj].L[2];  Mm[0][2] = mater[ii].M[2];  M2[0][2] = mater[jj].M[2];
  Lm[2][0] = mater[ii].L[2];  L2[2][0] = mater[jj].L[2];  Mm[2][0] = mater[ii].M[2];  M2[2][0] = mater[jj].M[2];
  Lm[1][1] = mater[ii].L[3];  L2[1][1] = mater[jj].L[3];  Mm[1][1] = mater[ii].M[3];  M2[1][1] = mater[jj].M[3];
  Lm[1][2] = mater[ii].L[4];  L2[1][2] = mater[jj].L[4];  Mm[1][2] = mater[ii].M[4];  M2[1][2] = mater[jj].M[4];
  Lm[2][1] = mater[ii].L[4];  L2[2][1] = mater[jj].L[4];  Mm[2][1] = mater[ii].M[4];  M2[2][1] = mater[jj].M[4];
  Lm[2][2] = mater[ii].L[5];  L2[2][2] = mater[jj].L[5];  Mm[2][2] = mater[ii].M[5];  M2[2][2] = mater[jj].M[5];
  Lm[3][3] = mater[ii].L[6];  L2[3][3] = mater[jj].L[6];  Mm[3][3] = mater[ii].M[6];  M2[3][3] = mater[jj].M[6];
  Lm[4][4] = mater[ii].L[7];  L2[4][4] = mater[jj].L[7];  Mm[4][4] = mater[ii].M[7];  M2[4][4] = mater[jj].M[7];
  Lm[5][5] = mater[ii].L[8];  L2[5][5] = mater[jj].L[8];  Mm[5][5] = mater[ii].M[8];  M2[5][5] = mater[jj].M[8];
  
  /* Comparison medium,  Lo = (a1*L1 + a2*L2)  */
  for (i=0;i<6;i++) for (j=0;j<6;j++) Lo[i][j] = 0.0;
  
  Lo[0][0] = matgeom->a1*mater[ii].L[0] + matgeom->a2*mater[jj].L[0];
  Lo[0][1] = matgeom->a1*mater[ii].L[1] + matgeom->a2*mater[jj].L[1];
  Lo[1][0] = matgeom->a1*mater[ii].L[1] + matgeom->a2*mater[jj].L[1];
  Lo[0][2] = matgeom->a1*mater[ii].L[2] + matgeom->a2*mater[jj].L[2];
  Lo[2][0] = matgeom->a1*mater[ii].L[2] + matgeom->a2*mater[jj].L[2];
  Lo[1][1] = matgeom->a1*mater[ii].L[3] + matgeom->a2*mater[jj].L[3];
  Lo[1][2] = matgeom->a1*mater[ii].L[4] + matgeom->a2*mater[jj].L[4];
  Lo[2][1] = matgeom->a1*mater[ii].L[4] + matgeom->a2*mater[jj].L[4];
  Lo[2][2] = matgeom->a1*mater[ii].L[5] + matgeom->a2*mater[jj].L[5];
  Lo[3][3] = matgeom->a1*mater[ii].L[6] + matgeom->a2*mater[jj].L[6];
  Lo[4][4] = matgeom->a1*mater[ii].L[7] + matgeom->a2*mater[jj].L[7];
  Lo[5][5] = matgeom->a1*mater[ii].L[8] + matgeom->a2*mater[jj].L[8];
  
  inv_I (Lo,Mo,6);
  
  /* Type of inclusion */
  switch (TYPE){
  case 0:{ /* Spherical inclusion in Isotropic medium */
    E = 1./Mo[0][0];  G = 1./Mo[3][3];  nu = -Mo[0][1]*E;  K = E/(3*(1 - 2*nu));  si = (1 + nu)/(3*(1 - nu));  th = 2*(4 - 5*nu)/(15*(1 - nu));
    
    P[0][0] = 1./3.*(si/3./K + th/G);    P[0][1] = 1./3.*(si/3./K - th/2./G); P[0][2] = 1./3.*(si/3./K - th/2./G);
    P[1][0] = 1./3.*(si/3./K - th/2./G); P[1][1] = 1./3.*(si/3./K + th/G);    P[1][2] = 1./3.*(si/3./K - th/2./G);
    P[2][0] = 1./3.*(si/3./K - th/2./G); P[2][1] = 1./3.*(si/3./K - th/2./G); P[2][2] = 1./3.*(si/3./K + th/G);   
    
    P[3][3] = th/G;  P[4][4] = th/G;  P[5][5] = th/G;
    
    /* Eshelby tensor */
    nas_AB (Lo,P,S,6,6,6);
    break;
  }
  default:{
    PGFEM_printf ("Non existing shape of inclusion\n");
    abort();
    break;
  }
  }
  
  inv_I (S,SI,6); for (i=0;i<6;i++) for (j=0;j<6;j++) B[i][j] = I[i][j] - S[i][j]; nas_AB (Lo,B,A,6,6,6); nas_AB (A,SI,Lst,6,6,6);
  
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A2[i][j] = Lst[i][j] + Lm[i][j]; Am[i][j] = Lst[i][j] + Lo[i][j];} inv_I (A2,C,6); nas_AB (C,Am,A,6,6,6);
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A2[i][j] = Lst[i][j] + L2[i][j]; Am[i][j] = Lst[i][j] + Lo[i][j];} inv_I (A2,C,6); nas_AB (C,Am,B,6,6,6);
  
  for (i=0;i<6;i++) for (j=0;j<6;j++) C[i][j] = matgeom->cm[kk]*A[i][j] + matgeom->cf[kk]*B[i][j]; inv_I (C,SI,6); nas_AB (A,SI,Am,6,6,6); nas_AB (B,SI,A2,6,6,6);
  
  /* Check concentration factors */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {C[i][j] = matgeom->cm[kk]*Am[i][j] + matgeom->cf[kk]*A2[i][j];
      if (i == j && (C[i][j] <= 0.9999 || C[i][j] >=  1.0001)) {PGFEM_printf("Incorrect Concentration tensors Am, A2\n"); abort();}
      if (i != j && (C[i][j] >= 1.e-15 || C[i][j] <= -1.e-15)) {PGFEM_printf("Incorrect Concentration tensors Am, A2\n"); abort();}
    }
  
  nas_AB (Lm,Am,A,6,6,6);  nas_AB (L2,A2,B,6,6,6);
  
  for (i=0;i<6;i++) for (j=0;j<6;j++) C[i][j] = matgeom->cm[kk]*A[i][j] + matgeom->cf[kk]*B[i][j];  inv_I (C,A,6);
  
  hommat[nn].L[0] = C[0][0];  hommat[nn].M[0] = A[0][0];
  hommat[nn].L[1] = C[0][1];  hommat[nn].M[1] = A[0][1];
  hommat[nn].L[2] = C[0][2];  hommat[nn].M[2] = A[0][2];
  hommat[nn].L[3] = C[1][1];  hommat[nn].M[3] = A[1][1];
  hommat[nn].L[4] = C[1][2];  hommat[nn].M[4] = A[1][2];
  hommat[nn].L[5] = C[2][2];  hommat[nn].M[5] = A[2][2];
  hommat[nn].L[6] = C[3][3];  hommat[nn].M[6] = A[3][3];
  hommat[nn].L[7] = C[4][4];  hommat[nn].M[7] = A[4][4];
  hommat[nn].L[8] = C[5][5];  hommat[nn].M[8] = A[5][5];
  
  dealoc2 (Mo,6); dealoc2 (Lo,6); dealoc2 (I,6); dealoc2 (S,6); dealoc2 (P,6); dealoc2 (SI,6); dealoc2 (Lst,6); dealoc2 (B,6); dealoc2 (A,6); dealoc2 (C,6);
  dealoc2 (L2,6); dealoc2 (M2,6); dealoc2 (Lm,6); dealoc2 (Mm,6); dealoc2 (Am,6); dealoc2 (A2,6);
}

void TFA_tensors (long ***a,long ne,long nmat,long nc,ELEMENT *elem,MATERIAL *mater,MATGEOM matgeom,HOMMAT *hommat,TFA *tfa,long SHAPE,long TYPE)
/*
  creates material matrices of the homogeneous medium
*/
{
  long i,j,k,nn;
  
  nn=0;
  for (i=0;i<nmat;i++){
    for (j=0;j<nmat;j++){
      for (k=0;k<nc;k++){
	if (a[i][j][k] == 1){
	  A_D_tensors (i,j,k,nn,mater,matgeom,hommat,tfa,SHAPE,TYPE);
	  a[i][j][k] = nn;
	  nn++;
	}
      }
    }
  }
}

void A_D_tensors (long ii,long jj,long kk,long nn,MATERIAL *mater,MATGEOM matgeom,HOMMAT *hommat,TFA *tfa,long SHAPE,long TYPE)
/*
       
 */
{
  long i,j;
  double **S,**Lo,**Mo,**I,**P,E,G,K,nu,si,th,**SI,**Lst,**Mm,**Lm,**M2,**L2,**Am,**A2,**L,**M;
  double **A,**B,**C;
  
  Lo = aloc2 (6,6); Mo = aloc2 (6,6); I = aloc2 (6,6); P = aloc2 (6,6); Lst = aloc2 (6,6); SI = aloc2 (6,6); A = aloc2 (6,6); B = aloc2 (6,6);
  M2 = aloc2 (6,6); L2 = aloc2 (6,6); Lm = aloc2 (6,6); Mm = aloc2 (6,6); A2 = aloc2 (6,6); Am = aloc2 (6,6); C = aloc2 (6,6); S = aloc2 (6,6);
  L = aloc2 (6,6); M = aloc2 (6,6);
  
  for (i=0;i<6;i++) for (j=0;j<6;j++){if (i == j) I[i][j] = 1.0; else I[i][j] = 0.0;}
  
  Lm[0][0] = mater[ii].L[0];  L2[0][0] = mater[jj].L[0];  Mm[0][0] = mater[ii].M[0];  M2[0][0] = mater[jj].M[0];
  Lm[0][1] = mater[ii].L[1];  L2[0][1] = mater[jj].L[1];  Mm[0][1] = mater[ii].M[1];  M2[0][1] = mater[jj].M[1];
  Lm[1][0] = mater[ii].L[1];  L2[1][0] = mater[jj].L[1];  Mm[1][0] = mater[ii].M[1];  M2[1][0] = mater[jj].M[1];
  Lm[0][2] = mater[ii].L[2];  L2[0][2] = mater[jj].L[2];  Mm[0][2] = mater[ii].M[2];  M2[0][2] = mater[jj].M[2];
  Lm[2][0] = mater[ii].L[2];  L2[2][0] = mater[jj].L[2];  Mm[2][0] = mater[ii].M[2];  M2[2][0] = mater[jj].M[2];
  Lm[1][1] = mater[ii].L[3];  L2[1][1] = mater[jj].L[3];  Mm[1][1] = mater[ii].M[3];  M2[1][1] = mater[jj].M[3];
  Lm[1][2] = mater[ii].L[4];  L2[1][2] = mater[jj].L[4];  Mm[1][2] = mater[ii].M[4];  M2[1][2] = mater[jj].M[4];
  Lm[2][1] = mater[ii].L[4];  L2[2][1] = mater[jj].L[4];  Mm[2][1] = mater[ii].M[4];  M2[2][1] = mater[jj].M[4];
  Lm[2][2] = mater[ii].L[5];  L2[2][2] = mater[jj].L[5];  Mm[2][2] = mater[ii].M[5];  M2[2][2] = mater[jj].M[5];
  Lm[3][3] = mater[ii].L[6];  L2[3][3] = mater[jj].L[6];  Mm[3][3] = mater[ii].M[6];  M2[3][3] = mater[jj].M[6];
  Lm[4][4] = mater[ii].L[7];  L2[4][4] = mater[jj].L[7];  Mm[4][4] = mater[ii].M[7];  M2[4][4] = mater[jj].M[7];
  Lm[5][5] = mater[ii].L[8];  L2[5][5] = mater[jj].L[8];  Mm[5][5] = mater[ii].M[8];  M2[5][5] = mater[jj].M[8];
  
  /* Comparison medium,  Lo = (a1*L1 + a2*L2)  */
  for (i=0;i<6;i++) for (j=0;j<6;j++) Lo[i][j] = 0.0;
  
  Lo[0][0] = matgeom->a1*mater[ii].L[0] + matgeom->a2*mater[jj].L[0];
  Lo[0][1] = matgeom->a1*mater[ii].L[1] + matgeom->a2*mater[jj].L[1];
  Lo[1][0] = matgeom->a1*mater[ii].L[1] + matgeom->a2*mater[jj].L[1];
  Lo[0][2] = matgeom->a1*mater[ii].L[2] + matgeom->a2*mater[jj].L[2];
  Lo[2][0] = matgeom->a1*mater[ii].L[2] + matgeom->a2*mater[jj].L[2];
  Lo[1][1] = matgeom->a1*mater[ii].L[3] + matgeom->a2*mater[jj].L[3];
  Lo[1][2] = matgeom->a1*mater[ii].L[4] + matgeom->a2*mater[jj].L[4];
  Lo[2][1] = matgeom->a1*mater[ii].L[4] + matgeom->a2*mater[jj].L[4];
  Lo[2][2] = matgeom->a1*mater[ii].L[5] + matgeom->a2*mater[jj].L[5];
  Lo[3][3] = matgeom->a1*mater[ii].L[6] + matgeom->a2*mater[jj].L[6];
  Lo[4][4] = matgeom->a1*mater[ii].L[7] + matgeom->a2*mater[jj].L[7];
  Lo[5][5] = matgeom->a1*mater[ii].L[8] + matgeom->a2*mater[jj].L[8];
  
  inv_I (Lo,Mo,6);
  
  /* Type of inclusion */
  switch (SHAPE){
  case 0:{ /* Spherical inclusion in Isotropic medium */
    E = 1./Mo[0][0];  G = 1./Mo[3][3];  nu = -Mo[0][1]*E;  K = E/(3*(1 - 2*nu));  si = (1 + nu)/(3*(1 - nu));  th = 2*(4 - 5*nu)/(15*(1 - nu));
    
    P[0][0] = 1./3.*(si/3./K + th/G);    P[0][1] = 1./3.*(si/3./K - th/2./G); P[0][2] = 1./3.*(si/3./K - th/2./G);
    P[1][0] = 1./3.*(si/3./K - th/2./G); P[1][1] = 1./3.*(si/3./K + th/G);    P[1][2] = 1./3.*(si/3./K - th/2./G);
    P[2][0] = 1./3.*(si/3./K - th/2./G); P[2][1] = 1./3.*(si/3./K - th/2./G); P[2][2] = 1./3.*(si/3./K + th/G);   
    
    P[3][3] = th/G;  P[4][4] = th/G;  P[5][5] = th/G;
    
    /* Eshelby tensor */
    nas_AB (Lo,P,S,6,6,6);
    break;
  }
  default:{
    PGFEM_printf ("Non existing shape of inclusion\n");
    abort();
    break;
  }
  }
  
  inv_I (S,SI,6); for (i=0;i<6;i++) for (j=0;j<6;j++) B[i][j] = I[i][j] - S[i][j]; nas_AB (Lo,B,A,6,6,6); nas_AB (A,SI,Lst,6,6,6);
  
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A2[i][j] = Lst[i][j] + Lm[i][j]; Am[i][j] = Lst[i][j] + Lo[i][j];} inv_I (A2,C,6); nas_AB (C,Am,A,6,6,6);
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A2[i][j] = Lst[i][j] + L2[i][j]; Am[i][j] = Lst[i][j] + Lo[i][j];} inv_I (A2,C,6); nas_AB (C,Am,B,6,6,6);
  
  for (i=0;i<6;i++) for (j=0;j<6;j++) C[i][j] = matgeom->cm[kk]*A[i][j] + matgeom->cf[kk]*B[i][j]; inv_I (C,SI,6); nas_AB (A,SI,Am,6,6,6); nas_AB (B,SI,A2,6,6,6);
  
  /* Check transformation factors */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {C[i][j] = matgeom->cm[kk]*Am[i][j] + matgeom->cf[kk]*A2[i][j];
      if (i == j && (C[i][j] <= 0.9999999 || C[i][j] >=  1.0000001)) {PGFEM_printf("Incorrect Transformation tensors Am, A2\n"); abort();}
      if (i != j && (C[i][j] >= 1.000e-10 || C[i][j] <= -1.000e-10)) {PGFEM_printf("Incorrect Transformation tensors Am, A2\n"); abort();}
    }
  
  nas_AB (Lm,Am,A,6,6,6);  nas_AB (L2,A2,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) L[i][j] = matgeom->cm[kk]*A[i][j] + matgeom->cf[kk]*B[i][j];  inv_I (L,M,6);
  
  /* Am, A2 */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {tfa[nn].Am[i*6+j] = Am[i][j];  tfa[nn].A2[i*6+j] = A2[i][j];}
  
  /* Dmm */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - Am[i][j]; B[i][j] = Lm[i][j] - L[i][j]; C[i][j] = I[i][j] - matgeom->cm[kk]*Am[j][i];}
  inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,Lm,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Dmm[i*6+j] = B[i][j];
  
  /* Dmb */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - Am[i][j]; B[i][j] = Lm[i][j] - L[i][j]; C[i][j] = -1.*matgeom->cf[kk]*A2[j][i];}
  inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,L2,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Dmb[i*6+j] = B[i][j];
  
  /* Dbb */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - A2[i][j]; B[i][j] = L2[i][j] - L[i][j]; C[i][j] = I[i][j] - matgeom->cf[kk]*A2[j][i];}
  inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,L2,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Dbb[i*6+j] = B[i][j];
  
  /* Dbm */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - A2[i][j]; B[i][j] = L2[i][j] - L[i][j]; C[i][j] = -1.*matgeom->cm[kk]*Am[j][i];}
  inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,Lm,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Dbm[i*6+j] = B[i][j];
  
  if (TYPE == 1){/* Three Phase Medium System */
    /* Dmd */
    for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - Am[i][j]; B[i][j] = Lm[i][j] - L[i][j]; C[i][j] = -1.*matgeom->cd[kk]*A2[j][i];}
    inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,L2,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Dmd[i*6+j] = B[i][j];
    
    /* Dbd */
    for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - A2[i][j]; B[i][j] = L2[i][j] - L[i][j]; C[i][j] = -1.*matgeom->cd[kk]*A2[j][i];}
    inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,L2,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Dbd[i*6+j] = B[i][j];
    
    /* Ddd */
    for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - A2[i][j]; B[i][j] = L2[i][j] - L[i][j]; C[i][j] = I[i][j] - matgeom->cd[kk]*A2[j][i];}
    inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,L2,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Ddd[i*6+j] = B[i][j];
    
    /* Ddm */
    for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - A2[i][j]; B[i][j] = L2[i][j] - L[i][j]; C[i][j] = -1.*matgeom->cm[kk]*Am[j][i];}
    inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,Lm,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Ddm[i*6+j] = B[i][j];
    
    /* Ddb */
    for (i=0;i<6;i++) for (j=0;j<6;j++) {A[i][j] = I[i][j] - A2[i][j]; B[i][j] = L2[i][j] - L[i][j]; C[i][j] = -1.*matgeom->cf[kk]*A2[j][i];}
    inv_I (B,SI,6);  nas_AB (A,SI,B,6,6,6); nas_AB (B,C,A,6,6,6); nas_AB (A,L2,B,6,6,6); for (i=0;i<6;i++) for (j=0;j<6;j++) tfa[nn].Ddb[i*6+j] = B[i][j];
  }
  
  /* Check concentration factors */
  
  /* SUM_s=1^N Drs = I - Ar */
  for (i=0;i<6;i++) for (j=0;j<6;j++){ A[i][j] = tfa[nn].Dmm[i*6+j] + tfa[nn].Dmb[i*6+j];  B[i][j] = I[i][j] - Am[i][j]; 
      C[i][j] = tfa[nn].Dbm[i*6+j] + tfa[nn].Dbb[i*6+j]; SI[i][j] = I[i][j] - A2[i][j];
      if (sqrt((A[i][j] - B[i][j])*(A[i][j] - B[i][j])) >= 1.e-10 || sqrt((C[i][j] - SI[i][j])*(C[i][j] - SI[i][j])) >= 1.e-10) {
	PGFEM_printf("Incorrect Concentration tensors | SUM_s=1^N Drs = I - Ar\n"); abort();}
    }
  
  /* SUM_s=1^N Drs*Ms = 0 */
  for (i=0;i<6;i++) for (j=0;j<6;j++){ A[i][j] = tfa[nn].Dmm[i*6+j]; B[i][j] = tfa[nn].Dmb[i*6+j];} nas_AB (A,Mm,C,6,6,6); nas_AB (B,M2,SI,6,6,6);
  for (i=0;i<6;i++) for (j=0;j<6;j++){ A[i][j] = C[i][j] + SI[i][j];
      if (A[i][j] >= 1.000e-10 || A[i][j] <= -1.000e-10) {PGFEM_printf("Incorrect Concentration tensors | SUM_s=1^N Drs*Ms = 0\n"); abort();}
    }
  for (i=0;i<6;i++) for (j=0;j<6;j++){ A[i][j] = tfa[nn].Dbm[i*6+j]; B[i][j] = tfa[nn].Dbb[i*6+j];} nas_AB (A,Mm,C,6,6,6); nas_AB (B,M2,SI,6,6,6);
  for (i=0;i<6;i++) for (j=0;j<6;j++){ A[i][j] = C[i][j] + SI[i][j];
      if (A[i][j] >= 1.000e-10 || A[i][j] <= -1.000e-10) {PGFEM_printf("Incorrect Concentration tensors | SUM_s=1^N Drs*Ms = 0\n"); abort();}
    }
  
  /* cr*Drs*Ms = cs*Mr*Dsr^T */
  for (i=0;i<6;i++) for (j=0;j<6;j++){ A[i][j] = matgeom->cm[kk]*tfa[nn].Dmb[i*6+j]; B[i][j] = matgeom->cf[kk]*tfa[nn].Dbm[i*6+j];} 
  nas_AB (A,M2,C,6,6,6); nas_ABT (Mm,B,SI,6,6,6);
  for (i=0;i<6;i++) for (j=0;j<6;j++) if (sqrt((C[i][j] - SI[i][j])*(C[i][j] - SI[i][j])) >= 1.e-10) {
	PGFEM_printf("Incorrect Concentration tensors | cr*Drs*Ms = cs*Mr*Dsr^T AAAAAAAA\n"); abort();}
  for (i=0;i<6;i++) for (j=0;j<6;j++){ A[i][j] = matgeom->cf[kk]*tfa[nn].Dbm[i*6+j]; B[i][j] = matgeom->cm[kk]*tfa[nn].Dmb[i*6+j];}
  nas_AB (A,Mm,C,6,6,6); nas_ABT (M2,B,SI,6,6,6);
  for (i=0;i<6;i++) for (j=0;j<6;j++) if (sqrt((C[i][j] - SI[i][j])*(C[i][j] - SI[i][j])) >= 1.e-10) {
	PGFEM_printf("Incorrect Concentration tensors | cr*Drs*Ms = cs*Mr*Dsr^T\n"); abort();}
  
  /* SUM_s=1^N cs*Dsr = 0 */
  for (i=0;i<6;i++) for (j=0;j<6;j++) {
      A[i][j] = matgeom->cm[kk]*tfa[nn].Dmm[i*6+j] + matgeom->cf[kk]*tfa[nn].Dbm[i*6+j];
      B[i][j] = matgeom->cm[kk]*tfa[nn].Dmb[i*6+j] + matgeom->cf[kk]*tfa[nn].Dbb[i*6+j];
      if (A[i][j] >= 1.000e-10 || A[i][j] <= -1.000e-10 || B[i][j] >= 1.000e-10 || B[i][j] <= -1.000e-10) 
	{PGFEM_printf("Incorrect Concentration tensors | SUM_s=1^N cs*Dsr = 0\n"); abort();}
    }
  
  if (TYPE == 1){/* Three Phase Medium System */
    
  }
  
  dealoc2 (Mo,6); dealoc2 (Lo,6); dealoc2 (I,6); dealoc2 (S,6); dealoc2 (P,6); dealoc2 (SI,6); dealoc2 (Lst,6); dealoc2 (B,6); dealoc2 (A,6); dealoc2 (C,6);
  dealoc2 (L2,6); dealoc2 (M2,6); dealoc2 (Lm,6); dealoc2 (Mm,6); dealoc2 (Am,6); dealoc2 (A2,6); dealoc2 (L,6); dealoc2 (M,6);
}
