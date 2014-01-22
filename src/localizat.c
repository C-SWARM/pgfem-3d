#include "localizat.h"

void localizat (double *k,double *lk,long *adr,long *cn,long ndofe)
{
  long i,j,ii,jj;
  
  for (i=0;i<ndofe;i++){
    ii = cn[i]-1;
    if (ii <= -1)  continue;
    for (j=0;j<ndofe;j++){
      jj = cn[j]-1;
      if (jj <= -1 || ii > jj)  continue;
      k[adr[jj]+jj-ii] += lk[i*ndofe+j];
    }
  }
}

void localizat_nonsym (double *kL,double *kU,double *lk,long *adr,long *cn,long ndofe)
{
  long i,j,ii,jj;
  
  for (i=0;i<ndofe;i++){/* row */
    ii = cn[i]-1;
    if (ii < 0)  continue;
    for (j=0;j<ndofe;j++){/* column */
      jj = cn[j]-1;
      if (jj < 0)  continue;
      if (ii > jj)
	kL[adr[ii]-(ii-jj)] += lk[i*ndofe+j];
      if (ii < jj)
	kU[adr[jj]-(jj-ii)] += lk[i*ndofe+j];
      if (ii == jj){
	kL[adr[jj]] += lk[i*ndofe+j];
	kU[adr[jj]] += lk[i*ndofe+j];
      }
    }
  }
}

void localizat_sparse (double *K,double *lk,long *Ai,long *Ap,long *cn,long ndofe)
{
  long i,j,k,ii,jj;
  
  for (i=0;i<ndofe;i++){/* column */
    ii = cn[i]-1; if (ii < 0)  continue;
    
    if ((Ap[ii+1]-Ap[ii]) > ndofe){
      for (k=Ap[ii];k<Ap[ii+1];k++){
	for (j=0;j<ndofe;j++){/* row */
	  jj = cn[j]-1;
	  if (Ai[k] != jj || jj < 0) continue;
	  else {K[k] += lk[j*ndofe+i]; break;}
	}
      }/* k < Ap[ii+1] */
      
    }
    else{
      for (j=0;j<ndofe;j++){/* row */
	jj = cn[j]-1;  if (jj < 0)  continue;
	for (k=Ap[ii];k<Ap[ii+1];k++){
	  if (Ai[k] != jj) continue;
	  else {K[k] += lk[j*ndofe+i]; break;}
	}
      }/* end j < ndofe */
      
    }
  }
}
