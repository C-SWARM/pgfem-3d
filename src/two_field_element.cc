/* HEADER */
#include "two_field_element.h"
#include <stdlib.h>
#include <string.h>
#include "mkl_cblas.h"
#include "PGFEM_mpi.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef INDEX_MACROS_H
#include "index_macros.h"
#endif

#ifndef UL_DEBUG
#define UL_DEBUG 0
#endif

static const int ndn = 3;

void UL_Kuu_at_ip(double *Kuu,
		  const int nne,
		  const int nne_t,
		  const double *ST,
		  const double *Fn,
		  const double *Fr,
		  const double *Fr_I,
		  const double Jn,
		  const double Jr,
		  const double pres,
		  const double *S,
		  const double *dSdF,
		  const double jj,
		  const double wt,
		  const int TYPE)
{
  double *AA_ab, *AA_wg, *B, *C, *D, *LAA_wg, *Fr_It;
  AA_ab = aloc1(9);
  AA_wg = aloc1(9);
  B = aloc1(9);
  C = aloc1(9);
  D = aloc1(9);
  double *E = aloc1(9);
  double *G = aloc1(9);
  LAA_wg = aloc1(9);
  Fr_It = aloc1(9);
  transpose(Fr_It,Fr_I,3,3);

  /* offsets for different matrices */
  int len_a, len_w, off_a, off_w;
 
  char fdebug[20];
  int err_rank;
  PGFEM_Error_rank(&err_rank);

  switch (TYPE){
  case 0: /* Kuu */
    len_a = nne;
    len_w = nne;
    off_a = 0;
    off_w = 0;
    sprintf(fdebug,"UL_Kuu_debug_%d.log",err_rank);
    break;

  case 1: /* Kbb */
    len_a = 1;
    len_w = 1;
    off_a = nne;
    off_w = nne;
    sprintf(fdebug,"UL_Kbb_debug_%d.log",err_rank);
    break;

  case 2: /* Kub */
    len_a = nne;
    len_w = 1;
    off_a = 0;
    off_w = nne;
    sprintf(fdebug,"UL_Kub_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
		   __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(UL_DEBUG) out = fopen(fdebug,"a");

  /* D = Fn S Fn' */
  for(int i=0; i<ndn; i++){
    for(int j=0; j<ndn; j++){
      D[idx_2(i,j)] = 0.0;
      for(int k=0; k<ndn; k++){
	for(int l=0; l<ndn; l++){
	  D[idx_2(i,j)] = (D[idx_2(i,j)] 
			   + Fn[idx_2(i,k)]*S[idx_2(k,l)]
			   *Fn[idx_2(j,l)]);
	}
      }
    }
  }

  for (int a=0; a<len_a; a++){
    for (int b=0; b<ndn; b++){
      if(UL_DEBUG){
      	PGFEM_fprintf(out,"=============== Row %d ===============\n",a*ndn+b);
      }

      const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
						   nne_t,ndn,ndn,ndn)];

      /*** compute AA{ab} = F^t ST{ab} ***/
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		  3,3,3,1.0,Fr,3,ptrST_ab,3,0.0,AA_ab,3);
 
      if(UL_DEBUG){
      	PGFEM_fprintf(out,"F'ST(a,b)\n");
      	print_array_d(out,AA_ab,9,3,3);
      }

      for (int w=0; w<len_w; w++){
	for (int g=0; g<ndn; g++){

	  const double* const ptrST_wg = &ST[idx_4_gen(off_w+w,g,0,0,
						       nne_t,ndn,ndn,ndn)];

	  /* compute LAA_wg = dSdF:ST_wg */
	  memset(LAA_wg,0,9*sizeof(double));
	  for(int j=0; j<9; j++){
	    if (ptrST_wg[j] == 0.0) continue;
	    for(int i=0; i<9; i++){
	      LAA_wg[i] = LAA_wg[i] +  dSdF[9*i +j] * ptrST_wg[j];
	    }
	  }


	  if(UL_DEBUG){
	    PGFEM_fprintf(out,"LAA = L:AA_wg\n");
	    print_array_d(out,LAA_wg,9,3,3);
	  }

	  /* compute triple matrix products (4 indices)*/
	  for(int i=0; i<ndn; i++){
	    for(int j=0; j<ndn; j++){
	      /* null matrices */
	      B[idx_2(i,j)] = C[idx_2(i,j)] = G[idx_2(i,j)] = 0.0;
	      E[idx_2(j,i)] = 0.0; /* null E' !!! */

	      for(int k=0; k<ndn; k++){ /* mat*mat */
		/* B = ST_wg' ST_ab */
		B[idx_2(i,j)] = (B[idx_2(i,j)]
				 + ptrST_wg[idx_2(k,i)]*ptrST_ab[idx_2(k,j)]);

		/* E' = Fr_I ST_wg */
		E[idx_2(j,i)] = (E[idx_2(j,i)] + Fr_I[idx_2(i,k)]*ptrST_wg[idx_2(k,j)]);

		/* G = Fr_I ST_ab */
		G[idx_2(i,j)] = (G[idx_2(i,j)] + Fr_I[idx_2(i,k)]*ptrST_ab[idx_2(k,j)]);

		for(int l=0; l<ndn; l++){ /* triple matrix product */
		  /* C = Fn LAA_wg Fn' */
		  C[idx_2(i,j)] = (C[idx_2(i,j)] 
				   + Fn[idx_2(i,k)]*LAA_wg[idx_2(k,l)]
				   *Fn[idx_2(j,l)]);
		}
	      }
	    }
	  }

	  /* compute dot products */
	  double CAA_ab = 0.0;
	  double DB = 0.0;
	  double PGS1 = 0.0;
	  double PGS2 = 0.0;
	  double PGS3 = 0.0;
	  for(int i=0; i<9; i++){
	    CAA_ab += C[i]*AA_ab[i];
	    DB += D[i]*B[i];
	    PGS1 += Fr_It[i]*ptrST_wg[i];
	    PGS2 += Fr_It[i]*ptrST_ab[i];
	    PGS3 += E[i]*G[i];
	  }


	  if(UL_DEBUG){
	    PGFEM_fprintf(out,"ST(w,g)'ST(a,b)\n");
	    print_array_d(out,B,9,3,3);
	  }

	  Kuu[idx_K_gen(a,b,w,g,
			len_a,ndn,
			len_w,ndn)] += ((CAA_ab+DB)/Jn + pres*Jr
					*(PGS1*PGS2-PGS3))*jj*wt;

	} /* g */
      } /* w */
    } /* b */
  } /* a */

  if(UL_DEBUG){
    PGFEM_fprintf(out,"K\n");
    print_array_d(out,Kuu,len_a*len_w*ndn*ndn,len_a*ndn,len_w*ndn);
  }

  free(AA_ab);
  free(AA_wg);
  free(B);
  free(C);
  free(D);
  free(E);
  free(G);
  free(LAA_wg);
  free(Fr_It);

  if(UL_DEBUG) fclose(out);
} /*** Kuu ***/

void UL_Kup_at_ip(double *Kup,
		  const int nne,
		  const int nne_t,
		  const double *Na,
		  const double *ST,
		  const double *Fr_I,
		  const double Jr, 
		  const double jj,
		  const double wt,
		  const int TYPE)
{
  double *B, *Fr_It;
  B = aloc1(9);
  Fr_It = aloc1(9);
  transpose(Fr_It,Fr_I,3,3);
  /* offsets for different matrices */
  int len_a, off_a;

  char fdebug[20];
  int err_rank;
  PGFEM_Error_rank(&err_rank);

  switch (TYPE){
  case 0: /* Kup */
    len_a = nne;
    off_a = 0;
    sprintf(fdebug,"UL_Kup_debug_%d.log",err_rank);
    break;

  case 1: /* Kbp */
    len_a = 1;
    off_a = nne;
    sprintf(fdebug,"UL_Kbp_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
		   __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(UL_DEBUG){
    out = fopen(fdebug,"a");
  }

  for (int a=0; a<len_a; a++){
    for(int b=0; b<ndn; b++){
      for(int w=0; w<nne; w++){
	const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
						     nne_t,ndn,ndn,ndn)];
	Kup[idx_K_gen(a,b,w,0,len_a,
		      ndn,nne,1)] += (jj*wt*Na[w]*Jr
				      *cblas_ddot(9,Fr_It,1,
						  ptrST_ab,1));
      }
    }
  }
 
  if(UL_DEBUG){
    PGFEM_fprintf(out,"K\n");
    print_array_d(out,Kup,len_a*ndn*nne,len_a*ndn,nne);

  }

  free(B);
  free(Fr_It);
  if(UL_DEBUG) fclose(out);
} /* Kup */

void UL_Kpu_at_ip(double *Kpu,
		  const int nne,
		  const int nne_t,
		  const double *Na,
		  const double *ST,
		  const double *Fr_I,
		  const double Jn,
		  const double Jr,
		  const double Upp, 
		  const double jj,
		  const double wt,
		  const int TYPE)
{
  double *B, *Fr_It;
  B = aloc1(9);
  Fr_It = aloc1(9);
  transpose(Fr_It,Fr_I,3,3);

  /* offsets for different matrices */
  int len_w, off_w;
  char fdebug[20];

  int err_rank;
  PGFEM_Error_rank(&err_rank);

  switch (TYPE){
  case 0: /* Kpu */
    len_w = nne;
    off_w = 0;
    sprintf(fdebug,"UL_Kup_debug_%d.log",err_rank);
    break;

  case 1: /* Kpb */
    len_w = 1;
    off_w = nne;
    sprintf(fdebug,"UL_Ktp_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
		   __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(UL_DEBUG){
    out = fopen(fdebug,"a");
  }

  for (int a=0; a<nne; a++){
    for(int w=0; w<len_w; w++){
      for(int g=0; g<ndn; g++){
	const double* const ptrST_wg = &ST[idx_4_gen(off_w+w,g,0,0,
						     nne_t,ndn,ndn,ndn)];
	Kpu[idx_K_gen(a,0,w,g,nne,
		      1,len_w,ndn)] += (jj*wt*Na[a]*Jr*Jn*Upp
					*cblas_ddot(9,Fr_It,1,
						    ptrST_wg,1));
      }
    }
  }
 
  if(UL_DEBUG){
    PGFEM_fprintf(out,"K\n");
    print_array_d(out,Kpu,len_w*ndn*nne,nne,len_w*ndn);

  }

  free(B);
  free(Fr_It);
  if(UL_DEBUG) fclose(out);
} /* Kpu */

void damage_UL_Kpu_at_ip(double *Kpu,
			 const int nne,
			 const int nne_t,
			 const double *Na,
			 const double *ST,
			 const double *Fr_I,
			 const double Jn,
			 const double Jr,
			 const double kappa,
			 const double Upp,
			 const double UP, /* n+1 */
			 const double P,
			 const double *SS,
			 const double omega, 
			 const double jj,
			 const double wt,
			 const int TYPE)
{
  /* offsets for different matrices */
  int len_w, off_w;
  char fdebug[20];

  int err_rank;
  PGFEM_Error_rank(&err_rank);

  switch (TYPE){
  case 0: /* Kpu */
    len_w = nne;
    off_w = 0;
    sprintf(fdebug,"dam_UL_Kup_debug_%d.log",err_rank);
    break;

  case 1: /* Kpb */
    len_w = 1;
    off_w = nne;
    sprintf(fdebug,"dam_UL_Ktp_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
		   __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(UL_DEBUG){
    out = fopen(fdebug,"a");
  }

  double  *Fr_It;
  Fr_It = aloc1(9);
  transpose(Fr_It,Fr_I,3,3);

  /* compute the kernel */
  double *kernel = aloc1(9);
  for(int i=0; i<9; i++){
    kernel[i] = (Jr*Jn*(1.-omega)*Upp*Fr_It[i]
		 - (UP-P/kappa)*SS[i]);
  }

  for (int a=0; a<nne; a++){
    for(int w=0; w<len_w; w++){
      for(int g=0; g<ndn; g++){
	const double* const ptrST_wg = &ST[idx_4_gen(off_w+w,g,0,0,
						     nne_t,ndn,ndn,ndn)];
	Kpu[idx_K_gen(a,0,w,g,nne,
		      1,len_w,ndn)] += jj*wt* Na[a] * cblas_ddot(9,kernel,1,
								 ptrST_wg,1);
      }
    }
  }
 
  if(UL_DEBUG){
    PGFEM_fprintf(out,"K\n");
    print_array_d(out,Kpu,len_w*ndn*nne,nne,len_w*ndn);

  }

  free(Fr_It);
  free(kernel);
  if(UL_DEBUG) fclose(out);
} /* Kpu */

void UL_Kpp_at_ip(double *Kpp,
		  const int nne,
		  const int nne_t,
		  const double *Na,
		  const double kappa,
		  const double jj,
		  const double wt)
{
  FILE *out;
  if(UL_DEBUG){
    int err_rank;
    char fname[50];
    PGFEM_Error_rank(&err_rank);
    sprintf(fname,"UL_Kpp_debug_%d.log",err_rank);
    out = fopen(fname,"a");
  }

  //PGFEM_printf("%12.12e\n",Upp);

  for (int a=0; a<nne; a++){
    for (int w=0; w<nne; w++){
      Kpp[idx_K(a,0,w,0,nne,1)] -= (Na[a]*Na[w]*jj*wt)/(kappa);
    }
  }

  if(UL_DEBUG){
    PGFEM_fprintf(out,"Kpp\n");
    print_array_d(out,Kpp,nne*nne,nne,nne);
    fclose(out);
  }
} /* Kpp */

void damage_UL_Kpp_at_ip(double *Kpp,
			 const int nne,
			 const int nne_t,
			 const double *Na,
			 const double kappa,
			 const double omega,
			 const double HH,
			 const double jj,
			 const double wt)
{
  FILE *out;
  if(UL_DEBUG){
    int err_rank;
    char fname[50];
    PGFEM_Error_rank(&err_rank);
    sprintf(fname,"dam_UL_Kpp_debug_%d.log",err_rank);
    out = fopen(fname,"a");
  }

  //PGFEM_printf("%12.12e\n",Upp);

  for (int a=0; a<nne; a++){
    for (int w=0; w<nne; w++){
      Kpp[idx_K(a,0,w,0,nne,1)] -= ((1.-omega)/kappa + HH)*(Na[a]*Na[w]*jj*wt);
    }
  }

  if(UL_DEBUG){
    PGFEM_fprintf(out,"Kpp\n");
    print_array_d(out,Kpp,nne*nne,nne,nne);
    fclose(out);
  }
} /* Kpp */

void UL_Ru_at_ip(double *Ru,
		 const int nne,
		 const int nne_t,
		 const double *ST,
		 const double *Fn,
		 const double *Fr,
		 const double *Fr_I,
		 const double Jn,
		 const double Jr,
		 const double *S,
		 const double pres,
		 const double jj,
		 const double wt,
		 const int TYPE)
{
  double *AA_ab, *BB, *CC, *Fr_It;
  AA_ab = aloc1(9);
  BB = aloc1(9);
  CC = aloc1(9);
  Fr_It = aloc1(9);
  transpose(Fr_It,Fr_I,3,3);

  /* offsets for different matrices */
  int len_a, off_a;
  char fdebug[20];

  int err_rank;
  PGFEM_Error_rank(&err_rank);

  switch (TYPE){
  case 0: /* Ru */
    len_a = nne;
    off_a = 0;
    sprintf(fdebug,"UL_Ru_debug_%d.log",err_rank);
    break;

  case 1: /* Rb */
    len_a = 1;
    off_a = nne;
    sprintf(fdebug,"UL_Rb_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
		   __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(UL_DEBUG){
    out = fopen(fdebug,"a");
    PGFEM_fprintf(out,"****************************************\n");
    PGFEM_fprintf(out,"Fr\n");
    print_array_d(out,Fr,9,1,9);
    PGFEM_fprintf(out,"Fr_I\n");
    print_array_d(out,Fr_I,9,1,9);
  }

  /***  Compute Fn S Fn' ***/
  for(int i=0; i<ndn; i++){
    for(int j=0; j<ndn; j++){
      CC[idx_2(i,j)] = 0.0;
      for(int k=0; k<ndn; k++){
	for(int l=0; l<ndn; l++){
	  CC[idx_2(i,j)] = (CC[idx_2(i,j)] 
			    + Fn[idx_2(i,k)]*S[idx_2(k,l)]
			    *Fn[idx_2(j,l)]);
	}
      }
    }
  }


  if(UL_DEBUG){
    PGFEM_fprintf(out,"Fn S Fn'\n");
    print_array_d(out,CC,9,1,9);
  }

  for (int a=0; a<len_a; a++){
    for (int b=0; b<ndn; b++){
      /*** compute AA{ab} = F^t ST{ab} ***/
      const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
						   nne_t,ndn,ndn,ndn)];
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		  3,3,3,1.0,Fr,3,ptrST_ab,3,0.0,AA_ab,3);

      if(UL_DEBUG){
	PGFEM_fprintf(out,"ST\n");
	print_array_d(out,ptrST_ab,9,1,9);
	PGFEM_fprintf(out,"Fr' ST\n");
	print_array_d(out,AA_ab,9,1,9);
      }

      for(int i=0; i<9; i++){
	Ru[a*ndn+b] += ((AA_ab[i]*CC[i]/Jn + pres*Jr
			 *Fr_It[i]*ptrST_ab[i])*jj*wt);
      }

    }
  }
  free(AA_ab);
  free(BB);
  free(CC);
  free(Fr_It);

  if(UL_DEBUG){
    print_array_d(out,Ru,len_a*ndn,len_a*ndn,1);
    fclose(out);
  }

}/* Ru */

void UL_Rp_at_ip(double *Rp,
		 const int nne,
		 const double *Na,
		 const double kappa,
		 const double Up,
		 const double pres,
		 const double jj,
		 const double wt)
{
  FILE *out;
  if(UL_DEBUG){
    int err_rank;
    char fname[50];
    PGFEM_Error_rank(&err_rank);
    sprintf(fname,"UL_Rp_debug_%d.log",err_rank);
    out = fopen(fname,"a");
  }

  for (int a=0; a<nne; a++){
    Rp[a] += (Up - pres/kappa)*Na[a]*jj*wt;
  }

  if(UL_DEBUG){
    print_array_d(out,Rp,nne,nne,1);
    fclose(out);
  }
}/* Rp */
