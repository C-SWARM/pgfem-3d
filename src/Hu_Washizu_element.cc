/* HEADER */
/**
 * AUTHORS:
 * Matthew Mosby
 */
#include "Hu_Washizu_element.h"
#include <string.h>
#include "mkl_cblas.h"

#include "PGFEM_mpi.h"
#include "PGFEM_io.h"
#include "allocation.h"
#include "utils.h"
#include "index_macros.h"

#ifndef HW_DEBUG
#define HW_DEBUG 0
#endif

/* This is TEMPORARY for testing */
#define JAC Jn

static const int ndn = 3;

void HW_Kuu_at_ip(double *Kuu,
          const int nne,
          const int nne_t,
          const double *ST,
          const double *Fn,
          const double *Fr,
          const double *Fr_I,
          const double Jn,
          const double Tn,
          const double Jr,
          const double pres,
          const double *S,
          const double *L,
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
  LAA_wg = aloc1(9);
  Fr_It = aloc1(9);
  transpose(Fr_It,Fr_I,3,3);

  double *zero;
  zero = aloc1(9);

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
    sprintf(fdebug,"HW_Kuu_debug_%d.log",err_rank);
    break;

  case 1: /* Kbb */
    len_a = 1;
    len_w = 1;
    off_a = nne;
    off_w = nne;
    sprintf(fdebug,"HW_Kbb_debug_%d.log",err_rank);
    break;

  case 2: /* Kub */
    len_a = nne;
    len_w = 1;
    off_a = 0;
    off_w = nne;
    sprintf(fdebug,"HW_Kub_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
        __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(HW_DEBUG) out = fopen(fdebug,"a");

  for (int a=0; a<len_a; a++){
    for (int b=0; b<ndn; b++){
      if(HW_DEBUG){
        PGFEM_fprintf(out,"=============== Row %d ===============\n",a*ndn+b);
      }

      const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
                           nne_t,ndn,ndn,ndn)];

      /*** compute AA{ab} = Fn' sym(Fr'ST{ab}) Fn ***/
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,Fr,3,ptrST_ab,3,0.0,B,3);
      symmetric_part(C,B,3);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,Fn,3,C,3,0.0,B,3);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
          3,3,3,1.0,B,3,Fn,3,0.0,AA_ab,3);

      if(HW_DEBUG){
        PGFEM_fprintf(out,"Fn' sym(Fr'ST{ab}) Fn\n");
        print_array_d(out,AA_ab,9,3,3);
      }

      for (int w=0; w<len_w; w++){
    for (int g=0; g<ndn; g++){

      const double* const ptrST_wg = &ST[idx_4_gen(off_w+w,g,0,0,
                               nne_t,ndn,ndn,ndn)];

      /*** compute AA{wg} = Fn' sym(Fr'ST{wg}) Fn ***/
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,Fr,3,ptrST_wg,3,0.0,B,3);
      symmetric_part(C,B,3);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,Fn,3,C,3,0.0,B,3);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              3,3,3,1.0,B,3,Fn,3,0.0,AA_wg,3);

      /* Compute L:AA{w,g} */
      for (int i=0; i<ndn; i++){
        for (int j=0; j<ndn; j++){
          LAA_wg[idx_2(i,j)] = cblas_ddot(9,&L[idx_4(i,j,0,0)],1,
                          AA_wg,1);
        }
      }

      if(HW_DEBUG){
        PGFEM_fprintf(out,"LAA = L:AA_wg\n");
        print_array_d(out,LAA_wg,9,3,3);
      }

      Kuu[idx_K_gen(a,b,w,g,
            len_a,ndn,
            len_w,ndn)] += (jj*wt/JAC
                    *cblas_ddot(9,AA_ab,1,LAA_wg,1));

      /*** Kuu Geometric Stiffness ***/
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,ptrST_wg,3,ptrST_ab,3,0.0,B,3);
      symmetric_part(C,B,3);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,Fn,3,C,3,0.0,D,3);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              3,3,3,1.0,D,3,Fn,3,0.0,B,3);

      if(HW_DEBUG){
        PGFEM_fprintf(out,"Fn' sym(ST(w,g)'ST(a,b)) Fn\n");
        print_array_d(out,B,9,3,3);
      }

      /* Push forward on S */
      Kuu[idx_K_gen(a,b,w,g,
            len_a,ndn,
            len_w,ndn)] += jj*wt/JAC*cblas_ddot(9,S,1,B,1);

      double PGS1, PGS2, PGS3;
      /*** Kuu Pressure Geometric Stiffness ***/
      /* PGS1 */
      PGS1 = cblas_ddot(9,Fr_It,1,ptrST_wg,1);

      /* PGS2 */
      PGS2 = cblas_ddot(9,Fr_It,1,ptrST_ab,1);


      /* PGS3 */

      cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,
              3,3,3,1.0,ptrST_wg,3,Fr_I,3,0.0,B,3);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              3,3,3,1.0,Fr_I,3,ptrST_ab,3,0.0,C,3);
      PGS3 = cblas_ddot(9,B,1,C,1);

      Kuu[idx_K_gen(a,b,w,g,
            len_a,ndn,
            len_w,ndn)] += Jn/JAC*pres*Jr*(PGS1*PGS2-PGS3)*jj*wt;

    } /* g */
      } /* w */
    } /* b */
  } /* a */

  if(HW_DEBUG){
    PGFEM_fprintf(out,"K\n");
    print_array_d(out,Kuu,len_a*len_w*ndn*ndn,len_a*ndn,len_w*ndn);
  }

  free(AA_ab);
  free(AA_wg);
  free(B);
  free(C);
  free(D);
  free(LAA_wg);
  free(Fr_It);
  free(zero);

  if(HW_DEBUG) fclose(out);
} /*** Kuu ***/

void HW_Kup_at_ip(double *Kup,
          const int nne,
          const int nne_t,
          const int nPres,
          const double *Np,
          const double *ST,
          const double *Fr_I,
          const double Jn,
          const double Tn,
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
    sprintf(fdebug,"HW_Kup_debug_%d.log",err_rank);
    break;

  case 1: /* Kbp */
    len_a = 1;
    off_a = nne;
    sprintf(fdebug,"HW_Kbp_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
        __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(HW_DEBUG){
    out = fopen(fdebug,"a");
  }

  for (int a=0; a<len_a; a++){
    for(int b=0; b<ndn; b++){
      for(int w=0; w<nPres; w++){
    const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
                             nne_t,ndn,ndn,ndn)];
    Kup[idx_K_gen(a,b,w,0,len_a,
              ndn,nPres,1)] += (jj*wt*Np[w]*Jr*Jn/JAC
                    *cblas_ddot(9,Fr_It,1,
                            ptrST_ab,1));
      }
    }
  }

  if(HW_DEBUG){
    PGFEM_fprintf(out,"K\n");
    print_array_d(out,Kup,len_a*ndn*nPres,len_a*ndn,nPres);

  }

  free(B);
  free(Fr_It);
  if(HW_DEBUG) fclose(out);
} /* Kup */

void HW_Ktt_at_ip(double *Ktt,
          const int nVol,
          const double *Nt,
          const double Tn,
          const double Jn,
          const double kappa,
          const double Upp,
          const double jj,
          const double wt)
{
  FILE *out;

  if(HW_DEBUG){
    char fname[50];
    int err_rank;
    PGFEM_Error_rank(&err_rank);
    sprintf(fname,"HW_Ktt_debug_%d.log",err_rank);
    out = fopen(fname,"a");
    PGFEM_fprintf(out,"*********************************************\n");
  }

  for(int a=0; a<nVol; a++){
    for(int w=0; w<nVol; w++){
      Ktt[idx_K(a,0,w,0,nVol,1)] += (jj*wt*Tn*Tn/JAC*kappa*Upp*Nt[a]*Nt[w]);
    }
  }

  if(HW_DEBUG){
    print_array_d(out,Ktt,nVol*nVol,nVol,nVol);
    fclose(out);
  }
}/* Ktt */

void HW_Kpt_at_ip(double *Kpt,
          const int nPres,
          const double *Np,
          const int nVol,
          const double *Nt,
          const double Tn,
          const double Jn,
          const double jj,
          const double wt)
{
  FILE *out;

  if(HW_DEBUG){
    char fname[50];
    int err_rank;
    PGFEM_Error_rank(&err_rank);
    sprintf(fname,"HW_Kpt_debug_%d.log",err_rank);
    out = fopen(fname,"a");
    PGFEM_fprintf(out,"*********************************************\n");
  }

  for(int a=0; a<nPres; a++){
    for(int w=0; w<nVol; w++){
      Kpt[idx_K_gen(a,0,w,0,nPres,1,nVol,1)] -= Tn/JAC*jj*wt*Np[a]*Nt[w];
    }
  }

  if(HW_DEBUG){
    print_array_d(out,Kpt,nPres*nVol,nPres,nVol);
    fclose(out);
  }
}/* Ktt */

void HW_Ru_at_ip(double *Ru,
         const int nne,
         const int nne_t,
         const double *ST,
         const double *Fn,
         const double *Fr,
         const double *Fr_I,
         const double Jn,
         const double Tn,
         const double Jr,
         const double Tr,
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
    sprintf(fdebug,"HW_Ru_debug_%d.log",err_rank);
    break;

  case 1: /* Rb */
    len_a = nne_t-nne;
    off_a = nne;
    sprintf(fdebug,"HW_Rb_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
        __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(HW_DEBUG){
    out = fopen(fdebug,"a");
    PGFEM_fprintf(out,"****************************************\n");
    PGFEM_fprintf(out,"Fr\n");
    print_array_d(out,Fr,9,1,9);
    PGFEM_fprintf(out,"Fr_I\n");
    print_array_d(out,Fr_I,9,1,9);
  }

  for (int a=0; a<len_a; a++){
    for (int b=0; b<ndn; b++){
      const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
                           nne_t,ndn,ndn,ndn)];

      /*** compute AA{ab} = Fn' sym(Fr'ST{ab}) Fn ***/
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,Fr,3,ptrST_ab,3,0.0,BB,3);
      symmetric_part(CC,BB,3);

      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,Fn,3,CC,3,0.0,BB,3);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
          3,3,3,1.0,BB,3,Fn,3,0.0,AA_ab,3);

      if(HW_DEBUG){
    PGFEM_fprintf(out,"ST\n");
    print_array_d(out,ptrST_ab,9,1,9);
    PGFEM_fprintf(out,"Fn' sym(Fr'ST{ab}) Fn\n");
    print_array_d(out,AA_ab,9,1,9);
      }

      Ru[a*ndn+b] += (cblas_ddot(9,AA_ab,1,S,1)*jj*wt/JAC + Jn*pres*Jr
                  *cblas_ddot(9,Fr_It,1,ptrST_ab,1)*jj*wt/JAC);

    }
  }
  free(AA_ab);
  free(BB);
  free(CC);
  free(Fr_It);

  if(HW_DEBUG){
    print_array_d(out,Ru,len_a*ndn,len_a*ndn,1);
    fclose(out);
  }

}/* Ru */

void HW_Ru1_at_ip(double *Ru,
          const int nne,
          const int nne_t,
          const double *ST,
          const double *Fn,
          const double *Fr,
          const double *Fr_I,
          const double Jn,
          const double Tn,
          const double Jr,
          const double Tr,
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
    sprintf(fdebug,"HW_Ru1_debug_%d.log",err_rank);
    break;

  case 1: /* Rb */
    len_a = nne_t-nne;
    off_a = nne;
    sprintf(fdebug,"HW_Rb1_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
        __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(HW_DEBUG){
    out = fopen(fdebug,"a");
    PGFEM_fprintf(out,"****************************************\n");
    PGFEM_fprintf(out,"Fr\n");
    print_array_d(out,Fr,9,1,9);
    PGFEM_fprintf(out,"Fr_I\n");
    print_array_d(out,Fr_I,9,1,9);
  }

  for (int a=0; a<len_a; a++){
    for (int b=0; b<ndn; b++){
      const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
                           nne_t,ndn,ndn,ndn)];

      /*** compute AA{ab} = Fn' sym(Fr'ST{ab}) Fn ***/
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,Fr,3,ptrST_ab,3,0.0,BB,3);
      symmetric_part(CC,BB,3);

      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,Fn,3,CC,3,0.0,BB,3);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
          3,3,3,1.0,BB,3,Fn,3,0.0,AA_ab,3);

      if(HW_DEBUG){
    PGFEM_fprintf(out,"ST\n");
    print_array_d(out,ptrST_ab,9,1,9);
    PGFEM_fprintf(out,"Fn' sym(Fr'ST{ab}) Fn\n");
    print_array_d(out,AA_ab,9,1,9);
      }

      Ru[a*ndn+b] += cblas_ddot(9,AA_ab,1,S,1)*jj*wt/JAC;
    }
  }
  free(AA_ab);
  free(BB);
  free(CC);
  free(Fr_It);

  if(HW_DEBUG){
    print_array_d(out,Ru,len_a*ndn,len_a*ndn,1);
    fclose(out);
  }

}/* Ru1 */


void HW_Ru2_at_ip(double *Ru,
          const int nne,
          const int nne_t,
          const double *ST,
          const double *Fn,
          const double *Fr,
          const double *Fr_I,
          const double Jn,
          const double Tn,
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
    sprintf(fdebug,"HW_Ru2_debug_%d.log",err_rank);
    break;

  case 1: /* Rb */
    len_a = nne_t-nne;
    off_a = nne;
    sprintf(fdebug,"HW_Rb2_debug_%d.log",err_rank);
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
        __func__);
    PGFEM_Abort();
    abort();
  }

  FILE *out;
  if(HW_DEBUG){
    out = fopen(fdebug,"a");
    PGFEM_fprintf(out,"****************************************\n");
    PGFEM_fprintf(out,"Fr\n");
    print_array_d(out,Fr,9,1,9);
    PGFEM_fprintf(out,"Fr_I\n");
    print_array_d(out,Fr_I,9,1,9);
  }

  for (int a=0; a<len_a; a++){
    for (int b=0; b<ndn; b++){
      const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
                           nne_t,ndn,ndn,ndn)];

      Ru[a*ndn+b] += (Jn*pres*Jr*cblas_ddot(9,Fr_It,1,ptrST_ab,1))*jj*wt/JAC;
    }
  }
  free(AA_ab);
  free(BB);
  free(CC);
  free(Fr_It);

  if(HW_DEBUG){
    print_array_d(out,Ru,len_a*ndn,len_a*ndn,1);
    fclose(out);
  }

}/* Ru2 */

void HW_Rp_at_ip(double *Rp,
         const int nPres,
         const double *Np,
         const double Jn,
         const double Jr,
         const double Tn,
         const double Tr,
         const double jj,
         const double wt)
{
  FILE *out;
  if(HW_DEBUG){
    int err_rank;
    char fname[50];
    PGFEM_Error_rank(&err_rank);
    sprintf(fname,"HW_Rp_debug_%d.log",err_rank);
    out = fopen(fname,"a");
  }

  for (int a=0; a<nPres; a++){
    Rp[a] += (Jn*Jr-Tn*Tr)*Np[a]*jj*wt/JAC;
  }

  if(HW_DEBUG){
    print_array_d(out,Rp,nPres,1,nPres);
    fclose(out);
  }
}/* Rp */

void HW_Rt_at_ip(double *Rt,
         const int nVol,
         const double *Nt,
         const double Tn,
         const double Jn,
         const double pres,
         const double kappa,
         const double Up,
         const double jj,
         const double wt)
{
  FILE *out;
  if(HW_DEBUG){
    int err_rank;
    char fname[50];
    PGFEM_Error_rank(&err_rank);
    sprintf(fname,"HW_Rt_debug_%d.log",err_rank);
    out = fopen(fname,"a");
  }

  for(int a=0; a<nVol; a++){
    Rt[a] += Tn/JAC*(kappa*Up-pres)*Nt[a]*jj*wt;
  }

  if(HW_DEBUG){
    print_array_d(out,Rt,nVol,1,nVol);
    fclose(out);
  }
}

void debug_NH_Ru_at_ip(double *Ru,
               const int nne,
               const int nne_t,
               const double *ST,
               const double *Fn,
               const double *Fr,
               const double *Fr_I,
               const double Jn,
               const double Tn,
               const double Jr,
               const double Tr,
               const double *S,
               const double pres,
               const double G,
               const double jj,
               const double wt,
               const int TYPE)
{
  /* offsets for different matrices */
  int len_a, off_a;
  int err_rank;
  PGFEM_Error_rank(&err_rank);

  switch (TYPE){
  case 0: /* Ru */
    len_a = nne;
    off_a = 0;
    break;

  case 1: /* Rb */
    len_a = nne_t-nne;
    off_a = nne;
    break;

  default: /* default error out */
    PGFEM_printerr("ERROR, unrecognised type in %s\n",
        __func__);
    PGFEM_Abort();
    abort();
  }

  double *AA_ab, *BB, *CC, *Fr_It;
  AA_ab = aloc1(9);
  BB = aloc1(9);
  CC = aloc1(9);
  Fr_It = aloc1(9);
  transpose(Fr_It,Fr_I,3,3);

  for (int a=0; a<len_a; a++){
    for (int b=0; b<ndn; b++){
      const double* const ptrST_ab = &ST[idx_4_gen(off_a+a,b,0,0,
                           nne_t,ndn,ndn,ndn)];

      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,Fr,3,ptrST_ab,3,0.0,BB,3);
      symmetric_part(CC,BB,3);

      /* double alpha = -1./3.*cblas_ddot(ndn*ndn,Fr_It,1,ptrST_ab,1); */
      /* cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, */
      /*          3,3,3,alpha,Fr,3,Fr,3,1.0,CC,3); */

      /* alpha = G*pow(det3x3(Fn)*Jr,-2./3.); */
      /* cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, */
      /*          3,3,3,alpha,CC,3,Fn,3,0.0,BB,3); */

      /* push forward of S */
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
          3,3,3,1.0,Fn,3,S,3,0.0,BB,3);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
          3,3,3,1.0,BB,3,Fn,3,0.0,AA_ab,3);

      Ru[a*ndn+b] += jj*wt/JAC;/* (cblas_ddot(9,AA_ab,1,CC,1)*jj*wt/JAC + Jn*pres*Jr */
                 /*  *cblas_ddot(9,Fr_It,1,ptrST_ab,1)*jj*wt/JAC); */

    }
  }
  free(AA_ab);
  free(BB);
  free(CC);
  free(Fr_It);
}
