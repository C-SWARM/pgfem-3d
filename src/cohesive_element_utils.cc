/* HEADER */
/** This file contains the utilities for cohesive elements */
#include "pgfem3d/Communication.hpp"
#include "cohesive_element_utils.h"
#include <math.h>

#include "PGFEM_io.h"
#include "allocation.h"
#include "elem3d.h"
#include "matice.h"
#include "tensors.h"

using namespace pgfem3d;

static const int ndn = 3;

/* ==== STATIC FUNCTION PROTO ==== */
static void int_tria_I (double *gk,
            double *ge,
            double *w);

/* ==== PUBLIC FUNCTIONS ==== */
long int_pointC (long nne)
{
  switch(nne){
  case 3: return 3;
  case 4: return 4;
  default:
    PGFEM_printerr("Unrecognized cohesive element! (%s:%s)\n",
        __FILE__,__func__);
    PGFEM_Abort();
    return -1;
  }
}

void integrateC (long nne,
         long *II,
         long *JJ,
         double *gk,
         double *ge,
         double *w)

{
  *II = *JJ = 0;
  switch(nne){
  case 3:
    { /* Linear triangle */
      int_tria_I (gk,ge,w);
      *II = 1;
      *JJ = 3;
      break;
    }
  case 4:
    { /* Quadrilateral */
      intpoints_2 (gk,w);
      *II = 2;
      *JJ = 2;
      break;
    }
  default:
    PGFEM_printerr("Unrecognized cohesive element! (%s:%s)\n",
        __FILE__,__func__);
    PGFEM_Abort();
    break;
  }
}

void shape_2DC (long nne,
        double ksi,
        double eta,
        double *N)
{
  /*
    ksi - integration point
    eta - integration point
    N - shape functions
  */
  switch(nne){
  case 1: /* constant */
    N[0] = 1.0;
    break;
  case 3:/* Triangle */
    N[0] = 1.-ksi-eta;
    N[1] = ksi;
    N[2] = eta;
    break;
  case 4:/* Quadrilateral */
    N[0] = 0.25*(1.-ksi)*(1.-eta);
    N[1] = 0.25*(1.+ksi)*(1.-eta);
    N[2] = 0.25*(1.+ksi)*(1.+eta);
    N[3] = 0.25*(1.-ksi)*(1.+eta);
    break;
  default:
    PGFEM_printerr("ERROR: unrecognized 2D element!\n");
    PGFEM_Abort();
  }
}

void dN_2D (long nne,
        double ksi,
        double eta,
        double *N_ksi,
        double *N_eta)
{
  switch(nne){
  case 3:
    {
      N_ksi[0] = -1.0; N_eta[0] = -1.0;
      N_ksi[1] =  1.0; N_eta[1] =  0.0;
      N_ksi[2] =  0.0; N_eta[2] =  1.0;
      break;
    }
  case 4:
    {
      N_ksi[0]  = 0.25*(eta - 1.0);  N_eta[0]  = 0.25*(ksi - 1.0);
      N_ksi[1]  = 0.25*(1.0 - eta);  N_eta[1]  = 0.25*(-1.0 - ksi);
      N_ksi[2]  = 0.25*(1.0 + eta);  N_eta[2]  = 0.25*(1.0 + ksi);
      N_ksi[3]  = 0.25*(-1.0 - eta); N_eta[3]  = 0.25*(1.0 - ksi);
      break;
    }
  case 6:
    {
      N_ksi[0] =     4.0*(eta+ksi)-3.0; N_eta[0] =     4.0*(eta + ksi)-3.0;
      N_ksi[1] =                   0.0; N_eta[1] =             4.0*ksi-1.0;
      N_ksi[2] =           4.0*eta-1.0; N_eta[2] =                     0.0;
      N_ksi[3] =              -4.0*ksi; N_eta[3] = 4.0*(1.0 - 2.0*ksi-eta);
      N_ksi[4] =               4.0*ksi; N_eta[4] =                 4.0*eta;
      N_ksi[5] = 4.0*(1.0-ksi-2.0*eta); N_eta[5] =                -4.0*eta;
      break;
    }
  default:
    PGFEM_printerr("Unrecognized cohesive element! (%s:%s)\n",
        __FILE__,__func__);
    PGFEM_Abort();
    break;
  }
}

void dN_dxb (long nne,
         double ksi,
         double eta,
         double *e1,
         double *e2,
         double *n,
         double ***Nxb)
{
  long Q,P,R,S;
  double *N_ksi,*N_eta,*a,***e,BB,***CC;

  N_ksi = aloc1 (nne);
  N_eta = aloc1 (nne);
  a = aloc1 (3);
  e = aloc3 (3,3,3);
  CC = aloc3 (3,ndn,nne);

  e[0][1][2] = e[1][2][0] = e[2][0][1] = +1.0;
  e[0][2][1] = e[2][1][0] = e[1][0][2] = -1.0;

  /* Derivative of the shape functions */
  dN_2D (nne,ksi,eta,N_ksi,N_eta);

  cross_product (e1,e2,a);
  BB = sqrt (ss (a,a,3));

  for (Q=0;Q<ndn;Q++){
    for (P=0;P<ndn;P++){
      for (R=0;R<nne;R++){
    CC[Q][P][R] = 0.0;
    for (S=0;S<3;S++){
      CC[Q][P][R] += ( e[Q][P][S]*N_ksi[R]*e2[S]
               + e[Q][S][P]*e1[S]*N_eta[R] );
    }
      }
    }
  }

  for (Q=0;Q<ndn;Q++){
    for (P=0;P<ndn;P++){
      for (R=0;R<nne;R++){
    Nxb[Q][P][R] = CC[Q][P][R]/BB;  /* aer */
    for (S=0;S<3;S++){
      Nxb[Q][P][R] -= 1./BB*n[Q]*n[S]*CC[S][P][R];
    }
      }
    }
  }

  /**************************************************************/
  /* ORTIZ : This is the same as I have only 1./2. is different */
  /**************************************************************/

  /*
    U = aloc2 (3,3);

    for (Q=0;Q<ndn;Q++){
    for (P=0;P<ndn;P++){
    U[Q][P] = 0.0;
    for (S=0;S<3;S++){
    U[Q][P] += n[S]*e[Q][P][S];
    }
    }
    }

    for (Q=0;Q<ndn;Q++){
    for (P=0;P<ndn;P++){
    for (R=0;R<nne;R++){
    Nxb[Q][P][R] = 0.0;
    for (S=0;S<3;S++){
    Nxb[Q][P][R] += (e[P][S][Q]/(2.*BB)*(e2[S]*N_ksi[R] - e1[S]*N_eta[R])
    - n[Q]*U[P][S]/(2.*BB)*(e2[S]*N_ksi[R] - e1[S]*N_eta[R]));
    }
    }
    }
    }

    dealoc2 (U,3);
  */

  dealoc1 (N_ksi);
  dealoc1 (N_eta);
  dealoc1 (a);
  dealoc3 (e,3,3);
  dealoc3 (CC,3,ndn);
}

double dN3_xy (double ksi,
           double eta,
           long nne,
           double *x,
           double *y,
           double *z,
           double *N_x,
           double *N_y)
{
  double *N_ksi,*N_eta,*dx,*dy,J;
  long i;

  N_ksi = aloc1 (nne);
  N_eta = aloc1 (nne);
  dx = aloc1 (2);
  dy = aloc1 (2);

  dN_2D (nne,ksi,eta,N_ksi,N_eta);

  for (i=0;i<nne;i++){
    dx[0] += N_ksi[i]*x[i];
    dx[1] += N_eta[i]*x[i];

    dy[0] += N_ksi[i]*y[i];
    dy[1] += N_eta[i]*y[i];
  }

  J = dx[0]*dy[1] - dx[1]*dy[0];
  if (J <= 0.0) {
    PGFEM_printf ("Negative determinant is isoparametric mapping"
        " (COHESIVE) : J = %12.12e || Bye Bye\n",J);
    PGFEM_Abort();
  }

  for (i=0;i<nne;i++){
    N_x[i] = 1./J * ( dy[1]*N_ksi[i] - dy[0]*N_eta[i]);
    N_y[i] = 1./J * (-dx[1]*N_ksi[i] + dx[0]*N_eta[i]);
  }

  dealoc1 (N_ksi);
  dealoc1 (N_eta);
  dealoc1 (dx);
  dealoc1 (dy);

  return (J);
}

void mean_map (long nne,
           double *x,
           double *y,
           double *z,
           double *r_u,
           double *xb,
           double *yb,
           double *zb)
{
  long M;

  /* Mean mapping */
  for (M=0;M<nne;M++){
    xb[M] = (0.5*(x[nne+M] + x[M])
         + 0.5*(r_u[nne*ndn + ndn*M + 0] + r_u[ndn*M + 0]));
    yb[M] = (0.5*(y[nne+M] + y[M])
         + 0.5*(r_u[nne*ndn + ndn*M + 1] + r_u[ndn*M + 1]));
    zb[M] = (0.5*(z[nne+M] + z[M])
         + 0.5*(r_u[nne*ndn + ndn*M + 2] + r_u[ndn*M + 2]));
  }
}

void base_vec (const long nne,
           const double ksi,
           const double eta,
           const double *x,
           const double *y,
           const double *z,
           double *e1,
           double *e2,
           double *e2h,
           double *n,
           const int myrank)
{
  long i;
  double *N_ksi,*N_eta;

  N_ksi = aloc1 (nne);
  N_eta = aloc1 (nne);

  /* Derivative of the shape functions */
  dN_2D (nne,ksi,eta,N_ksi,N_eta);

  for (i=0;i<3;i++){e1[i] = 0.0; e2h[i] = 0.0;}

  for (i=0;i<nne;i++){
    e1[0] += N_ksi[i]*x[i];  e2h[0] += N_eta[i]*x[i];
    e1[1] += N_ksi[i]*y[i];  e2h[1] += N_eta[i]*y[i];
    e1[2] += N_ksi[i]*z[i];  e2h[2] += N_eta[i]*z[i];
  }
  /* Compute normal */
  cross_product (e1,e2h,n);
  /* Norm base vectors */
  nor_vec_serial (e1,3,myrank);
  nor_vec_serial (e2h,3,myrank);
  nor_vec_serial (n,3,myrank);
  /* Compute e2 */
  cross_product (n,e1,e2);

  dealoc1 (N_ksi); dealoc1 (N_eta);
}

void tran_coord (long nne,
         double *x,
         double *y,
         double *z,
         double *e1,
         double *e2,
         double *n,
         double *xl,
         double *yl,
         double *zl,
         long TO)
{
  /*

    TO = 0 || L -> G :: vg = T   . vl
    TO = 1 || G -> L :: vl = T^T . vg

    T is assembled by columns (stored in 1D array by rows)

  */

  long i;
  double *T,*XL,*XG;

  T = aloc1 (9);
  XL = aloc1 (3);
  XG = aloc1 (3);

  T[0] = e1[0];
  T[1] = e2[0];
  T[2] = n[0];

  T[3] = e1[1];
  T[4] = e2[1];
  T[5] = n[1];

  T[6] = e1[2];
  T[7] = e2[2];
  T[8] = n[2];

  for (i=0;i<nne;i++){
    if (TO == 0){
      XL[0] = xl[i];
      XL[1] = yl[i];
      XL[2] = zl[i];
      mv (T,XL,XG,3,3);
      x[i] = XG[0];
      y[i] = XG[1];
      z[i] = XG[2];
    }
    if (TO == 1){
      XG[0] =  x[i];
      XG[1] =  y[i];
      XG[2] =  z[i];
      mtv(T,XG,XL,3,3);
      xl[i] = XL[0];
      yl[i] = XL[1];
      zl[i] = XL[2];
    }
  }

  dealoc1 (T);
  dealoc1 (XL);
  dealoc1 (XG);
}

/* compute the TOTAL jump across the interface at an integration
   point */
int get_jump (const long nne,
          const double *x,
          const double *y,
          const double *z,
          const double *r_u,
          const double *N,
          double *jump)
{
  int err = 0;
  long i,j;
  double xij{};

  for (i=0;i<ndn;i++){
    jump[i] = 0.0;
    for (j=0;j<nne;j++){
      if (i == 0){
    xij = x[nne+j] - x[j];
      }
      if (i == 1){
    xij = y[nne+j] - y[j];
      }
      if (i == 2){
    xij = z[nne+j] - z[j];
      }

      jump[i] += xij*N[j]
    + (r_u[nne*ndn + ndn*j + i] - r_u[ndn*j + i])*N[j];
    }
  }
  return err;
}


int get_jump_nt(double *jump_n,
        double *jump_t,
        const double *jump,
        const double *normal)
{
  int err = 0;

  *jump_n = *jump_t = 0.0;

  /* compute normal component */
  for(int i=0; i<ndn; i++){
    *jump_n += jump[i]*normal[i];
  }

  /* compute shear component */
  for(int i=0; i<ndn; i++){
    *jump_t += pow(jump[i]- *jump_n*normal[i],2);
  }
  *jump_t = sqrt(*jump_t);

  return err;
}

int get_eff_jump(double *eff_jump,
         const double jump_n,
         const double jump_t,
         const double beta)
{
  int err = 0;
  *eff_jump = sqrt(beta*beta*jump_t*jump_t
           + jump_n*jump_n);
  return err;
}

double get_Xxi (const double bet,
        const double *Xi,
        const double *n)
{
  double Xxi = 0.0;
  const double Xn = ss(Xi,n,ndn);
  double Xt2 = 0.0;
  for(int i=0; i<ndn; i++){
    Xt2 += (Xi[i]-Xn*n[i])*(Xi[i]-Xn*n[i]);
  }
  Xxi = sqrt(bet*bet*Xt2 + Xn*Xn);
  return (Xxi);
}

/* ==== STATIC FUNCTIONS ==== */
static void int_tria_I (double *gk,
            double *ge,
            double *w)

{
  static const int TYPE = 3;

  switch(TYPE){
  case 0:
    {
      /* 3-point formula || degree of precision 2 */
      w[0] = 0.333333333333333/2.;
      gk[0] = 0.666666666666667;
      ge[0] = 0.166666666666667;

      w[1] = 0.333333333333333/2.;
      gk[1] = 0.166666666666667;
      ge[1] = 0.166666666666667;

      w[2] = 0.333333333333333/2.;
      gk[2] = 0.166666666666667;
      ge[2] = 0.666666666666667;
      break;
    }
  case 1:
    {
      /* 4-point formula || degree of precision 3 */
      w[0] = -0.562500000000000/2.;
      gk[0] = 0.333333333333333;
      ge[0] = 0.333333333333333;

      w[1] =  0.520833333333333/2.;
      gk[1] = 0.600000000000000;
      ge[1] = 0.200000000000000;

      w[2] =  0.520833333333333/2.;
      gk[2] = 0.200000000000000;
      ge[2] = 0.200000000000000;

      w[3] =  0.520833333333333/2.;
      gk[3] = 0.200000000000000;
      ge[3] = 0.600000000000000;
      break;
    }
  case 2:
    {
      /* Newton-Cotes quadrature || degree of precision 1 */
      w[0] = 0.333333333333333/2.;
      gk[0] = 0.000000000000000;
      ge[0] = 0.000000000000000;

      w[1] = 0.333333333333333/2.;
      gk[1] = 0.000000000000000;
      ge[1] = 1.000000000000000;

      w[2] = 0.333333333333333/2.;
      gk[2] = 1.000000000000000;
      ge[2] = 0.000000000000000;
      break;
    }
  case 3:
    {
      /* Newton-Cotes quadrature || degree of precision 2 */
      w[0] = 0.333333333333333/2.;
      gk[0] = 0.500000000000000;
      ge[0] = 0.500000000000000;

      w[1] = 0.333333333333333/2.;
      gk[1] = 0.500000000000000;
      ge[1] = 0.000000000000000;

      w[2] = 0.333333333333333/2.;
      gk[2] = 0.000000000000000;
      ge[2] = 0.500000000000000;
      break;
    }
  default:
    {
      PGFEM_printerr("Unrecognized cohesive integration type! (%s:%s)\n",
          __FILE__,__func__);
      PGFEM_Abort();
      break;
    }
  }/* end switch */
}


/* ==== OLD/REPLACED FUNCTIONS ==== */
/* double get_t (long TYPE, */
/*        double Sc, */
/*        double Xxi, */
/*        double Xc, */
/*        double bet, */
/*        double kap, */
/*        double Xn, */
/*        double Xmax, */
/*        double tmax, */
/*        double nor_min) */
/* { */
/*   double txi,n; */

/*   /\* nor_min *= 1.e-3; *\/ */

/*   if (fil == 1) {/\* Elastic instability of rubber films *\/ */
/*     txi = Sc/(6.0*PII*bet*bet*bet)*(pow((bet/Xxi),3.0) - pow((bet/Xxi),9.0)); */
/*   } */
/*   else{ */
/*     if (Xn >= 0.0){/\* Tension *\/ */
/*       if ((Xxi + 1.e-10) >= Xmax || Xmax < Xc) {/\* Loading + elastic regime *\/ */
/*  if (TYPE == 0) txi = Sc*Xxi/Xc*exp(-1.*(Xxi-Xc)/Xc); /\* Needleman *\/ */
/*  if (TYPE == 1){ /\* Matous-Arciniega *\/ */
/*    if (Xxi <= kap*Xc) txi = Sc*Xxi/Xc*exp(-1.*(Xxi-Xc)/Xc); */
/*    if (Xc < Xxi && Xxi < kap*Xc) txi = Sc; */
/*    if (Xxi >= kap*Xc) { */
/*      n = kap*kap*PII/(4.0*(3.0-kap)*(3.0-kap)); */
/*      txi = Sc*exp(-n*(Xxi/(kap*Xc)-1.0)*(Xxi/(kap*Xc)-1.0)); */
/*    } */
/*  }/\* end TYPE 1 *\/ */
/*       }/\* end tension *\/ */
/*       else txi = Xxi*tmax/Xmax; /\* Unloading *\/ */
/*     } */
/*     else{ /\* Compression *\/ */
/*       /\* exponential penalty *\/ */
/*       /\* txi = Xxi*(Sc*(Xxi+Xc)/(Xc*Xc)*exp((Xxi+Xc)/Xc)); *\/ */

/*       /\* linear penalty *\/ */
/*       /\* txi = Xxi*Sc/Xc*exp(1.0)*MULT; *\/ */

/*       /\* quadradic penalty *\/ */
/*       txi = Xxi*Sc/Xc*exp(1.0) + MULT*Xxi*Xxi; */

/*     }/\* end compression *\/ */
/*   }/\* end fill == 0 *\/ */

/*   return (txi); */
/* } */

/* double get_dt (long TYPE, */
/*         double Sc, */
/*         double Xxi, */
/*         double Xc, */
/*         double bet, */
/*         double kap, */
/*         double Xn, */
/*         double Xmax, */
/*         double tmax, */
/*         double nor_min) */
/* { */
/*   double H1,n; */

/*   /\* nor_min *= 1.e-3; *\/ */

/*   if (fil == 1){/\* Elastic instability of rubber films *\/ */
/*     H1 = Sc/(2.0*PII*pow(Xxi,10.0))*(3.0*pow(bet,6.0) - pow(Xxi,6.0)); */
/*   } */
/*   else{ */
/*     if (Xn >= 0.0) {/\* Tension *\/ */
/*       if ((Xxi + 1.e-10) >= Xmax || Xmax < Xc){/\* Loading + elastic regime *\/ */
/*  if (TYPE == 0) */
/*    H1 = Sc*(Xc-Xxi)/(Xc*Xc)*exp(-1.*(Xxi-Xc)/Xc); /\* Needleman *\/ */
/*  if (TYPE == 1){ /\* Matous-Arciniega *\/ */
/*    if (Xxi <= kap*Xc) H1 = Sc*(Xc-Xxi)/(Xc*Xc)*exp(-1.*(Xxi-Xc)/Xc); */
/*    if (Xc < Xxi && Xxi < kap*Xc) H1 = 0.0; */
/*    if (Xxi >= kap*Xc){ */
/*      n = kap*kap*PII/(4.0*(3.0-kap)*(3.0-kap)); */
/*      H1 = ( -2.0*Sc*n/(kap*Xc)*(Xxi/(kap*Xc)-1.0) */
/*         *exp(-n*(Xxi/(kap*Xc)-1.0)*(Xxi/(kap*Xc)-1.0)) ); */
/*    } */
/*  }/\* end TYPE  1 *\/ */
/*       }/\* end tension *\/ */
/*       else H1 = tmax/Xmax; /\* Unloading *\/ */
/*     } */
/*     else {/\* Compression *\/ */
/*       /\* exponential penalty *\/ */
/*       /\* H1 = Sc*(3.*Xxi*Xc + Xc*Xc + Xxi*Xxi)/(Xc*Xc*Xc)*exp((Xxi+Xc)/Xc); *\/ */

/*       /\* linear penalty *\/ */
/*       /\* H1 = Sc/Xc*exp(1.0)*MULT; *\/ */

/*       /\* quadradic penalty *\/ */
/*       H1 = Sc/Xc*exp(1.0) + 2*MULT*Xxi; */
/*     }/\* end compression *\/ */
/*   } */

/*   return (H1); */
/* } */

