#include "PLC.h"
#include <math.h>
#include <string.h>

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef INCL_H
#include "incl.h"
#endif

double PLC_diff_g_ga (long ii,
              long ip,
              long nss,
              long mat,
              double dt,
              CRPL *crpl,
              double *GA,
              EPS *eps)
{
  long k,M,N;
  double diff=0.0,N1,N3,Gama,ga,Rof,Rom,Xi;
  const char *err[]={"inf","-inf","nan"};
  char str[500];

  ga = 0.0; for (k=0;k<nss;k++) ga += fabs(GA[k]);

  Gama = eps[ii].il[ip].GAMA + dt*ga;

  Rof = (crpl[mat].l1/crpl[mat].l2 + crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama))*(crpl[mat].l1/crpl[mat].l2 + crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama));
  Rom = 2.0*crpl[mat].c1*crpl[mat].c2/(crpl[mat].b*crpl[mat].l2)*(exp(-1./2.*crpl[mat].l2*Gama) - crpl[mat].c3);
  Xi = pow ((crpl[mat].b*Rom/(sqrt(Rof)*crpl[mat].to*ga)),1./3.);

  N1 = crpl[mat].c1*(crpl[mat].l1 + crpl[mat].l2*crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama))*exp(-1./2.*crpl[mat].l2*Gama);
  N3 = (dt*ga*crpl[mat].c1*(crpl[mat].b*Rom*(crpl[mat].l1 + crpl[mat].l2*crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama)) -
                2.*crpl[mat].c2*Rof)*exp(-1./2.*crpl[mat].l2*Gama) -
    2.0*crpl[mat].b*Rom*Rof)/(6.*pow((crpl[mat].b*Rom/(sqrt(Rof)*crpl[mat].to*ga)),2./3.)*pow(Rof,3./2.)*crpl[mat].to*ga*ga);

  diff = -1.0*dt*crpl[mat].nu*crpl[mat].b*N1/(6.0*sqrt(Rof)) + crpl[mat].fo*exp(-1.*Xi)*N3;

  sprintf (str,"%f",diff); for (N=0;N<3;N++){ M = 10; M = strcmp(err[N],str); if (M == 0) return (0.0);}

  return (diff);
}

long PLC_Inaccessible (long ii,
               long ip,
               double nor_min,
               long nss,
               long mat,
               double dt,
               CRPL *crpl,
               double *GA,
               double *GA1,
               EPS *eps,
               double **UU,
               double HAR)
{
  long j,k,M,N,INFO = 0;
  double N4,ga,GAGA{},nor,dga,R,Gama,Rof,Rom,Xi,PP[3][3],*GA2;
  const char *err[]={"inf","-inf","nan"};
  char str[500];

  GA2 = aloc1 (nss);

  ga = 0.0; for (k=0;k<nss;k++) {ga += fabs(GA[k]); GA2[k] = GA[k];}
  if (ga < 1.e-20 || ga > 100) return (0);

  N4 = PLC_diff_g_ga (ii,ip,nss,mat,dt,crpl,GA,eps);
  if (fabs(N4) < 1.e-20 || N4 >= 0.0) return (0);

  for (k=0;k<nss;k++){

    if (fabs(GA[k]) < 1.e-20) continue;

    // @todo The following two branches leave GAGA uninitialized when GA[k] ==
    //       GA1[k]. I assume that this condition is not possible, and have
    //       added this assert to verify. @cp should review and make sure that
    //       this is correct. In order to suppress uninitialized warnings I have
    //       added the constructor call at its definition.
    assert(GA[k] != GA1[k]);

    if (GA[k] < GA1[k]) GAGA = GA2[k] = 1.e-10;
    if (GA[k] > GA1[k]) GAGA = GA2[k] = 1.0;

    nor = 10.0;
    while (nor > nor_min){

      ga = 0.0; for (j=0;j<nss;j++) ga += fabs(GA2[j]);

      N4 = PLC_diff_g_ga (ii,ip,nss,mat,dt,crpl,GA2,eps);
      Gama = eps[ii].il[ip].GAMA + dt*ga;

      Rof = (crpl[mat].l1/crpl[mat].l2 + crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama))*(crpl[mat].l1/crpl[mat].l2 + crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama));
      Rom = 2.*crpl[mat].c1*crpl[mat].c2/(crpl[mat].b*crpl[mat].l2)*(exp(-1./2.*crpl[mat].l2*Gama) - crpl[mat].c3);
      Xi = pow ((crpl[mat].b*Rom/(sqrt(Rof)*crpl[mat].to*ga)),1./3.);

      /* Compute residual */
      R = crpl[mat].To + 1./3.*crpl[mat].nu*crpl[mat].b*sqrt(Rof) + crpl[mat].fo*(1. - exp(-1.*Xi)) - HAR;

      dga = -R/(N4*GA2[k]/fabs(GA2[k])); GAGA += dga; GA2[k] = GAGA;

      nor = sqrt (R*R)/HAR;
      sprintf (str,"%f",nor);
      for (N=0;N<3;N++){ M = 10; M = strcmp(err[N],str);
      if (M == 0) {PGFEM_printf("%ld %ld || Error in PLC_Inaccessible : nor = %s\n",ii,ip,err[N]); return (1);}}
    }/* end while */
  }/* end k < nss */

  N4 = PLC_diff_g_ga (ii,ip,nss,mat,dt,crpl,GA1,eps); if (N4 < 0.0) PGFEM_printf ("Pozor pozor N4 < 0.0\n");

  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      if (M == N) nor = 1.0; else nor = 0.0;
      UU[M][N] = nor;
    }
  }

  for (k=0;k<nss;k++){
    GA[k] = GA2[k];

    PP[0][0] = crpl[mat].P[k][3]*crpl[mat].P[k][0];  PP[0][1] = crpl[mat].P[k][3]*crpl[mat].P[k][1];  PP[0][2] = crpl[mat].P[k][3]*crpl[mat].P[k][2];
    PP[1][0] = crpl[mat].P[k][4]*crpl[mat].P[k][0];  PP[1][1] = crpl[mat].P[k][4]*crpl[mat].P[k][1];  PP[1][2] = crpl[mat].P[k][4]*crpl[mat].P[k][2];
    PP[2][0] = crpl[mat].P[k][5]*crpl[mat].P[k][0];  PP[2][1] = crpl[mat].P[k][5]*crpl[mat].P[k][1];  PP[2][2] = crpl[mat].P[k][5]*crpl[mat].P[k][2];

    for (M=0;M<3;M++) for (N=0;N<3;N++) UU[M][N] -= dt*GA2[k]*PP[M][N];
  }

  dealoc1 (GA2);

  return (INFO);
}
