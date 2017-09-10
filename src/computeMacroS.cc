#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "computeMacroS.h"
#include "allocation.h"
#include "elem3d.h"
#include "index_macros.h"
#include "utils.h"
#include <mkl_cblas.h>
#include <cstring>

double* computeMacroS(Element *elem,
                      long ne,
                      NODE *node,
                      long nn,
                      SIG *sig,
                      double oVolume,
                      MPI_Comm mpi_comm)
{
  long *nod;
  double  *S, *gk, *ge, *gz, *w;
  double  *x, *y, *z, *N_x, *N_y, *N_z, *Np;

  /* alocate everything that will not change size */
  S = aloc1(6);
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);

  S[0] = S[1] = S[2] = S[3] = S[4] = S[5] = 0.0;

  for(int i=0; i<ne; i++){
    long nne = elem[i].toe;

    if(nne != 4){
      continue;
    }

    /* allocate */
    nod = aloc1l(nne);

    x = aloc1(nne);
    y = aloc1(nne);
    z = aloc1(nne);
    N_x = aloc1 (nne);
    N_y = aloc1 (nne);
    N_z = aloc1 (nne);
    Np  = aloc1 (nne);


    elemnodes(i,nne,nod,elem);
    nodecoord_total(nne,nod,node,x,y,z);

    /* integration */
    long JJ, KK, MM;
    integrate(nne, &JJ, &KK, &MM, gk, ge, gz, w);

    int ip = 0;
    for(int j=0; j<JJ; j++){
      for(int k=0; k<KK; k++){
        for(int m=0; m<MM; m++){
          double ksi, eta, zet, ai, aj, ak, J;

          ksi = gk[m];
          eta = ge[m];
          zet = gz[m];
          ai = w[m];
          aj = ak = 1.0;

          shape_func(ksi,eta,zet,nne,Np);
          J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);

          S[0] += ai*aj*ak*J* sig[i].il[ip].o[0];
          S[1] += ai*aj*ak*J* sig[i].il[ip].o[1];
          S[2] += ai*aj*ak*J* sig[i].il[ip].o[2];
          S[3] += ai*aj*ak*J* sig[i].il[ip].o[3];
          S[4] += ai*aj*ak*J* sig[i].il[ip].o[4];
          S[5] += ai*aj*ak*J* sig[i].il[ip].o[5];

          ip++;
        }
      }
    }

    /* free */
    free(nod);
    free(x);
    free(y);
    free(z);
    free(N_x);
    free(N_y);
    free(N_z);
    free(Np);
  }/* for each element (i) */

  /* Now the SPK is integrated over the domain.  Divide by the
     reference volume and communicate. */

  S[0] /= oVolume;
  S[1] /= oVolume;
  S[2] /= oVolume;
  S[3] /= oVolume;
  S[4] /= oVolume;
  S[5] /= oVolume;

  double *GS;
  GS = aloc1(6);

  GS[0] = GS[1] = GS[2] = GS[3] = GS[4] = GS[5] = 0.0;

  MPI_Allreduce(S,GS,6,MPI_DOUBLE,MPI_SUM,mpi_comm);

  free(S);
  free(gk);
  free(ge);
  free(gz);
  free(w);

  return(GS);
}

double* computeMacroP(Element *elem,
                      long ne,
                      NODE *node,
                      long nn,
                      SIG *sig,
                      EPS *eps,
                      double oVolume,
                      MPI_Comm mpi_comm)
{
  long *nod;
  double  *S, *P, *gk, *ge, *gz, *w;
  double  *x, *y, *z, *N_x, *N_y, *N_z, *Np;

  /* alocate everything that will not change size */
  S = aloc1(9);
  P = aloc1(9);
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);

  /* null P */
  memset(P,0,9*sizeof(double));

  for(int i=0; i<ne; i++){
    long nne = elem[i].toe;

    if(nne != 4){
      continue;
    }

    /* allocate */
    nod = aloc1l(nne);

    x = aloc1(nne);
    y = aloc1(nne);
    z = aloc1(nne);
    N_x = aloc1 (nne);
    N_y = aloc1 (nne);
    N_z = aloc1 (nne);
    Np  = aloc1 (nne);


    elemnodes(i,nne,nod,elem);
    nodecoord_total(nne,nod,node,x,y,z);

    /* integration */
    long JJ, KK, MM;
    integrate(nne, &JJ, &KK, &MM, gk, ge, gz, w);

    int ip = 0;
    for(int j=0; j<JJ; j++){
      for(int k=0; k<KK; k++){
        for(int m=0; m<MM; m++){
          double ksi, eta, zet, ai, aj, ak, J;

          ksi = gk[m];
          eta = ge[m];
          zet = gz[m];
          ai = w[m];
          aj = ak = 1.0;

          shape_func(ksi,eta,zet,nne,Np);
          J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);

          S[idx_2(0,0)] = sig[i].il[ip].o[0];
          S[idx_2(1,1)] = sig[i].il[ip].o[1];
          S[idx_2(2,2)] = sig[i].il[ip].o[2];
          S[idx_2(1,2)] = S[idx_2(2,1)] = sig[i].il[ip].o[3];
          S[idx_2(0,2)] = S[idx_2(2,0)] = sig[i].il[ip].o[4];
          S[idx_2(0,1)] = S[idx_2(1,0)] = sig[i].il[ip].o[5];

          cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                      3,3,3,(ai*aj*ak*J),eps[i].il[ip].F,3,
                      S,3,1.0,P,3);
          ip++;
        }
      }
    }

    /* free */
    free(nod);
    free(x);
    free(y);
    free(z);
    free(N_x);
    free(N_y);
    free(N_z);
    free(Np);
  }/* for each element (i) */

  /* Now the 1PK is integrated over the domain.  Divide by the
     reference volume and communicate. */
  cblas_dscal(9,1./oVolume,P,1);

  double *GS;
  GS = aloc1(9);

  MPI_Allreduce(P,GS,9,MPI_DOUBLE,MPI_SUM,mpi_comm);

  free(S);
  free(P);
  free(gk);
  free(ge);
  free(gz);
  free(w);

  return(GS);
}

