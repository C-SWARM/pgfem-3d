#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "computeMacroF.h"
#include "elem3d.h"
#include "incl.h"
#include "utils.h"

double* computeMacroF(Element *elem,
                      long ne,
                      NODE *node,
                      long nn,
                      EPS *eps,
                      double oVolume,
                      MPI_Comm mpi_comm)
{
  long *nod;
  double *F, *gk, *ge, *gz, *w;
  double  *x, *y, *z, *N_x, *N_y, *N_z, *Np;


  /* alocate everything that will not change size */
  F = aloc1(9);
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);

  F[0] = F[1] = F[2] = 0.0;
  F[3] = F[4] = F[5] = 0.0;
  F[6] = F[7] = F[8] = 0.0;


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

    /* 4 point integration */
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

          F[0] += ai*aj*ak*J* eps[i].il[ip].F[0];
          F[1] += ai*aj*ak*J* eps[i].il[ip].F[1];
          F[2] += ai*aj*ak*J* eps[i].il[ip].F[2];

          F[3] += ai*aj*ak*J* eps[i].il[ip].F[3];
          F[4] += ai*aj*ak*J* eps[i].il[ip].F[4];
          F[5] += ai*aj*ak*J* eps[i].il[ip].F[5];

          F[6] += ai*aj*ak*J* eps[i].il[ip].F[6];
          F[7] += ai*aj*ak*J* eps[i].il[ip].F[7];
          F[8] += ai*aj*ak*J* eps[i].il[ip].F[8];

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

  /*Now the deformation gradient is integrated over the whole
    domain. Divide by the reference volume and communicate.*/
  F[0] /= oVolume;
  F[1] /= oVolume;
  F[2] /= oVolume;

  F[3] /= oVolume;
  F[4] /= oVolume;
  F[5] /= oVolume;

  F[6] /= oVolume;
  F[7] /= oVolume;
  F[8] /= oVolume;

  double *GF;
  GF = aloc1(9);
  GF[0] = GF[1] = GF[2] = 0.0;
  GF[3] = GF[4] = GF[5] = 0.0;
  GF[6] = GF[7] = GF[8] = 0.0;

  MPI_Allreduce(F,GF,9,MPI_DOUBLE,MPI_SUM,mpi_comm);

  free(F);
  free(gk);
  free(ge);
  free(gz);
  free(w);

  return GF;

}
