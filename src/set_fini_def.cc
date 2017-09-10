#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "set_fini_def.h"
#include "elem3d.h"
#include "enumerations.h"
#include <cstring>

static constexpr int periodic = 0;

void set_fini_def (long ne,
                   long npres,
                   Element *elem,
                   EPS *eps,
                   SIG *sig,
                   const int analysis)
{
  long ii,ip,i,nne,II,P,R;
  double dij;

  if (periodic == 1){
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
        if (P == R)
          dij = 1.0;
        else dij = 0.0;
        eps[0].P[P][R] = eps[0].S[P][R] = 0.0;
        eps[0].F[P][R] = eps[0].Fn[P][R] =
       eps[0].FB[P][R] = eps[0].Fe[P][R] = dij;
        if (analysis == FS_CRPL){
          eps[0].Fp[P][R] = dij;
        }
      }
    }
  }/* end PERIODIC */

  for (ii=0;ii<ne;ii++){

    nne = elem[ii].toe;

    /* Integration */
    int_point (nne,&II);

    for (ip=0;ip<II;ip++){
      eps[ii].il[ip].F[0] = 1.0;
      eps[ii].il[ip].F[1] = 0.0;
      eps[ii].il[ip].F[2] = 0.0;

      eps[ii].il[ip].F[3] = 0.0;
      eps[ii].il[ip].F[4] = 1.0;
      eps[ii].il[ip].F[5] = 0.0;

      eps[ii].il[ip].F[6] = 0.0;
      eps[ii].il[ip].F[7] = 0.0;
      eps[ii].il[ip].F[8] = 1.0;

      eps[ii].il[ip].Un_1 = eps[ii].il[ip].Un = 0.0;
      eps[ii].il[ip].Jn_1 = 1.0;

      eps[ii].il[ip].eff = eps[0].eff = 0.0;

      if (periodic == 1) {
        switch(analysis){
         case FS_CRPL:
         case FINITE_STRAIN:
         case STABILIZED:
          eps[ii].il[ip].Fe[0] = eps[ii].il[ip].Fe1[0] = 1.0;
          eps[ii].il[ip].Fe[1] = eps[ii].il[ip].Fe1[1] = 0.0;
          eps[ii].il[ip].Fe[2] = eps[ii].il[ip].Fe1[2] = 0.0;
          eps[ii].il[ip].Fe[3] = eps[ii].il[ip].Fe1[3] = 0.0;
          eps[ii].il[ip].Fe[4] = eps[ii].il[ip].Fe1[4] = 1.0;
          eps[ii].il[ip].Fe[5] = eps[ii].il[ip].Fe1[5] = 0.0;
          eps[ii].il[ip].Fe[6] = eps[ii].il[ip].Fe1[6] = 0.0;
          eps[ii].il[ip].Fe[7] = eps[ii].il[ip].Fe1[7] = 0.0;
          eps[ii].il[ip].Fe[8] = eps[ii].il[ip].Fe1[8] = 1.0;
          break;
         default: break;
        }
        eps[ii].il[ip].F[0] = eps[ii].il[ip].F[4] = eps[ii].il[ip].F[8] = 0.0;
      }
    }/*end ip < II */

    /* PRESSURE PART */
    if (analysis == STABILIZED || analysis == MINI /* || analysis == MINI_3F */){
      /* Integration */
      if(analysis == STABILIZED)
        int_point (10,&II); /*if (gem == 1) II *= elem[ii].nse;*/
      if(analysis == MINI || analysis == MINI_3F)
        int_point (5,&II);

      for (ip=0;ip<II;ip++){
        eps[ii].st[ip].Fpp[0] = 1.0;
        eps[ii].st[ip].Fpp[1] = 0.0;
        eps[ii].st[ip].Fpp[2] = 0.0;

        eps[ii].st[ip].Fpp[3] = 0.0;
        eps[ii].st[ip].Fpp[4] = 1.0;
        eps[ii].st[ip].Fpp[5] = 0.0;

        eps[ii].st[ip].Fpp[6] = 0.0;
        eps[ii].st[ip].Fpp[7] = 0.0;
        eps[ii].st[ip].Fpp[8] = 1.0;

        eps[ii].st[ip].Un_1 = eps[ii].st[ip].Un = 0.0;
        eps[ii].st[ip].Jn_1 = 1.0;

        if (periodic == 1/* || gem == 1*/) {
          eps[ii].st[ip].Fpp[0] = eps[ii].st[ip].Fpp[4]
                                = eps[ii].st[ip].Fpp[8] = 0.0;
        }
      }
    }/* end analysis  == STAB || MINI */

    for (i=0;i<npres;i++){
      if (analysis == FS_CRPL || analysis == FINITE_STRAIN) {
        eps[ii].T[i] = 1.0;
        eps[ii].d_T[i] = 1.0;
      }
      sig[ii].p[i] = 0.0;
      sig[ii].d_p[i] = 0.0;
      sig[ii].pn_1[i] = 0.0;
    }

    if(analysis == MINI_3F){
      for (i=0;i<4;i++){
        eps[ii].T[i] = 1.0;
        eps[ii].d_T[i] = 1.0;
      }
    }


    if (elem[ii].n_bub*elem[ii].n_bub_dofs > 0){
      memset(elem[ii].bub_dofs,0,
             elem[ii].n_bub*elem[ii].n_bub_dofs*sizeof(double));
      memset(elem[ii].d_bub_dofs,0,
             elem[ii].n_bub*elem[ii].n_bub_dofs*sizeof(double));
    }

  }/* ii < ne */
}

void set_fini_def_pl (long ne,
                      long npres,
                      Element *elem,
                      EPS *eps,
                      SIG *sig,
                      CRPL *crpl,
                      const int analysis,
                      const int plc)
{
  long ii,ip,i,nne,II,mat;

  for (ii=0;ii<ne;ii++){

    nne = elem[ii].toe; mat = elem[ii].mat[2];

    /* Integration */
    int_point (nne,&II);

    for (ip=0;ip<II;ip++){
      eps[ii].il[ip].UU[0] = 1.0; eps[ii].il[ip].UU[1] = 0.0; eps[ii].il[ip].UU[2] = 0.0;
      eps[ii].il[ip].UU[3] = 0.0; eps[ii].il[ip].UU[4] = 1.0; eps[ii].il[ip].UU[5] = 0.0;
      eps[ii].il[ip].UU[6] = 0.0; eps[ii].il[ip].UU[7] = 0.0; eps[ii].il[ip].UU[8] = 1.0;

      eps[ii].il[ip].Fp[0] = 1.0; eps[ii].il[ip].Fp[1] = 0.0; eps[ii].il[ip].Fp[2] = 0.0;
      eps[ii].il[ip].Fp[3] = 0.0; eps[ii].il[ip].Fp[4] = 1.0; eps[ii].il[ip].Fp[5] = 0.0;
      eps[ii].il[ip].Fp[6] = 0.0; eps[ii].il[ip].Fp[7] = 0.0; eps[ii].il[ip].Fp[8] = 1.0;

      /* Stress on slip systems */
      for (i=0;i<crpl[mat].nss;i++) {sig[ii].il[ip].Tau[i] = eps[ii].il[ip].GA[i] = eps[ii].il[ip].GA1[i] = 0.0;}

      /* Hardening */
      sig[ii].il[ip].Har = sig[ii].il[ip].Har1 = crpl[mat].To;

      /* Lagrange Multiplier */
      eps[ii].il[ip].lam = 0.0;

      /* Acumulative plastic slip : PLC */
      eps[ii].il[ip].GAMA = 0.0; if (plc == 1) for (i=0;i<crpl[mat].nss;i++) eps[ii].il[ip].PLC_B[i] = 0;
    }
  }
}
