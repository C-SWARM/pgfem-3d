#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "res_fini_def.h"

#include "elem3d.h"
#include "enumerations.h"

static constexpr int periodic = 0;

void res_fini_def (long ne,
                   long npres,
                   Element *elem,
                   EPS *eps,
                   SIG *sig,
                   CRPL *crpl,
                   const int analysis)
{
  long ii, ip, i, nne, II, P, R, U, W, mat;

  for (ii=0;ii<ne;ii++){

    nne = elem[ii].toe;
    mat = elem[ii].mat[2];

    /* Integration */
    int_point (nne,&II);

    for (ip=0;ip<II;ip++){
      if (analysis == FS_CRPL){
        sig[ii].il[ip].Har1 = sig[ii].il[ip].Har;

        for (P=0;P<crpl[mat].nss;P++)
          eps[ii].il[ip].GA1[P] = eps[ii].il[ip].GA[P];

        eps[ii].il[ip].UU[0] = 1.0;
        eps[ii].il[ip].UU[1] = 0.0;
        eps[ii].il[ip].UU[2] = 0.0;

        eps[ii].il[ip].UU[3] = 0.0;
        eps[ii].il[ip].UU[4] = 1.0;
        eps[ii].il[ip].UU[5] = 0.0;

        eps[ii].il[ip].UU[6] = 0.0;
        eps[ii].il[ip].UU[7] = 0.0;
        eps[ii].il[ip].UU[8] = 1.0;

        for (P=0;P<3;P++){
          for (R=0;R<3;R++){
            eps[ii].il[ip].dUU_Tr[P][R] = eps[ii].il[ip].dUU_Tr_n[P][R];
            for (U=0;U<3;U++){
              for (W=0;W<3;W++){
                eps[ii].il[ip].dUU_Fr[P][R][U][W] =
                eps[ii].il[ip].dUU_Fr_n[P][R][U][W];
              }
            }/* end U */
          }
        }/* end  P < 3 */
      }/* analysis == FS_CRPL */

      if (periodic == 1){
        eps[ii].il[ip].Fe1[0] = eps[ii].il[ip].Fe[0];
        eps[ii].il[ip].Fe1[1] = eps[ii].il[ip].Fe[1];
        eps[ii].il[ip].Fe1[2] = eps[ii].il[ip].Fe[2];

        eps[ii].il[ip].Fe1[3] = eps[ii].il[ip].Fe[3];
        eps[ii].il[ip].Fe1[4] = eps[ii].il[ip].Fe[4];
        eps[ii].il[ip].Fe1[5] = eps[ii].il[ip].Fe[5];

        eps[ii].il[ip].Fe1[6] = eps[ii].il[ip].Fe[6];
        eps[ii].il[ip].Fe1[7] = eps[ii].il[ip].Fe[7];
        eps[ii].il[ip].Fe1[8] = eps[ii].il[ip].Fe[8];
      }
    }/* ip < II */
    for (i=0;i<npres;i++){
      sig[ii].d_p[i] = 0.0;
      eps[ii].d_T[i] = 1.0;
      if (periodic == 1)
        eps[ii].d_T[i] = eps[ii].T[i];
    }
  }/* end ii < ne */
}
