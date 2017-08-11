/* HEADER */
#include "press_theta.h"
#include <math.h>
#include "enumerations.h"
#include "get_dof_ids_on_elem.h"
#include "cast_macros.h"
#include "def_grad.h"
#include "elem3d.h"
#include "allocation.h"
#include "integration.h"
#include "pressu_shape.h"
#include "Re1_Re2_Re3.h"
#include "stress_strain.h"
#include "tensors.h"
#include "utils.h"

static const int periodic = 0;

void press_theta (long ne,
          long ndofn,
          long npres,
          ELEMENT *elem,
          NODE *node,
          double *d_r,
          double *rr,
          SUPP sup,
          MATGEOM matgeom,
          HOMMAT *hommat,
          EPS *eps,
          SIG *sig,
          long iter,
          double nor_min,
          double dt,
          CRPL *crpl,
          const PGFem3D_opt *opts,
          const int mp_id)
{

  long ii,i,j,k,II,JJ,KK,ip,M,N,X,P,Q,R,U,W,*nod,nne,ndn,mat,*cn;
  double *Psi,*gk,*ge,*gz,*w,*N_x,*N_y,*N_z,*x,*y,*z,*r_e,ksi;
  double eta,zet,ai,aj,ak,J,**Fn,**Fr,**Fr_I,Jr,Jn,Tr,Tn;
  double L[3][3][3][3],****FF,**f,***AA,***aa,***BB;
  double ***bb,**DD,**dd,**MM,**mm,****ST,**u,*r_r,pp,**Re1;
  double **re1,*Re2,*re2,*Re3,*re3,**E,*D_T,*D_P,*sup_def,**UU;
  double ****pFFp,**pfp,***BB1,***bb1,***BB2,***bb2;
  double **MMp,**mmp,**S,**FnB,****pGp,dij;

  /* 3D */ ndn = 3;

  Psi = aloc1 (npres);
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);
  x = aloc1 (10);
  y = aloc1 (10);
  z = aloc1 (10);
  Fn = aloc2 (3,3);
  Fr = aloc2 (3,3);
  f = aloc2 (3,3);
  DD = aloc2 (npres,npres);
  S = aloc2 (3,3);
  MM = aloc2 (npres,npres);
  nod = aloc1l (10);
  Fr_I = aloc2 (3,3);
  E = aloc2 (3,3);
  Re2 = aloc1(npres);
  Re3 = aloc1(npres);
  D_T = aloc1 (npres);
  D_P = aloc1 (npres);
  if(sup->npd > 0){
    sup_def = aloc1 (sup->npd);
  } else {
    sup_def = NULL;
  }
  FnB = aloc2 (3,3);


  if (opts->analysis_type == FS_CRPL){
    UU = aloc2 (3,3);
    pfp = aloc2 (3,3);
  }

  for (ii=0;ii<ne;ii++){

    nne = elem[ii].toe;
    mat = elem[ii].mat[2];
    elemnodes (ii,nne,nod,elem);
    cn = aloc1l (nne*ndofn);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    if (opts->analysis_type == FS_CRPL){
      pFFp = aloc4 (3,3,ndn,nne);
      MMp = aloc2 (npres,npres);
      mmp = aloc2 (npres,npres);
      BB1 = aloc3 (ndn,npres,nne);
      bb1 = aloc3 (ndn,npres,nne);
      BB2 = aloc3 (ndn,npres,nne);
      bb2 = aloc3 (ndn,npres,nne);
      pGp = aloc4 (3,3,ndn,nne);
    }

    r_e = aloc1 (nne*ndofn);
    r_r = aloc1 (nne*ndofn);
    N_x = aloc1 (nne);
    N_y = aloc1 (nne);
    N_z = aloc1 (nne);
    FF = aloc4 (3,3,ndn,nne);
    AA = aloc3 (ndn,npres,nne);
    aa = aloc3 (ndn,npres,nne);
    u = aloc2 (ndn,nne);
    BB = aloc3 (ndn,npres,nne);
    bb = aloc3 (ndn,npres,nne);
    ST = aloc4 (3,3,ndn,nne);
    Re1 = aloc2 (ndn,nne);
    re1 = aloc2 (ndn,nne);
    dd = aloc2 (npres,npres);
    mm = aloc2 (npres,npres);
    re2 = aloc1 (npres);
    re3 = aloc1 (npres);

    nodecoord_updated (nne,nod,node,x,y,z);

    /* Integration */
    integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);

    /* deformation on element */
    if (iter != 0){
      for (j=0;j<sup->npd;j++){
    sup_def[j] = sup->defl_d[j];
    sup->defl_d[j] = 0.0;
      }
    }
    def_elem (cn,nne*ndofn,rr,elem,node,r_r,sup,0);
    if (iter != 0){
      for (i=0;i<sup->npd;i++)
    sup->defl_d[i] = sup_def[i];
    }

    if (iter == 0){
      for (j=0;j<sup->npd;j++){
    sup_def[j] = sup->defl_d[j];
    sup->defl_d[j] = 0.0;
      }
    }
    def_elem (cn,nne*ndofn,d_r,elem,node,r_e,sup,0);
    if (iter == 0){
      for (i=0;i<sup->npd;i++)
    sup->defl_d[i] = sup_def[i];
    }

    ip = 0;
    for(i=0;i<II;i++){
      for(j=0;j<JJ;j++){
    for (k=0;k<KK;k++){

      if (nne == 4)  {
        ksi = *(gk+k);
        eta = *(ge+k);
        zet = *(gz+k);
        ai = *(w+k);
        aj = 1.0;
        ak = 1.0;
      }
      if (nne == 10) {
        ksi = *(gk+k);
        eta = *(ge+k);
        zet = *(gz+k);
        ai = *(w+k);
        aj = 1.0;
        ak = 1.0;
      }
      if (nne == 8)  {
        ksi = *(gk+i);
        eta = *(gk+j);
        zet = *(gk+k);
        ai = *(w+i);
        aj = *(w+j);
        ak = *(w+k);
      }

      /* Derivatives of shape functions */
      J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);

      Fn[0][0] = eps[ii].il[ip].F[0];
      Fn[0][1] = eps[ii].il[ip].F[1];
      Fn[0][2] = eps[ii].il[ip].F[2];

      Fn[1][0] = eps[ii].il[ip].F[3];
      Fn[1][1] = eps[ii].il[ip].F[4];
      Fn[1][2] = eps[ii].il[ip].F[5];

      Fn[2][0] = eps[ii].il[ip].F[6];
      Fn[2][1] = eps[ii].il[ip].F[7];
      Fn[2][2] = eps[ii].il[ip].F[8];

      shape_tensor (nne,ndofn,N_x,N_y,N_z,ST);
      def_grad_get (nne,ndofn,CONST_4(double) ST,r_e,Fr);
      def_grad_inv (CCONST_2(double) Fr,Fr_I);
      Jr = def_grad_det (CCONST_2(double) Fr);
      Jn = def_grad_det (CCONST_2(double) Fn);

      /* Pressure shape functions */
      pressu_shape (npres,ksi,eta,zet,Psi);

      Tn = 0.0; Tr = 0.0; pp = 0.0;
      for (M=0;M<npres;M++){
        Tn += Psi[M]*eps[ii].T[M];
        Tr += Psi[M]*eps[ii].d_T[M];
        pp += Psi[M]*(sig[ii].p[M] + sig[ii].d_p[M]);
      }

      /* Set Fn-BAR */
      for (Q=0;Q<3;Q++){
        for (R=0;R<3;R++){
          FnB[Q][R] = pow(Tn,1./3.)*pow(Jn,-1./3.)*Fn[Q][R];
        }
      }

      /****** UNIT CELL APPROACH *********/
      if (periodic == 1) {
        if (opts->analysis_type == FS_CRPL){
          S[0][0] = eps[ii].il[ip].Fp[0];
          S[0][1] = eps[ii].il[ip].Fp[1];
          S[0][2] = eps[ii].il[ip].Fp[2];

          S[1][0] = eps[ii].il[ip].Fp[3];
          S[1][1] = eps[ii].il[ip].Fp[4];
          S[1][2] = eps[ii].il[ip].Fp[5];

          S[2][0] = eps[ii].il[ip].Fp[6];
          S[2][1] = eps[ii].il[ip].Fp[7];
          S[2][2] = eps[ii].il[ip].Fp[8];

          def_grad_inv (CCONST_2(double) S,FnB);
        }/* opts->analysis_type == FS_CRPL */
        else{
          for (N=0;N<3;N++){
        for (P=0;P<3;P++){
          if (N == P) dij = 1.0;
          else dij = 0.0;
          FnB[N][P] = dij;
        }
          }
        }

        for (N=0;N<3;N++){
          for (P=0;P<3;P++){
        Fr[N][P] += eps[0].F[N][P] + Fn[N][P];
          }
        }

        Jr = def_grad_det (CCONST_2(double) Fr);
        def_grad_inv (CCONST_2(double) Fr,Fr_I);
        Jn = Tn = 1.;

      }/* end PERIODIC */

      /* Material stiffness matrix */
      matrix_tensor_3D (elem[ii].mat[2],hommat,L);

      /***** CRYSTAL PLASTICITY ********/
      if (opts->analysis_type == FS_CRPL){
        UU[0][0] = eps[ii].il[ip].UU[0];
        UU[0][1] = eps[ii].il[ip].UU[1];
        UU[0][2] = eps[ii].il[ip].UU[2];

        UU[1][0] = eps[ii].il[ip].UU[3];
        UU[1][1] = eps[ii].il[ip].UU[4];
        UU[1][2] = eps[ii].il[ip].UU[5];

        UU[2][0] = eps[ii].il[ip].UU[6];
        UU[2][1] = eps[ii].il[ip].UU[7];
        UU[2][2] = eps[ii].il[ip].UU[8];

        for (M=0;M<3;M++){
          for (N=0;N<3;N++){
        S[M][N] = 0.0;
        for (P=0;P<3;P++){
          S[M][N] += FnB[M][P]*UU[P][N];
        }
          }
        }

        /* Strain */
        get_GL_strain (S,Fr,Jr,Tr,E);
      }/* opts->analysis_type == FS_CRPL */
      else {
        /* Strain */
        get_GL_strain (FnB,Fr,Jr,Tr,E);
      }

      /* Stress */
      get_SPK_stress (L,E,S);

      /* Tensors used in opts->analysis_type */
      tensors_FF_f (nne,ndn,FnB,Fr,Fr_I,Jr,ST,FF,f);

      /************ CRYSTAL PLASTICITY ***********************/
      if (opts->analysis_type == FS_CRPL){

        for (P=0;P<3;P++){
          for (R=0;R<3;R++){
        /* eFn+1 and F* */
        pfp[P][R] = 0.0;

        for (U=0;U<3;U++){
          for (W=0;W<3;W++){
            pfp[P][R] += UU[U][P] * f[U][W] * UU[W][R];
          }
        }

        for (M=0;M<ndn;M++){
          for (N=0;N<nne;N++){
            pFFp[P][R][M][N] = 0.0;

            for (U=0;U<3;U++){
              for (W=0;W<3;W++){
            pFFp[P][R][M][N] += (UU[U][P] * FF[U][W][M][N]
                         * UU[W][R]);
              }
            }
          }
        }
          }/* end R < 3 */
        }/* end P < 3 */

        /* Elasto-Plastic tensors */
        tensors_aa_bb_dd_mm_plast (ii,ip,nne,ndn,npres,ai,aj,ak,J,
                       Psi,L,Fn,Fr,Fr_I,Jn,Jr,Tn,Tr,ST,
                       FF,f,AA,aa,BB,bb,DD,dd,MM,mm,S,
                       pFFp,pfp,UU,eps,BB1,bb1,
                       BB2,bb2,MMp,mmp,pGp);
      }/* opts->analysis_type == FS_CRPL */
      else
        tensors_aa_bb_dd_mm (nne,ndn,npres,ai,aj,ak,J,Psi,L,FnB,
                 Fr,Fr_I,Jn,Jr,Tn,Tr,ST,FF,f,AA,aa,BB,
                 bb,DD,dd,MM,mm);

      /* Residuals */
      if (opts->analysis_type == FS_CRPL){
        Re1_Re2_Re3 (nne,ndn,npres,ai,aj,ak,J,Psi,
             Jr,Jn,Tr,Tn,Fr_I,ST,pFFp,S,pfp,
             pp,Re1,re1,Re2,re2,Re3,re3);
      } else {
        Re1_Re2_Re3 (nne,ndn,npres,ai,aj,ak,J,Psi,
             Jr,Jn,Tr,Tn,Fr_I,ST,FF,S,f,pp,
             Re1,re1,Re2,re2,Re3,re3);
      }

      ip++;
    }/* end k < KK */
      }/* end j < JJ */
    }/* end i < II */

    /* for (Q=0;Q<nne;Q++){u[0][Q] = r_e[ndofn*Q+0] + r_r[ndofn*Q+0];
       u[1][Q] = r_e[ndofn*Q+1] + r_r[ndofn*Q+1]; u[2][Q] =
       r_e[ndofn*Q+2] + r_r[ndofn*Q+2];} */

    for (Q=0;Q<nne;Q++){
      u[0][Q] = r_r[ndofn*Q+0];
      u[1][Q] = r_r[ndofn*Q+1];
      u[2][Q] = r_r[ndofn*Q+2];
    }

    /***** CRYSTAL PLASTICITY ******/
    if (opts->analysis_type == FS_CRPL){
      /* Tensors BB and bb */
      for (P=0;P<ndn;P++){
    for (M=0;M<npres;M++){
      for (N=0;N<nne;N++){
        bb[P][M][N] += bb2[P][M][N];
      }
    }
      }
      for (P=0;P<npres;P++){
    for (M=0;M<npres;M++){
      mm[P][M] += mmp[P][M];
    }
      }
    }/* end opts->analysis_type == FS_CRPL */

    for (N=0;N<npres;N++){
      D_T[N] = 0.0; D_P[N] = 0.0;
      for (X=0;X<npres;X++){
    D_T[N] += 1./dd[N][X]*re2[X];
    D_P[N] += 1./dd[N][X]*re3[X];

    for (R=0;R<ndn;R++){
      for (P=0;P<nne;P++){
        D_T[N] += 1./dd[N][X]*aa[R][X][P]*u[R][P];
        D_P[N] += 1./dd[N][X]*bb[R][X][P]*u[R][P];

        for (U=0;U<npres;U++){
          for (W=0;W<npres;W++){
        D_P[N] += 1./dd[N][X]*mm[X][U]*1./dd[U][W]*aa[R][W][P]*u[R][P];
          }
        }
      }
    }/* end R < ndn */

    for (R=0;R<npres;R++){
      for (U=0;U<npres;U++){
        D_P[N] += 1./dd[N][X]*mm[X][R]*1./dd[R][U]*re2[U];
      }
    }

      }
    }

    for (N=0;N<npres;N++){
      eps[ii].d_T[N] += D_T[N];
      sig[ii].d_p[N] += D_P[N];
    }

    /* PGFEM_printf("dT = %12.20f || dP = %12.20f :: ddT = %12.20f || ddP =
       %12.20f\n",D_T[0],D_P[0],eps[ii].d_T[0],sig[ii].d_p[0]);  */

    dealoc1 (r_e);
    dealoc1 (r_r);
    dealoc1 (N_x);
    dealoc1 (N_y);
    dealoc1 (N_z);
    dealoc4 (FF,3,3,ndn);
    dealoc3 (AA,ndn,npres);
    dealoc3 (aa,ndn,npres);
    dealoc2 (u,ndn);
    dealoc3 (BB,ndn,npres);
    dealoc3 (bb,ndn,npres);
    dealoc4 (ST,3,3,ndn);
    dealoc2 (Re1,ndn);
    dealoc2 (re1,ndn);
    dealoc2 (dd,npres);
    dealoc2 (mm,npres);
    dealoc1 (re2);
    dealoc1 (re3);
    dealoc1l (cn);


    if (opts->analysis_type == FS_CRPL){
      dealoc4 (pFFp,3,3,ndn);
      dealoc2 (MMp,npres);
      dealoc2 (mmp,npres);
      dealoc3 (BB1,ndn,npres);
      dealoc3 (bb1,ndn,npres);
      dealoc3 (BB2,ndn,npres);
      dealoc3 (bb2,ndn,npres);
      dealoc4 (pGp,3,3,ndn);
    }

  }/* end ii < ne */

  dealoc1 (Psi);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (x);
  dealoc1 (y);
  dealoc1 (z);
  dealoc2 (Fn,3);
  dealoc2 (Fr,3);
  dealoc2 (f,3);
  dealoc2 (DD,npres);
  dealoc2 (S,3);
  dealoc2 (MM,npres);
  dealoc1l (nod);
  dealoc2 (Fr_I,3);
  dealoc2 (E,3);
  dealoc1 (Re2);
  dealoc1 (Re3);
  dealoc1 (D_T);
  dealoc1 (D_P);
  dealoc1 (sup_def);
  dealoc2 (FnB,3);

  if (opts->analysis_type == FS_CRPL) {
    dealoc2 (UU,3);
    dealoc2 (pfp,3);
  }
}
