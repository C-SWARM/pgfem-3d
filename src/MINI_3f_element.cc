/**
 * AUTHORS:
 * Matthew Mosby
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "MINI_3f_element.h"
#include "Hu_Washizu_element.h"
#include "PGFEM_io.h"
#include "allocation.h"
#include "cast_macros.h"
#include "def_grad.h"
#include "elem3d.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "index_macros.h"
#include "new_potentials.h"
#include "tensors.h"
#include "utils.h"
#include <mkl_cblas.h>
#include <errno.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace pgfem3d;
using namespace pgfem3d::net;

#ifndef MINI_3F_DEBUG
#define MINI_3F_DEBUG 0
#endif

#ifndef MINI_3F_P1_P1
#define MINI_3F_P1_P1 0
#endif

/* For testing */
static double compute_potential(const double *C,
                                const double Jn,
                                const double Tn,
                                const double Jr,
                                const double Tr,
                                const double p,
                                const double G,
                                const double K)
{
  const double trC = C[0] + C[4] + C[8];
  return (0.5*(
              G*(pow(Jr*Jn,-2./3.)*trC-3)
              +K*pow((Tr*Tn-1),2))
          +p*(Jr*Jn-Tr*Tn));
}

static void integration_help(const int ip,
                             const int ii,
                             const int nne,
                             const int total_nne,
                             const int ndn,
                             const int nPres,
                             const double *p,
                             const int nVol,
                             const double *int_pt_ksi,
                             const double *int_pt_eta,
                             const double *int_pt_zet,
                             const double *weights,
                             const double * x,
                             const double * y,
                             const double * z,
                             const SIG *sig,
                             const EPS *eps,
                             double *wt,
                             double *J,
                             double *pressure,
                             double *Tn,
                             double *Tr,
                             double *Na,
                             double *N_x,
                             double *N_y,
                             double *N_z,
                             double *Np,
                             double *Nt,
                             double ****ST_tensor,
                             double *ST)
{
  double ksi,eta,zet;
  ksi = int_pt_ksi[ip];
  eta = int_pt_eta[ip];
  zet = int_pt_zet[ip];
  *wt = weights[ip];

  shape_func(ksi,eta,zet,total_nne,Na);
  shape_func(ksi,eta,zet,nPres,Np);
  shape_func(ksi,eta,zet,nVol,Nt);
  *J = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
  get_bubble_grad(total_nne,ksi,eta,zet,x,y,z,N_x,N_y,N_z);
  shape_tensor (total_nne,ndn,N_x,N_y,N_z,ST_tensor);
  shapeTensor2array(ST,CONST_4(double) ST_tensor,total_nne);

  /* compute pressure */
  *pressure = *Tn = *Tr = 0.0;
  for(int i=0; i<nPres; i++){
    *pressure += Np[i]*(sig[ii].p[i] + p[i]);
  }

  for(int i=0; i<nVol; i++){
    *Tn += Nt[i]*eps[ii].T[i];
    *Tr += Nt[i]*eps[ii].d_T[i];
  }
}

void MINI_3f_reset(Element *elem,
                   const int nelem,
                   const int npres,
                   const int nvol,
                   SIG *sig,
                   EPS *eps)
{
  for(int i=0; i<nelem; i++){
    for(int j=0; j<npres; j++){
      sig[i].d_p[j] = 0.0;
    }
    for(int j=0; j<nvol; j++){
      eps[i].d_T[j] = 1.0;
    }

    memset(elem[i].d_bub_dofs,0,
           elem[i].n_bub*elem[i].n_bub_dofs*sizeof(double));
  }
}

int MINI_3f_stiffmat_el(double *Ks,            /**< Element stiffmat */
                        const int ii,          /**< element id */
                        const int ndofn,
                        const int nne,
                        const double *x,
                        const double *y,
                        const double *z,
                        const Element *elem,
                        const HOMMAT *hommat,
                        const long *nod,
                        const Node *node,
                        const EPS *eps,
                        const SIG *sig,
                        const double *r_e)    /**< dof values on elem */
{
  int count; /* ALWAYS reset before use */
  int err = 0;

  /* 3D */
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  const int nPres = nne;
  const int nVol = 4;

  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += ndofn;
  }

  const int n_bub = elem[ii].n_bub;
  const int n_bub_dofs = elem[ii].n_bub_dofs;
  const int total_nne = nne + n_bub;

  double *disp, *p;
  disp = aloc1(total_nne*ndn);
  p = aloc1(nPres);

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
        disp[i*ndn + j] = r_e[count+j];
      } else if (j==ndn){                       /* pressure dof */
        p[i] = r_e[count+j];
      } else {                                         /* ERROR */
        PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
        PGFEM_Abort();
      }
    } /* ndofn */
    count += ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] =
      elem[ii].bub_dofs[i*n_bub_dofs + j]
      + elem[ii].d_bub_dofs[i*n_bub_dofs + j];
    }
  }

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi = PGFEM_calloc (double, npt_z);
  double *int_pt_eta = PGFEM_calloc (double, npt_z);
  double *int_pt_zet = PGFEM_calloc (double, npt_z);
  double *weights = PGFEM_calloc (double, npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  Np = aloc1 (nPres);
  Nt = aloc1 (nVol);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne);

  /* allocate space for the stiffness tensors */
  double *Kuu, *Kbb, *Kub, *Kup, *Kbp, *Kpt, *Ktt, *Kpp;
  Kuu = aloc1(nne*ndn*nne*ndn);
  Kbb = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Kub = aloc1(nne*ndn*n_bub*n_bub_dofs);
  Kup = aloc1(nne*ndn*nne);
  Kbp = aloc1(n_bub*n_bub_dofs*nPres);
  Kpt = aloc1(nPres*nVol);
  Ktt = aloc1(nVol*nVol);
  Kpp = aloc1(nPres*nPres);

  /* allocate space for deformation gradients etc */
  double *Fn, **Fr_mat, *Fr, *Fr_lin, *Fr_I;
  double *C, *C_I, *AA;
  double *L, *S;
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_lin = aloc1(9);
  Fr_I = aloc1(9);
  C = aloc1(9);
  C_I = aloc1(9);
  AA = aloc1(9);
  L = aloc1(81);
  S = aloc1(9);

  /* element integration */
  double wt = 0.0;

  /* pressure */
  double pressure, Upp, Tn, Tr;

  memcpy(Fn,eps[ii].il[0].F,9*sizeof(double));
  const double Jn = getJacobian(Fn,ii,&err);

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
             int_pt_ksi,int_pt_eta,int_pt_zet,weights);
  for (int ip=0; ip<npt_z; ip++){
    integration_help(ip,ii,nne,total_nne,ndn,nPres,p,nVol,
                     int_pt_ksi,int_pt_eta,int_pt_zet,weights,
                     x,y,z,sig,eps,&wt,&J,&pressure,&Tn,&Tr,Na,
                     N_x,N_y,N_z,Np,Nt,ST_tensor,ST);

    /* get the linear portion of Fr */
    if(ip == 0){
      def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
      mat2array(Fr_lin,CONST_2(double) Fr_mat,3,3);
    }

    /* add the bubble portion of Fr */
    compute_disp_grad(total_nne,ST,disp,Fr,nne);
    cblas_daxpy(ndn*ndn,1.0,Fr_lin,1,Fr,1);

    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    /* total F */
    double alpha = 1.0;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                3,3,3,alpha,Fr,3,Fn,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                3,3,3,1.0,AA,3,AA,3,0.0,C,3);

    // Get Deviatoric 2 P-K stress
    devStressFuncPtr Stress;
    matStiffFuncPtr Stiffness;
    d2UdJ2FuncPtr UPP;
    Stress = getDevStressFunc(1,&hommat[mat]);
    Stiffness = getMatStiffFunc(1,&hommat[mat]);
    UPP = getD2UdJ2Func(1,&hommat[mat]);

    Stress(C,&hommat[mat],S);
    Stiffness(C,&hommat[mat],L);
    UPP(Tn*Tr,&hommat[mat],&Upp);

    /* Begin computing matrices */

    /* Kuu (12 x 12) */
    HW_Kuu_at_ip(Kuu,nne,total_nne,ST,
                 Fn,Fr,Fr_I,Jn,Tn,Jr,pressure,S,L,J,wt,0);

    /* Kbb (3 x 3) */
    HW_Kuu_at_ip(Kbb,nne,total_nne,ST,
                 Fn,Fr,Fr_I,Jn,Tn,Jr,pressure,S,L,J,wt,1);

    /* Kub (12 x 3) */
    HW_Kuu_at_ip(Kub,nne,total_nne,ST,
                 Fn,Fr,Fr_I,Jn,Tn,Jr,pressure,S,L,J,wt,2);

    /* Kup (12 x 4) */
    HW_Kup_at_ip(Kup,nne,total_nne,nPres,Np,ST,Fr_I,Jn,Tn,Jr,J,wt,0);

    /* Kbp (3 x 4) */
    HW_Kup_at_ip(Kbp,nne,total_nne,nPres,Np,ST,Fr_I,Jn,Tn,Jr,J,wt,1);

    /* Ktt (4 x 4) */
    HW_Ktt_at_ip(Ktt,nVol,Nt,Tn,Jn,kappa,Upp,J,wt);

    /* Kpt (4 x 4) */
    HW_Kpt_at_ip(Kpt,nPres,Np,nVol,Nt,Tn,Jn,J,wt);

  } /* int z-dir */

  /* free un-needed */
  free(disp);
  free(p);
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  free(Fn);
  free(Fr);
  free(Fr_lin);
  dealoc2(Fr_mat,3);
  free(Fr_I);
  free(AA);
  free(C);
  free(C_I);
  free(L);
  free(S);

  /* stuff for debugging */
  char filename[50];
  FILE *debug_log;
  int err_rank = 0;
  if(MINI_3F_DEBUG){
    PGFEM_Error_rank(&err_rank);
    sprintf(filename,"MINI_3f_stiff_%d_%d.log",ii,err_rank);
    debug_log = fopen(filename,"a"); /* APPEND TO FILE */
    PGFEM_fprintf(debug_log,"========================================\n");
    PGFEM_fprintf(debug_log,"============= New Computation ==========\n");
    PGFEM_fprintf(debug_log,"========================================\n\n");
    PGFEM_fprintf(debug_log,"Kuu\n");
    print_array_d(debug_log,Kuu,nne*ndn*nne*ndn,nne*ndn,nne*ndn);
    PGFEM_fprintf(debug_log,"Kbb\n");
    print_array_d(debug_log,Kbb,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
                  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kub\n");
    print_array_d(debug_log,Kub,nne*ndn*n_bub*n_bub_dofs,nne*ndn,
                  n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kup\n");
    print_array_d(debug_log,Kup,nne*ndn*nPres,nne*ndn,nPres);
    PGFEM_fprintf(debug_log,"Kbp\n");
    print_array_d(debug_log,Kbp,nPres*n_bub*n_bub_dofs,n_bub*n_bub_dofs,nPres);
    PGFEM_fprintf(debug_log,"Ktt\n");
    print_array_d(debug_log,Ktt,nVol,nVol,nVol);
    PGFEM_fprintf(debug_log,"Kpt\n");
    print_array_d(debug_log,Kpt,nPres*nVol,nPres,nVol);
  }

  /*** OPTIMIZE Unroll matrix multiplications */

  /*** Construct Condensed Matrices ***/
  double *Kbb_I, *Ktt_I, *KubKbbI, *KpbKbbI, *KptKttI;
  Kbb_I = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Ktt_I = aloc1(nVol*nVol);
  KubKbbI = aloc1(nne*ndn*n_bub*n_bub_dofs);
  KpbKbbI = aloc1(nPres*n_bub*n_bub_dofs);
  KptKttI = aloc1(nPres*nVol);

  if(!MINI_3F_P1_P1){
    if(inverse(Kbb,n_bub*n_bub_dofs,Kbb_I) != 0){
      PGFEM_printerr("ERROR computing Kbb_I in %s!\n",__func__);
      PGFEM_Abort();
    }
  }

  if(inverse(Ktt,nVol,Ktt_I) != 0){
    PGFEM_printerr("ERROR computing Ktt_I in %s!\n",__func__);
    PGFEM_Abort();
  }

  /* KubKbbI (12 x 3) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*ndn,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
              Kub,n_bub*n_bub_dofs,Kbb_I,n_bub*n_bub_dofs,
              0.0,KubKbbI,n_bub*n_bub_dofs);

  /* KpbKbbI (4 x 3) */
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              nPres,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
              Kbp,nPres,Kbb_I,n_bub*n_bub_dofs,
              0.0,KpbKbbI,n_bub*n_bub_dofs);

  /* KptKttI (4 x 1) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nPres,nVol,nVol,1.0,Kpt,nVol,Ktt_I,nVol,
              0.0,KptKttI,nVol);

  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Kbb_I\n");
    print_array_d(debug_log,Kbb_I,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
                  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kub Kbb_I\n");
    print_array_d(debug_log,KubKbbI,n_bub*n_bub_dofs*nne*ndn,
                  nne*ndn,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kpb Kbb_I\n");
    print_array_d(debug_log,KpbKbbI,nPres*n_bub*n_bub_dofs,nPres,
                  n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kpt Ktt_I\n");
    print_array_d(debug_log,KptKttI,nPres*nVol,nPres,nVol);
  }

  /*** Kuu  ***/
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
              nne*ndn,nne*ndn,n_bub*n_bub_dofs,-1.0,
              KubKbbI,n_bub*n_bub_dofs,Kub,n_bub*n_bub_dofs,
              1.0,Kuu,nne*ndn);

  /*** Kup ***/
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*ndn,nPres,n_bub*n_bub_dofs,-1.0,
              KubKbbI,n_bub*n_bub_dofs,Kbp,nPres,
              1.0,Kup,nne);

  /*** Kpp ***/
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nPres,nPres,n_bub*n_bub_dofs,-1.0,
              KpbKbbI,n_bub*n_bub_dofs,Kbp,nPres,
              0.0,Kpp,nPres);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
              nPres,nPres,nVol,-1.0,
              KptKttI,nVol,Kpt,nVol,
              1.0,Kpp,nPres);


  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Kuu BAR\n");
    print_array_d(debug_log,Kuu,nne*nne*ndn*ndn,nne*ndn,nne*ndn);
    PGFEM_fprintf(debug_log,"Kup BAR\n");
    print_array_d(debug_log,Kup,nne*ndn*nPres,nne*ndn,nPres);
    PGFEM_fprintf(debug_log,"Kpp BAR\n");
    print_array_d(debug_log,Kpp,nPres*nPres,nPres,nPres);
  }

  /* free un-needed */
  free(Kbb_I);
  free(Ktt_I);
  free(KubKbbI);
  free(KpbKbbI);
  free(KptKttI);

  /*** Scatter the matrix into the proper form ***/
  for(int a=0; a<nne; a++){
    for(int b=0; b<ndofn; b++){
      for(int w=0; w<nne; w++){
        for(int g=0; g<ndofn; g++){
          if ((b<ndn) && (g<ndn)){ /* from Kuu */
            Ks[idx_K(a,b,w,g,nne,ndofn)] = Kuu[idx_K(a,b,w,g,nne,ndn)];
          } else if ((b<ndn) && (g==ndn)){ /* from Kup */
            Ks[idx_K(a,b,w,g,nne,ndofn)] =
            Kup[idx_K_gen(a,b,w,0,nne,ndn,nPres,1)];
          } else if ((b==ndn) && (g<ndn)){ /* from Kpu = (Kup)' */
            Ks[idx_K(a,b,w,g,nne,ndofn)] =
            Kup[idx_K_gen(w,g,a,0,nne,ndn,nPres,1)];
          } else if ((b==ndn) && (g==ndn)){ /* from Kpp */
            Ks[idx_K(a,b,w,g,nne,ndofn)] = Kpp[idx_K(a,0,w,0,nPres,1)];
          } else if ((b==ndn) && (g>ndn)){ /* Kpt */
            Ks[idx_K(a,b,w,g,nne,ndofn)] =
            Kpt[idx_K_gen(a,0,w,0,nPres,1,nVol,1)];
          } else if ((b>ndn) && (g==ndn)){ /* Ktp = (Kpt)' */
            Ks[idx_K(a,b,w,g,nne,ndofn)] =
            Kpt[idx_K_gen(w,0,a,0,nPres,1,nVol,1)];
          } else {
            PGFEM_printerr("ERROR constructing stiffness matrix!\n");
            PGFEM_Abort();
            abort();
          }
        }
      }
    }
  }

  /* free remaining */
  free(Kuu);
  free(Kbb);
  free(Kpp);
  free(Ktt);
  free(Kub);
  free(Kup);
  free(Kbp);
  free(Kpt);


  /* print out the element stiffness */
  if (MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Element Stiffness\n");
    print_array_d(debug_log,Ks,nne*ndofn*nne*ndofn,nne*ndofn,nne*ndofn);
    fclose(debug_log);
  }

  return err;
} /* MINI_3f_stiffmat_el */


int MINI_3f_resid_el(double *Res,         /**< Element residual */
                     const int ii,        /**< element id */
                     const int ndofn,
                     const int nne,
                     const double *x,
                     const double *y,
                     const double *z,
                     const Element *elem,
                     const long *nod,
                     const Node *node,
                     const HOMMAT *hommat,
                     const EPS *eps,
                     const SIG *sig,
                     const double *r_e)    /**< dof values on elem */
{
  int count; /* ALWAYS reset before use */
  int err = 0;

  /* 3D */
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  const int nPres = nne;
  const int nVol = 4;

  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += ndofn;
  }

  int n_bub = elem[ii].n_bub;
  int n_bub_dofs = elem[ii].n_bub_dofs;
  int total_nne = nne + n_bub;

  double *disp, *p;
  disp = aloc1(total_nne*ndn);
  p = aloc1(nPres);

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
        disp[i*ndn + j] = r_e[count+j];
      } else if (j==ndn){                       /* pressure dof */
        p[i] = r_e[count+j];
      } else {                                         /* ERROR */
        PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
        PGFEM_Abort();
      }
    } /* ndofn */
    count += ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] =
      elem[ii].bub_dofs[i*n_bub_dofs + j]
      + elem[ii].d_bub_dofs[i*n_bub_dofs + j];
    }
  }

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi = PGFEM_calloc (double, npt_z);
  double *int_pt_eta = PGFEM_calloc (double, npt_z);
  double *int_pt_zet = PGFEM_calloc (double, npt_z);
  double *weights = PGFEM_calloc (double, npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  Np = aloc1 (nPres);
  Nt = aloc1 (nVol);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

  /* allocate space for deformation gradients etc */
  double *Fn, **Fr_mat, *Fr, *Fr_lin, *Fr_I;
  double *C, *AA;
  double *L, *S;
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_lin = aloc1(9);
  Fr_I = aloc1(9);
  AA = aloc1(9);
  C = aloc1(9);
  L = aloc1(81);
  S = aloc1(9);

  /* Allocate space for the residuals */
  double *Ru, *Rb, *Rp, *Rt;
  Ru = aloc1(nne*ndn);
  Rb = aloc1(n_bub*n_bub_dofs);
  Rp = aloc1(nne);
  Rt = aloc1(nne);

  /* allocate space for the stiffness tensors */
  double *Kbb, *Kub, *Kbp, *Kpt, *Ktt;
  Kbb = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Kub = aloc1(nne*ndn*n_bub*n_bub_dofs);
  Kbp = aloc1(n_bub*n_bub_dofs*nne);
  Kpt = aloc1(nPres*nVol);
  Ktt = aloc1(nVol*nVol);

  /* integration point in natural coords */
  double wt;

  /* pressure */
  double pressure, Up, Upp, Tn, Tr;

  memcpy(Fn,eps[ii].il[0].F,9*sizeof(double));
  const double Jn = getJacobian(Fn,ii,&err);

  /* stuff for debugging */
  char filename[50];
  FILE *debug_log;
  int err_rank = 0;
  if(MINI_3F_DEBUG){
    PGFEM_Error_rank(&err_rank);
    sprintf(filename,"MINI_3f_resid_%d_%d.log",ii,err_rank);
    debug_log = fopen(filename,"a"); /* APPEND TO FILE */
    PGFEM_fprintf(debug_log,"========================================\n");
    PGFEM_fprintf(debug_log,"============= New Computation ==========\n");
    PGFEM_fprintf(debug_log,"========================================\n\n");
    PGFEM_fprintf(debug_log,"Coords:\n");
    print_coords(debug_log,total_nne,x,y,z);
    PGFEM_fprintf(debug_log,"Disp\n");
    print_array_d(debug_log,disp,total_nne*ndn,total_nne,ndn);
    PGFEM_fprintf(debug_log,"pressure: ");
    print_array_d(debug_log,sig[ii].p,nPres,1,nPres);
    PGFEM_fprintf(debug_log,"d_pressure: ");
    print_array_d(debug_log,p,nPres,1,nPres);
    PGFEM_fprintf(debug_log,"Tn: ");
    print_array_d(debug_log,eps[ii].T,nVol,1,nVol);
    PGFEM_fprintf(debug_log,"Tr: ");
    print_array_d(debug_log,eps[ii].d_T,nVol,1,nVol);
  }

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
             int_pt_ksi,int_pt_eta,int_pt_zet,weights);
  for (int ip=0; ip<npt_z; ip++){
    integration_help(ip,ii,nne,total_nne,ndn,nPres,p,nVol,
                     int_pt_ksi,int_pt_eta,int_pt_zet,weights,
                     x,y,z,sig,eps,&wt,&J,&pressure,&Tn,&Tr,Na,
                     N_x,N_y,N_z,Np,Nt,ST_tensor,ST);

    /* get the linear portion of Fr */
    if(ip == 0){
      def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
      mat2array(Fr_lin,CONST_2(double) Fr_mat,3,3);
    }

    /* add the bubble portion of Fr */
    compute_disp_grad(total_nne,ST,disp,Fr,nne);
    cblas_daxpy(ndn*ndn,1.0,Fr_lin,1,Fr,1);

    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    /* Compute total F; F=FrFn */
    double alpha = 1.0;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                3,3,3,alpha,Fr,3,Fn,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                3,3,3,1.0,AA,3,AA,3,0.0,C,3);

    // Get Deviatoric 2 P-K stress
    devStressFuncPtr Stress;
    matStiffFuncPtr Stiffness;
    dUdJFuncPtr UP;
    d2UdJ2FuncPtr UPP;
    Stress = getDevStressFunc(1,&hommat[mat]);
    Stiffness = getMatStiffFunc(1,&hommat[mat]);
    UP = getDUdJFunc(1,&hommat[mat]);
    UPP = getD2UdJ2Func(1,&hommat[mat]);

    Stress(C,&hommat[mat],S);
    Stiffness(C,&hommat[mat],L);
    UP(Tn*Tr,&hommat[mat],&Up);
    UPP(Tn*Tr,&hommat[mat],&Upp);

    if(MINI_3F_DEBUG){
      PGFEM_fprintf(debug_log,"=== ip: %d ===\n",ip);
      PGFEM_fprintf(debug_log,"Jn: %22.15e Jr: %22.15e Tn: %22.15e Tr: %22.15e "
                    "p: %22.15e\n",Jn,Jr,Tn,Tr,pressure);
      PGFEM_fprintf(debug_log,"(JrJn-TrTn): %22.15e\n",(Jr*Jn-Tr*Tn));
      PGFEM_fprintf(debug_log,"(K Up - p):  %22.15e\n",(kappa*Up-pressure));
      PGFEM_fprintf(debug_log,"W(C,T)+p(J-T): %22.15e\n",
                    compute_potential(C,Jn,Tn,Jr,Tr,pressure,hommat[mat].G,kappa));
      PGFEM_fprintf(debug_log,"S:");
      print_array_d(debug_log,S,9,1,9);
      PGFEM_fprintf(debug_log,"Fr:");
      print_array_d(debug_log,Fr,9,1,9);
      PGFEM_fprintf(debug_log,"Fn:");
      print_array_d(debug_log,Fn,9,1,9);
      PGFEM_fprintf(debug_log,"C:");
      print_array_d(debug_log,C,9,1,9);
    }

    /* Compute residuals */
    /* Ru (12 x 1) */
    HW_Ru_at_ip(Ru,nne,total_nne,ST,Fn,Fr,Fr_I,
                Jn,Tn,Jr,Tr,S,pressure,J,wt,0);

    /* Rb (3 x 1) */
    HW_Ru_at_ip(Rb,nne,total_nne,ST,Fn,Fr,Fr_I,
                Jn,Tn,Jr,Tr,S,pressure,J,wt,1);

    /* Rp (4 x 1) */
    HW_Rp_at_ip(Rp,nPres,Np,Jn,Jr,Tn,Tr,J,wt);

    /* Rt (4 x 1) */
    HW_Rt_at_ip(Rt,nVol,Nt,Tn,Jn,pressure,kappa,Up,J,wt);

    /* Compute needed tangents */
    /* Kbb (3 x 3) */
    HW_Kuu_at_ip(Kbb,nne,total_nne,ST,Fn,Fr,Fr_I,
                 Jn,Tn,Jr,pressure,S,L,J,wt,1);

    /* Kub (12 x 3) */
    HW_Kuu_at_ip(Kub,nne,total_nne,ST,Fn,Fr,Fr_I,
                 Jn,Tn,Jr,pressure,S,L,J,wt,2);

    /* Kbp (3 x 4) */
    HW_Kup_at_ip(Kbp,nne,total_nne,nPres,Np,ST,Fr_I,Jn,Tn,Jr,J,wt,1);

    /* Kpt (4 x 1) */
    HW_Kpt_at_ip(Kpt,nPres,Np,nVol,Nt,Tn,Jn,J,wt);

    /* Ktt (1 x 1) */
    HW_Ktt_at_ip(Ktt,nVol,Nt,Tn,Jn,kappa,Upp,J,wt);

  }/* end integration */

  /* free un-needed */
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  free(Fn);
  dealoc2(Fr_mat,3);
  free(Fr);
  free(Fr_lin);
  free(Fr_I);
  free(AA);
  free(C);
  free(L);
  free(S);

  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Ru\n");
    print_array_d(debug_log,Ru,nne*ndn,1,nne*ndn);
    PGFEM_fprintf(debug_log,"Rb\n");
    print_array_d(debug_log,Rb,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Rp\n");
    print_array_d(debug_log,Rp,nPres,1,nPres);
    PGFEM_fprintf(debug_log,"Rt\n");
    print_array_d(debug_log,Rt,nVol,1,nVol);
  }

  free(disp);
  free(p);

  /*** OPTIMIZE Unroll mat mults ***/

  /*** Construct condensed residuals ***/
  double *Kbb_I, *Ktt_I, *KubKbbI, *KpbKbbI, *KptKttI;
  Kbb_I = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Ktt_I = aloc1(nVol*nVol);
  KubKbbI = aloc1(nne*ndn*n_bub*n_bub_dofs);
  KpbKbbI = aloc1(nne*n_bub*n_bub_dofs);
  KptKttI = aloc1(nPres*nVol);

  if(!MINI_3F_P1_P1){
    if(inverse(Kbb,n_bub*n_bub_dofs,Kbb_I) != 0){
      PGFEM_printerr("ERROR computing Kbb_I in %s!\n",__func__);
      PGFEM_Abort();
    }
  }

  if(inverse(Ktt,nVol,Ktt_I) != 0){
    PGFEM_printerr("ERROR computing Ktt_I in %s!\n",__func__);
    PGFEM_Abort();
  }

  /* KubKbbI (12 x 3) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*ndn,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
              Kub,n_bub*n_bub_dofs,Kbb_I,n_bub*n_bub_dofs,
              0.0,KubKbbI,n_bub*n_bub_dofs);

  /* KpbKbbI (4 x 3) */
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              nne,n_bub*n_bub_dofs,n_bub*n_bub_dofs,1.0,
              Kbp,nne,Kbb_I,n_bub*n_bub_dofs,
              0.0,KpbKbbI,n_bub*n_bub_dofs);

  /* KptKttI (4 x 1) */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nPres,nVol,nVol,1.0,Kpt,nVol,Ktt_I,nVol,
              0.0,KptKttI,nVol);

  /*** Ru (12 x 1) ***/
  cblas_dgemv(CblasRowMajor,CblasNoTrans,
              nne*ndn,n_bub*n_bub_dofs,1.0,
              KubKbbI,n_bub*n_bub_dofs,Rb,1,
              -1.0,Ru,1);

  /*** Rp (4 x 1) ***/
  cblas_dgemv(CblasRowMajor,CblasNoTrans,
              nPres,n_bub*n_bub_dofs,1.0,
              KpbKbbI,n_bub*n_bub_dofs,Rb,1,
              -1.0,Rp,1);

  cblas_dgemv(CblasRowMajor,CblasNoTrans,
              nPres,nVol,1.0,
              KptKttI,nVol,Rt,1,
              1.0,Rp,1);

  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Ru BAR\n");
    print_array_d(debug_log,Ru,nne*ndn,1,nne*ndn);
    PGFEM_fprintf(debug_log,"Rp BAR\n");
    print_array_d(debug_log,Rp,nPres,1,nPres);

    /* print matrices */
    PGFEM_fprintf(debug_log,"Kub\n");
    print_array_d(debug_log,Kub,nne*ndn*n_bub*n_bub_dofs,
                  nne*ndn,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kbb\n");
    print_array_d(debug_log,Kbb,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
                  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kbp\n");
    print_array_d(debug_log,Kbp,n_bub*n_bub_dofs*nPres,
                  n_bub*n_bub_dofs,nPres);
    PGFEM_fprintf(debug_log,"Kpt\n");
    print_array_d(debug_log,Kpt,nPres*nVol,nPres,nVol);
    PGFEM_fprintf(debug_log,"Ktt\n");
    print_array_d(debug_log,Ktt,nVol*nVol,nVol,nVol);
  }
  /* free un-needed */
  free(Kbb);
  free(Kbp);
  free(Kub);
  free(Ktt);
  free(Kpt);
  free(Kbb_I);
  free(Ktt_I);
  free(KubKbbI);
  free(KpbKbbI);
  free(KptKttI);
  free(Rb);
  free(Rt);

  int sign = -1;
  /*** Scatter the residuals (u,v,w,p,t) ***/
  /* NOTE residuals are negative because subtracted from RHS later */
  for (int a=0; a<nne; a++){
    for (int b=0; b<ndofn; b++){
      if (b<ndn){ /* from Ru */
        Res[a*ndofn+b] = sign*Ru[a*ndn+b];
      } else if (b == ndn){ /* from Rp */
        Res[a*ndofn+b] = sign*Rp[a];
      }
    }
  }

  /* print out the element residual */
  if (MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"========== Element Residual ==========\n");
    print_array_d(debug_log,Res,nne*ndofn,1,nne*ndofn);
    fclose(debug_log);
  }

  /* free remaining */
  free(Ru);
  free(Rp);

  return err;

} /* MINI_3f_resid_el */

int MINI_3f_update_bubble_el(Element *elem,
                             const int ii, /* id of element working on */
                             const int nne,
                             const Node *node,
                             const int ndofn,
                             const double *x,
                             const double *y,
                             const double *z,
                             const long *nod, /* list of node ids on elem */
                             const EPS *eps,
                             const SIG *sig,
                             const HOMMAT *hommat,
                             const double *sol_e, /* accum. incr on elem */
                             const double *dsol_e) /* increment on element */
{
  int err = 0;
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  const int nPres = nne;
  const int nVol = 4;

  int count; /* ALWAYS reset before use */

  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += ndofn;
  }

  const int n_bub = elem[ii].n_bub;
  const int n_bub_dofs = elem[ii].n_bub_dofs;
  const int total_nne = nne + n_bub;

  /* stuff for debugging purposes */
  FILE *debug_log;
  int err_rank = 0;
  errno = 0;

  if(MINI_3F_DEBUG){
    PGFEM_Error_rank(&err_rank);
    char fname[256];
    sprintf(fname,"MINI_3f_up_%d_%d.log",ii,err_rank);
    debug_log = fopen(fname,"a");
    PGFEM_fprintf(debug_log,"***********************************\n");
    PGFEM_fprintf(debug_log,"Current bub dofs\n");
    print_array_d(debug_log,elem[ii].bub_dofs,
                  n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  /* disp contains displacements from nodes + bub */
  double *disp, *ddisp, *p, *dp;
  disp = aloc1(total_nne*ndn);
  ddisp = aloc1(nne*ndn);
  p = aloc1(nne);
  dp = aloc1(nne);

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
        disp[i*ndn + j] = sol_e[count+j];   /* accum */
        ddisp[i*ndn + j] = dsol_e[count+j]; /* increment */
      } else if (j==ndn){                   /* pressure dof */
        p[i] = sol_e[count+j];              /* accum */
        dp[i] = dsol_e[count+j];            /* increment */
      } else {                                         /* ERROR */
        PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
        PGFEM_Abort();
      }
    } /* ndofn */
    count += ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] =
      elem[ii].bub_dofs[i*n_bub_dofs + j]
      + elem[ii].d_bub_dofs[i*n_bub_dofs + j];
    }
  }

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi = PGFEM_calloc (double, npt_z);
  double *int_pt_eta = PGFEM_calloc (double, npt_z);
  double *int_pt_zet = PGFEM_calloc (double, npt_z);
  double *weights = PGFEM_calloc (double, npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  Np = aloc1 (nPres);
  Nt = aloc1 (nVol);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne);

  /* allocate space for the stiffness tensors */
  double *Kbb, *Kub, *Kbp, *Ktt, *Kpt, *Rb, *Rt;
  Kbb = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Kub = aloc1(nne*ndn*n_bub*n_bub_dofs);
  Kbp = aloc1(n_bub*n_bub_dofs*nne);
  Ktt = aloc1(nVol*nVol);
  Kpt = aloc1(nPres*nVol);
  Rb = aloc1(n_bub*n_bub_dofs);
  Rt = aloc1(nVol);

  /* allocate space for deformation gradients etc */
  double *Fn, **Fr_mat, *Fr, *Fr_lin, *Fr_I;
  double *C, *C_I, *AA;
  double *L, *S;
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_lin = aloc1(9);
  Fr_I = aloc1(9);
  C = aloc1(9);
  C_I = aloc1(9);
  AA = aloc1(9);
  L = aloc1(81);
  S = aloc1(9);

  /* integration point in natural coords */
  double wt;

  /* pressure */
  double pressure, Tr, Tn, Up, Upp;

  memcpy(Fn,eps[ii].il[0].F,9*sizeof(double));
  const double Jn = getJacobian(Fn,ii,&err);

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
             int_pt_ksi,int_pt_eta,int_pt_zet,weights);
  for (int ip=0; ip<npt_z; ip++){
    integration_help(ip,ii,nne,total_nne,ndn,nPres,p,nVol,
                     int_pt_ksi,int_pt_eta,int_pt_zet,weights,
                     x,y,z,sig,eps,&wt,&J,&pressure,&Tn,&Tr,Na,
                     N_x,N_y,N_z,Np,Nt,ST_tensor,ST);

    /* get the linear portion of Fr */
    if(ip == 0){
      def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
      mat2array(Fr_lin,CONST_2(double) Fr_mat,3,3);
    }

    /* add the bubble portion of Fr */
    compute_disp_grad(total_nne,ST,disp,Fr,nne);
    cblas_daxpy(ndn*ndn,1.0,Fr_lin,1,Fr,1);

    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    double alpha = 1.0;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                3,3,3,alpha,Fr,3,Fn,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                3,3,3,1.0,AA,3,AA,3,0.0,C,3);

    // Get Deviatoric 2 P-K stress
    devStressFuncPtr Stress;
    matStiffFuncPtr Stiffness;
    dUdJFuncPtr UP;
    d2UdJ2FuncPtr UPP;
    Stress = getDevStressFunc(1,&hommat[mat]);
    Stiffness = getMatStiffFunc(1,&hommat[mat]);
    UP = getDUdJFunc(1,&hommat[mat]);
    UPP = getD2UdJ2Func(1,&hommat[mat]);

    Stress(C,&hommat[mat],S);
    Stiffness(C,&hommat[mat],L);
    UP(Tn*Tr,&hommat[mat],&Up);
    UPP(Tn*Tr,&hommat[mat],&Upp);

    /* Begin computing matrices */
    /* Kbb (3 x 3) */
    HW_Kuu_at_ip(Kbb,nne,total_nne,ST,Fn,Fr,Fr_I,
                 Jn,Tn,Jr,pressure,S,L,J,wt,1);

    /* Kub (12 x 3) */
    HW_Kuu_at_ip(Kub,nne,total_nne,ST,Fn,Fr,Fr_I,
                 Jn,Tn,Jr,pressure,S,L,J,wt,2);

    /* Kbp (3 x 4) */
    HW_Kup_at_ip(Kbp,nne,total_nne,nPres,Np,ST,Fr_I,Jn,Tn,Jr,J,wt,1);

    /* Ktt (1 x 1) */
    HW_Ktt_at_ip(Ktt,nVol,Nt,Tn,Jn,kappa,Upp,J,wt);

    /* Kpt (4 x 1) */
    HW_Kpt_at_ip(Kpt,nPres,Np,nVol,Nt,Tn,Jn,J,wt);

    /* Rb (3 x 1) */
    HW_Ru_at_ip(Rb,nne,total_nne,ST,Fn,Fr,Fr_I,
                Jn,Tn,Jr,Tr,S,pressure,J,wt,1);

    /* Rt (1 x 1) */
    HW_Rt_at_ip(Rt,nVol,Nt,Tn,Jn,pressure,kappa,Up,J,wt);

  } /* int z-dir */

  /* free un-needed */
  free(disp);
  free(p);
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  free(Fn);
  dealoc2(Fr_mat,3);
  free(Fr);
  free(Fr_lin);
  free(Fr_I);
  free(AA);
  free(C);
  free(C_I);
  free(L);
  free(S);

  /*** Compute the bubble update ***/
  double *Kbb_I, *Ktt_I, *ddb, *ddT;
  ddb = aloc1(n_bub*n_bub_dofs);
  ddT = aloc1(nVol);
  Kbb_I = aloc1(n_bub*n_bub*n_bub_dofs*n_bub_dofs);
  Ktt_I = aloc1(nVol*nVol);

  if(!MINI_3F_P1_P1){
    if(inverse(Kbb,n_bub*n_bub_dofs,Kbb_I) != 0){
      PGFEM_printerr("ERROR computing Kbb_I in %s!\n",__func__);
      PGFEM_Abort();
    }
  }

  if(inverse(Ktt,nVol,Ktt_I) != 0){
    PGFEM_printerr("ERROR computing Ktt_I in %s!\n",__func__);
    PGFEM_Abort();
  }

  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Kbb\n");
    print_array_d(debug_log,Kbb,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
                  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kbb_I\n");
    print_array_d(debug_log,Kbb_I,n_bub*n_bub_dofs*n_bub*n_bub_dofs,
                  n_bub*n_bub_dofs,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kub\n");
    print_array_d(debug_log,Kub,n_bub*n_bub_dofs*nne*ndn,
                  nne*ndn,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Kbp\n");
    print_array_d(debug_log,Kbp,n_bub*n_bub_dofs*nPres,
                  n_bub*n_bub_dofs,nPres);

    PGFEM_fprintf(debug_log,"Ktt\n");
    print_array_d(debug_log,Ktt,nVol*nVol,nVol,nVol);
    PGFEM_fprintf(debug_log,"Ktt_I\n");
    print_array_d(debug_log,Ktt_I,nVol*nVol,nVol,nVol);
    PGFEM_fprintf(debug_log,"Kpt\n");
    print_array_d(debug_log,Kpt,nPres*nVol,nPres,nVol);

    PGFEM_fprintf(debug_log,"ddisp\n");
    print_array_d(debug_log,ddisp,nne*ndn,1,nne*ndn);
    PGFEM_fprintf(debug_log,"dp\n");
    print_array_d(debug_log,dp,nPres,1,nPres);
    PGFEM_fprintf(debug_log,"Volume\n");
    print_array_d(debug_log,eps[ii].d_T,nVol,1,nVol);

    PGFEM_fprintf(debug_log,"Rb\n");
    print_array_d(debug_log,Rb,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Rt\n");
    print_array_d(debug_log,Rt,nVol,1,nVol);
  }

  /*** OPTIMIZE Unroll mat ops ***/

  /* Rb += Kub'ddu */
  cblas_dgemv(CblasRowMajor,CblasTrans,nne*ndn,n_bub_dofs*n_bub,
              1.0,Kub,n_bub_dofs*n_bub,ddisp,1,1.0,Rb,1);

  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Rb += Kub' ddisp\n");
    print_array_d(debug_log,Rb,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  /* Rb += Kbp ddp */
  cblas_dgemv(CblasRowMajor,CblasNoTrans,n_bub*n_bub_dofs,nPres,
              1.0,Kbp,nPres,dp,1,1.0,Rb,1);

  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Rb += Kbp dp\n");
    print_array_d(debug_log,Rb,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  /* ddb = -Kbb_I Rb */
  cblas_dgemv(CblasRowMajor,CblasNoTrans,n_bub*n_bub_dofs,n_bub*n_bub_dofs,
              -1.0,Kbb_I,3,Rb,1,0.0,ddb,1);

  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"db\n");
    print_array_d(debug_log,ddb,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
  }

  for (int i=0; i<n_bub*n_bub_dofs; i++){
    elem[ii].d_bub_dofs[i] += ddb[i];
  }

  /* Rt += Kpt' dp */
  cblas_dgemv(CblasRowMajor,CblasTrans,nPres,nVol,
              1.0,Kpt,nVol,dp,1,1.0,Rt,1);
  if(MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Rt += Kpt' dp\n");
    print_array_d(debug_log,Rt,nVol,1,nVol);
  }

  /* ddT = -Ktt_I Rt */
  cblas_dgemv(CblasRowMajor,CblasNoTrans,nVol,nVol,
              -1.0,Ktt_I,nVol,Rt,1,0.0,ddT,1);

  for(int i=0; i<nVol; i++){
    eps[ii].d_T[i] += ddT[i];
  }

  if (MINI_3F_DEBUG){
    PGFEM_fprintf(debug_log,"Updated bubble dofs\n");
    print_array_d(debug_log,elem[ii].bub_dofs,
                  n_bub*n_bub_dofs,
                  1,n_bub*n_bub_dofs);
    PGFEM_fprintf(debug_log,"Updated bubble vol\n");
    print_array_d(debug_log,eps[ii].d_T,nVol,1,nVol);
    fclose(debug_log);
  }

  free(Kbb);
  free(Kbb_I);
  free(Kub);
  free(Kbp);
  free(Rb);
  free(ddb);
  free(ddisp);
  free(dp);

  return err;
}/* MINI_3f_update_bubble_el */

int MINI_3f_update_bubble(Element *elem,
                          const int nelem,
                          const Node *node,
                          const int ndofn,
                          const SUPP sup,
                          const EPS *eps,
                          const SIG *sig,
                          const HOMMAT *hommat,
                          const double *sol, /* accum. solution  on incr */
                          const double *dsol, /* sol from current iter */
                          const int iter,
                          const int mp_id)
{
  int err = 0;

  int nne, nne_t;

  long *nod, *cn;
  double *x, *y, *z, *sup_def;
  double *sol_e, *dsol_e;

  sup_def = aloc1(sup->npd);

  /* for each element */
  for (int i=0; i<nelem; i++){
    nne = elem[i].toe;
    nne_t = nne + elem[i].n_bub;

    // @todo Commented as dead code. @cp please review. LD
    // int mat = elem[i].mat[2];

    /* allocate */
    nod = aloc1l(nne);
    cn = aloc1l(nne*ndofn);
    x = aloc1(nne_t);
    y = aloc1(nne_t);
    z = aloc1(nne_t);
    sol_e = aloc1(nne*ndofn);
    dsol_e = aloc1(nne*ndofn);

    /* get node ids on element */
    elemnodes(i,nne,nod,elem);

    /* get local dof ids on elemnt */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* Get the node and bubble coordinates */
    nodecoord_updated(nne,nod,node,x,y,z);
    element_center(nne,x,y,z);

    /*** compute the iteration and accumlated solution vectors on the
         element ***/
    /* solution for this iteration */
    if (iter != 0){
      for (int j=0;j<sup->npd;j++){
        sup_def[j] = sup->defl_d[j];
        sup->defl_d[j] = 0.0;
      }
    }
    def_elem(cn,nne*ndofn,dsol,elem,node,dsol_e,sup,0);
    if (iter != 0){
      for (int j=0;j<sup->npd;j++)
        sup->defl_d[j] = sup_def[j];
    }
    /* accumulated solution for this load increment */
    if (iter == 0){
      for (int j=0;j<sup->npd;j++){
        sup_def[j] = sup->defl_d[j];
        sup->defl_d[j] = 0.0;
      }
    }
    def_elem(cn,nne*ndofn,sol,elem,node,sol_e,sup,0);
    if (iter == 0){
      for (int j=0;j<sup->npd;j++)
        sup->defl_d[j] = sup_def[j];
    }

    /* compute the bubble update on the element */
    err = MINI_3f_update_bubble_el(elem,i,nne,node,ndofn,x,y,z,nod,
                                   eps,sig,hommat,sol_e,dsol_e);


    free(nod);
    free(cn);
    free(x);
    free(y);
    free(z);
    free(sol_e);
    free(dsol_e);
  }
  free(sup_def);

  return err;

}/* MINI_3f_update_bubble */

void MINI_3f_increment_el(Element *elem,
                          const int ii, /* id of element working on */
                          const int nne,
                          const Node *node,
                          const long *nod,
                          const int ndofn,
                          const double *x,
                          const double *y,
                          const double *z,
                          EPS *eps,
                          SIG *sig,
                          const HOMMAT *hommat,
                          const double *sol_e)
{
  int err = 0;
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  const int nPres = nne;
  const int nVol = 4;

  int count; /* ALWAYS reset before use */

  /* Element Dof */
  int ndofe = 0; /* ndofs on element without bubble */
  for (int i=0; i<nne; i++){
    ndofe += ndofn;
  }

  const int n_bub = elem[ii].n_bub;
  const int n_bub_dofs = elem[ii].n_bub_dofs;
  const int total_nne = nne + n_bub;

  double *disp, *p;
  disp = aloc1(total_nne*ndn);
  p = aloc1(nne);

  count = 0;
  for (int i=0; i<nne; i++){
    for (int j=0; j<ndofn; j++){

      /* filter displacement and pressure from element unknowns */
      if (j<ndn){                           /* displacement dof */
        disp[i*ndn + j] = sol_e[count+j];
      } else if (j==ndn){                       /* pressure dof */
        p[i] = sol_e[count+j];
      } else {                                         /* ERROR */
        PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
        PGFEM_Abort();
      }
    } /* ndofn */
    count += ndofn;
  } /* nne */

  /* displacements from bubble */
  for (int i=0; i<n_bub; i++){
    for (int j=0; j<n_bub_dofs; j++){
      disp[nne*ndn + i*n_bub_dofs + j] =
      elem[ii].bub_dofs[i*n_bub_dofs + j]
      + elem[ii].d_bub_dofs[i*n_bub_dofs + j];
    }
  }

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(total_nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);

  double *weights;
  weights = aloc1(npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  Na = aloc1 (total_nne);
  Np = aloc1 (nPres);
  Nt = aloc1 (nVol);
  N_x = aloc1 (total_nne);
  N_y = aloc1 (total_nne);
  N_z = aloc1 (total_nne);
  ST_tensor = aloc4 (3,3,ndn,total_nne);
  ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

  /* allocate space for deformation gradients etc */
  double *Fn, **Fr_mat, *Fr, *Fr_lin, *Fr_I;
  double *C, *C_I, *E, *AA, *F_total, *F_total_I;
  double *S;
  Fn = aloc1(9);
  Fr_mat = aloc2(3,3);
  Fr = aloc1(9);
  Fr_lin = aloc1(9);
  Fr_I = aloc1(9);
  C = aloc1(9);
  C_I = aloc1(9);
  E = aloc1(9);
  AA = aloc1(9);
  F_total = aloc1(9);
  F_total_I = aloc1(9);
  S = aloc1(9);

  /* identity */
  double *ident;
  ident = aloc1(9);
  ident[0] = ident[4] = ident[8] = 1.0;

  /* integration point in natural coords */
  double wt;

  /* pressure */
  double pressure,Tr,Tn;

  /* element volume */
  const double volume = Tetra_V(x,y,z);

  /* zero the stress & strain */
  memset(sig[ii].el.o,0,6*sizeof(double));
  memset(eps[ii].el.o,0,6*sizeof(double));

  /* linear Fn */
  memcpy(Fn,eps[ii].il[0].F,9*sizeof(double));
  const double Jn = getJacobian(Fn,ii,&err);

  /* get the integration points and weights */
  integrate (total_nne,&npt_x,&npt_y,&npt_z,
             int_pt_ksi,int_pt_eta,int_pt_zet,weights);
  for (int ip=0; ip<npt_z; ip++){
    integration_help(ip,ii,nne,total_nne,ndn,nPres,p,nVol,
                     int_pt_ksi,int_pt_eta,int_pt_zet,weights,
                     x,y,z,sig,eps,&wt,&J,&pressure,&Tn,&Tr,Na,
                     N_x,N_y,N_z,Np,Nt,ST_tensor,ST);

    /* get the linear portion of Fr */
    if(ip == 0){
      def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
      mat2array(Fr_lin,CONST_2(double) Fr_mat,3,3);
    }

    /* add the bubble portion of Fr */
    compute_disp_grad(total_nne,ST,disp,Fr,nne);
    cblas_daxpy(ndn*ndn,1.0,Fr_lin,1,Fr,1);

    double Jr = getJacobian(Fr,ii,&err);
    inverse(Fr,3,Fr_I);

    /* CORRECTED TOTAL deformation gradient */
    double alpha = pow(Tn*Tr/(Jn*Jr),1./3.);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                3,3,3,alpha,Fr,3,Fn,3,0.0,F_total,3);

    /* Green-Lagrange deformation tensor */
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                3,3,3,1.0,F_total,3,F_total,3,0.0,C,3);
    inverse(C,3,C_I);


    // Get Deviatoric 2 P-K stress
    devStressFuncPtr Stress;
    Stress = getDevStressFunc(1,&hommat[mat]);
    Stress(C,&hommat[mat],S);

    /* Compute total stress S = dev(S) + p*Tn*Tr*C_I */
    cblas_daxpy(9,pressure*Tn*Tr,C_I,1,S,1);

    /* Compute G-L strain E = 0.5(C - 1) */
    for(int i=0; i<9; i++){
      E[i] = 0.5*(C[i] - ident[i]);
    }

    /* Compute Cauchy stress (Tn*Tr)^-1 F S F' */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                3,3,3,1.0,F_total,3,S,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                3,3,3,1./(Tn*Tr),AA,3,F_total,3,0.0,S,3);

    sig[ii].el.o[0] += wt*J/volume*S[idx_2(0,0)];
    sig[ii].el.o[1] += wt*J/volume*S[idx_2(1,1)];
    sig[ii].el.o[2] += wt*J/volume*S[idx_2(2,2)];

    sig[ii].el.o[3] += wt*J/volume*S[idx_2(1,2)];
    sig[ii].el.o[4] += wt*J/volume*S[idx_2(0,2)];
    sig[ii].el.o[5] += wt*J/volume*S[idx_2(0,1)];

    /* Compute logarithmic strain */
    inverse(F_total,3,F_total_I);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                3,3,3,1.0,F_total_I,3,E,3,0.0,AA,3);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                3,3,3,1.0,AA,3,F_total_I,3,0.0,E,3);

    eps[ii].el.o[0] += wt*J/volume*E[idx_2(0,0)];
    eps[ii].el.o[1] += wt*J/volume*E[idx_2(1,1)];
    eps[ii].el.o[2] += wt*J/volume*E[idx_2(2,2)];

    eps[ii].el.o[3] += 2.*wt*J/volume*E[idx_2(1,2)];
    eps[ii].el.o[4] += 2.*wt*J/volume*E[idx_2(0,2)];
    eps[ii].el.o[5] += 2.*wt*J/volume*E[idx_2(0,1)];

    /* Store COMPUTATIONAL Fn at the ip */
    alpha = 1.0;
    def_grad_get(nne,ndn,CONST_4(double)ST_tensor,disp,Fr_mat);
    mat2array(Fr,CONST_2(double) Fr_mat,3,3);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                3,3,3,alpha,Fr,3,Fn,3,0.0,F_total,3);
    if(ip == 0){
      memcpy(eps[ii].il[0].F,F_total,9*sizeof(double));
    }
  } /* integration */

  /* update pressure and volume */
  for (int i=0; i<nPres; i++){
    sig[ii].p[i] += p[i];
  }
  for (int i=0; i<nVol; i++){
    eps[ii].T[i] *= eps[ii].d_T[i];
    eps[ii].d_T[i] = 1.0;
  }

  /* update bubble dofs */
  for(int i=0; i<elem[ii].n_bub*elem[ii].n_bub_dofs; i++){
    elem[ii].bub_dofs[i] += elem[ii].d_bub_dofs[i];
    elem[ii].d_bub_dofs[i] = 0;
  }

  free(disp);
  free(p);
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,3,3,ndn);
  free(ST);
  free(Fn);
  dealoc2(Fr_mat,3);
  free(Fr);
  free(Fr_lin);
  free(Fr_I);
  free(AA);
  free(F_total);
  free(F_total_I);
  free(C);
  free(C_I);
  free(S);
  free(E);
  free(ident);
}/* MINI_3f_increment_el */


void MINI_3f_increment(Element *elem,
                       const int nelem,
                       Node *node,
                       const int nnodes,
                       const int ndofn,
                       const SUPP sup,
                       EPS *eps,
                       SIG *sig,
                       const HOMMAT *hommat,
                       const double *sol,
		       const CommunicationStructure *com,
                       const int mp_id)
{
  const int ndn = 3;
  int nne, nne_t, II;

  long *nod, *cn;
  double *x, *y, *z;
  double *sol_e;

  /* for each element */
  for (int i=0; i<nelem; i++){
    nne = elem[i].toe;
    nne_t = nne + elem[i].n_bub;

    // @todo Commented as dead code. @cp please review. LD
    // int mat = elem[i].mat[2];

    /* allocate */
    nod = aloc1l(nne);
    cn = aloc1l(nne*ndofn);
    x = aloc1(nne_t);
    y = aloc1(nne_t);
    z = aloc1(nne_t);
    sol_e = aloc1(nne*ndofn);

    /* get node ids on element */
    elemnodes(i,nne,nod,elem);

    /* get local dof ids on elemnt */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* Get the node and bubble coordinates */
    nodecoord_updated(nne,nod,node,x,y,z);
    element_center(nne,x,y,z);

    def_elem(cn,nne*ndofn,sol,elem,node,sol_e,sup,0);

    /* increment the element quantities */
    MINI_3f_increment_el(elem,i,nne,node,nod,ndofn,
                         x,y,z,eps,sig,hommat,sol_e);

    free(nod);
    free(cn);
    free(x);
    free(y);
    free(z);
    free(sol_e);
  } /* For each element */

  /*** Update Coordinates ***/
  for (int ii=0;ii<nnodes;ii++){
    for (int i=0;i<ndn;i++){
      II = node[ii].id_map[mp_id].id[i];
      if (II > 0){
        if (i == 0) node[ii].x1 += sol[II-1];
        if (i == 1) node[ii].x2 += sol[II-1];
        if (i == 2) node[ii].x3 += sol[II-1];
      }
      if (II < 0){
        if (i == 0) node[ii].x1 += sup->defl_d[abs(II)-1];
        if (i == 1) node[ii].x2 += sup->defl_d[abs(II)-1];
        if (i == 2) node[ii].x3 += sup->defl_d[abs(II)-1];
      }
    }
  }/* end ii < nn */

  int myrank = com->rank;
  double PL, GPL;
  PL = T_VOLUME (nelem,ndn,elem,node);
  /* Gather Volume from all domains */
  com->net->allreduce(&PL,&GPL,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  if (myrank == 0) {
    PGFEM_printf ("AFTER DEF - VOLUME = %12.12f\n",GPL);
  }
}/* MINI_3f_increment */

void MINI_3f_check_resid(const int ndofn,
                         const int ne,
                         const Element *elem,
                         const Node *node,
                         const HOMMAT *hommat,
                         const EPS *eps,
                         const SIG *sig,
                         const double *d_r,
                         const SUPP sup,
                         const double *RR,
                         const long *DomDof,
                         const int ndofd,
                         const int GDof,
			 const CommunicationStructure *com,
                         const int mp_id)
{
  /* compute the norm of the residauls for each variable */
  int myrank = com->rank;
  const int ndn = 3;
  const int nVol = 4;

  int err = 0;

  int count; /* ALWAYS reset before use */
  long *nod, *cn;
  double *r_e, *x, *y, *z, *fe, *rue1, *rue2;
  double nrT, nrB, nrB1, nrB2, nrBd;
  nrT = nrB = nrB1 = nrB2 = nrBd = 0.0;

  double *f, *f_u, *f_u1, *f_u2, *f_p;
  f = aloc1(ndofd);
  f_u = aloc1(ndofd);
  f_u1 = aloc1(ndofd);
  f_u2 = aloc1(ndofd);
  f_p = aloc1(ndofd);

  FILE *debug_log;

  for(int ii=0; ii<ne; ii++){

    /* Number of element nodes */
    const int nne = elem[ii].toe;
    const int nPres = nne;
    const int total_nne = nne + elem[ii].n_bub;
    /* Nodes on element */
    nod = aloc1l (nne);
    elemnodes (ii,nne,nod,elem);
    /* Element Dof */
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

    const int n_bub = elem[ii].n_bub;
    const int n_bub_dofs = elem[ii].n_bub_dofs;
    const int mat = elem[ii].mat[2];
    const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

    /* allocation */
    r_e = aloc1 (ndofe); /* container for accum solution */
    x = aloc1 (total_nne);
    y = aloc1 (total_nne);
    z = aloc1 (total_nne);
    fe = aloc1 (ndofe); /* container for element global res */
    rue1 = aloc1 (ndofe);
    rue2 = aloc1 (ndofe);
    cn = aloc1l (ndofe);

    /* coordinates */
    nodecoord_updated (nne,nod,node,x,y,z);
    element_center(nne,x,y,z);

    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* deformation on element */
    def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);

    /*============== COMPUTE ELEMENT RESIDUALS ================*/
    double *disp, *p;
    disp = aloc1(total_nne*ndn);
    p = aloc1(nPres);

    count = 0;
    for (int i=0; i<nne; i++){
      for (int j=0; j<ndofn; j++){

        /* filter displacement and pressure from element unknowns */
        if (j<ndn){                           /* displacement dof */
          disp[i*ndn + j] = r_e[count+j];
        } else if (j==ndn){                       /* pressure dof */
          p[i] = r_e[count+j];
        } else {                                         /* ERROR */
          PGFEM_printerr("Too many DOFs on node in %s!\n",__func__);
          PGFEM_Abort();
        }
      } /* ndofn */
      count += ndofn;
    } /* nne */

    /* displacements from bubble */
    for (int i=0; i<n_bub; i++){
      for (int j=0; j<n_bub_dofs; j++){
        disp[nne*ndn + i*n_bub_dofs + j] =
        elem[ii].bub_dofs[i*n_bub_dofs + j];
      }
    }

    if(MINI_3F_DEBUG){
      char fname[50];
      double vol = Tetra_V(x,y,z);
      sprintf(fname,"MINI_3f_check_%d_%d.log",ii,myrank);
      debug_log = fopen(fname,"a");
      PGFEM_fprintf(debug_log,"***********************************\n");
      PGFEM_fprintf(debug_log,"Volume = %25.12e\n",vol);
      PGFEM_fprintf(debug_log,"Displacements:\n");
      print_array_d(debug_log,disp,total_nne*ndn,1,total_nne*ndn);
      PGFEM_fprintf(debug_log,"Pressures:\n");
      print_array_d(debug_log,p,nPres,1,nPres);
      PGFEM_fprintf(debug_log,"Volume:\n");
      print_array_d(debug_log,eps[ii].d_T,nVol,1,nVol);
    }

    /*=== BEGIN INTEGRATION ===*/
    long npt_x, npt_y, npt_z;
    int_point(total_nne,&npt_z);

    double *int_pt_ksi = PGFEM_calloc (double, npt_z);
    double *int_pt_eta = PGFEM_calloc (double, npt_z);
    double *int_pt_zet = PGFEM_calloc (double, npt_z);
    double *weights = PGFEM_calloc (double, npt_z);

    /* allocate space for the shape functions, derivatives etc */
    double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
    Na = aloc1 (total_nne);
    Np = aloc1 (nPres);
    Nt = aloc1 (nVol);
    N_x = aloc1 (total_nne);
    N_y = aloc1 (total_nne);
    N_z = aloc1 (total_nne);
    ST_tensor = aloc4 (3,3,ndn,total_nne);
    ST = aloc1(3*3*ndn*total_nne); /* index space i,j,node,dof */

    /* allocate space for deformation gradients etc */
    double *Fn, **Fr_mat, *Fr, *Fr_lin, *Fr_I;
    double *C, *AA;
    double *S;
    Fn = aloc1(9);
    Fr_mat = aloc2(3,3);
    Fr = aloc1(9);
    Fr_lin = aloc1(9);
    Fr_I = aloc1(9);
    AA = aloc1(9);
    C = aloc1(9);
    S = aloc1(9);

    /* Allocate space for the residuals */
    double *Ru,*Ru1,*Ru2, *dRu, *Rb,*Rb1,*Rb2, *dRb, *Rp, *Rt;
    Ru = aloc1(nne*ndn);
    Ru1 = aloc1(nne*ndn);
    Ru2 = aloc1(nne*ndn);
    dRu = aloc1(nne*ndn);
    Rb = aloc1(n_bub*n_bub_dofs);
    Rb1 = aloc1(n_bub*n_bub_dofs);
    Rb2 = aloc1(n_bub*n_bub_dofs);
    dRb = aloc1(n_bub*n_bub_dofs);
    Rp = aloc1(nne);
    Rt = aloc1(nne);

    /* integration point in natural coords */
    double wt;

    /* pressure */
    double pressure, Up, Tn, Tr;

    memcpy(Fn,eps[ii].il[0].F,9*sizeof(double));
    const double Jn = getJacobian(Fn,ii,&err);

    /* get the integration points and weights */
    integrate (total_nne,&npt_x,&npt_y,&npt_z,
               int_pt_ksi,int_pt_eta,int_pt_zet,weights);
    for (int ip=0; ip<npt_z; ip++){
      integration_help(ip,ii,nne,total_nne,ndn,nPres,p,nVol,
                       int_pt_ksi,int_pt_eta,int_pt_zet,weights,
                       x,y,z,sig,eps,&wt,&J,&pressure,&Tn,&Tr,Na,
                       N_x,N_y,N_z,Np,Nt,ST_tensor,ST);

      /* get the linear portion of Fr */
      if(ip == 0){
        def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
        mat2array(Fr_lin,CONST_2(double) Fr_mat,3,3);
      }

      /* add the bubble portion of Fr */
      compute_disp_grad(total_nne,ST,disp,Fr,nne);
      cblas_daxpy(ndn*ndn,1.0,Fr_lin,1,Fr,1);

      double Jr = getJacobian(Fr,ii,&err);
      inverse(Fr,3,Fr_I);


      /* Compute total F; F=FrFn */
      double alpha = 1.0;
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                  3,3,3,alpha,Fr,3,Fn,3,0.0,AA,3);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                  3,3,3,1.0,AA,3,AA,3,0.0,C,3);

      // Get Deviatoric 2 P-K stress
      devStressFuncPtr Stress;
      dUdJFuncPtr UP;
      Stress = getDevStressFunc(1,&hommat[mat]);
      UP = getDUdJFunc(1,&hommat[mat]);

      Stress(C,&hommat[mat],S);
      UP(Tn*Tr,&hommat[mat],&Up);

      /* Compute residuals */
      /* Ru (12 x 1) */
      HW_Ru_at_ip(Ru,nne,total_nne,ST,Fn,Fr,Fr_I,
                  Jn,Tn,Jr,Tr,S,pressure,J,wt,0);
      HW_Ru1_at_ip(Ru1,nne,total_nne,ST,Fn,Fr,Fr_I,
                   Jn,Tn,Jr,Tr,S,pressure,J,wt,0);
      HW_Ru2_at_ip(Ru2,nne,total_nne,ST,Fn,Fr,Fr_I,
                   Jn,Tn,Jr,S,pressure,J,wt,0);
      debug_NH_Ru_at_ip(dRu,nne,total_nne,ST,Fn,Fr,Fr_I,
                        Jn,Tn,Jr,Tr,S,pressure,hommat[mat].G,J,wt,0);

      /* Rb (3 x 1) */
      HW_Ru_at_ip(Rb,nne,total_nne,ST,Fn,Fr,Fr_I,
                  Jn,Tn,Jr,Tr,S,pressure,J,wt,1);
      HW_Ru1_at_ip(Rb1,nne,total_nne,ST,Fn,Fr,Fr_I,
                   Jn,Tn,Jr,Tr,S,pressure,J,wt,1);
      HW_Ru2_at_ip(Rb2,nne,total_nne,ST,Fn,Fr,Fr_I,
                   Jn,Tn,Jr,S,pressure,J,wt,1);
      debug_NH_Ru_at_ip(dRb,nne,total_nne,ST,Fn,Fr,Fr_I,
                        Jn,Tn,Jr,Tr,S,pressure,hommat[mat].G,J,wt,1);

      /* Rp (4 x 1) */
      HW_Rp_at_ip(Rp,nPres,Np,Jn,Jr,Tn,Tr,J,wt);

      /* Rt (4 x 1) */
      HW_Rt_at_ip(Rt,nVol,Nt,Tn,Jn,pressure,kappa,Up,J,wt);

    }/* end integration */

    /* free un-needed */
    free(int_pt_ksi);
    free(int_pt_eta);
    free(int_pt_zet);
    free(weights);
    free(Na);
    free(Np);
    free(Nt);
    free(N_x);
    free(N_y);
    free(N_z);
    dealoc4(ST_tensor,3,3,ndn);
    free(ST);
    free(Fn);
    dealoc2(Fr_mat,3);
    free(Fr);
    free(Fr_lin);
    free(Fr_I);
    free(AA);
    free(C);
    free(S);
    /*=== END INTEGRATION ===*/

    int sign = 1;
    /*** Scatter the GLOBAL residuals (u,v,w,p,t) ***/
    /* NOTE residuals are positive because subtracted from RHS later */
    for (int a=0; a<nne; a++){
      for (int b=0; b<ndofn; b++){
        if (b<ndn){ /* from Ru */
          fe[a*ndofn+b] = sign*Ru[a*ndn+b];
          rue1[a*ndofn+b] = sign*dRu[a*ndn+b];
          rue2[a*ndofn+b] = sign*Ru2[a*ndn+b];
        } else if (b == ndn){ /* from Rp */
          fe[a*ndofn+b] = sign*Rp[a];
        }
      }
    }

    if(MINI_3F_DEBUG){
      PGFEM_fprintf(debug_log,"Ru:\n");
      print_array_d(debug_log,Ru,nne*ndn,1,ndn*nne);
      PGFEM_fprintf(debug_log,"Ru1:\n");
      print_array_d(debug_log,Ru1,nne*ndn,1,ndn*nne);
      PGFEM_fprintf(debug_log,"Ru2:\n");
      print_array_d(debug_log,Ru2,nne*ndn,1,ndn*nne);
      PGFEM_fprintf(debug_log,"Rb:\n");
      print_array_d(debug_log,Rb,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
      PGFEM_fprintf(debug_log,"Rb1:\n");
      print_array_d(debug_log,Rb1,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
      PGFEM_fprintf(debug_log,"Rb2:\n");
      print_array_d(debug_log,Rb2,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
      PGFEM_fprintf(debug_log,"dRb:\n");
      print_array_d(debug_log,dRb,n_bub*n_bub_dofs,1,n_bub*n_bub_dofs);
      PGFEM_fprintf(debug_log,"Rt:\n");
      print_array_d(debug_log,Rt,nVol,1,nVol);
      PGFEM_fprintf(debug_log,"Rp:\n");
      print_array_d(debug_log,Rp,nPres,1,nPres);
      fclose(debug_log);
    }

    /* Compute norm of LOCAL residuals */
    nrT += cblas_ddot(nVol,Rt,1,Rt,1);
    nrB += cblas_ddot(n_bub*n_bub_dofs,Rb,1,Rb,1);
    nrB1 += cblas_ddot(n_bub*n_bub_dofs,Rb1,1,Rb1,1);
    nrB2 += cblas_ddot(n_bub*n_bub_dofs,Rb2,1,Rb2,1);
    nrBd += cblas_ddot(n_bub*n_bub_dofs,dRb,1,dRb,1);

    free(Ru);
    free(Ru1);
    free(Ru2);
    free(dRu);
    free(Rb);
    free(Rb1);
    free(Rb2);
    free(dRb);
    free(Rt);
    free(Rp);
    /*============== END COMPUTE ELEMENT RESIDUALS ============*/

    /* Assemble global residual for U and P, but separate into
       separate arrays */
    int j = 0;
    for (int k=0;k<nne;k++){
      for (int kk=0;kk<ndofn;kk++){
        int II = node[nod[k]].id_map[mp_id].id[kk]-1;
        if (II < 0) continue;
        if (kk<ndofn-1){
          f_u[II] += fe[j+kk];
          f_u1[II] += rue1[j+kk];
          f_u2[II] += rue2[j+kk];
        } else {
          f_p[II] += fe[j+kk];
        }
      }/*end kk*/
      j += ndofn;
    }/*end k*/

    free(r_e);
    free(x);
    free(y);
    free(z);
    free(fe);
    free(rue1);
    free(rue2);
    free(cn);

  }/* end ii < ne */

  /* Compute norms */
  double normal, tmp;
  double *BS_f;
  BS_f = aloc1(DomDof[myrank]);

  /*=== Norm of Ru ===*/
  for (int i=0;i<ndofd;i++) f[i] = RR[i] - f_u[i];
  LToG(f,BS_f,ndofd,com);
  normal = cblas_ddot(DomDof[myrank],BS_f,1,BS_f,1);
  com->net->allreduce(&normal,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Ru = %12.12e\n",normal);
  FILE *BS_f_out;
  char fname[50];
  sprintf(fname,"Ru_%d.out",myrank);
  BS_f_out = fopen(fname,"a");
  print_array_d(BS_f_out,BS_f,DomDof[myrank],1,DomDof[myrank]);
  fclose(BS_f_out);
  /* memset(BS_f,0,DomDof[myrank]*sizeof(double)); */

  /*=== Norm of dRu ===*/
  for (int i=0;i<ndofd;i++) f[i] = RR[i] - f_u1[i];
  LToG(f,BS_f,ndofd,com);
  normal = cblas_ddot(DomDof[myrank],BS_f,1,BS_f,1);
  com->net->allreduce(&normal,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM dRu = %12.12e\n",normal);
  /* memset(BS_f,0,DomDof[myrank]*sizeof(double)); */

  /*=== Norm of Rp ===*/
  LToG(f_p,BS_f,ndofd,com);
  normal = cblas_ddot(DomDof[myrank],BS_f,1,BS_f,1);
  com->net->allreduce(&normal,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Rp = %12.12e\n",normal);

  /*=== Norm of Rb ===*/
  com->net->allreduce(&nrB,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Rb = %12.12e\n",normal);

  com->net->allreduce(&nrBd,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM dRb = %12.12e\n",normal);

  com->net->allreduce(&nrB1,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Rb1 = %12.12e\n",normal);

  com->net->allreduce(&nrB2,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Rb2 = %12.12e\n",normal);

  /*=== Norm of Rt ===*/
  com->net->allreduce(&nrT,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  normal = sqrt (tmp);
  if (myrank == 0) PGFEM_printf("NORM Rt = %12.12e\n",normal);

  free(f_u);
  free(f_u1);
  free(f_u2);
  free(f_p);
  free(f);

}/* MINI_3f_check_resid */
