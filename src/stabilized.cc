/* HEADER */
/**
 * Authors:
 * Karel Matous
 * Matthew Mosby
 */
#include "stabilized.h"
#include <time.h>
#include <string.h>
#include <math.h>
#include "mkl_cblas.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "cast_macros.h"
#include "def_grad.h"
#include "allocation.h"
#include "tensors.h"
#include "new_potentials.h"
#include "utils.h"
#include "index_macros.h"
#include "two_field_element.h"
#include "elem3d.h"

#ifndef STAB_DEBUG
#define STAB_DEBUG 0
#endif

#define DELTA(X,Y) ((double)((X)==(Y)))

static const int periodic = 0;
static const double PII = 3.141592653589793238462643;
static const double identity[9]={1,0,0,0,1,0,0,0,1};

/*=== Declare the static helper functions ===*/
/** Compute common quantities used in integration */
static int integration_help(const int ip,
                const int elem_id,
                const int nne,
                const int i,
                const int j,
                const int k,
                const double *x,
                const double *y,
                const double *z,
                const double *int_pt_ksi,
                const double *int_pt_eta,
                const double *int_pt_zet,
                const double *weights,
                const double *disp,
                const double *pn_1,
                const double *pn,
                const double *p,
                const double *Fn,
                const SUPP sup,
                double *wt,
                double *jj,
                double *Na,
                double *N_x,
                double *N_y,
                double *N_z,
                double *Np_x,
                double *ST,
                double *Fr,
                double *F,
                double *C,
                double *Fr_I,
                double *Jr,
                double *Pn_1,
                double *Pn,
                double *pressure,
                double *grad_delP,
                double *grad_P);

/** Compute common quantities used in pressure term integration */
static int quadradic_integration_help(const int ip,
                      const int elem_id,
                      const int nne,
                      const int i,
                      const int j,
                      const int k,
                      const double *x,
                      const double *y,
                      const double *z,
                      const double *int_pt_ksi,
                      const double *int_pt_eta,
                      const double *int_pt_zet,
                      const double *weights,
                      const double *disp,
                      const double *pn,
                      const double *p,
                      const SUPP sup,
                      double *wt,
                      double *jj,
                      double *Na,
                      double *N_x,
                      double *N_y,
                      double *N_z,
                      double *ST,
                      double *Fr,
                      double *Jr,
                      double *delP,
                      double *Pn);

/** Compute the material deviatoric stress at an integration point */
static int get_material_stress(double *devSbar,
                   const double *C,
                   const HOMMAT *mat);

/** Compute the material deviatoric tangent (dSdFr) at an integration
    point NOTE: This is not purely the deviatoric part if there is
    damage */
static int get_material_tangent(double *L,
                double *SS,
                double *Ybar,
                const double kappa,
                const double Un_1, /* U(n-1) */
                const double Jn_1, /* J(n-1) */
                const double Jn,
                const double Jr,
                const double *Fn,
                const double *Fr,
                const double *C,
                const double *C_I,
                const double Pn_1, /* p(n-1) */
                const double Pn,
                const double pressure,
                const double *devS,
                const damage *dam,
                const HOMMAT *mat);

/** Compute the pressure at an integration point */
static int get_material_pres(double *Up,
                 const double Jn,
                 const double Jr,
                 const HOMMAT *mat,
                 const damage *dam);

/** Compute the linearized pressure term at an integration point */
static int get_material_lin_pres(double *Upp,
                 const double Jn,
                 const double Jr,
                 const HOMMAT *mat);

/** Compute the displacement residual at an integration point */
static int get_Ru_at_ip(double *Ru,
            const int nne,
            const double *ST,
            const double *Fn,
            const double *Fr,
            const double *Fr_I,
            const double Jr,
            const double Jn,
            const double *devSbar,
            const double pres,
            const damage *dam,
            const double jj,
            const double wt);

/** Compute the typical pressure residual at an integration
    point. NOTE: This formulation requires quadradic integration on
    this term */
static int get_Rp_at_ip(double *Rp,
            const int TL,
            const int nne,
            const double *Np,
            const double kappa,
            const double deltaUp,
            const double deltaPres,
            const double P,
            const damage *dam,
            const double jj,
            const double wt);

/** Compute the stabilization term for the residual */
static int get_Rp_stab_at_ip(double *Rp,
                 const int TL,
                 const int nne,
                 const double *Np_x,/* Np_x [ndn x nne] */
                 const double *grad_delP,
                 const double *grad_P, /* n+1 */
                 const double Jr,
                 const double *Fr_I,
                 const double stab,
                 const damage *dam,
                 const double jj,
                 const double wt);

/** Compute Kuu at an integration point */
static int get_Kuu_at_ip(double *Kuu,
             const int nne,
             const double *ST,
             const double *Fn,
             const double *Fr,
             const double *Fr_I,
             const double Jn,
             const double Jr,
             const damage *dam,
             const double totalP,
             const double *devSbar,
             const double *L,
             const double jj,
             const double wt);

/** Compute Kup at an integration point */
static int get_Kup_at_ip(double *Kup,
             const int nne,
             const double *Np,
             const double *ST,
             const double *Fr_I,
             const double Ybar,
             const double Jn_1,
             const double Jn,
             const double Jr,
             const damage *dam,
             const double jj,
             const double wt);

/** Compute Kpu at an integration point */
static int get_Kpu_at_ip(double *Kpu,
             const int nne,
             const double *Np,
             const double *ST,
             const double *Fr_I,
             const double Jn,
             const double Jr,
             const double kappa,
             const double Upp,
             const double UP, /* n+1 only */
             const double P, /* total pressure */
             const double *SS,
             const damage *dam,
             const double jj,
             const double wt);

/** Compute Kpp at an integration point. NOTE: This formulation
    requires quadratic integration for this term. */
static int get_Kpp_at_ip(double *Kpp,
             const int nne,
             const double *Np,
             const double kappa,
             const double UP, /* n+1 only */
             const double P, /* total pressure */
             const double Jn_1,
             const double Jn,
             const double Jr,
             const double H,
             const damage *dam,
             const double jj,
             const double wt);

/** Compute the stabilization term for the Kpu tangent */
static int get_Kpu_stab_at_ip(double *Kpu,
                  const int TL,
                  const int nne,
                  const double *ST,
                  const double *Np_x,
                  const double *Fr,
                  const double *Fr_I,
                  const double Jr,
                  const double *grad_delP,
                  const double *grad_P, /*n+1*/
                  const double stab,
                  const damage *dam,
                  const double *SS,
                  const double jj,
                  const double wt);

/** Compute the stabilization term for the Kpp tangent */
static int get_Kpp_stab_at_ip(double *Kpp,
                  const int nne,
                  const double *Np,
                  const double *Np_x,
                  const double *Fr_I,
                  const double Jn_1,
                  const double Jn,
                  const double Jr,
                  const double stab,
                  const double *grad_P, /*n+1*/
                  const double Ybar,
                  const damage *dam,
                  const double jj,
                  const double wt);

/* increment an element */
static double st_incr_elem (const long ii,
                const long nne,
                const long ndofn,
                const double *x,
                const double *y,
                const double *z,
                const double *r_u,
                const double *a_a,
                const double *p,
                const double dt,
                const HOMMAT *hommat,
                const ELEMENT *elem,
                EPS *eps,
                SIG *sig,
                const SUPP sup);

/*=======================================*/
/*   Main API function definitions       */
/*=======================================*/

int stab_get_material_potential(double *Wbar,
                const double kappa,
                const double Un_1, /* U(n-1) */
                const double Jn_1, /* J(n-1) */
                const double Jnp1, /* J(n+1) */
                const double Pn_1, /* p(n-1) */
                const double Pn,
                const double P,    /* p(n+1) */
                const double *C,
                const HOMMAT *mat)
{
  int err = 0;

  /* Integrate volumetric energy by Simpson's rule */
  /* double U = Un_1 + (Jnp1-Jn_1)/(6.*kappa)*(P + 4.*Pn + Pn_1); */

  /* Use Volumetric potential from J */
  UFuncPtr VolPot = getUFunc(0,mat);
  double U = 0.0;
  VolPot(Jnp1,mat,&U);

  devPotentialFuncPtr dev_potential = getDevPotentialFunc(0,mat);
  *Wbar = 0.0;
  dev_potential(C,mat,Wbar);

  /* *Wbar += kappa*U; */
  *Wbar = (*Wbar + kappa*U);

  return err;
}/* static int get_material_potential() */

void res_stab_def (long ne,
           long npres,
           ELEMENT *elem,
           EPS *eps,
           SIG *sig,
           double stab)
{
  long ii,ip,i,nne,II,mat;
  for (ii=0;ii<ne;ii++){

    nne = elem[ii].toe;
    mat = elem[ii].mat[2];

    /* Integration */
    int_point (nne,&II);

    if (periodic == 1){
      for (ip=0;ip<II;ip++){
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
    }

    for (i=0;i<npres;i++){
      sig[ii].d_p[i] = 0.0;
    }
  }/* end ii<ne */
}

int resid_st_elem (long ii,
           long ndofn,
           long nne,
           ELEMENT *elem,
           long *nod,
           NODE *node,
           HOMMAT *hommat,
           double *x,
           double *y,
           double *z,
           EPS *eps,
           SIG *sig,
           SUPP sup,
           double *r_e,
           double nor_min,
           double *fe,
           double dt,
           double stab)
{
  /* compute constants */
  static const int ndn = 3;
  const int mat = elem[ii].mat[2];
  const double  kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  const HOMMAT *ptrMat = &hommat[mat];


  const double Vol = Tetra_V (x,y,z);
  const double he = 2.*pow(3.*Vol/4./PII,1./3.)/1.7320508;
  const double SStab = stab*he*he/(2.*hommat[mat].G);

  int err = 0;

  /*=== ALLOCATION/declaration ===*/
  /* integration */
  long npt_x, npt_y, npt_z;
  /* alloce enough space for quadradic integration */
  int_point(10,&npt_z);
  double *int_pt_ksi = aloc1(npt_z);
  double *int_pt_eta = aloc1(npt_z);
  double *int_pt_zet = aloc1(npt_z);
  double *weights = aloc1(npt_z);
  double *Na = aloc1(nne);
  double *N_x = aloc1(nne);
  double *N_y = aloc1(nne);
  double *N_z = aloc1(nne);
  double *Np_x = aloc1(3*nne);
  double *ST = aloc1(3*3*ndn*nne);

  /* kinematics */
  /* Fn will be a temporary pointer within the integration loop and
     will not need allocation/dealocation */
  double *Fr = aloc1(9);
  double *Fr_I = aloc1(9);
  double *F = aloc1(9);
  double *C = aloc1(9);

  /* physical quantities */
  double *disp = aloc1(ndn*nne);
  double *p = aloc1(nne);
  double *grad_delP = aloc1(3);
  double *grad_P = aloc1(3);
  double *devSbar = aloc1(9);

  /* residuals */
  double *Ru = aloc1(ndn*nne);
  double *Rp = aloc1(nne);


  /*=== get nodal values ===*/
  {
    int k = 0;
    for(int i=0; i<nne; i++){
      for(int j=0; j<ndofn; j++){
    if(j<ndn){
      disp[i*ndn+j] = r_e[k+j];
    } else if(j==ndn){
      p[i] = r_e[k+j];
    }
      }
      k += ndofn;
    }
  }

  /*=== main integration loop ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
        int_pt_ksi,int_pt_eta,int_pt_zet,
        weights);
  int ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){

    /* declare new variables for integration point */
    double Jr,Jn,pressure,wt,jj,Pn_1,Pn;

    /* Get pointer to Fn and compute Jn */
    /* DO NOT free Fn!!! */
    const double *Fn = NULL;
    if (sup->multi_scale){
      Fn = identity;
      Jn = 1.0;
    } else {
      Fn = eps[ii].il[ip].F;
      Jn = getJacobian(Fn,ii,&err);
    }

    /* get pointer to damage object */
    const damage *ptrDam = &eps[ii].dam[ip];

    /* compute various objects needed for integration */
    err += integration_help(ip,ii,nne,i,j,k,x,y,z,
                int_pt_ksi,int_pt_eta,int_pt_zet,
                weights,disp,sig[ii].pn_1,sig[ii].p,p,
                Fn,sup,&wt,&jj,Na,N_x,N_y,N_z,Np_x,ST,Fr,F,
                C,Fr_I,&Jr,&Pn_1,&Pn,&pressure,grad_delP,
                grad_P);

    /* compute the material stress */
    err += get_material_stress(devSbar,C,ptrMat);

    /* compute Ru contribution at integration point */
    err += get_Ru_at_ip(Ru,nne,ST,Fn,Fr,Fr_I,Jr,Jn,devSbar,
                /*n+1*/ pressure,ptrDam,jj,wt);

    /* compute Rp stabilization contribution at integration point */
    err += get_Rp_stab_at_ip(Rp,sup->multi_scale,nne,Np_x,grad_delP,grad_P,Jr,
                 Fr_I,SStab,ptrDam,jj,wt);

    /* increment ip counter */
    ip++;

    /* error detected, deallocate and exit */
    if(err != 0) goto dealocate;
      }
    }
  }

  /*=== quadradic integration ===*/
  integrate(10,&npt_x,&npt_y,&npt_z,
        int_pt_ksi,int_pt_eta,int_pt_zet,
        weights);
  ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){

    /* declare fresh variables for integration */
    double Jr,Jn,delP,Pn,Up,wt,jj;

    /* get pointer to Fn and compute Jn */
    /* DO NOT free Fn!!! */
    const double *Fn = NULL;
    if (sup->multi_scale){
      Fn = identity;
      Jn = 1.0;
    } else {
      Fn = eps[ii].st[ip].Fpp;
      Jn = getJacobian(Fn,ii,&err);
    }

    const damage *ptrDam = &eps[ii].st[ip].dam;

    /* compute various quantities for integration */
    err += quadradic_integration_help(ip,ii,nne,i,j,k,x,y,z,
                      int_pt_ksi,int_pt_eta,int_pt_zet,
                      weights,disp,sig[ii].p,p,sup,&wt,&jj,
                      Na,N_x,N_y,N_z,ST,Fr,&Jr,&delP,&Pn);

    /* compute Up = (1-w)Up(n+1) - (1-wn)Upn(n) */
    err += get_material_pres(&Up,Jn,Jr,ptrMat,ptrDam);

    /* compute contribution to Rp from integration point */
    err += get_Rp_at_ip(Rp,sup->multi_scale,nne,Na,kappa,Up,delP,Pn+delP,ptrDam,jj,wt);

    /* increment ip counter */
    ip++;

    /* error detected, deallocate and exit */
    if(err != 0) goto dealocate;
      }
    }
  }

 dealocate:
  /*=== DEALLOCATION ===*/
  /* integration */
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  free(Np_x);
  free(ST);

  /* kinematics */
  free(Fr);
  free(Fr_I);
  free(F);
  free(C);

  /* if there is an error, exit; do not assemble */
  if(err != 0) goto exit_function;

  /*=== Residual Assembly ===*/
  {
    int k = 0;
    for (int i=0;i<nne;i++){
      for (int j=0;j<ndofn;j++){
    if (j  < ndn){
      fe[k+j] = Ru[i*ndn+j];
    } else if (j == ndn){
      fe[k+j] = Rp[i];
    } else {
      fe[k+j] = 0.0;
    }
      }
      k += ndofn;
    }
  }

  if(STAB_DEBUG){
    int err_rank = 0;
    char filename[50];
    FILE *debug_log;
    PGFEM_Error_rank(&err_rank);
    sprintf(filename,"new_stab_resid_%ld_%d.log",ii,err_rank);
    debug_log = fopen(filename,"a");
    PGFEM_fprintf(debug_log,
        "=========================================\n"
        "=========== New Computation =============\n"
        "=========================================\n");
    PGFEM_fprintf(debug_log,"Disp\n");
    print_array_d(debug_log,disp,nne*ndn,1,nne*ndn);
    PGFEM_fprintf(debug_log,"Pressure\n");
    print_array_d(debug_log,p,nne,1,nne);
    PGFEM_fprintf(debug_log,"Ru\n");
    print_array_d(debug_log,Ru,nne*ndn,1,nne*ndn);
    PGFEM_fprintf(debug_log,"Rp\n");
    print_array_d(debug_log,Rp,nne,1,nne);
    PGFEM_fprintf(debug_log,"Assembled\n");
    print_array_d(debug_log,fe,nne*ndofn,1,nne*ndofn);
    fclose(debug_log);
  }

 exit_function:
  /* physical quantities */
  free(disp);
  free(p);
  free(grad_delP);
  free(grad_P);
  free(devSbar);

  /* deallocate residuals */
  free(Ru);
  free(Rp);

  if(STAB_DEBUG && err != 0){
    PGFEM_printf("error in %s [elem %ld]\n",__func__,ii);
  }

  return err;
}/* int resid_st_elem () */

int stiffmatel_st (long ii,
           long ndofn,
           long nne,
           double *x,
           double *y,
           double *z,
           ELEMENT *elem,
           HOMMAT *hommat,
           long *nod,
           NODE *node,
           SIG *sig,
           EPS *eps,
           SUPP sup,
           double *r_e,
           long npres,
           double nor_min,
           double *Ks,
           double dt,
           double stab,
           long FNR,
           double lm,
           double *fe)
{
  /* compute constants */
  static const int ndn = 3;
  const int mat = elem[ii].mat[2];
  const double  kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  const HOMMAT *ptrMat = &hommat[mat];


  const double Vol = Tetra_V (x,y,z);
  const double he = 2.*pow(3.*Vol/4./PII,1./3.)/1.7320508;
  const double SStab = stab*he*he/(2.*hommat[mat].G);

  int err = 0;

  /*=== ALLOCATION/declaration ===*/
  /* integration */
  long npt_x, npt_y, npt_z;
  /* alloce enough space for quadradic integration */
  int_point(10,&npt_z);
  double *int_pt_ksi = aloc1(npt_z);
  double *int_pt_eta = aloc1(npt_z);
  double *int_pt_zet = aloc1(npt_z);
  double *weights = aloc1(npt_z);
  double *Na = aloc1(nne);
  double *N_x = aloc1(nne);
  double *N_y = aloc1(nne);
  double *N_z = aloc1(nne);
  double *Np_x = aloc1(3*nne);
  double *ST = aloc1(3*3*ndn*nne);

  /* kinematics */
  /* Fn will be a temporary pointer within the integration loop and
     will not need allocation/deallocation */
  double *Fr = aloc1(9);
  double *Fr_I = aloc1(9);
  double *F = aloc1(9);
  double *C = aloc1(9);
  double *C_I = aloc1(9);

  /* physical quantities */
  double *disp = aloc1(ndn*nne);
  double *p = aloc1(nne);
  double *grad_delP = aloc1(3);
  double *grad_P = aloc1(3);
  double *devSbar = aloc1(9);
  double *L = aloc1(81);
  double *SS = aloc1(9);

  /* tangents */
  double *Kuu = aloc1(ndn*nne*ndn*nne);
  double *Kup = aloc1(ndn*nne*nne);
  double *Kpu = aloc1(ndn*nne*nne);
  double *Kpp = aloc1(nne*nne);

  /*=== get nodal values ===*/
  {
    int k = 0;
    for(int i=0; i<nne; i++){
      for(int j=0; j<ndofn; j++){
    if(j<ndn){
      disp[i*ndn+j] = r_e[k+j];
    } else if(j==ndn){
      p[i] = r_e[k+j];
    }
      }
      k += ndofn;
    }
  }

  /*=== main integration loop ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
        int_pt_ksi,int_pt_eta,int_pt_zet,
        weights);
  int ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){

    /* declare new variables for integration point */
    double Jr,Jn,pressure,wt,jj,Upp,Ybar,Pn,Pn_1,UP;

    /* Get pointer to Fn and compute Jn */
    /* DO NOT free Fn!!! */
    const double *Fn = NULL;
    if (sup->multi_scale){
      Fn = identity;
      Jn = 1.0;
    } else {
      Fn = eps[ii].il[ip].F;
      Jn = getJacobian(Fn,ii,&err);
    }

    /* get pointer to damage object */
    const damage *ptrDam = &eps[ii].dam[ip];


    const double Un_1 = eps[ii].il[ip].Un_1;
    const double Jn_1 = eps[ii].il[ip].Jn_1;

    /* compute various objects needed for integration */
    err += integration_help(ip,ii,nne,i,j,k,x,y,z,
                int_pt_ksi,int_pt_eta,int_pt_zet,
                weights,disp,sig[ii].pn_1,sig[ii].p,p,
                Fn,sup,&wt,&jj,Na,N_x,N_y,N_z,Np_x,ST,Fr,F,
                C,Fr_I,&Jr,&Pn_1,&Pn,&pressure,grad_delP,
                grad_P);

    /* compute C_I */
    err += inverse(C,3,C_I);

    /* get devSbar */
    err += get_material_stress(devSbar,C,ptrMat);

    /* get material stiffness */
    err += get_material_tangent(L,SS,&Ybar,kappa,Un_1,Jn_1,Jn,Jr,
                    Fn,Fr,C,C_I,Pn_1,Pn,pressure,
                    devSbar,ptrDam,ptrMat);

    /* get linearized pressure term */
    err += get_material_lin_pres(&Upp,Jn,Jr,ptrMat);

    /* get undamaged dUdJ|(n+1) */
    {
      dUdJFuncPtr dudJ = getDUdJFunc(0,ptrMat);
      dudJ(Jn*Jr,ptrMat,&UP);
    }

    /* get contributions to tangents from integration point */
    err += get_Kuu_at_ip(Kuu,nne,ST,Fn,Fr,Fr_I,Jn,Jr,ptrDam,
                 pressure,devSbar,L,jj,wt);

    err += get_Kup_at_ip(Kup,nne,Na,ST,Fr_I,Ybar,Jn_1,Jn,Jr,ptrDam,jj,wt);

    err += get_Kpu_at_ip(Kpu,nne,Na,ST,Fr_I,Jn,Jr,kappa,Upp,UP,
                 pressure,SS,ptrDam,jj,wt);

    err += get_Kpu_stab_at_ip(Kpu,sup->multi_scale,nne,ST,Np_x,Fr,Fr_I,Jr,grad_delP,
                  grad_P,SStab,ptrDam,SS,jj,wt);

    err += get_Kpp_stab_at_ip(Kpp,nne,Na,Np_x,Fr_I,Jn_1,Jn,Jr,SStab,
                  grad_P,Ybar,ptrDam,jj,wt);

    /* increment ip counter */
    ip++;

    /* error detected, deallocate and exit */
    if(err != 0) goto dealocate;
      }
    }
  }

  /*=== quadradic integration ===*/
  integrate(10,&npt_x,&npt_y,&npt_z,
        int_pt_ksi,int_pt_eta,int_pt_zet,
        weights);
  ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){

    /* declare fresh variables for integration */
    double Jr,delP,jj,wt,Pn_1,Pn,UP,H;

    double Jn = 0.0;
    const double *Fn = NULL;
    if (sup->multi_scale){
      Fn = identity;
      Jn = 1.0;
    } else {
      Fn =  eps[ii].st[ip].Fpp;
      Jn = getJacobian(Fn,ii,&err);
    }

    /* get pointer to damage object */
    const damage *ptrDam = &eps[ii].st[ip].dam;
    const double Un_1 = eps[ii].st[ip].Un_1;
    const double Jn_1 = eps[ii].st[ip].Jn_1;

    /* compute various quantities for integration */
    err += quadradic_integration_help(ip,ii,nne,i,j,k,x,y,z,
                      int_pt_ksi,int_pt_eta,int_pt_zet,
                      weights,disp,sig[ii].p,p,sup,&wt,&jj,
                      Na,N_x,N_y,N_z,ST,Fr,&Jr,&delP,&Pn);

    /* get Pn_1 */
    Pn_1 = 0.0;
    for(int l=0; l<nne; l++){
      Pn_1 += Na[l]*sig[ii].pn_1[l];
    }

    /* compute C */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
            3,3,3,1.0,Fr,3,Fn,3,0.0,F,3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
            3,3,3,1.0,F,3,F,3,0.0,C,3);


    /* get undamaged dUdJ|(n+1) and H*/
    {
      dUdJFuncPtr dudJ = getDUdJFunc(0,ptrMat);
      dudJ(Jn*Jr,ptrMat,&UP);
      if(ptrDam->damaged){
        double Ybar = 0.0;
        err += stab_get_material_potential(&Ybar,kappa,Un_1,Jn_1,
                           Jn*Jr,Pn_1,Pn,delP+Pn,
                           C,ptrMat);
        H = (ptrDam->dmu/(1.+ptrDam->dmu)
         *ptrDam->evolution(Ybar,&ptrDam->params));
      } else H = 0.0;
    }

    /* get contribution to Kpp at integration point */
    err += get_Kpp_at_ip(Kpp,nne,Na,kappa,UP,delP+Pn,
                 Jn_1,Jn,Jr,H,ptrDam,jj,wt);

    /* increment ip counter */
    ip++;

    /* error detected, deallocate and exit */
    if(err != 0) goto dealocate;
      }
    }
  }

 dealocate:
  /*=== DEALLOCATION ===*/
  /* integration */
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  free(Np_x);
  free(ST);

  /* kinematics */
  free(Fr);
  free(Fr_I);
  free(F);
  free(C);
  free(C_I);

  /* physical quantities */
  free(disp);
  free(p);
  free(grad_delP);
  free(grad_P);
  free(devSbar);
  free(L);
  free(SS);

  /* if there is an error, exit; do not assemble */
  if(err != 0) goto exit_function;

  /*=== ASSEMBLY ===*/
  for(int a=0; a<nne; a++){
    for(int b=0; b<ndofn; b++){
      for(int w=0; w<nne; w++){
    for(int g=0; g<ndofn; g++){
      if(b < ndn && g< ndn){
        Ks[idx_K(a,b,w,g,nne,ndofn)] = Kuu[idx_K(a,b,w,g,nne,ndn)];
      } else if(b<ndn && g==ndn){
        Ks[idx_K(a,b,w,g,nne,ndofn)] = Kup[idx_K_gen(a,b,w,0,
                             nne,ndn,nne,1)];
      } else if(b==ndn && g<ndn){
        Ks[idx_K(a,b,w,g,nne,ndofn)] = Kpu[idx_K_gen(a,0,w,g,
                             nne,1,nne,ndn)];
      } else if(b==ndn && g==ndn){
        Ks[idx_K(a,b,w,g,nne,ndofn)] = Kpp[idx_K(a,0,w,0,nne,1)];
      } else {
        PGFEM_printf("Error in assembly: %s [rank: %d, elen: %ld]\n",
           __func__,0,ii);
        Ks[idx_K(a,b,w,g,nne,ndofn)] = 0.0;
      }
    }
      }
    }
  }

  if(STAB_DEBUG){
    int err_rank = 0;
    char filename[50];
    FILE *debug_log;
    PGFEM_Error_rank(&err_rank);
    sprintf(filename,"new_stab_stiff_%ld_%d.log",ii,err_rank);
    debug_log = fopen(filename,"a");
    PGFEM_fprintf(debug_log,
        "=========================================\n"
        "=========== New Computation =============\n"
        "=========================================\n");
    PGFEM_fprintf(debug_log,"Kuu\n");
    print_array_d(debug_log,Kuu,nne*ndn*nne*ndn,nne*ndn,nne*ndn);
    PGFEM_fprintf(debug_log,"Kup\n");
    print_array_d(debug_log,Kup,nne*ndn*nne,nne*ndn,nne);
    PGFEM_fprintf(debug_log,"Kpu\n");
    print_array_d(debug_log,Kpu,nne*ndn*nne,nne,nne*ndn);
    PGFEM_fprintf(debug_log,"Kpp\n");
    print_array_d(debug_log,Kpp,nne*nne,nne,nne);
    PGFEM_fprintf(debug_log,"Assembled\n");
    print_array_d(debug_log,Ks,nne*ndofn*nne*ndofn,nne*ndofn,nne*ndofn);
    fclose(debug_log);
  }

 exit_function:
  /* deallocate tangents */
  free(Kuu);
  free(Kup);
  free(Kpu);
  free(Kpp);

  if(STAB_DEBUG && err != 0){
    PGFEM_printf("error in %s [elem %ld]\n",__func__,ii);
  }
  return err;
}/* int stiffmatel_st ()  */

int st_increment (long ne,
          long nn,
          long ndofn,
          long ndofd,
          MATGEOM matgeom,
          HOMMAT *hommat,
          ELEMENT *elem,
          NODE *node,
          SUPP sup,
          EPS *eps,
          SIG *sig,
          double *d_r,
          double *r,
          double nor_min,
          double stab,
          double dt,
          long nce,
          COEL *coel,
          double *pores,
          MPI_Comm mpi_comm,
          const int coh,
          const int mp_id)
{
  static const int ndn = 3;
  int err = 0;
  long ii,nne,*nod,i,j,k,II,M,P,R,ndofe,U,*cn;
  double *x,*y,*z,*r_e,*r_u,*p,AA[3][3],PL,EL_e,
    FoN[3][3],pom,BB[3][3],*a,*r_a,*a_a;
  double GEL_e,GPL,Gpores,*LPf,*GPf;

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  /************ UNIT CELL APPROACH + GFEM - TOTAL LAGRANGIAN *******/
  if (periodic == 1){
    EL_e = 0.0;
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
    FoN[P][R] = eps[0].FB[P][R];
    eps[0].P[P][R] = eps[0].FB[P][R] = eps[0].Fe[P][R] = 0.0;
      }
    }
  }/* end PERIODIC */

  /* Deformation rate */
  for (P=0;P<3;P++){
    for (R=0;R<3;R++){
      eps[0].Dp[P][R] = 0.0;
    }
  }

  for (ii=0;ii<ne;ii++){

    /* Number of element nodes */
    nne = elem[ii].toe;
    /* Nodes on element */
    nod = aloc1l (nne);
    elemnodes (ii,nne,nod,elem);
    /* Element Dof */
    ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

    /* allocation */
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);
    r_e = aloc1 (ndofe);
    r_u = aloc1 (nne*ndn);
    p = aloc1 (nne);
    cn = aloc1l (ndofe);
    r_a = aloc1 (nne*ndn);
    a_a = aloc1 (ndn*nne);
    a = aloc1 (ndofe);

    /* Id numbers */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* Coordinates/deformation on element */
    if(sup->multi_scale){
      nodecoord_total (nne,nod,node,x,y,z);
      def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
    } else {
      nodecoord_updated (nne,nod,node,x,y,z);
      def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);
    }

    /* Null fields */
    for (k=0;k<6;k++){
      eps[ii].el.o[k] = 0.0;
      sig[ii].el.o[k] = 0.0;
    }

    /* Displacement and pressure unknowns */
    k=0;
    for (i=0;i<nne;i++){
      for (j=0;j<ndofn;j++){
    if (j  < ndn){
      r_u[i*ndn+j] = r_e[k+j];
    }
    if (j == ndn){
      p[i] = r_e[k+j];
    }
    if (j  > ndn) {
      r_a[i*ndn+j-ndofn] = r_e[k+j];
      a_a[i*ndn+j-ndofn] = a[k+j];
    }
      }
      k += j;
    }

    /* Update the element quantities */
    EL_e += st_incr_elem (ii,nne,ndofn,x,y,z,r_u,a_a,
              p,dt,hommat,elem,eps,sig,sup);

    /* PRESSURE CHANGES */
    for (M=0;M<nne;M++){
      sig[ii].pn_1[M] = sig[ii].p[M];
      sig[ii].p[M] += p[M];
    }

    dealoc1l (nod);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1 (r_e);
    dealoc1 (r_u);
    dealoc1 (p);
    dealoc1l (cn);
    dealoc1 (r_a);
    dealoc1 (a_a);
    dealoc1 (a);
  }/* end ii < ne */

  if (periodic == 1){
    /* Averaging over the domains */
    GPf = aloc1 (nproc*27);
    LPf = aloc1 (27);

    U = 0;
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
    LPf[U+0]  = eps[0].P[P][R];
    LPf[U+9]  = eps[0].FB[P][R];
    LPf[U+18] = eps[0].Fn[P][R];
    U++;
      }
    }

    MPI_Allgather (LPf,27,MPI_DOUBLE,GPf,27,MPI_DOUBLE,mpi_comm);

    for (M=0;M<nproc;M++){
      if (M == myrank) continue;
      U = 0;
      for (P=0;P<3;P++){
    for (R=0;R<3;R++){
      eps[0].P[P][R] += GPf[M*27+U+0];
      eps[0].FB[P][R] += GPf[M*27+U+9];
      eps[0].Fn[P][R] += GPf[M*27+U+18];
      U++;
    }
      }
    }

    dealoc1 (GPf);
    dealoc1 (LPf);

    /* MACRO ENERGY || P:dF */
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
    AA[P][R] = 1./dt * (eps[0].F[P][R] - FoN[P][R]);
      }
    }

    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
    BB[P][R] = 0.0;
    for (U=0;U<3;U++){
      BB[P][R] += 1./2. *(AA[U][P]*eps[0].FB[U][R]
                  + eps[0].FB[U][P]*AA[U][R]);
    }
      }
    }

    pom = 0.0;
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
    pom  += eps[0].P[P][R] * AA[P][R];
    eps[0].Fn[P][R] = eps[0].F[P][R];
      }
    }

    MPI_Allreduce(&EL_e,&GEL_e,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    if ( myrank == 0){
      PGFEM_printf("P:dF = %12.12f || Aver. of Mic. S:dE = %12.12f\n",pom,GEL_e);
    }
  }/* end periodic */

  GPf = aloc1(nproc*9);
  LPf = aloc1(9);

  U = 0;
  for (P=0;P<3;P++){
    for (R=0;R<3;R++){
      LPf[U] = eps[0].Dp[P][R];
      U++;
    }
  }

  MPI_Allgather (LPf,9,MPI_DOUBLE,GPf,9,MPI_DOUBLE,mpi_comm);

  for (M=0;M<nproc;M++){
    if (M == myrank){
      continue;
    }
    U = 0;
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
    eps[0].Dp[P][R] += GPf[M*9+U];
    U++;
      }
    }
  }

  dealoc1 (GPf);
  dealoc1 (LPf);

  pom = 0.0;
  for (P=0;P<3;P++){
    for (R=0;R<3;R++){
      pom += eps[0].Dp[P][R]*eps[0].Dp[P][R];
    }
  }

  /* Macro-scale effective plastic or elastic strain */
  eps[0].eff += dt*sqrt(2./3.*pom); /* int_0^t (sqrt (2/3*d.d))dt */


  /*********************/
  /* Coordinate update */
  /*********************/
   for (ii=0;ii<nn;ii++){
    for (i=0;i<ndn;i++){
      II = node[ii].id_map[mp_id].id[i];
      if (II > 0){
    if (i == 0) node[ii].x1 = node[ii].x1_fd + r[II-1] + d_r[II-1];
    else if (i == 1) node[ii].x2 = node[ii].x2_fd + r[II-1] + d_r[II-1];
    else if (i == 2) node[ii].x3 = node[ii].x3_fd + r[II-1] + d_r[II-1];
      }
      else if (II < 0){
    if (i == 0) node[ii].x1 = (node[ii].x1_fd
                   + sup->defl[abs(II)-1]
                   + sup->defl_d[abs(II)-1]);

    else if (i == 1) node[ii].x2 = (node[ii].x2_fd
                    + sup->defl[abs(II)-1]
                    + sup->defl_d[abs(II)-1]);

    else if (i == 2) node[ii].x3 = (node[ii].x3_fd
                    + sup->defl[abs(II)-1]
                    + sup->defl_d[abs(II)-1]);
      }
    }
  }/* end ii < nn */

   /* update the coordinates including macroscale deformations */
   if(sup->multi_scale){
     const double *F = sup->F0;
     double X[ndn];
     double Y[ndn];
     for(int i=0; i<nn; i++){
       X[0] = node[i].x1_fd;
       X[1] = node[i].x2_fd;
       X[2] = node[i].x3_fd;
       cblas_dgemv(CblasRowMajor,CblasNoTrans,ndn,ndn,1.0,
           F,ndn,X,1,0.0,Y,1);
       node[i].x1 += Y[0];
       node[i].x2 += Y[1];
       node[i].x3 += Y[2];
     }
   }

  PL = T_VOLUME (ne,ndofn-1,elem,node);
  /* Gather Volume from all domains */
  MPI_Allreduce(&PL,&GPL,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

  if (coh == 1) {
    MPI_Allreduce(pores,&Gpores,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    *pores = Gpores;
  }

  if (myrank == 0) {
    if (coh == 0){
      PGFEM_printf ("AFTER DEF - VOLUME = %12.12f\n",GPL);
    } else {
      PGFEM_printf ("AFTER DEF - VOLUME = %12.12f || Volume of voids %12.12f || Vv/V = %12.12f\n",
            GPL,*pores,*pores/GPL);
      /*VVolume --> GPL*/ /* VVolume doesn't make sense...*/
    }
  }

    return err;
}

/*=====================================================================*/
/*                  STATIC HELPER FUNCTIONS                            */
/*=====================================================================*/

/**** NOTE: Potential/stress functions use new potentials!!! ****/

static int integration_help(const int ip,
                const int elem_id,
                const int nne,
                const int i,
                const int j,
                const int k,
                const double *x,
                const double *y,
                const double *z,
                const double *int_pt_ksi,
                const double *int_pt_eta,
                const double *int_pt_zet,
                const double *weights,
                const double *disp,
                const double *pn_1,
                const double *pn,
                const double *p,
                const double *Fn,
                const SUPP sup,
                double *wt,
                double *jj,
                double *Na,
                double *N_x,
                double *N_y,
                double *N_z,
                double *Np_x,
                double *ST,
                double *Fr,
                double *F,
                double *C,
                double *Fr_I,
                double *Jr,
                double *Pn_1,
                double *Pn,
                double *pressure,
                double *grad_delP,
                double *grad_P)
{
  static const int ndn = 3;
  int err = 0;
  double ksi,eta,zet;

  if(nne == 8){/* hexahedron */
    ksi = int_pt_ksi[i];
    eta = int_pt_eta[j];
    zet = int_pt_zet[k];
    (*wt) = weights[i]*weights[j]*weights[k];
  } else { /* tetrahedron type */
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
    (*wt) = weights[ip];
  }

  double ****ST_tensor, **Fr_mat;
  ST_tensor = aloc4 (3,3,ndn,nne);
  Fr_mat = aloc2(3,3);

  /* compute the shape functions and their derivatives */
  shape_func(ksi,eta,zet,nne,Na);
  (*jj) = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);

  /* compute the pressure gradient operator and pressure */
  *Pn_1 = *Pn = *pressure = 0.0;
  memset(grad_delP,0,3*sizeof(double));
  memset(grad_P,0,3*sizeof(double));
  for(int i=0; i<nne; i++){
    Np_x[idx_2_gen(0,i,3,nne)] = N_x[i];
    Np_x[idx_2_gen(1,i,3,nne)] = N_y[i];
    Np_x[idx_2_gen(2,i,3,nne)] = N_z[i];

    *Pn_1 += Na[i]*(pn_1[i]);
    *Pn += Na[i]*(pn[i]);
    *pressure += Na[i]*(pn[i] + p[i]);

    grad_delP[0] += p[i]*Np_x[idx_2_gen(0,i,3,nne)];
    grad_delP[1] += p[i]*Np_x[idx_2_gen(1,i,3,nne)];
    grad_delP[2] += p[i]*Np_x[idx_2_gen(2,i,3,nne)];

    grad_P[0] += (pn[i]+p[i])*Np_x[idx_2_gen(0,i,3,nne)];
    grad_P[1] += (pn[i]+p[i])*Np_x[idx_2_gen(1,i,3,nne)];
    grad_P[2] += (pn[i]+p[i])*Np_x[idx_2_gen(2,i,3,nne)];
  }

  /* build the shape tensor */
  shape_tensor (nne,ndn,N_x,N_y,N_z,ST_tensor);
  shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);

  /* compute the current deformation gradient */
  def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
  mat2array(Fr,CONST_2(double) Fr_mat,3,3);

  /*=== compute F and C ===*/
  /*
   * Add multiscale contribution.
   * For multi-scale analysis, Fn = Identity, F = Fr
   */
  if(sup->multi_scale){
    for(int i=0; i<ndn*ndn; i++){
      Fr[i] +=  sup->F0[i];
    }
    memcpy(F,Fr,ndn*ndn*sizeof(double));
  } else {
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
        3,3,3,1.0,Fr,3,Fn,3,0.0,F,3);
  }

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,F,3,F,3,0.0,C,3);

  /* compute the Jacobian determinant and get the error flag */
  *Jr = getJacobian(Fr,elem_id,&err);

  /* compute Fr_I */
  err += inverse(Fr,3,Fr_I);

  /* dealocate matrix/tensor containers */
  dealoc4(ST_tensor,3,3,ndn);
  dealoc2(Fr_mat,3);

  return err;
}/* static int integration_help() */

static int quadradic_integration_help(const int ip,
                      const int elem_id,
                      const int nne,
                      const int i,
                      const int j,
                      const int k,
                      const double *x,
                      const double *y,
                      const double *z,
                      const double *int_pt_ksi,
                      const double *int_pt_eta,
                      const double *int_pt_zet,
                      const double *weights,
                      const double *disp,
                      const double *pn,
                      const double *p,
                      const SUPP sup,
                      double *wt,
                      double *jj,
                      double *Na,
                      double *N_x,
                      double *N_y,
                      double *N_z,
                      double *ST,
                      double *Fr,
                      double *Jr,
                      double *delP,
                      double *Pn)
{
  static const int ndn = 3;
  int err = 0;
  double ksi,eta,zet;

  if(nne == 8){/* hexahedron */
    ksi = int_pt_ksi[i];
    eta = int_pt_eta[j];
    zet = int_pt_zet[k];
    (*wt) = weights[i]*weights[j]*weights[k];
  } else { /* tetrahedron type */
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
    (*wt) = weights[ip];
  }

  double ****ST_tensor, **Fr_mat;
  ST_tensor = aloc4 (3,3,ndn,nne);
  Fr_mat = aloc2(3,3);

  /* compute the shape functions and their derivatives */
  shape_func(ksi,eta,zet,nne,Na);
  (*jj) = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);

  /* build the shape tensor */
  shape_tensor (nne,ndn,N_x,N_y,N_z,ST_tensor);
  shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);

  /* compute the current deformation gradient */
  def_grad_get(nne,ndn,CONST_4(double) ST_tensor,disp,Fr_mat);
  mat2array(Fr,CONST_2(double) Fr_mat,3,3);

  /* Add multiscale contribution */
  if(sup->multi_scale){
    for(int i=0; i<ndn*ndn; i++){
      Fr[i] += sup->F0[i];
    }
  }

  /* compute delP */
  *Pn = *delP = 0.0;
  for(int i=0; i<nne; i++){
    *Pn += Na[i]*pn[i];
    *delP += Na[i]*p[i];
  }

 /* compute the Jacobian determinant and get the error flag */
  *Jr = getJacobian(Fr,elem_id,&err);

  /* dealocate matrix/tensor containers */
  dealoc4(ST_tensor,3,3,ndn);
  dealoc2(Fr_mat,3);

  return err;
}/* static int quadradic_integration_help() */

static int get_material_stress(double *devSbar,
                   const double *C,
                   const HOMMAT *mat)
{
  int err = 0;
  devStressFuncPtr dev_stress = getDevStressFunc(0,mat);
  dev_stress(C,mat,devSbar);
  return err;
}/* static int get_material_stress() */

static int get_material_tangent(double *L,
                double *SS,
                double *Ybar,
                const double kappa,
                const double Un_1, /* U(n-1) */
                const double Jn_1, /* J(n-1) */
                const double Jn,
                const double Jr,
                const double *Fn,
                const double *Fr,
                const double *C,
                const double *C_I,
                const double Pn_1, /* p(n-1) */
                const double Pn,
                const double pressure,
                const double *devS,
                const damage *dam,
                const HOMMAT *mat)
{
  /*** NOTE: stabilized formulation needs dSdF ****/
  int err = 0;

  /* compute damage */
  double H = 0.0;
  *Ybar = 0.0;
  if(dam->damaged){
    err = stab_get_material_potential(Ybar,kappa,Un_1,Jn_1,Jn*Jr,
                      Pn_1,Pn,pressure,C,mat);
    H = dam->dmu/(1+dam->dmu)*dam->evolution(*Ybar,&(dam->params));
  }

  /* compute dSdC */
  /*=== NOTE:
    dSdC has been multiplied by 2 for symmetry in other
    formulations and needs to be accounted for here
    ===*/
  double *dSdC = aloc1(81);
  matStiffFuncPtr dev_stiff = getMatStiffFunc(0,mat);
  dev_stiff(C,mat,dSdC);

  double UP = 0.0;
  UFuncPtr VolPot = getUFunc(0,mat);
  VolPot(Jn*Jr,mat,&UP);

  /* compute compSbar = devS + Jr Jn p C_I */
  double *compSbar = aloc1(9);
  double *Sbar = aloc1(9);
  double *SoxS = aloc1(81);

  /* compute effective comp and continuum stress and clear SS */
  for(int i=0; i<9; i++){
    SS[i] = 0.0;
    /* Sbar[i] = devS[i] + Jr*Jn*C_I[i]*(pressure + 4*Pn + Pn_1)/6.; */
    Sbar[i] = devS[i] + Jr*Jn*C_I[i]*kappa*UP;
    compSbar[i] = devS[i] + Jr*Jn*C_I[i]*pressure;
  }

  /* compute 0.5*dCdFr & compSbar x Sbar */
  double *FrFn = aloc1(9);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
          3,3,3,1.0,Fr,3,Fn,3,0.0,FrFn,3);
  double *dCdFr = aloc1(81);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
    for(int l=0; l<3; l++){
      dCdFr[idx_4(i,j,k,l)] = 0.5*(Fn[idx_2(l,i)]*FrFn[idx_2(k,j)]
                       + Fn[idx_2(l,j)]*FrFn[idx_2(k,i)]);
      SoxS[idx_4(i,j,k,l)] = compSbar[idx_2(i,j)]*Sbar[idx_2(k,l)];
    }
      }
    }
  }

  for(int i=0; i<81; i++){
    dSdC[i] = (1-dam->w)*dSdC[i] - H*SoxS[i];
    L[i] = 0.0;
  }

  /* compute material tangent */
  for(int m=0; m<3; m++){
    for(int n=0; n<3; n++){
      for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      /* Using the idx_2 function and setting the cblas offset to
         9 allows contracion of the 4th order tensors by
         cblas_ddot */
      L[idx_4(i,j,m,n)] = cblas_ddot(9,&dSdC[idx_4(i,j,0,0)],1,
                     &dCdFr[idx_2(m,n)],9);
      SS[idx_2(m,n)] += H*Sbar[idx_2(i,j)]*dCdFr[idx_4(i,j,m,n)];
    }
      }
    }
  }

  free(dSdC);
  free(Sbar);
  free(compSbar);
  free(SoxS);
  free(FrFn);
  free(dCdFr);

  return err;
}/* static int get_material_tangent() */

static int get_material_pres(double *Up,
                 const double Jn,
                 const double Jr,
                 const HOMMAT *mat,
                 const damage *dam)
{
  int err = 0;
  /* stabilized formulation uses the change in pressure, not total
     pressure */
  dUdJFuncPtr dudJ = getDUdJFunc(0,mat);
  double Upn;
  *Up = Upn = 0.0;

  /* compute dU/dJ|n+1 */
  dudJ(Jn*Jr,mat,Up);
  /* compute dU/dJ|n */
  dudJ(Jn,mat,&Upn);

  *Up = (1-dam->w)*(*Up) - (1-dam->wn)*Upn;
  return err;
}/* static int get_materaial_pres() */

static int get_material_lin_pres(double *Upp,
                 const double Jn,
                 const double Jr,
                 const HOMMAT *mat)
{
  int err = 0;
  d2UdJ2FuncPtr upp_func = getD2UdJ2Func(0,mat);
  upp_func(Jn*Jr,mat,Upp);

  return err;
}/* static int get_materail_lin_pres() */

static int get_Ru_at_ip(double *Ru,
            const int nne,
            const double *ST,
            const double *Fn,
            const double *Fr,
            const double *Fr_I,
            const double Jr,
            const double Jn,
            const double *devSbar,
            const double pres,
            const damage *dam,
            const double jj,
            const double wt)
{
  int err = 0;
  /* compute devS from devSbar */
  double *devS = aloc1(9);
  for(int i=0; i<9; i++){
    devS[i] = (1.-dam->w)*devSbar[i];
  }

  /* add contribution of integration point to Ru */
  UL_Ru_at_ip(Ru,nne,nne,ST,Fn,Fr,Fr_I,Jn,Jr,devS,
          (1. - dam->w)*pres,jj,wt,0);

  free(devS);
  return err;
}/* static int get_Ru_at_ip() */


static int get_Rp_at_ip(double *Rp,
            const int TL,
            const int nne,
            const double *Np,
            const double kappa,
            const double deltaUp,
            const double deltaPres,
            const double P,
            const damage *dam,
            const double jj,
            const double wt)
{
  int err = 0;

  /* compute dP for residual */
  double dP = (1.-dam->wn)*deltaPres + (dam->wn - dam->w)*P;
  if(TL) dP = (1.-dam->w)*P;

  /* add contribution of integration point to Rp */
  UL_Rp_at_ip(Rp,nne,Np,kappa,deltaUp,dP,jj,wt);

  return err;
}/* static int get_Rp_at_ip() */


static int get_Rp_stab_at_ip(double *Rp,
                 const int TL,
                 const int nne,
                 const double *Np_x,/* Np_x [ndn x nne] */
                 const double *grad_delP,
                 const double *grad_P,
                 const double Jr,
                 const double *Fr_I,
                 const double stab,
                 const damage *dam,
                 const double jj,
                 const double wt)
{
  int err = 0;
  if(stab == 0.0) {/* zero stabilization term */
    return err;
  }

  double *FrFr = aloc1(9);
  double *gPgP = aloc1(9);
  double *gradP = aloc1(3);

  /* precompute stab*Jr*Fr_I*Fr_I' */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
          3,3,3,stab*Jr,Fr_I,3,Fr_I,3,0.0,FrFr,3);

  /* compute damage gradP */
  if(TL){
    for(int i=0; i<3; i++){
      gradP[i] = (1.-dam->w)*grad_P[i];
    }
  } else {
    for(int i=0; i<3; i++){
      gradP[i] = (1.-dam->wn)*grad_delP[i] + (dam->wn - dam->w)*grad_P[i];
    }
  }

  /* add contribution of stabilization term to Rp at the integration
     point */
  for(int i=0; i<nne; i++){
    memset(gPgP,0,9*sizeof(double));
    cblas_dger(CblasRowMajor,3,3,1.0,gradP,1,&Np_x[i],nne,gPgP,3);
    Rp[i] -= cblas_ddot(9,FrFr,1,gPgP,1)*jj*wt;
  }

  free(FrFr);
  free(gPgP);
  free(gradP);
  return err;
}/* static int get_Rp_stab_at_ip() */


static int get_Kuu_at_ip(double *Kuu,
             const int nne,
             const double *ST,
             const double *Fn,
             const double *Fr,
             const double *Fr_I,
             const double Jn,
             const double Jr,
             const damage *dam,
             const double totalP,
             const double *devSbar,
             const double *L,
             const double jj,
             const double wt)
{
  int err = 0;

  /* Compute devS from devSbar */
  double *devS = aloc1(9);
  for(int i=0; i<9; i++){
    devS[i] = (1. - dam->w)*devSbar[i];
  }

  /* add contribution of integration point to Kuu */
  UL_Kuu_at_ip(Kuu,nne,nne,ST,Fn,Fr,Fr_I,Jn,Jr,
           (1.-dam->w)*totalP,devS,L,jj,wt,0);

  free(devS);
  return err;
}/* static int get_Kuu_at_ip() */


static int get_Kup_at_ip(double *Kup,
             const int nne,
             const double *Np,
             const double *ST,
             const double *Fr_I,
             const double Ybar,
             const double Jn_1,
             const double Jn,
             const double Jr,
             const damage *dam,
             const double jj,
             const double wt)
{
  int err = 0;

  double mult = (1.-dam->w);
  /* if(dam->damaged){ */
  /*   mult -= (dam->dmu/(1. + dam->dmu)*dam->evolution(Ybar,&dam->params) */
  /*         *(Jn*Jr-Jn_1)/6.); */
  /* } */

  /* add contribution of integration point to Kup */
  UL_Kup_at_ip(Kup,nne,nne,Np,ST,Fr_I,mult*Jr,jj,wt,0);

  return err;
}/* static int get_Kup_at_ip() */


static int get_Kpu_at_ip(double *Kpu,
             const int nne,
             const double *Np,
             const double *ST,
             const double *Fr_I,
             const double Jn,
             const double Jr,
             const double kappa,
             const double Upp,
             const double UP, /* n+1 only */
             const double P, /* total pressure */
             const double *SS,
             const damage *dam,
             const double jj,
             const double wt)
{
  int err = 0;

  /* add contribution of integration point to Kpu */
  /* UL_Kpu_at_ip(Kpu,nne,nne,Np,ST,Fr_I,Jn,Jr,Upp,jj,wt,0); */

  damage_UL_Kpu_at_ip(Kpu,nne,nne,Np,ST,Fr_I,Jn,Jr,
              kappa,Upp,UP,P,SS,dam->w,jj,wt,0);

  return err;
}/* static int get_Kpu_at_ip() */


static int get_Kpp_at_ip(double *Kpp,
             const int nne,
             const double *Np,
             const double kappa,
             const double UP, /* n+1 only */
             const double P, /* total pressure */
             const double Jn_1,
             const double Jn,
             const double Jr,
             const double H,
             const damage *dam,
             const double jj,
             const double wt)
{
  int err = 0;

  /* add contribution of integration point to Kpu */
  /* UL_Kpp_at_ip(Kpp,nne,nne,Np,kappa,jj,wt); */

  double HH = 0.0; /* (UP-P/kappa)*H*(Jn*Jr-Jn_1)/6.; */
  damage_UL_Kpp_at_ip(Kpp,nne,nne,Np,kappa,dam->w,HH,jj,wt);

  return err;
}/* static int get_Kpp_at_ip() */


static int get_Kpu_stab_at_ip(double *Kpu,
                  const int TL,
                  const int nne,
                  const double *ST,
                  const double *Np_x,
                  const double *Fr,
                  const double *Fr_I,
                  const double Jr,
                  const double *grad_delP,
                  const double *grad_P,
                  const double stab,
                  const damage *dam,
                  const double *SS,
                  const double jj,
                  const double wt)
{
  const int ndn = 3;
  int err = 0;
  if(stab == 0.0) {/* zero stabilization term */
    return err;
  }

  double *Cr_I = aloc1(9);
  double *dam_gPgP = aloc1(9);
  double *gPgP = aloc1(9);
  double *AA = aloc1(9);
  double *dam_gradP = aloc1(3);

  /* compute Cr_I */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
          3,3,3,1.0,Fr_I,3,Fr_I,3,0.0,Cr_I,3);

  /* compute damage grad P */
  if(TL){
    for(int i=0; i<3; i++){
      dam_gradP[i] = (1.-dam->w)*grad_P[i];
    }
  } else {
    for(int i=0; i<3; i++){
      dam_gradP[i] = (1.-dam->wn)*grad_delP[i] + (dam->wn-dam->w)*grad_P[i];
    }
  }

  for(int a=0; a<nne; a++){

    /* compute gradP ox gradP terms */
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
    int c_idx = idx_2(i,j);
    dam_gPgP[c_idx] = dam_gradP[i]*Np_x[j*nne+a];
    gPgP[c_idx] = grad_P[i]*Np_x[j*nne+a];
      }
    }

    for(int w=0; w<nne; w++){
      for(int g=0; g<3 /*ndn*/; g++){
    const double* const ptrST_wg = &ST[idx_4_gen(w,g,0,0,
                             nne,3,3,3)];

    /* compute help2 = Fr_I':ST_wg */
    /* compute help3 = SS : ST_wg */
    double help2 = 0.0;
    double help3 = 0.0;
    for(int i=0; i<ndn; i++){
      for(int j=0; j<ndn; j++){
        /* AA = help2 Cr_I - Fr_I ST_wg Cr_I - Cr_I ST_wg' Fr_I' */
        const int c_idx = idx_2(i,j);
        AA[c_idx] = 0.0;
        if(i==0 && j==0){/* compute help2 and help3 as well */
          for(int k=0; k<ndn; k++){
        for(int l=0; l<ndn; l++){
          help2 += Fr_I[idx_2(l,k)]*ptrST_wg[idx_2(k,l)];
          help3 += SS[idx_2(k,l)]*ptrST_wg[idx_2(k,l)];
          AA[c_idx] += (-(Fr_I[idx_2(i,k)]*ptrST_wg[idx_2(k,l)]
                    *Cr_I[idx_2(l,j)])
                  -(Cr_I[idx_2(i,k)]*ptrST_wg[idx_2(l,k)]
                    *Fr_I[idx_2(j,l)]));
        }
          }
        } else {
          for(int k=0; k<ndn; k++){
        for(int l=0; l<ndn; l++){
          AA[c_idx] += (-(Fr_I[idx_2(i,k)]*ptrST_wg[idx_2(k,l)]
                    *Cr_I[idx_2(l,j)])
                  -(Cr_I[idx_2(i,k)]*ptrST_wg[idx_2(l,k)]
                    *Fr_I[idx_2(j,l)]));
        }
          }
        }
          AA[c_idx] += help2*Cr_I[c_idx];
      }
    }

    const int k_idx = idx_K_gen(a,0,w,g,nne,1,nne,3);
    for(int i=0; i<9; i++){
      Kpu[k_idx] += (stab*Jr*jj*wt*((help3 * Cr_I[i] * gPgP[i])
                    -(AA[i] * dam_gPgP[i])) );
    }

      }
    }
  }

  free(Cr_I);
  free(AA);
  free(dam_gPgP);
  free(gPgP);
  free(dam_gradP);
  return err;
}/* static int get_Kpu_stab_tan_at_ip() */


static int get_Kpp_stab_at_ip(double *Kpp,
                  const int nne,
                  const double *Np,
                  const double *Np_x,
                  const double *Fr_I,
                  const double Jn_1,
                  const double Jn,
                  const double Jr,
                  const double stab,
                  const double *grad_P,
                  const double Ybar,
                  const damage *dam,
                  const double jj,
                  const double wt)
{
  int err = 0;
  if(stab == 0.0) {/* zero stabilization term */
    return err;
  }

  double *gPgP = aloc1(9);

  /* pre compute AA = stab*Jr* Fr_I Fr_I' */
  double *AA = aloc1(9);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
          3,3,3,stab*Jr,Fr_I,3,Fr_I,3,0.0,AA,3);

  /* double H = 0.0; */
  /* if(dam->damaged){ */
  /*   H = (dam->dmu/(1.+dam->dmu)*dam->evolution(Ybar,&dam->params) */
  /*     *(Jr*Jn-Jn_1)/6.); */
  /* } */

  for(int a=0; a<nne; a++){
    double damage_term = 0.0;
    /* if(dam->damaged){ */
    /*   for(int i=0; i<3; i++){ */
    /*  for(int j=0; j<3; j++){ */
    /*    damage_term += (H*AA[idx_2(i,j)]*grad_P[i]*Np_x[j*nne+a]); */
    /*  } */
    /*   } */
    /* } */
    for(int w=0; w<nne; w++){
      memset(gPgP,0,9*sizeof(double));
      cblas_dger(CblasRowMajor,3,3,1.0,&Np_x[w],nne,&Np_x[a],nne,gPgP,3);
      Kpp[idx_K(a,0,w,0,nne,1)] += jj*wt*(damage_term*Np[w]
                      -(1-dam->w)*cblas_ddot(9,AA,1,
                                 gPgP,1) );
    }
  }

  free(gPgP);
  free(AA);
  return err;
}/* static int get_Kpp_stab_tan_at_ip() */

static double st_incr_elem (const long ii,
                const long nne,
                const long ndofn,
                const double *x,
                const double *y,
                const double *z,
                const double *r_u,
                const double *a_a,
                const double *p,
                const double dt,
                const HOMMAT *hommat,
                const ELEMENT *elem,
                EPS *eps,
                SIG *sig,
                const SUPP sup)
{
  double result = 0.0;
  /* compute constants */
  const int ndn = 3;
  const int mat = elem[ii].mat[2];
  const double  kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  const HOMMAT *ptrMat = &hommat[mat];

  const double volume = Tetra_V (x,y,z);

  int err = 0;

  /*=== ALLOCATION/declaration ===*/
  /* integration */
  long npt_x, npt_y, npt_z;
  /* alloce enough space for quadradic integration */
  int_point(10,&npt_z);
  double *int_pt_ksi = aloc1(npt_z);
  double *int_pt_eta = aloc1(npt_z);
  double *int_pt_zet = aloc1(npt_z);
  double *weights = aloc1(npt_z);
  double *Na = aloc1(nne);
  double *N_x = aloc1(nne);
  double *N_y = aloc1(nne);
  double *N_z = aloc1(nne);
  double *Np_x = aloc1(3*nne);
  double *ST = aloc1(3*3*ndn*nne);

  /* kinematics */
  /* Fn will be a temporary pointer within the integration loop and
     will not need allocation/dealocation */
  double *Fr = aloc1(9);
  double *Fr_I = aloc1(9);
  double *F = aloc1(9);
  double *F_I = aloc1(9);
  double *C = aloc1(9);
  double *C_I = aloc1(9);
  double *E = aloc1(9);

  /* physical quantities */
  double *grad_delP = aloc1(3);
  double *grad_P = aloc1(3);
  double *devSbar = aloc1(9);
  double *sigma = aloc1(9);
  double *almansi = aloc1(9);

  /*=== null output quantities ===*/
  memset(eps[ii].el.o,0,6*sizeof(double));
  memset(sig[ii].el.o,0,6*sizeof(double));

  /*=== main integration loop ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
        int_pt_ksi,int_pt_eta,int_pt_zet,
        weights);
  int ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){

    /* declare new variables for integration point */
    double Jr,Jn,pressure,wt,jj,Pn_1, Pn;

    /* Get pointer to Fn and compute Jn */
    /* DO NOT free Fn!!! */
    const double *Fn = NULL;
    if (sup->multi_scale){
      Fn = identity;
      Jn = 1.0;
    } else {
      Fn = eps[ii].il[ip].F;
      Jn = getJacobian(Fn,ii,&err);
    }

    /* get pointer to damage object */
    const damage *ptrDam = &eps[ii].dam[ip];

    /* compute various objects needed for integration */
    err += integration_help(ip,ii,nne,i,j,k,x,y,z,
                int_pt_ksi,int_pt_eta,int_pt_zet,
                weights,r_u,sig[ii].pn_1,sig[ii].p,p,
                Fn,sup,&wt,&jj,Na,N_x,N_y,N_z,Np_x,ST,Fr,
                F,C,Fr_I,&Jr,&Pn_1,&Pn,&pressure,
                grad_delP,grad_P);

    /* record deformation gradient at integration point */
    memcpy(eps[ii].il[ip].F,F,9*sizeof(double));

    /* compute C_I & F_I */
    err += inverse(C,3,C_I);
    err += inverse(F,3,F_I);

    /* compute the PK2 stress and E */
    err += get_material_stress(devSbar,C,ptrMat);
    for(int I=0; I<3; I++){
      for(int J=0; J<3; J++){
        devSbar[idx_2(I,J)] += pressure*Jr*Jn*C_I[idx_2(I,J)];
        E[idx_2(I,J)] = 0.5*(C[idx_2(I,J)] - DELTA(I,J));
      }
    }

    /*=== Store inelastic quantities ===*/
    /* Second Piola Kirchoff stress */
    sig[ii].il[ip].o[0] = (1.-ptrDam->w)*devSbar[0];
    sig[ii].il[ip].o[1] = (1.-ptrDam->w)*devSbar[4];
    sig[ii].il[ip].o[2] = (1.-ptrDam->w)*devSbar[8];

    sig[ii].il[ip].o[3] = (1.-ptrDam->w)*devSbar[5]; /* yz */
    sig[ii].il[ip].o[4] = (1.-ptrDam->w)*devSbar[2]; /* xz */
    sig[ii].il[ip].o[5] = (1.-ptrDam->w)*devSbar[1]; /* xy */

    /* Elastic Green Lagrange strain */
    eps[ii].il[ip].o[0] = E[0];
    eps[ii].il[ip].o[1] = E[4];
    eps[ii].il[ip].o[2] = E[8];

    eps[ii].il[ip].o[3] = 2.*E[5];
    eps[ii].il[ip].o[4] = 2.*E[2];
    eps[ii].il[ip].o[5] = 2.*E[1];

    /*=== integrate output quantities ===*/
    /* compute sigma and almansi */
    for(int I=0; I<3; I++){
      for(int J=0; J<3; J++){
        sigma[idx_2(I,J)] = almansi[idx_2(I,J)] = 0.0;
        for(int K=0; K<3; K++){
          for(int L=0; L<3; L++){
        sigma[idx_2(I,J)] += ((1. - ptrDam->w)/(Jr*Jn)
                      *F[idx_2(I,K)]*devSbar[idx_2(K,L)]
                      *F[idx_2(K,J)]);
        almansi[idx_2(I,J)] += (F_I[idx_2(K,I)]*E[idx_2(K,L)]
                    *F_I[idx_2(L,J)]);
          }
        }
      }
    }

    /* Almansi and/or Logarithmic strain */
    eps[ii].el.o[0] += wt*jj/volume*almansi[0];
    eps[ii].el.o[1] += wt*jj/volume*almansi[4];
    eps[ii].el.o[2] += wt*jj/volume*almansi[8];

    eps[ii].el.o[3] += 2.*wt*jj/volume*almansi[5];
    eps[ii].el.o[4] += 2.*wt*jj/volume*almansi[2];
    eps[ii].el.o[5] += 2.*wt*jj/volume*almansi[1];

    /* Cauchy Stress */
    sig[ii].el.o[0] += wt*jj/volume*sigma[0];
    sig[ii].el.o[1] += wt*jj/volume*sigma[4];
    sig[ii].el.o[2] += wt*jj/volume*sigma[8];

    sig[ii].el.o[3] += wt*jj/volume*sigma[5];
    sig[ii].el.o[4] += wt*jj/volume*sigma[2];
    sig[ii].el.o[5] += wt*jj/volume*sigma[1];

    /* update volumetric energy, J(n-1) */
    const double Jn_1 = eps[ii].il[ip].Jn_1;
    const double Un_1 = eps[ii].il[ip].Un_1;

    eps[ii].il[ip].Un_1 = eps[ii].st[ip].Un;
    eps[ii].il[ip].Jn_1 = Jn;
    eps[ii].il[ip].Un = (Un_1 + (Jr*Jn-Jn_1)/(6.*kappa)
                 *(pressure + 4*Pn + Pn_1));

    /* update damage variables */
    update_damage(&eps[ii].dam[ip]);

    /* increment ip counter */
    ip++;

    /* error detected, deallocate and exit */
    if(err != 0) goto dealocate;
      }
    }
  }

  /*=== quadradic integration ===*/
  integrate(10,&npt_x,&npt_y,&npt_z,
        int_pt_ksi,int_pt_eta,int_pt_zet,
        weights);
  ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){

    /* declare fresh variables for integration */
    double Jr,delP,jj,wt,Pn;

    double Jn = 0.0;
    const double *Fn = NULL;
    if (sup->multi_scale){
      Fn = identity;
      Jn = 1.0;
    } else {
      Fn =  eps[ii].st[ip].Fpp;
      Jn = getJacobian(Fn,ii,&err);
    }

    double Un_1 = eps[ii].st[ip].Un_1;
    double Jn_1 = eps[ii].st[ip].Jn_1;

    /* compute various quantities for integration */
    err += quadradic_integration_help(ip,ii,nne,i,j,k,x,y,z,
                      int_pt_ksi,int_pt_eta,int_pt_zet,
                      weights,r_u,sig[ii].p,p,sup,&wt,&jj,
                      Na,N_x,N_y,N_z,ST,Fr,&Jr,&delP,&Pn);

    /* get Pn_1 */
    double Pn_1 = 0.0;
    for(int l=0; l<nne; l++){
      Pn_1 += Na[i]*sig[ii].pn_1[i];
    }

    /* compute and store F(n+1) at integration point */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
            3,3,3,1.0,Fr,3,Fn,3,0.0,F,3);
    memcpy(eps[ii].st[ip].Fpp,F,9*sizeof(double));

    /* update volumetric energy, J(n-1) */
    eps[ii].st[ip].Un_1 = eps[ii].st[ip].Un;
    eps[ii].st[ip].Jn_1 = Jn;
    eps[ii].st[ip].Un = (Un_1 + (Jr*Jn-Jn_1)/(6.*kappa)
                 *((delP + Pn) + 4*Pn + Pn_1));

    /* update damage variables */
    update_damage(&eps[ii].st[ip].dam);

    /* increment ip counter */
    ip++;

    /* error detected, deallocate and exit */
    if(err != 0) goto dealocate;
      }
    }
  }

 dealocate:
  /*=== DEALLOCATION ===*/
  /* integration */
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  free(Np_x);
  free(ST);

  /* kinematics */
  free(Fr);
  free(Fr_I);
  free(F);
  free(F_I);
  free(C);
  free(C_I);
  free(E);

  /* physical quantities */
  free(grad_delP);
  free(grad_P);
  free(devSbar);
  free(sigma);
  free(almansi);

  if(STAB_DEBUG && err != 0){
    PGFEM_printf("error in %s [elem %ld]\n",__func__,ii);
  }

  return result;
}/* static double st_incr_elem () */
