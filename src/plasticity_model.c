/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 *  Sangmin Lee, University of Notre Dame, Notre Dame, IN, <slee43@nd.edu>
 */

#include "plasticity_model.h"
#include "constitutive_model.h"
#include "new_potentials.h"
#include "data_structure_c.h"
#include <string.h>
#include <math.h>

#include "CM.h"

#define DIM 3
#define TENSOR_LEN 9
#define TENSOR4_LEN 81


enum variable_names {
  VAR_L_n,
  VAR_L_np1,
  VAR_g_n,
  VAR_g_np1,
  VAR_gamma_dot_0, /* this and the rest are actually material parameters */
  VAR_gamma_dot_s,
  VAR_m,
  VAR_g0,
  VAR_G0,
  VAR_gs_0,
  VAR_w
};

enum tensor_names {
  TENSOR_Fn,
  TENSOR_pFn,
  TENSOR_Fnp1,
  TENSOR_pFnp1,
  TENSOR_tau,
  TENSOR_gamma_dot
};

Define_Matrix(double);

/**
 * Private structure for use exclusively with this model and
 * associated functions.
 */
typedef struct plasticity_ctx {
  double F[TENSOR_LEN];
  double dt; /* time increment */
} plasticity_ctx;

static double compute_bulk_mod(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static int plasticity_int_alg(Constitutive_model *m,
                              const void *ctx)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  Matrix_double Fnp1, pFnp1, Fe_n;
  Matrix_construct(double, Fnp1);
  Matrix_init_w_array(Fnp1, DIM, DIM, CTX->F);
  Matrix_construct_redim(double, pFnp1, DIM, DIM);
  Matrix_construct_redim(double, Fe_n, DIM, DIM);
  err += m->param->get_eFn(m,&Fe_n);

  /* NOTE: m->...pFnp1 is set in this function. The returned value is
     not needed. */
  err += plasticity_model_integration_ip(&pFnp1, m, &Fnp1, &Fe_n, CTX->dt);

  Matrix_cleanup(Fnp1);
  Matrix_cleanup(pFnp1);
  Matrix_cleanup(Fe_n);
  return err;
}

static int cp_compute_eC(const double * restrict eF,
                         double * restrict eC)
{
  int err = 0;
  memset(eC, 0, TENSOR_LEN * sizeof(*eC));
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      for (int k = 0; k < DIM; k++) {
        eC[idx_2(i,j)] += eF[idx_2(k,i)] * eF[idx_2(k,j)];
      }
    }
  }
  return err;
}

static int cp_dev_stress(const double *eC,
                         const HOMMAT *p_hmat,
                         double *Sdev)
{
  int err = 0;
  devStressFuncPtr Stress = getDevStressFunc(-1, p_hmat);
  Stress(eC, p_hmat, Sdev);
  return err;
}

static int cp_dev_tangent(const double *eC,
                          const HOMMAT *p_hmat,
                          double *Ldev)
{
  int err = 0;
  matStiffFuncPtr Tangent = getMatStiffFunc(-1, p_hmat);
  Tangent(eC, p_hmat, Ldev);
  return err;
}

static int plasticity_update(Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = (m->vars).Fs;
  double *state_var = (m->vars).state_vars[0].m_pdata;
  Matrix_AeqB(Fs[TENSOR_Fn], 1.0,Fs[TENSOR_Fnp1]);
  Matrix_AeqB(Fs[TENSOR_pFn],1.0,Fs[TENSOR_pFnp1]);
  state_var[VAR_g_n] = state_var[VAR_g_np1];
  state_var[VAR_L_n] = state_var[VAR_L_np1];
  return err;
}

static int plasticity_reset(Constitutive_model *m)
{
  int err = 0;
  Matrix(double) *Fs = (m->vars).Fs;
  double *state_var = (m->vars).state_vars[0].m_pdata;
  Matrix_AeqB(Fs[TENSOR_Fnp1], 1.0,Fs[TENSOR_Fn]);
  Matrix_AeqB(Fs[TENSOR_pFnp1],1.0,Fs[TENSOR_pFn]);
  state_var[VAR_g_np1] = state_var[VAR_g_n];
  state_var[VAR_L_np1] = state_var[VAR_L_n];
  return err;
}

static int plasticity_info(Model_var_info **info)
{
  *info = malloc(sizeof(**info));
  int Fno = 6;
  
  (*info)->n_Fs = Fno;
  (*info)->F_names = (char **)malloc(sizeof(char*)*Fno);
  for(int a=0; a<Fno; a++)
    (*info)->F_names[a] = (char *)malloc(sizeof(char)*1024);

  sprintf((*info)->F_names[TENSOR_Fn],        "Fn");
  sprintf((*info)->F_names[TENSOR_pFn],       "pFn");
  sprintf((*info)->F_names[TENSOR_Fnp1],      "Fnp1");
  sprintf((*info)->F_names[TENSOR_pFnp1],     "pFnp1");
  sprintf((*info)->F_names[TENSOR_tau],       "tau");
  sprintf((*info)->F_names[TENSOR_gamma_dot], "gamma_dot");
  
  int varno = 11;
  (*info)->n_vars = varno;
  (*info)->var_names = (char **)malloc(sizeof(char*)*varno);
  for(int a=0; a<varno; a++)
    (*info)->var_names[a] = (char *)malloc(sizeof(char)*1024);
  
  sprintf((*info)->var_names[VAR_L_n],        "L_n");
  sprintf((*info)->var_names[VAR_L_np1],      "L_np1");
  sprintf((*info)->var_names[VAR_g_n],        "g_n");
  sprintf((*info)->var_names[VAR_g_np1],      "g_np1");

  /* Parameters, temporarily stored as variables */
  sprintf((*info)->var_names[VAR_gamma_dot_0],"gamma_dot_0");
  sprintf((*info)->var_names[VAR_gamma_dot_s],"gamma_dot_s"); 
  sprintf((*info)->var_names[VAR_m],          "m");          
  sprintf((*info)->var_names[VAR_g0],         "g0");         
  sprintf((*info)->var_names[VAR_G0],         "G0");         
  sprintf((*info)->var_names[VAR_gs_0],       "gs_0");       
  sprintf((*info)->var_names[VAR_w],           "w");           

  return 0;
}

static int plasticity_get_pF(const Constitutive_model *m,
                             Matrix_double *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars.Fs[TENSOR_pFnp1]);
  return err;
}

static int plasticity_get_Fn(const Constitutive_model *m,
                             Matrix_double *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars.Fs[TENSOR_Fn]);
  return err;
}

static int plasticity_get_pFn(const Constitutive_model *m,
                              Matrix_double *F)
{
  int err = 0;
  Matrix_AeqB(*F,1.0,m->vars.Fs[TENSOR_pFn]);
  return err;
}

static int plasticity_get_eF(const Constitutive_model *m,
                             Matrix_double *F)
{
  int err = 0;
  Matrix_double invFp;
  Matrix_construct_redim(double,invFp,DIM,DIM);
  Matrix_inv(m->vars.Fs[TENSOR_pFnp1],invFp);
  Matrix_AxB(*F, 1.0, 0.0, m->vars.Fs[TENSOR_Fnp1], 0, invFp, 0);
  Matrix_cleanup(invFp);
  return err;
}

static int plasticity_get_eFn(const Constitutive_model *m,
                              Matrix_double *F)
{
  int err = 0;
  Matrix_double invFp;
  Matrix_construct_redim(double,invFp,DIM,DIM);
  Matrix_inv(m->vars.Fs[TENSOR_pFn],invFp);
  Matrix_AxB(*F, 1.0, 0.0, m->vars.Fs[TENSOR_Fn], 0, invFp, 0);
  Matrix_cleanup(invFp);
  return err;
}

static int plasticity_get_hardening(const Constitutive_model *m,
                                    double *var)
{
  int err = 0;
  *var = m->vars.state_vars->m_pdata[VAR_g_n];
  return err;
}

static int plasticity_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 Matrix_double *stress)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  double eC[TENSOR_LEN] = {};
  Matrix_double eFnp1;
  Matrix_construct_redim(double, eFnp1, DIM, DIM);
  err += plasticity_get_eF(m,&eFnp1);
  err += cp_compute_eC(eFnp1.m_pdata,eC);
  err += cp_dev_stress(eC,m->param->p_hmat,stress->m_pdata);
  return err;
}

static int plasticity_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *dudj)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
  Matrix_double eFnp1;
  Matrix_construct_redim(double, eFnp1, DIM, DIM);
  err += plasticity_get_eF(m,&eFnp1);
  double J = 0.0;
  Matrix_det(eFnp1, J);
  Pressure(J,m->param->p_hmat,dudj);
  return err;
}

static int plasticity_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  Matrix_double *tangent)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  double eC[TENSOR_LEN] = {};
  Matrix_double eFnp1;
  Matrix_construct_redim(double, eFnp1, DIM, DIM);
  err += plasticity_get_eF(m,&eFnp1);
  err += cp_compute_eC(eFnp1.m_pdata,eC);
  err += cp_dev_tangent(eC,m->param->p_hmat,tangent->m_pdata);
  return err;
}

static int plasticity_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *d2udj2)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  Matrix_double eFnp1;
  Matrix_construct_redim(double, eFnp1, DIM, DIM);
  err += plasticity_get_eF(m,&eFnp1);
  double J = 0.0;
  Matrix_det(eFnp1, J);
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
  D_Pressure(J,m->param->p_hmat,d2udj2);
  return err;
}

static int plasticity_compute_dMdu(const Constitutive_model *m,
                                   const void *ctx,
                                   const double *Grad_op,
                                   const int nne,
                                   const int ndofn,
                                   double *dM_du)
{
  int err = 0;
  /* The existing function compute_dMdu takes Grad_op for a given
     alpha,beta pair and computes the corresponding dMdu. We will
     generate the required inputs for this function and call it
     alpha*beta times. This should be improved */
  const plasticity_ctx *CTX = ctx;
  Matrix_double eFn, eFnp1;

  /* extract the elastic deformations at (n) and (n + 1) */
  Matrix_construct_redim(double, eFn, DIM, DIM);
  Matrix_construct_redim(double, eFnp1, DIM, DIM);
  err += m->param->get_eFn(m,&eFn);
  err += m->param->get_eF(m,&eFnp1);

  /* compute M at n+1, again, this information is in CM */
  Matrix_double M;
  {
    Matrix_double pFnp1_I;
    Matrix_construct_redim(double, pFnp1_I, DIM, DIM);
    Matrix_construct_redim(double, M, DIM, DIM);
    Matrix_inv(m->vars.Fs[TENSOR_pFnp1], pFnp1_I);
    Matrix_AxB(M, 1.0, 0.0, m->vars.Fs[TENSOR_pFn], 0, pFnp1_I, 0);
    Matrix_cleanup(pFnp1_I);
  }

  /* compute the elastic stress/tangent (S, L) using
     constitutive_model_update_elasticity. This needs to become part
     of the CM interface */
  Matrix_double L,S;
  Matrix_construct_redim(double, L, TENSOR4_LEN, 1);
  Matrix_construct_redim(double, S, DIM, DIM);
  err += constitutive_model_update_elasticity(m, &eFnp1, CTX->dt, &L, &S, 1);

  /* make successive calls to compute_dMdu for each node/dof. To avoid
     copying, I am abusing access to the internal data structure of
     the Matrix structure. */
  Matrix_double dMdu_ab, Grad_op_ab;
  Matrix_construct(double, dMdu_ab);
  Matrix_construct(double, Grad_op_ab);
  for (int a = 0; a < nne; a++) {
    for(int b = 0; b < ndofn; b++) {
      int idx_ab = idx_4_gen(a,b,0,0,nne,ndofn,DIM,DIM);
      /* reset dimensions of the matrix objects and set pointer */
      dMdu_ab.m_row = DIM;
      dMdu_ab.m_col = DIM;
      dMdu_ab.m_pdata = dM_du + idx_ab;

      /* need to copy Grad_op due to const qualifier */
      Matrix_init_w_array(Grad_op_ab, DIM, DIM, Grad_op + idx_ab);

      /* call to compute_dMdu */
      err += compute_dMdu(m, &dMdu_ab, &Grad_op_ab, &eFn, &eFnp1, &M,
                          &S, &L, CTX->dt);
    }
  }

  /* clean up */
  Matrix_cleanup(M);
  Matrix_cleanup(S);
  Matrix_cleanup(L);
  Matrix_cleanup(eFn);
  Matrix_cleanup(eFnp1);

  return err;
}

int plasticity_model_initialize(Model_parameters *p)
{
  int err = 0;

  /* set functions */
  p->integration_algorithm = plasticity_int_alg;
  p->compute_dev_stress = plasticity_dev_stress;
  p->compute_dudj = plasticity_dudj;
  p->compute_dev_tangent = plasticity_dev_tangent;
  p->compute_d2udj2 = plasticity_d2udj2;
  p->update_state_vars = plasticity_update;
  p->reset_state_vars = plasticity_reset;
  p->get_var_info = plasticity_info;
  p->get_Fn = plasticity_get_Fn;
  p->get_pF = plasticity_get_pF;
  p->get_pFn = plasticity_get_pFn;
  p->get_eF = plasticity_get_eF;
  p->get_eFn = plasticity_get_eFn;
  p->get_hardening = plasticity_get_hardening;
  p->destroy_ctx = plasticity_model_ctx_destroy;
  p->compute_dMdu = plasticity_compute_dMdu;

  return err;
}

int plasticity_model_ctx_build(void **ctx,
                               const double *F,
                               const double dt)
{
  int err = 0;
  plasticity_ctx *t_ctx = malloc(sizeof(plasticity_ctx));

  /* copy data into context */
  memcpy(t_ctx->F, F, TENSOR_LEN * sizeof(*F));
  t_ctx->dt = dt;

  /* assign handle */
  *ctx = t_ctx;
  return err;
}

int plasticity_model_ctx_destroy(void **ctx)
{
  int err = 0;
  plasticity_ctx *t_ctx = *ctx;
  /* invalidate handle */
  *ctx = NULL;

  /* there are no internal pointers */

  /* free object memory */
  free(t_ctx);
  return err;
}

static int compute_P_alpha(const Constitutive_model *m,
                           const int alpha,
                           Matrix(double) *Pa)
{
  const double * restrict P_sys = ((m->param)->Psys)->m_pdata;
  Matrix_redim(*Pa, DIM,DIM);
  for(int a=0; a<DIM; a++)
  {
    for(int b=0;b<DIM;b++)
      Mat_v(*Pa,a+1,b+1) = P_sys[Index_3D(alpha,a,b,DIM,DIM)];
  }
  return 0;
}

static int compute_C_D_alpha(const Constitutive_model *m,
                             Matrix(double) *aC,
                             Matrix(double) *aD,
                             const Matrix(double) *eFn,
                             const Matrix(double) *eFnp1,
                             const Matrix(double) *M,
                             const Matrix(double) *Pa,
                             const Matrix(double) *S,
                             const Matrix(double) *L,
                             const Matrix(double) *C)
{
  Matrix_init(*aC, 0.0);  
  Matrix_init(*aD, 0.0);
  
  Matrix(double) LC, AA, CAA, eFnp1AA, eFnp1AAMT, MI;
  Matrix_construct_redim(double,LC,       DIM,DIM);
  Matrix_construct_redim(double,AA,       DIM,DIM);
  Matrix_construct_redim(double,CAA,      DIM,DIM);
  Matrix_construct_redim(double,eFnp1AA,  DIM,DIM);
  Matrix_construct_redim(double,eFnp1AAMT,DIM,DIM);    
  Matrix_construct_redim(double,MI,       DIM,DIM);      

  Matrix_AxB(AA,1.0,0.0,*Pa,0,*S,0);   // AA = Pa*S
  Matrix_AxB(AA,1.0,1.0,*S,0,*Pa,1);   // AA = AA + S*Pa' 
       
  Matrix_Tns4_dd_Tns2(LC, *L, *C);     // LC = L:C
  Matrix_AxB(AA,1.0,1.0,LC,0,*Pa,0);   // AA = AA + L:C*Pa

  Matrix_inv(*M,MI);
  Matrix_AxB(CAA,1.0,0.0,*C,0,AA,0);   
  Matrix_AxB(*aC,1.0,0.0,MI,1,CAA,0);

  Matrix_AxB(eFnp1AA,1.0,0.0,*eFnp1,0,AA,0);
  Matrix_AxB(eFnp1AAMT,1.0,0.0,eFnp1AA,0,*M,1);
  Matrix_AxB(*aD,1.0,0.0,eFnp1AAMT,0,*eFn,1);
  
  Matrix_cleanup(LC);
  Matrix_cleanup(AA);
  Matrix_cleanup(CAA);
  Matrix_cleanup(eFnp1AA);
  Matrix_cleanup(eFnp1AAMT);       
  Matrix_cleanup(MI);          
  return 0;
}

int compute_dMdu(const Constitutive_model *m,
                 Matrix(double) *dMdu,
                 const Matrix(double) *Grad_du,
                 const Matrix(double) *eFn,
                 const Matrix(double) *eFnp1,
                 const Matrix(double) *M,
                 const Matrix(double) *S,
                 const Matrix(double) *L,
                 const double dt)
{
  // compute dMdu:U = -grad(du):B
  // Grad_du = Grad(du)

  const double *state_var = (m->vars).state_vars[0].m_pdata;
  
  const int N_SYS          = (m->param)->N_SYS;
  const double g_n         = state_var[VAR_g_n];
  const double g_np1       = state_var[VAR_g_np1];
  const double gamma_dot_0 = state_var[VAR_gamma_dot_0];
  const double gamma_dot_s = state_var[VAR_gamma_dot_s];
  const double mm          = state_var[VAR_m];
  const double g0          = state_var[VAR_g0];
  const double G0          = state_var[VAR_G0];
  const double gs_0        = state_var[VAR_gs_0];
  const double w           = state_var[VAR_w];

  const double *tau_np1 = (m->vars).Fs[TENSOR_tau].m_pdata;
  const double *gamma_dots = (m->vars).Fs[TENSOR_gamma_dot].m_pdata;

  // --------------> define variables
  Matrix(double) C;
  Matrix_construct_redim(double, C, DIM,DIM);
  Matrix_AxB(C,1.0,0.0,*eFnp1,1,*eFnp1,0);
  
  Matrix(double) U,UI,II,B,aCxPa,CxP,aDxPa,DxP;
  
  Matrix_construct_redim(double, U,     TENSOR4_LEN,1); // 3x3x3x3 tensor
  Matrix_construct_redim(double, UI,    TENSOR_LEN,TENSOR_LEN);
  Matrix_construct_redim(double, II,    TENSOR4_LEN,1);
  Matrix_construct_init( double, B,     TENSOR4_LEN,1,0.0);
  Matrix_construct_redim(double, aCxPa, TENSOR4_LEN,1);
  Matrix_construct_redim(double, CxP,   TENSOR4_LEN,1);
  Matrix_construct_redim(double, aDxPa, TENSOR4_LEN,1);
  Matrix_construct_redim(double, DxP,   TENSOR4_LEN,1);  

  Matrix(double) aC,Pa,aD,sum_aC,sum_Pa,sum_aD;
  Matrix_construct_redim(double, aC, DIM,DIM);  
  Matrix_construct_redim(double, Pa, DIM,DIM);
  Matrix_construct_redim(double, aD, DIM,DIM);
        
  Matrix_construct_init( double, sum_aC, DIM,DIM,0.0);
  Matrix_construct_init( double, sum_Pa, DIM,DIM,0.0);
  Matrix_construct_init( double, sum_aD, DIM,DIM,0.0);  
  
  Matrix_Tns4_eye(II);
  Matrix_AeqB(U, 1.0/dt, II);

  // <-------------- define variables  
  
  double gamma_dot = 0.0;  
  for(int a = 0; a<N_SYS; a++)
    gamma_dot += fabs(gamma_dots[a]);

  double gm_gms   = gamma_dot/gamma_dot_s;
  double sign_gm_gms = (gm_gms < 0) ? -1.0 : 1.0;
//vvvvvvv need to verify
  double R3 = 0.0;
  double gs_np1 = 0.0;  
  if(fabs(gm_gms)>1.0e-15)
  {
    R3 = gs_0*w/gamma_dot_s*sign_gm_gms*pow(fabs(gm_gms), w-1.0);
    gs_np1 = gs_0*pow(fabs(gm_gms),w);
  }   
//^^^^^^^ need to verify

  double AA = R3*gamma_dot*(g_n - g0 + dt*G0*gamma_dot) + gs_np1*(gs_np1 - g0 - g_n) + g0*g_n;
  double BB = gs_np1 - g0  - dt*G0*gamma_dot;
  double R4 = dt*G0*AA/BB/BB;
  
  double sum_1gm1gm = 0.0;
  
  for(int a = 0; a<N_SYS; a++)
  {
    double drdtau = gamma_dot_0/mm/g_np1*pow(fabs(tau_np1[a]/g_np1), 1.0/mm - 1.0);
    double drdg   = -drdtau*tau_np1[a]/g_np1;

    double R2_a = ((gamma_dots[a] < 0) ? -1.0 : 1.0)*drdtau;
    sum_1gm1gm += ((gamma_dots[a] < 0) ? -1.0 : 1.0)*drdg;  
        
    compute_P_alpha(m,a,&Pa);
    compute_C_D_alpha(m,&aC, &aD,eFn,eFnp1,M,&Pa,S,L,&C);
    
    Matrix_AOxB(aCxPa, aC, Pa);
    Matrix_AOxB(aDxPa, aD, Pa);    
    
    Matrix_AplusB(sum_aC, 1.0, sum_aC, R2_a, aC);
    Matrix_AplusB(sum_Pa, 1.0, sum_Pa, drdg, Pa);
    Matrix_AplusB(sum_aD, 1.0, sum_aD, R2_a, aD);
        
    Matrix_AOxB(aCxPa, aC, Pa);
    Matrix_AOxB(aDxPa, aD, Pa);
    
    
    Matrix_AplusB(U, 1.0, U, drdtau, aCxPa);    
    Matrix_AplusB(B, 1.0, B, drdtau, aDxPa);            
  } 
  
  double R1 = R4/(1.0-R4*sum_1gm1gm);
  
  Matrix_AOxB(CxP, sum_aC, sum_Pa);
  Matrix_AOxB(DxP, sum_aD, sum_Pa);
  Matrix_AplusB(U, 1.0, U, R1, CxP);
  Matrix_AplusB(B, 1.0, B, R1, DxP);

/*  
  Matrix(double) V;
  Matrix_construct_init(double, V, DIM,DIM,0.0);
  Matrix_Tns2_dd_Tns4(V,*Grad_du,B);
  Matrix_AeqB(V, -1.0, V);
  
  Matrix_Tns4_mat_9x9(U);
  Matrix_inv(U, UI);
  
  dMdu->m_row = 1;
  dMdu->m_col = TENSOR_LEN;
  Matrix_Mat2Vec(V);        
  Matrix_AxB(*dMdu,1.0,0.0,V,1,UI,0);   
  Matrix_Vec2Mat(*dMdu,DIM,DIM);*/
  
  Matrix(double) V;
  Matrix_construct_redim(double, V, DIM,DIM);
  Matrix_Tns4_dd_Tns2(V,B,*Grad_du);
  Matrix_AeqB(V, -1.0, V);
  
  Matrix_Tns4_mat_9x9(U);
  Matrix_inv(U, UI);

  Matrix_Mat2Vec(V);
  Matrix_Mat2Vec(*dMdu);
  Matrix_AxB(*dMdu,1.0,0.0,UI,0,V,0); 
  Matrix_Vec2Mat(*dMdu,DIM,DIM);

  // clear variables
  Matrix_cleanup(U);
  Matrix_cleanup(UI);
  Matrix_cleanup(II);
  Matrix_cleanup(B);
  Matrix_cleanup(aCxPa);
  Matrix_cleanup(CxP);
  Matrix_cleanup(aDxPa);
  Matrix_cleanup(DxP);  
  Matrix_cleanup(aC);  
  Matrix_cleanup(Pa);
  Matrix_cleanup(aD);
  Matrix_cleanup(sum_aC);
  Matrix_cleanup(sum_Pa);
  Matrix_cleanup(sum_aD);   
  Matrix_cleanup(V);
  
  return 0;
}

int plasticity_model_slip_system(Matrix_double *P)
{
  int err = 0;
  int N_SYS = 12; // depended on slip system

  int j_max = DIM;  
  Matrix_redim(*P, N_SYS*j_max*j_max, 1);
  double *P_sys = P->m_pdata;
  
  for (int k = 0; k<N_SYS; k++)
  {
    for (int j = 0; j<j_max; j++)
    {
      for (int i = 0; i<j_max; i++)
        P_sys[Index_3D(k,j,i,j_max,j_max)] = 0.0;
    }
    P_sys_sp(k, P_sys);
  }
  
  return N_SYS;  
}

/**
 * Function to set the values of Fnp1 under the stand-alone CM
 * integrator. This is a band-aid since the integrator does not
 * directly access/modify the Constitutive_model data.
 */
static void CM_set_Fnp1(MaterialProperties *Props,
                        double **Fnp1)
{
  assert(Props->cm != NULL);
  double *F = Props->cm->vars.Fs[TENSOR_Fnp1].m_pdata;
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      F[idx_2(i,j)] = Fnp1[i][j];
    }
  }
}

/**
 * Function to set the values of pFnp1 under the stand-alone CM
 * integrator. This is a band-aid since the integrator does not
 * directly access/modify the Constitutive_model data.
 */
static void CM_set_pFnp1(MaterialProperties *Props,
                         double **pFnp1)
{
  assert(Props->cm != NULL);
  double *F = Props->cm->vars.Fs[TENSOR_pFnp1].m_pdata;
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      F[idx_2(i,j)] = pFnp1[i][j];
    }
  }
}

void elastic_stress(MaterialProperties *Props, double *F, double *S)
{
  Matrix(double) Se, Fe;
  Matrix_construct_redim(double,Se,DIM,DIM);
  Matrix_construct(double, Fe);   
  Matrix_init_w_array(Fe,DIM,DIM,F);  
  constitutive_model_update_elasticity(Props->cm, &Fe, 0.0, NULL, &Se, 0);

  for(int a=0;a<TENSOR_LEN;a++)
    S[a] = Se.m_pdata[a];
    
  Matrix_cleanup(Se);
  Matrix_cleanup(Fe);   
}

static void elastic_tangent(MaterialProperties *Props, double *F, double *L)
{
  Matrix(double) Se, Fe, Le;
  Matrix_construct_redim(double,Se,DIM,DIM); 
  Matrix_construct(double, Fe);  
  Matrix_construct_redim(double,Le,TENSOR4_LEN,1);     
  Matrix_init_w_array(Fe,DIM,DIM,F);  
  constitutive_model_update_elasticity(Props->cm, &Fe, 0.0, &Le, &Se, 1);

  memcpy(L, Le.m_pdata, TENSOR4_LEN * sizeof(*L));
    
  Matrix_cleanup(Se);
  Matrix_cleanup(Le);
  Matrix_cleanup(Fe);
}

static int plasticity_model_staggered_NR(Matrix(double) *pFnp1,
                                         double *g_np1,
                                         double *L_np1,
                                         const Matrix(double) *Fnp1,
                                         const Matrix(double) *eFn,
                                         const Matrix(double) *pFn,
                                         MaterialProperties *Props,
                                         MaterialParameters *Param,
                                         MaterialStructure *Struc,
                                         SolverInformation *Solver,
                                         Constitutive_model *m,
                                         const double dt,
                                         const double g_n)
{
  const int print_step_info = 0;   
  const int max_itr = 10;   
  const double computer_zero = 1.0e-15;  
  const double tol_g = 1.0e-6;   
  
  int err = 0;
  Matrix(double) *Fs;  
  enum {eFn_I, pFnp1_k, pFnp1_I, eFnp1, PK2, Fend};
  
  const int N_SYS = (m->param)->N_SYS;
  double *P_sys = ((m->param)->Psys)->m_pdata;
  Matrix_double *CM_pFnp1 = &(m->vars.Fs[TENSOR_pFnp1]);
  Fs = (Matrix(double) *) malloc(Fend*sizeof(Matrix(double)));
  
  for(int a=0; a<Fend; a++)
    Matrix_construct_redim(double, Fs[a],DIM,DIM);  
  
  Matrix(double) Tau_Array, gamma_RateArray;
  Matrix_construct_redim(double,Tau_Array      , N_SYS,1);
  Matrix_construct_redim(double,gamma_RateArray, N_SYS,1);  
      
  double L_np1_k   = 0.0;
  double L_np1_kp1 = 0.0;
  
  double g_np1_k   = g_n;
  double g_np1_kp1 = g_n;
  Matrix_inv(*eFn, Fs[eFn_I]);
  Matrix_AxB(Fs[pFnp1_k],1.0,0.0,Fs[eFn_I],0,*Fnp1,0);  



  double err_g_0 = computer_zero;   

  for(int k=0; k<max_itr; k++)
  {
    Matrix_init(*pFnp1, 0.0);
        *L_np1 = 0.0;

        /* update the internal datastructure for pF */
        Matrix_AeqB(*CM_pFnp1,1.0,Fs[pFnp1_k]);

        Staggered_NewtonRapson_Testing_sp(Props,Param,Struc,Solver,
                                          0.0,dt,P_sys,g_np1_kp1,pFn->m_pdata,Fnp1->m_pdata,
                                          L_np1_k,&L_np1_kp1,Fs[pFnp1_k].m_pdata,pFnp1->m_pdata);

        /* update the internal datastructure for pF */
        Matrix_AeqB(*CM_pFnp1,1.0,*pFnp1);

        Matrix_inv(*pFnp1, Fs[pFnp1_I]);
     Matrix_AxB(Fs[eFnp1],1.0,0.0,*Fnp1,0,Fs[pFnp1_I],0);         
    elastic_stress(Props,Fs[eFnp1].m_pdata,Fs[PK2].m_pdata);
    Matrix_init(Tau_Array, 0.0);
    Matrix_init(gamma_RateArray, 0.0);
    double gmdot_tmp = 0.0;
    for (int q = 0; q<N_SYS; q++)
    {
      double tau_q = Tau_Rhs_sp(q, P_sys, Fs[eFnp1].m_pdata, Fs[PK2].m_pdata);
      Vec_v(Tau_Array,q+1) = tau_q;
      Vec_v(gamma_RateArray,q+1) = gamma_Rate_PL(Param,g_np1_kp1,tau_q);
      gmdot_tmp += fabs(Vec_v(gamma_RateArray,q+1));
    }

    g_np1_kp1 = g_Rate_VK_implicit(Param,Struc,g_n,dt,gamma_RateArray.m_pdata);
    double err_g = sqrt((g_np1_kp1-g_np1_k)*(g_np1_kp1-g_np1_k));

    if(k==0 && err_g>computer_zero)  
      err_g_0 = err_g;
    

    g_np1_k = g_np1_kp1; 
    L_np1_k = L_np1_kp1;
    if(print_step_info)
      printf("staggered iteration: %d, err of g: %e/%e = %e\n", k, err_g, err_g_0, err_g/err_g_0);

    if((err_g/err_g_0)<tol_g)
      break;
    
  }
  *g_np1 = g_np1_kp1;
  *L_np1 = L_np1_kp1;
  
  for(int a=0; a<Fend; a++)
    Matrix_cleanup(Fs[a]);        
  
  free(Fs);
  Matrix_cleanup(Tau_Array);
  Matrix_cleanup(gamma_RateArray);  
  
  return err;    
}

int plasticity_model_integration_ip(Matrix_double *pFnp1,
                                    Constitutive_model *m,
                                    const Matrix_double *Fnp1,
                                    const Matrix_double *Fe_n,
                                    const double dt)
{
  const HOMMAT *hmat = m->param->p_hmat;
  double *state_var = (m->vars).state_vars[0].m_pdata;
  Matrix(double) *Fs = (m->vars).Fs;
  
  int err = 0;
  int N_SYS = (m->param)->N_SYS;
  double *P_sys = ((m->param)->Psys)->m_pdata;
  double *pFn    = Fs[TENSOR_pFn].m_pdata;
    
  /*--------Simulation_Settings--------*/
  SolverInformation Solver;
  Solver.Solver_Type = IMPLICIT;/*ImplicitWithIncompresibility_or_IC*/
  Solver.AMatrix_Size = 100;
  Solver.BVector_Size = 10;
  Solver.NR_ML_MAX = 100;
  Solver.NR_G_MAX = 100;
  Solver.SNR_MAX = 100;
  Solver.Fp_TOL = 10e-12;
  Solver.L_TOL = 10e-12;
  Solver.g_TOL = 10e-12;
  
  /*--------MaterialStructure_Settings--------*/
  MaterialStructure Struc;
  Struc.NUM_GRAIN = 1;
  Struc.N_SYS = N_SYS;
  Struc.Structure_Type = FCC;

  /*--------MaterialProperties_Settings--------*/
  MaterialProperties Props;
  Props.Lame_I =hmat->nu*hmat->E/(1.0+hmat->nu)/(1.0-2.0*hmat->nu);; 
  Props.Lame_II =hmat->G;
  Props.Modulus_Elastic = hmat->E;
  Props.Modulus_Shear = hmat->G;
  Props.Poissons_Ratio= hmat->nu;
  Props.Modulus_Bulk = (2.0*hmat->G*(1.0+hmat->nu))/(3.0*(1.0 - 2.0*hmat->nu));
  Props.use_hyperelastic = 1;
  Props.cm = m;
  Props.compute_elastic_stress = elastic_stress;
  Props.compute_elastic_tangent = elastic_tangent;
  Props.set_cm_Fnp1 = CM_set_Fnp1;
  Props.set_cm_pFnp1 = CM_set_pFnp1;

  /*--------MaterialParameters_Settings--------*/
  MaterialParameters Params;
  Params.Model_Type = PL_VK;
  Params.Parameters_Count = 7;
  Params.gam0dot = state_var[VAR_gamma_dot_0];
  Params.m_matl  = state_var[VAR_m];
  Params.gs0     = state_var[VAR_gs_0];
  Params.gamsdot = state_var[VAR_gamma_dot_s];
  Params.w       = state_var[VAR_w];
  Params.G0      = state_var[VAR_G0];
  Params.g0      = state_var[VAR_g0];  
  
  int j_max = DIM;

  Matrix(double) Fe_I;
  Matrix_construct_redim(double,Fe_I,DIM,DIM);
  Matrix_inv(*Fe_n, Fe_I);

  double g_n = state_var[VAR_g_n];
  double g_np1 = g_n;
  double L_np1 = state_var[VAR_L_n];

  /* Update the internally stored total F, pF */
  Matrix_AeqB(Fs[TENSOR_Fnp1], 1.0, *Fnp1);
  Matrix_AeqB(Fs[TENSOR_pFnp1], 1.0, *pFnp1);

  plasticity_model_staggered_NR(pFnp1, &g_np1, &L_np1,
                                Fnp1, Fe_n, (Fs+TENSOR_pFn),
                                &Props, &Params, &Struc, &Solver, m, dt, g_n);

  /* update pF after the integration algorithm */
  Matrix_AeqB(Fs[TENSOR_pFnp1], 1.0, *pFnp1);

  state_var[VAR_g_np1] =  g_np1;
  state_var[VAR_L_np1] =  L_np1;

  Matrix(double) S_n, pFnp1_I, eFnp1;
  Matrix_construct_redim(double,S_n,DIM,DIM);       
  Matrix_construct_redim(double,pFnp1_I,DIM,DIM);       
  Matrix_construct_redim(double,eFnp1,DIM,DIM);           
  
  Matrix_inv(*pFnp1,pFnp1_I);
  Matrix_AxB(eFnp1,1.0,0.0,*Fnp1,0,pFnp1_I,0);
  
  elastic_stress(&Props, eFnp1.m_pdata, S_n.m_pdata);

  double gamma_dot = 0.0;    //////////////////////// compute g_np1
  for (int k = 0; k<N_SYS; k++)
  {
    double tau_k = Tau_Rhs_sp(k, P_sys, eFnp1.m_pdata, S_n.m_pdata);
    Vec_v(Fs[TENSOR_tau],k+1) = tau_k;
    Vec_v(Fs[TENSOR_gamma_dot],k+1) = gamma_Rate_PL(&Params,g_n,tau_k);
    gamma_dot += fabs(Vec_v(Fs[TENSOR_gamma_dot],k+1)); //////////////////////// compute g_np1
  }

  Matrix_cleanup(Fe_I);
  Matrix_cleanup(S_n);
  Matrix_cleanup(eFnp1);
  Matrix_cleanup(pFnp1_I);
  
  return err;
}


int plasticity_model_integration_ip_(Matrix_double *pFnp1, Constitutive_model *m, Matrix_double *Fnp1, Matrix_double *Fe_n, double dt)
{
  const HOMMAT *hmat = m->param->p_hmat;
  double *state_var = (m->vars).state_vars[0].m_pdata;
  Matrix(double) *Fs = (m->vars).Fs;
  
  int err = 0;
  int N_SYS = (m->param)->N_SYS;
  double *P_sys = ((m->param)->Psys)->m_pdata;
  double *pFn    = Fs[TENSOR_pFn].m_pdata;
    
  /*--------Simulation_Settings--------*/
  SolverInformation Solver;
  Solver.Solver_Type = IMPLICIT;/*ImplicitWithIncompresibility_or_IC*/
  Solver.AMatrix_Size = 100;
  Solver.BVector_Size = 10;
  Solver.NR_ML_MAX = 100;
  Solver.NR_G_MAX = 100;
  Solver.SNR_MAX = 100;
  Solver.Fp_TOL = 10e-12;
  Solver.L_TOL = 10e-12;
  Solver.g_TOL = 10e-12;
  
  /*--------MaterialStructure_Settings--------*/
  MaterialStructure Struc;
  Struc.NUM_GRAIN = 1;
  Struc.N_SYS = N_SYS;
  Struc.Structure_Type = FCC;

  /*--------MaterialProperties_Settings--------*/
  MaterialProperties Props;
  Props.Lame_I =hmat->nu*hmat->E/(1.0+hmat->nu)/(1.0-2.0*hmat->nu);; 
  Props.Lame_II =hmat->G;
  Props.Modulus_Elastic = hmat->E;
  Props.Modulus_Shear = hmat->G;
  Props.Poissons_Ratio= hmat->nu;
  Props.Modulus_Bulk = (2.0*hmat->G*(1.0+hmat->nu))/(3.0*(1.0 - 2.0*hmat->nu));
  Props.use_hyperelastic = 1;
  Props.cm = m;
  Props.compute_elastic_stress = elastic_stress;
  Props.compute_elastic_tangent = elastic_tangent;
  Props.set_cm_Fnp1 = CM_set_Fnp1;
  Props.set_cm_pFnp1 = CM_set_pFnp1;

  /*--------MaterialParameters_Settings--------*/
  MaterialParameters Params;
  Params.Model_Type = PL_VK;
  Params.Parameters_Count = 7;
  Params.gam0dot = state_var[VAR_gamma_dot_0];
  Params.m_matl  = state_var[VAR_m];
  Params.gs0     = state_var[VAR_gs_0];
  Params.gamsdot = state_var[VAR_gamma_dot_s];
  Params.w       = state_var[VAR_w];
  Params.G0      = state_var[VAR_G0];
  Params.g0      = state_var[VAR_g0];  
  
  int j_max = DIM;

  Matrix(double) Fe_I,Fp_np1k;
  
  Matrix_construct_redim(double,Fe_I,DIM,DIM);
  Matrix_construct_redim(double,Fp_np1k,DIM,DIM);

  double g_n = state_var[VAR_g_n];
  
  double L_np1k = state_var[VAR_L_n];
  if(fabs(L_np1k)<1.0e-12)
    L_np1k = 0.001;
    
  double L_np1 = 0.0;
  Matrix_inv(*Fe_n, Fe_I);
  Matrix_AxB(Fp_np1k,1.0,0.0,Fe_I,0,*Fnp1,0);

  Staggered_NewtonRapson_Testing_sp(&Props,&Params,&Struc,&Solver,0.0,dt,
                                    P_sys,g_n,pFn,Fnp1->m_pdata,L_np1k,&L_np1,Fp_np1k.m_pdata,pFnp1->m_pdata);

  Matrix(double) S_n, pFnp1_I, eFnp1;
  Matrix_construct_redim(double,S_n,DIM,DIM);       
  Matrix_construct_redim(double,pFnp1_I,DIM,DIM);       
  Matrix_construct_redim(double,eFnp1,DIM,DIM);           

  Matrix_inv(*pFnp1,pFnp1_I);
  Matrix_AxB(eFnp1,1.0,0.0,*Fnp1,0,pFnp1_I,0);

  elastic_stress(&Props, eFnp1.m_pdata, S_n.m_pdata);

  double gamma_dot = 0.0;    //////////////////////// compute g_np1
  for (int k = 0; k<N_SYS; k++)
    {
      double tau_k = Tau_Rhs_sp(k, P_sys, eFnp1.m_pdata, S_n.m_pdata);
      Vec_v(Fs[TENSOR_tau],k+1) = tau_k;
      Vec_v(Fs[TENSOR_gamma_dot],k+1) = gamma_Rate_PL(&Params,g_n,tau_k);
      gamma_dot += fabs(Vec_v(Fs[TENSOR_gamma_dot],k+1)); //////////////////////// compute g_np1
    }
  double g_Rhs = g_Rate_VK(&Params,&Struc,g_n, Fs[TENSOR_gamma_dot].m_pdata);
  state_var[VAR_g_np1] =  g_n + dt*g_Rhs;
  state_var[VAR_L_np1] =  L_np1;
  
  
  //////////////////////// compute g_np1
  double gm_gms   = gamma_dot/Params.gamsdot;
  double gs_np1 = 0.0;   
  if(fabs(gm_gms)>1.0e-15)
    gs_np1 = Params.gs0*pow(fabs(gm_gms),Params.w);
  
    
  double gg = ((gs_np1-Params.g0)*g_n + dt*Params.G0*gs_np1*gamma_dot)/(gs_np1 - Params.g0 + dt*Params.G0*gamma_dot);
//  printf("%e %e %e\n", g_n + dt*g_Rhs, gg, g_n + fabs(dt*g_Rhs - gg));  
  state_var[VAR_g_np1] = gg;
  //////////////////////// compute g_np1  

  Matrix_cleanup(Fe_I);
  Matrix_cleanup(Fp_np1k);
  Matrix_cleanup(S_n);
  Matrix_cleanup(eFnp1);
  Matrix_cleanup(pFnp1_I);    
  
  Matrix_AeqB(Fs[TENSOR_pFnp1], 1.0, *pFnp1);
  Matrix_AeqB(Fs[TENSOR_Fnp1], 1.0, *Fnp1);
  
  return err;    
}


int plasticity_model_read_parameters(Constitutive_model *m)
{
  int err = 0;
  double *state_var = (m->vars).state_vars[0].m_pdata;

  /* set parameter values -- these remain unchanged */
  state_var[VAR_gamma_dot_0] = 1.0;
  state_var[VAR_gamma_dot_s] = 50.0e+9;  
  state_var[VAR_m]           = 0.05;
  state_var[VAR_g0]          = 210.0;  
  state_var[VAR_G0]          = 200.0;    
  state_var[VAR_gs_0]        = 330.0;
  state_var[VAR_w]           = 0.005;

  /* set state variables to initial values */
  state_var[VAR_g_n] = state_var[VAR_g0];
  state_var[VAR_g_np1] = state_var[VAR_g0];    
  state_var[VAR_L_n] = 0.0;  
  state_var[VAR_L_np1] = 0.0;

  /* set the dimensions for the tau and gamma_dot varables */
  Matrix(double) *Fs = (m->vars).Fs;
  const int N_SYS = (m->param)->N_SYS;
  Matrix_redim(Fs[TENSOR_tau],N_SYS, 1);
  Matrix_redim(Fs[TENSOR_gamma_dot],N_SYS, 1);

  /* intitialize to zeros */
  Matrix_init(Fs[TENSOR_tau], 0.0);
  Matrix_init(Fs[TENSOR_gamma_dot], 0.0);
  return err;
}

int test_set_CM_interface_values(MaterialProperties *Props, MaterialParameters *Param, 
                                 MaterialStructure *Struc, SolverInformation *Solver,
                                 Constitutive_model *m, const HOMMAT *hmat, int N_SYS)
{
  int err = 0;
  /*--------Simulation_Settings--------*/
  Solver->Solver_Type = IMPLICIT;/*ImplicitWithIncompresibility_or_IC*/
  Solver->AMatrix_Size = 100;
  Solver->BVector_Size = 10;
  Solver->NR_ML_MAX = 100;
  Solver->NR_G_MAX = 100;
  Solver->SNR_MAX = 100;
  Solver->Fp_TOL = 10e-12;
  Solver->L_TOL = 10e-12;
  Solver->g_TOL = 10e-12;

  /*--------MaterialStructure_Settings--------*/
  Struc->NUM_GRAIN = 1;
  Struc->N_SYS = N_SYS;
  Struc->Structure_Type = FCC;

  // 7.160177e+04 13050 0 26100.0 0 0 3.716814e-01 0 0 1.0 1.0 1.0 1.7e+03 1 2
  /*--------MaterialProperties_Settings--------*/
  Props->Lame_I =hmat->nu*hmat->E/(1.0+hmat->nu)/(1.0-2.0*hmat->nu);; 
  Props->Lame_II =hmat->G;
  Props->Modulus_Elastic = hmat->E;
  Props->Modulus_Shear = hmat->G;
  Props->Poissons_Ratio= hmat->nu;
  Props->Modulus_Bulk = (2.0*hmat->G*(1.0+hmat->nu))/(3.0*(1.0 - 2.0*hmat->nu));

  Props->use_hyperelastic = 1;
  Props->cm = m;
  Props->compute_elastic_stress = elastic_stress;
  Props->compute_elastic_tangent = elastic_tangent;
  Props->set_cm_Fnp1 = CM_set_Fnp1;
  Props->set_cm_pFnp1 = CM_set_pFnp1;    

  /*--------MaterialParameters_Settings--------*/
  Param->Model_Type = PL_VK;
  Param->Parameters_Count = 7;
  Param->gam0dot = 1.0;
  Param->m_matl = 0.05;
  Param->gs0 = 330.0;
  Param->gamsdot = 50000000000.0;
  Param->w = 0.005;
  Param->G0 = 200.0;
  Param->g0 = 210.0;
  return err;
}

void test_load_type(Matrix(double) *L, int type, double Load_History, double t)
{
  Matrix_init(*L,0.0);
  switch(type)
  {
    case UNIAXIAL_TENSION:
      // Tension      
      Mat_v(*L,1,1) = -(0.5)*Load_History;
      Mat_v(*L,2,2) = -(0.5)*Load_History;
      Mat_v(*L,3,3) = +(1.0)*Load_History;
      break;
    case SIMPLE_SHEAR: 
      // shear
      Mat_v(*L, 1,2) = Load_History;
      break;
    case PLAIN_STRAIN_COMPRESSION:     
      // plain_strain_compression
      Mat_v(*L,1,1) = (1.0)*Load_History;
      Mat_v(*L,2,2) = (0.0)*Load_History;
      Mat_v(*L,3,3) = -(1.0)*Load_History;
      break;      
    case UNIAXIAL_COMPRESSION:
      Mat_v(*L,1,1) = (0.5)*Load_History;
      Mat_v(*L,2,2) = (0.5)*Load_History;
      Mat_v(*L,3,3) = -(1.0)*Load_History;
      break;
    case CYCLIC_LOADING:
      if(t<=0.5)
      {
        Mat_v(*L,1,1) = (0.5)*Load_History;
        Mat_v(*L,2,2) = (0.5)*Load_History;
        Mat_v(*L,3,3) = -(1.0)*Load_History;
      }      
      if((0.5 < t) && (t <= 0.9))
      {
        double new_rate = 1.0;
        Mat_v(*L,1,1)=-(0.5)*new_rate;
        Mat_v(*L,2,2)=-(0.5)*new_rate;
        Mat_v(*L,3,3)=+(1.0)*new_rate;
      }
      if(0.9<t && (t <= 1.6))
      {
        double new_rate = -1.0;
        Mat_v(*L,1,1)=-(0.5)*new_rate;
        Mat_v(*L,2,2)=-(0.5)*new_rate;
        Mat_v(*L,3,3)=+(1.0)*new_rate;
      }
      break;      
    case STRESS_RELAXATION:
      if(t <= 0.5)
      {
        Mat_v(*L,1,1) = (0.5)*Load_History;
        Mat_v(*L,2,2) = (0.5)*Load_History;
        Mat_v(*L,3,3) = -(1.0)*Load_History;        
      }    
      break;      
    default:
      break;
  }      
}

int plasticity_model_test_staggered(const HOMMAT *hmat, Matrix(double) *L_in, int Load_Type)
{ 
  int err = 0;  
  double tol_g = 1.0e-6;
  double computer_zero = 1.0e-15;
   
  Constitutive_model m;
  Model_parameters p;
  
  constitutive_model_construct(&m);
  model_parameters_construct(&p);  
  model_parameters_initialize(&p, NULL, NULL, hmat, CRYSTAL_PLASTICITY);
  constitutive_model_initialize(&m, &p);
  
  p.N_SYS = plasticity_model_slip_system(p.Psys);
  plasticity_model_read_parameters(&m);  
  
  int N_SYS = p.N_SYS;
  
  SolverInformation Solver; // Simulation_Settings
  MaterialStructure Struc;  // MaterialStructure_Settings
  MaterialProperties Props; // MaterialProperties_Settings
  MaterialParameters Param; // MaterialParameters_Settings
  err += test_set_CM_interface_values(&Props,&Param,&Struc,&Solver,&m,hmat,N_SYS);  

  double T_Initial = 0.0;
  double T_Final = 1.0;
  double dt = 0.001;
  double Load_History = 1.0;
  
  char fn_out[1024];
  
  if(L_in==NULL)
  {  
    switch(Load_Type)
    {
      case UNIAXIAL_TENSION:
        sprintf(fn_out, "uniaxial_tension.txt"); 
        break;
      case UNIAXIAL_COMPRESSION:
        sprintf(fn_out, "uniaxial_compression.txt");
        break;
      case SIMPLE_SHEAR:
        sprintf(fn_out, "simple_shear.txt");
        break;
      case PLAIN_STRAIN_COMPRESSION:
        sprintf(fn_out, "plain_strain_compression.txt");
        break;
      case CYCLIC_LOADING:
        sprintf(fn_out, "cyclic_loading.txt");
        break;
      case STRESS_RELAXATION:
        sprintf(fn_out, "stress_relaxation.txt");
        break;
      default:
        Load_Type = UNIAXIAL_TENSION;
        sprintf(fn_out, "uniaxial_tension.txt");  
        break;
    }
  }
  else
  {  
    sprintf(fn_out, "user_define_L.txt");
    Load_Type = -1;
  }
  
  FILE *fp = fopen(fn_out, "w");
  
  int Num_Steps=(int)(ceil(T_Final/dt));  
  
  int j_max = DIM;
  Matrix(double) *P_sys = p.Psys;  
  Matrix(double) *Fs;  
  
  enum {Fnp1,Fn,eFn,eFnp1,eFnp1_I,pFn,pFnp1,pFnp1_I,PK2,sigma,L,eFnp1PK2,
        PK2dev,sigma_dev,C,E,Fend};
  
  Fs = (Matrix(double) *) malloc(Fend*sizeof(Matrix(double)));  
  
  for(int a=0; a<Fend; a++)
    Matrix_construct_redim(double, Fs[a],j_max,j_max); 

  double g0=Param.g0;
  double g_n = g0;    
    
  Matrix_eye(Fs[Fn],   j_max);
  Matrix_eye(Fs[Fnp1], j_max);
  Matrix_eye(Fs[pFn],  j_max);
  Matrix_eye(Fs[eFn],  j_max);

  double L_np1 = 0.0;
  double L_n = 0.0;

  double det_Fe;
  double g_np1 = g_n;
  
  double err_g_0 = computer_zero;  
  
  int print_step_info = 0;
  for (int n=1; n<Num_Steps+1; n++)
  {
     double t=((double) n)*dt;
     if(print_step_info)
     {  
       printf("----------------------------------------------------------\n");          
       printf("time = %e\n", t);
       printf("----------------------------------------------------------\n");     
     }

    if(L_in != NULL)
    {
      if(n==1)  
        Matrix_AeqB(Fs[L], 1.0, *L_in);
    }
    else     
       test_load_type((Fs+L), Load_Type, Load_History, t);
     
    //compute Fnp1    
    Matrix_init(Fs[Fnp1], 0.0);
    F_Implicit_sp(dt, Fs[Fn].m_pdata, Fs[L].m_pdata, Fs[Fnp1].m_pdata);

    plasticity_model_staggered_NR((Fs+pFnp1), &g_np1, &L_np1,
                     (Fs+Fnp1), (Fs+eFn), (Fs+pFn),                    
                     &Props, &Param, &Struc, &Solver,
                     &m, dt, g_n);   

    // compute stress (PK2)
    Matrix_AxB(Fs[C],1.0,0.0,Fs[Fnp1],1,Fs[Fnp1],0);
    Matrix_eye(Fs[E], DIM);
    Matrix_AplusB(Fs[E],0.5,Fs[C],-0.5,Fs[E]);
    
//    printf("%e, %e\n", t, Mat_v(Fs[E], 1,1));                                       
    Matrix_inv(Fs[pFnp1], Fs[pFnp1_I]);  
    Matrix_AxB(Fs[eFnp1],1.0,0.0,Fs[Fnp1],0,Fs[pFnp1_I],0);
    elastic_stress(&Props,Fs[eFnp1].m_pdata,Fs[PK2].m_pdata);
    
    // compute Causy stress sigma = 1/det(Fe)*Fe*PK2*Fe'
    Matrix_det(Fs[eFnp1], det_Fe);
    Matrix_AxB(Fs[eFnp1PK2],1.0,0.0,Fs[eFnp1],0,Fs[PK2],0);
    Matrix_AxB(Fs[sigma],1.0/det_Fe,0.0,Fs[eFnp1PK2],0,Fs[eFnp1],1);        

    // update for next time step values
    L_n = L_np1;
    g_n = g_np1;    
    Matrix_AeqB(Fs[Fn], 1.0,Fs[Fnp1]);      
    Matrix_AeqB(Fs[pFn],1.0,Fs[pFnp1]);
    Matrix_AeqB(Fs[eFn],1.0,Fs[eFnp1]);
    
    // print result at time t
    double trPK2, tr_sigma;
    
    Matrix_trace(Fs[PK2],trPK2);
    Matrix_trace(Fs[sigma],tr_sigma);
    Matrix_eye(Fs[PK2dev], DIM);
    Matrix_eye(Fs[sigma_dev], DIM);
        
    Matrix_AplusB(Fs[PK2dev],    1.0, Fs[PK2],      -trPK2/3.0, Fs[PK2dev]);
    Matrix_AplusB(Fs[sigma_dev], 1.0, Fs[sigma], -tr_sigma/3.0, Fs[sigma_dev]);    
    
    double norm_sigma, norm_PK2;
    Matrix_ddot(Fs[PK2dev],Fs[PK2dev],norm_PK2);    
    Matrix_ddot(Fs[sigma_dev],Fs[sigma_dev],norm_sigma);
    
    double sigma_eff=sqrt(3.0/2.0*norm_sigma);
    double PK2_eff = sqrt(3.0/2.0*norm_PK2);    

    fprintf(fp, "%e %e %e %e %e %e\n",t,sigma_eff,PK2_eff, g_n, Mat_v(Fs[E],1,1), Mat_v(Fs[PK2],1,1));    
  }

  constitutive_model_destroy(&m);
  model_parameters_destroy(&p); 
  
  for(int a=0; a<Fend; a++)
    Matrix_cleanup(Fs[a]);  

  fclose(fp);           
  return err;
}

int plasticity_model_test_staggered_F_of_t(const HOMMAT *hmat)
{ 
  // test for defined F
  // F = [1 - t,       0,       0
  //          0, 1 + t/2,       0
  //          0,       0, 1 + t/2];
  int err = 0;  
  double tol_g = 1.0e-6;
  double computer_zero = 1.0e-15;
   
  Constitutive_model m;
  Model_parameters p;
  
  constitutive_model_construct(&m);
  model_parameters_construct(&p);  
  model_parameters_initialize(&p, NULL, NULL, hmat, CRYSTAL_PLASTICITY);
  constitutive_model_initialize(&m, &p);
  
  p.N_SYS = plasticity_model_slip_system(p.Psys);
  plasticity_model_read_parameters(&m);  
  
  int N_SYS = p.N_SYS;
  
  SolverInformation Solver; // Simulation_Settings
  MaterialStructure Struc;  // MaterialStructure_Settings
  MaterialProperties Props; // MaterialProperties_Settings
  MaterialParameters Param; // MaterialParameters_Settings
  err += test_set_CM_interface_values(&Props,&Param,&Struc,&Solver,&m,hmat,N_SYS);  

  double T_Initial = 0.0;
  double T_Final = 0.5;
  double dt = 0.001;
  double Load_History = 1.0;
  
  char fn_out[1024];
  sprintf(fn_out, "uniaxial_compression_F_of_t.txt");
  
  FILE *fp = fopen(fn_out, "w");
  
  int Num_Steps=(int)(ceil(T_Final/dt));  
  
  int j_max = 3;
  Matrix(double) *P_sys = p.Psys;  
  Matrix(double) *Fs;  
  
  enum {Fnp1,Fn,eFn,eFnp1,eFnp1_I,pFn,pFnp1,pFnp1_I,PK2,sigma,eFnp1PK2,
        PK2dev,sigma_dev,C,E,Fend};
  
  Fs = (Matrix(double) *) malloc(Fend*sizeof(Matrix(double)));  
  
  for(int a=0; a<Fend; a++)
    Matrix_construct_redim(double, Fs[a],j_max,j_max); 

  double g0=Param.g0;
  double g_n = g0;    
    
  Matrix_eye(Fs[Fn],   j_max);
  Matrix_eye(Fs[Fnp1], j_max);
  Matrix_eye(Fs[pFn],  j_max);
  Matrix_eye(Fs[eFn],  j_max);

  double L_np1 = 0.0;
  double L_n = 0.0;

  double det_Fe;
  double g_np1 = g_n;
  
  double err_g_0 = computer_zero;  
  
  int print_step_info = 0;
  for (int n=1; n<Num_Steps+1; n++)
  {
     double t=((double) n)*dt;
     if(print_step_info)
     {  
       printf("----------------------------------------------------------\n");          
       printf("time = %e\n", t);
       printf("----------------------------------------------------------\n");     
     }

     
    //define Fnp1    
    Matrix_init(Fs[Fnp1], 0.0);
    Mat_v(Fs[Fnp1],1,1) = 1.0 - t;
    Mat_v(Fs[Fnp1],2,2) = Mat_v(Fs[Fnp1],3,3) = 1.0 + t*0.5;

    plasticity_model_staggered_NR((Fs+pFnp1), &g_np1, &L_np1,
                     (Fs+Fnp1), (Fs+eFn), (Fs+pFn),                    
                     &Props, &Param, &Struc, &Solver,
                     &m, dt, g_n);   

    // compute stress (PK2)
    Matrix_AxB(Fs[C],1.0,0.0,Fs[Fnp1],1,Fs[Fnp1],0);
    Matrix_eye(Fs[E], 3);
    Matrix_AplusB(Fs[E],0.5,Fs[C],-0.5,Fs[E]);
    
//    printf("%e, %e\n", t, Mat_v(Fs[E], 1,1));                                       
    Matrix_inv(Fs[pFnp1], Fs[pFnp1_I]);  
    Matrix_AxB(Fs[eFnp1],1.0,0.0,Fs[Fnp1],0,Fs[pFnp1_I],0);
    elastic_stress(&Props,Fs[eFnp1].m_pdata,Fs[PK2].m_pdata);
    
    // compute Causy stress sigma = 1/det(Fe)*Fe*PK2*Fe'
    Matrix_det(Fs[eFnp1], det_Fe);
    Matrix_AxB(Fs[eFnp1PK2],1.0,0.0,Fs[eFnp1],0,Fs[PK2],0);
    Matrix_AxB(Fs[sigma],1.0/det_Fe,0.0,Fs[eFnp1PK2],0,Fs[eFnp1],1);        

    // update for next time step values
    L_n = L_np1;
    g_n = g_np1;    
    Matrix_AeqB(Fs[Fn], 1.0,Fs[Fnp1]);      
    Matrix_AeqB(Fs[pFn],1.0,Fs[pFnp1]);
    Matrix_AeqB(Fs[eFn],1.0,Fs[eFnp1]);
    
    // print result at time t
    double trPK2, tr_sigma;
    
    Matrix_trace(Fs[PK2],trPK2);
    Matrix_trace(Fs[sigma],tr_sigma);
    Matrix_eye(Fs[PK2dev], DIM);
    Matrix_eye(Fs[sigma_dev], DIM);
        
    Matrix_AplusB(Fs[PK2dev],    1.0, Fs[PK2],      -trPK2/3.0, Fs[PK2dev]);
    Matrix_AplusB(Fs[sigma_dev], 1.0, Fs[sigma], -tr_sigma/3.0, Fs[sigma_dev]);    
    
    double norm_sigma, norm_PK2;
    Matrix_ddot(Fs[PK2dev],Fs[PK2dev],norm_PK2);    
    Matrix_ddot(Fs[sigma_dev],Fs[sigma_dev],norm_sigma);
    
    double sigma_eff=sqrt(3.0/2.0*norm_sigma);
    double PK2_eff = sqrt(3.0/2.0*norm_PK2);    

    fprintf(fp, "%e %e %e %e %e %e\n",t,sigma_eff,PK2_eff, g_n, Mat_v(Fs[E],1,1), Mat_v(Fs[PK2],1,1));    
  }

  constitutive_model_destroy(&m);
  model_parameters_destroy(&p); 
  
  for(int a=0; a<Fend; a++)
    Matrix_cleanup(Fs[a]);  

  fclose(fp);           
  return err;
}


int plasticity_model_test_no_staggered(const HOMMAT *hmat, Matrix(double) *L_in, int Load_Type)
{ 
  int err = 0;  
  double computer_zero = 1.0e-15;
   
  Constitutive_model m;
  Model_parameters p;
  
  constitutive_model_construct(&m);
  model_parameters_construct(&p);  
  model_parameters_initialize(&p, NULL, NULL, hmat, CRYSTAL_PLASTICITY);
  constitutive_model_initialize(&m, &p);
  
  p.N_SYS = plasticity_model_slip_system(p.Psys);
  plasticity_model_read_parameters(&m);  
  
  int N_SYS = p.N_SYS;
  
  SolverInformation Solver; // Simulation_Settings
  MaterialStructure Struc;  // MaterialStructure_Settings
  MaterialProperties Props; // MaterialProperties_Settings
  MaterialParameters Param; // MaterialParameters_Settings
  err += test_set_CM_interface_values(&Props,&Param,&Struc,&Solver,&m,hmat,N_SYS);  

  double T_Initial = 0.0;
  double T_Final = 1.0;
  double dt = 0.001;
  double Load_History = 1.0;

  char fn_out[1024];
  
  if(L_in==NULL)
  {  
    switch(Load_Type)
    {
      case UNIAXIAL_TENSION:
        sprintf(fn_out, "uniaxial_tension.txt"); 
        break;
      case UNIAXIAL_COMPRESSION:
        sprintf(fn_out, "uniaxial_compression.txt");
        break;
      case SIMPLE_SHEAR:
        sprintf(fn_out, "simple_shear.txt");
        break;
      case PLAIN_STRAIN_COMPRESSION:
        sprintf(fn_out, "plain_strain_compression.txt");
        break;
      case CYCLIC_LOADING:
        sprintf(fn_out, "cyclic_loading.txt");
        break;
      case STRESS_RELAXATION:
        sprintf(fn_out, "stress_relaxation.txt");
        break;
      default:
        Load_Type = UNIAXIAL_TENSION;
        sprintf(fn_out, "uniaxial_tension.txt");  
        break;
    }
  }
  else
  {  
    sprintf(fn_out, "user_define_L.txt");
    Load_Type = -1;
  }        
   
  FILE *fp = fopen(fn_out, "w");
  
  int Num_Steps=(int)(ceil(T_Final/dt));  
  
  int j_max = DIM;
  Matrix(double) *P_sys = p.Psys;
  Matrix(double) *Fs;
  
  enum {Fnp1,Fn,eFnp1,eFnp1_I,pFn,pFnp1,pFnp1_I,pFnp1_k,PK2,sigma,L,eFnp1PK2,
        PK2dev,sigma_dev,Fend};

  Fs = (Matrix(double) *) malloc(Fend*sizeof(Matrix(double)));  

  for(int a=0; a<Fend; a++)
    Matrix_construct_redim(double, Fs[a],j_max,j_max); 
  
  Matrix(double) Tau_Array, gamma_RateArray;  
  Matrix_construct_redim(double,Tau_Array      , N_SYS,1);
  Matrix_construct_redim(double,gamma_RateArray, N_SYS,1);

  double g0=Param.g0;
  double g_n = g0;    
    
  Matrix_eye(Fs[Fn],      j_max);
  Matrix_eye(Fs[Fnp1],    j_max);
  Matrix_eye(Fs[eFnp1_I], j_max);
  Matrix_eye(Fs[pFn],     j_max);

  double L_np1 = 0.0;
  double L_n = 0.0;

  double det_Fe;
 
  double g_np1;
  g_np1 = g_n;
  
  int max_itr = 10;
  double err_g_0 = computer_zero;  
  
  int print_step_info = 0;
  for (int n=0; n<Num_Steps+1; n++)
  {
     double t=((double) n)*dt;
     if(print_step_info)
     {  
       printf("----------------------------------------------------------\n");          
       printf("time = %e\n", t);
       printf("----------------------------------------------------------\n");     
     }

    if(L_in != NULL)
    {
      if(n==1)  
        Matrix_AeqB(Fs[L], 1.0, *L_in);
    }
    else     
       test_load_type((Fs+L), Load_Type, Load_History, t);
     
    //compute Fnp1
    Matrix_init(Fs[Fnp1], 0.0);
    F_Implicit_sp(dt, Fs[Fn].m_pdata, Fs[L].m_pdata, Fs[Fnp1].m_pdata);
         
    Matrix_AxB(Fs[pFnp1_k],1.0,0.0,Fs[eFnp1_I],0,Fs[Fnp1],0);
    
     Matrix_init(Fs[pFnp1], 0.0);
    L_np1 = 0.0;

    Staggered_NewtonRapson_Testing_sp(&Props,&Param,&Struc,&Solver,
                             t,dt,P_sys->m_pdata,g_np1,Fs[pFn].m_pdata,Fs[Fnp1].m_pdata,
                             L_n,&L_np1,Fs[pFnp1_k].m_pdata,Fs[pFnp1].m_pdata);

    // compute stress (PK2)
     Matrix_inv(Fs[pFnp1], Fs[pFnp1_I]);
    Matrix_AxB(Fs[eFnp1],1.0,0.0,Fs[Fnp1],0,Fs[pFnp1_I],0);
    elastic_stress(&Props,Fs[eFnp1].m_pdata,Fs[PK2].m_pdata);
    
    // compute Causy stress sigma = 1/det(Fe)*Fe*PK2*Fe' 
    Matrix_det(Fs[eFnp1], det_Fe);
    Matrix_AxB(Fs[eFnp1PK2],1.0,0.0,Fs[eFnp1],0,Fs[PK2],0);
    Matrix_AxB(Fs[sigma],1.0/det_Fe,0.0,Fs[eFnp1PK2],0,Fs[eFnp1],1);
    
    // update hardness explicitly 
    Matrix_init(Tau_Array, 0.0);
    Matrix_init(gamma_RateArray, 0.0);
    double gmdot_tmp = 0.0;
    for (int q = 0; q<N_SYS; q++)
    {
      double tau_q = Tau_Rhs_sp(q, P_sys->m_pdata, Fs[eFnp1].m_pdata, Fs[PK2].m_pdata);
      Vec_v(Tau_Array,q+1) = tau_q;
      Vec_v(gamma_RateArray,q+1) = gamma_Rate_PL(&Param,g_np1,tau_q);
      gmdot_tmp += fabs(Vec_v(gamma_RateArray,q+1));
    }

    double g_Rhs = g_Rate_VK(&Param,&Struc,g_n, gamma_RateArray.m_pdata);
    g_np1 = g_n + dt*(g_Rhs);
    g_Rhs = 0.0;

    // update for next time step
    L_n = L_np1;
    g_n = g_np1;      
    Matrix_AeqB(Fs[Fn], 1.0,Fs[Fnp1]);      
    Matrix_AeqB(Fs[pFn],1.0,Fs[pFnp1]);

    // print result at time t    
    double trPK2, tr_sigma;
    
    Matrix_trace(Fs[PK2],trPK2);
    Matrix_trace(Fs[sigma],tr_sigma);
    Matrix_eye(Fs[PK2dev], 3);
    Matrix_eye(Fs[sigma_dev], 3);
        
    Matrix_AplusB(Fs[PK2dev],    1.0, Fs[PK2],      -trPK2/3.0, Fs[PK2dev]);
    Matrix_AplusB(Fs[sigma_dev], 1.0, Fs[sigma], -tr_sigma/3.0, Fs[sigma_dev]);    
    
    double norm_sigma, norm_PK2;
    Matrix_ddot(Fs[PK2dev],Fs[PK2dev],norm_PK2);    
    Matrix_ddot(Fs[sigma_dev],Fs[sigma_dev],norm_sigma);
    
    double sigma_eff=sqrt(3.0/2.0*norm_sigma);
    double PK2_eff = sqrt(3.0/2.0*norm_PK2);    

    fprintf(fp, "%e \t %e %e %e\n",t,sigma_eff,PK2_eff, g_n);    
  }
  
  constitutive_model_destroy(&m);
  model_parameters_destroy(&p);     
  
  for(int a=0; a<Fend; a++)
    Matrix_cleanup(Fs[a]);

  Matrix_cleanup(Tau_Array);
  Matrix_cleanup(gamma_RateArray);
  fclose(fp);       
  return err;
}

int plasticity_model_test(const HOMMAT *hmat, Matrix(double) *L_in, int Load_Type)
{ 
  int err = 0;
   
  if(1)
  { 
    if(Load_Type<0) 
      err += plasticity_model_test_staggered_F_of_t(hmat);
    else
      err += plasticity_model_test_staggered(hmat, L_in, Load_Type);
  }
  else
    err += plasticity_model_test_no_staggered(hmat, L_in, Load_Type);    
   
  return err;
}

