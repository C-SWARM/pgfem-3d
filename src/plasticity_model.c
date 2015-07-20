/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */

#include "plasticity_model.h"
#include "constitutive_model.h"
#include "new_potentials.h"
#include "data_structure_c.h"
#include <string.h>

Define_Matrix(double);

/**
 * Private structure for use exclusively with this model and
 * associated functions.
 */
typedef struct plasticity_ctx {
  const double *C; /*< pointer to the left Cauchy-Green deformation tensor */
  const double *J; /*< pointer to the Jacobian -OR- Volume term */
} plasticity_ctx;

static double compute_bulk_mod(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static int plasticity_int_alg(Constitutive_model *m,
                                   const void *ctx)
{
  int err = 0;
  /* hyperelastic, no integration algorithm */
  return err;
}

static int plasticity_dev_stress(const Constitutive_model *m,
                                      const void *ctx,
                                      Matrix_double *stress)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  devStressFuncPtr Stress = getDevStressFunc(-1,m->param->p_hmat);
  Stress(CTX->C,m->param->p_hmat,stress->m_pdata);
  return err;
}

static int plasticity_dudj(const Constitutive_model *m,
                                const void *ctx,
                                double *dudj)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  dUdJFuncPtr Pressure = getDUdJFunc(-1,m->param->p_hmat);
  Pressure(*(CTX->J),m->param->p_hmat,dudj);
  return err;
}

static int plasticity_dev_tangent(const Constitutive_model *m,
                                       const void *ctx,
                                       Matrix_double *tangent)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,m->param->p_hmat);
  Tangent(CTX->C,m->param->p_hmat,tangent->m_pdata);
  return err;
}

static int plasticity_d2udj2(const Constitutive_model *m,
                                  const void *ctx,
                                  double *d2udj2)
{
  int err = 0;
  const plasticity_ctx *CTX = ctx;
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,m->param->p_hmat);
  D_Pressure(*(CTX->J),m->param->p_hmat,d2udj2);
  return err;
}

static int plasticity_update(Constitutive_model *m)
{
  int err = 0;
  return err;
}

static int plasticity_reset(Constitutive_model *m)
{
  int err = 0;
  return err;
}

static int plasticity_info(Model_var_info **info)
{
  *info = malloc(sizeof(**info));
  int Fno = 3;
  
  (*info)->n_Fs = Fno;
  (*info)->F_names = (char **)malloc(sizeof(char*)*Fno);
  for(int a=0; a<Fno; a++)
    (*info)->F_names[a] = (char *)malloc(sizeof(char)*1024);

  sprintf((*info)->F_names[TENSOR_F],  "F");
  sprintf((*info)->F_names[TENSOR_Fp], "Fp");
  sprintf((*info)->F_names[TENSOR_Fpnp1],  "Fpnp1");

  int varno = 3;
  (*info)->n_vars = varno;
  (*info)->var_names = (char **)malloc(sizeof(char*)*varno);  
  for(int a=0; a<varno; a++)
    (*info)->var_names[a] = (char *)malloc(sizeof(char)*1024);
  
  sprintf((*info)->F_names[VAR_Ns],  "Ns");
  sprintf((*info)->F_names[VAR_L_n], "L_n");
  sprintf((*info)->F_names[VAR_g_n], "g_n");  
  

  return 0;
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

  return err;
}

int plasticity_model_ctx_build(void **ctx,
                                    const double *C,
                                    const double *J_or_Theta)
{
  int err = 0;
  plasticity_ctx *t_ctx = malloc(sizeof(plasticity_ctx));

  /* assign internal pointers. NOTE: We are copying the pointer NOT
     the value. No additional memory is allocated. */
  t_ctx->C = C;
  t_ctx->J = J_or_Theta;

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

  /* we do not control memory for internal pointers */

  /* free object memory */
  free(t_ctx);
  return err;
}

int compute_P_alpha(Matrix(double) Pa)
{
  Matrix_eye(Pa, 3);
  return 0;
}

int compute_C_D_alpha(Constitutive_model *m, Matrix(double) *aC, Matrix(double) *Da, 
                    Matrix(double) *Fe, Matrix(double) *Pa, Matrix(double) *S,
                    Matrix(double) *L, Matrix(double) *C)
{
  Matrix_init(*aC, 0.0);  
  Matrix_init(*Da, 0.0);
  
  Matrix(double) LC, AA, CAA, FeAA, FeAAMT;
  Matrix_construct_redim(double,LC,    3,3);
  Matrix_construct_redim(double,AA,    3,3);
  Matrix_construct_redim(double,CAA,   3,3);
  Matrix_construct_redim(double,FeAA,  3,3);
  Matrix_construct_redim(double,FeAAMT,3,3);    

  Matrix_AxB(AA,1.0,0.0,*Pa,0,*S,0);   // AA = Pa*S
  Matrix_AxB(AA,1.0,1.0,*S,0,*Pa,1);   // AA = AA + S*Pa' 
       
  Matrix_Tns4_dd_Tns2(LC, *L, *C);     // LC = L:C
  Matrix_AxB(AA,1.0,1.0,LC,0,*Pa,0);   // AA = AA + L:C*Pa

  Matrix_AxB(CAA,1.0,0.0,*C,0,AA,0);   
  Matrix_AxB(*aC,1.0,0.0,m->vars.Fs[1],1,CAA,0);     // Fp: m->vars.Fs[1]

  Matrix_AxB(FeAA,1.0,0.0,*Fe,0,AA,0);
  Matrix_AxB(FeAAMT,1.0,0.0,FeAA,0,m->vars.Fs[2],1); // M: m->vars.Fs[2]
  Matrix_AxB(*Da,1.0,0.0,FeAAMT,0,*Fe,1);
  
  Matrix_cleanup(LC);
  Matrix_cleanup(AA);
  Matrix_cleanup(CAA);
  Matrix_cleanup(FeAA);
  Matrix_cleanup(FeAAMT);       
         
  return 0;
}

int compute_dMdu(Constitutive_model *m, Matrix(double) *dMdu, Matrix(double) *Grad_du, Matrix_double *Fe, 
                 Matrix(double) *S, Matrix(double) *L, double dt)
{
  // compute dMdu:U = -grad(du):B
  // Grad_du = Grad(du)

  // --------------> define variables
  Matrix(double) C;
  Matrix_construct_redim(double, C, 3,3);
  Matrix_AxB(C,1.0,0.0,*Fe,1,*Fe,0);
  
  Matrix(double) U,UI,II,B,aCxPa,CxP,aDxPa,DxP;
  
  Matrix_construct_redim(double, U,     81,1); // 3x3x3x3 tensor
  Matrix_construct_redim(double, UI,    81,1);
  Matrix_construct_redim(double, II,    81,1);
  Matrix_construct_init( double, B,     81,1,0.0);
  Matrix_construct_redim(double, aCxPa, 81,1);
  Matrix_construct_redim(double, CxP,   81,1);
  Matrix_construct_redim(double, aDxPa, 81,1);
  Matrix_construct_redim(double, DxP,   81,1);  

  Matrix(double) aC,Pa,aD,sum_aC,sum_Pa,sum_aD;
  Matrix_construct_redim(double, aC, 3,3);  
  Matrix_construct_redim(double, Pa, 3,3);
  Matrix_construct_redim(double, aD, 3,3);
        
  Matrix_construct_init( double, sum_aC, 3,3,0.0);
  Matrix_construct_init( double, sum_Pa, 3,3,0.0);
  Matrix_construct_init( double, sum_aD, 3,3,0.0);  
  
  Matrix_Tns4_eye(II);
  Matrix_AeqB(U, 1.0/dt, II);
  double R1, R2_a, drdg, drdtau;
  R1 = R2_a = drdg = drdtau = 0.0;
  // <-------------- define variables  
  
  int Ns = (int) Vec_v((m->vars).state_vars[0], VAR_Ns+1);
  for(int a = 1; a<=Ns; a++)
  {
    compute_P_alpha(Pa);
    compute_C_D_alpha(m,&aC, &aD,Fe,&Pa,S,L,&C);
    
    Matrix_AOxB(aCxPa, aC, Pa);
    Matrix_AOxB(aDxPa, aD, Pa);    
    
    Matrix_AplusB(sum_aC, 1.0, sum_aC, R2_a, aC);
    Matrix_AplusB(sum_Pa, 1.0, sum_Pa, drdtau, Pa);
    Matrix_AplusB(sum_aD, 1.0, sum_aD, R2_a, aC);
        
    Matrix_AOxB(aCxPa, aC, Pa);
    Matrix_AOxB(aDxPa, aD, Pa);    
    Matrix_AplusB(U, drdtau, aCxPa, 1.0, U);    
    Matrix_AplusB(B, drdtau, aDxPa, 1.0, B);            
  } 
  
  Matrix_AOxB(CxP, sum_aC, sum_Pa);
  Matrix_AOxB(DxP, sum_aD, sum_Pa);
  Matrix_AplusB(U, R1, CxP, 1.0, U);
  Matrix_AplusB(B, R1, DxP, 1.0, B);
  
  Matrix(double) V;
  Matrix_construct_redim(double, V, 3,3);
  Matrix_Tns4_dd_Tns2(V,B,*Grad_du);
  Matrix_AeqB(V, -1.0, V);
  Matrix_Mat2Vec(V);
  
  Matrix_Tns4_mat_9x9(U);
///////////////////////////////
  Matrix_eye(UI, 9);
///////////////////////////////  

  Matrix_inv(U, UI);

  Matrix_AxB(*dMdu,1.0,0.0,UI,0,V,0);   
  Matrix_Vec2Mat(*dMdu,3,3);

///////////////////////////////  
  Matrix_init(*dMdu, 0.0);
///////////////////////////////
  
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

