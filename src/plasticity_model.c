/**
 * Authors:
 *  Matt Mosby, University of Notre Dame, Notre Dame, IN, <mmosby1@nd.edu>
 */

#include "plasticity_model.h"
#include "constitutive_model.h"
#include "new_potentials.h"
#include "data_structure_c.h"
#include <string.h>
#include <math.h>

#include "CM.h"

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
  int Fno = 6;
  
  (*info)->n_Fs = Fno;
  (*info)->F_names = (char **)malloc(sizeof(char*)*Fno);
  for(int a=0; a<Fno; a++)
    (*info)->F_names[a] = (char *)malloc(sizeof(char)*1024);

  sprintf((*info)->F_names[TENSOR_F],         "F");
  sprintf((*info)->F_names[TENSOR_Fp],        "Fp");
  sprintf((*info)->F_names[TENSOR_M],         "M");
  sprintf((*info)->F_names[TENSOR_P_sys],     "P_sys");  
  sprintf((*info)->F_names[TENSOR_tau],       "tau"); 
  sprintf((*info)->F_names[TENSOR_gamma_dot], "gamma_dot");     
  
  int varno = 11;
  (*info)->n_vars = varno;
  (*info)->var_names = (char **)malloc(sizeof(char*)*varno);  
  for(int a=0; a<varno; a++)
    (*info)->var_names[a] = (char *)malloc(sizeof(char)*1024);
  
  sprintf((*info)->var_names[VAR_Ns],         "Ns");         
  sprintf((*info)->var_names[VAR_L_np1],      "L_n");        
  sprintf((*info)->var_names[VAR_g_n],        "g_n");        
  sprintf((*info)->var_names[VAR_g_np1],      "g_np1");      
  sprintf((*info)->var_names[VAR_gamma_dot_0],"gamma_dot_0");
  sprintf((*info)->var_names[VAR_gamma_dot_s],"gamma_dot_s"); 
  sprintf((*info)->var_names[VAR_m],          "m");          
  sprintf((*info)->var_names[VAR_g0],         "g0");         
  sprintf((*info)->var_names[VAR_G0],         "G0");         
  sprintf((*info)->var_names[VAR_gs_0],       "gs_0");       
  sprintf((*info)->var_names[VAR_w],           "w");           

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

int compute_P_alpha(Constitutive_model *m, int alpha, Matrix(double) *Pa)
{
  int Ns = (int)(m->vars).state_vars[0].m_pdata[VAR_Ns];
  double *P_sys = ((m->vars).Fs[TENSOR_P_sys]).m_pdata;
  Matrix_redim(*Pa, 3,3);
  for(int a=0; a<3; a++)
  {
    for(int b=0;b<3;b++)
      Mat_v(*Pa,a+1,b+1) = P_sys[Index_3D(alpha,a,b,Ns,3)];
  }
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
  Matrix_AxB(*aC,1.0,0.0,m->vars.Fs[TENSOR_Fp],1,CAA,0);     // Fp: m->vars.Fs[TENSOR_Fp]

  Matrix_AxB(FeAA,1.0,0.0,*Fe,0,AA,0);
  Matrix_AxB(FeAAMT,1.0,0.0,FeAA,0,m->vars.Fs[TENSOR_M],1); // M: m->vars.Fs[TENSOR_M]
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

  double *state_var = (m->vars).state_vars[0].m_pdata;
  
  int Ns        = (int)state_var[VAR_Ns];
  double g_n         = state_var[VAR_g_n];
  double g_np1       = state_var[VAR_g_np1];  
  double gamma_dot_0 = state_var[VAR_gamma_dot_0];
  double gamma_dot_s = state_var[VAR_gamma_dot_s];  
  double mm          = state_var[VAR_m];
  double g0          = state_var[VAR_g0];  
  double G0          = state_var[VAR_G0];    
  double gs_0        = state_var[VAR_gs_0];
  double w           = state_var[VAR_w];
  
  
  double *tau_np1 = (m->vars).Fs[TENSOR_tau].m_pdata;
  double *gamma_dots = (m->vars).Fs[TENSOR_gamma_dot].m_pdata;
    

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

  // <-------------- define variables  
  
  double gamma_dot = 0.0;
  double sum_1gm1gm = 0.0;
  for(int a = 0; a<Ns; a++)
  {
    gamma_dot += fabs(gamma_dots[a]);
    sum_1gm1gm += fabs(gamma_dots[a])/g_np1;
  }

  double gm_gms   = gamma_dot/gamma_dot_s;
  double sign_gm_gms = (gm_gms < 0) ? -1.0 : 1.0;
  double R3 = gs_0*w/gamma_dot_s*sign_gm_gms*pow(fabs(sign_gm_gms), w-1.0);
  double gs_np1 = gs_0*pow(fabs(gm_gms),w); 
  double AA = R3*gamma_dot*(g_n - g0 + dt*G0*gamma_dot) + gs_np1*(gs_np1 - g0 - g_n) + g0*g_n;
  double BB = gs_np1 - g0  - dt*G0*gamma_dot;
  double R4 = dt*G0*AA/BB/BB;

  double R1 = R4/(1.0-R4*sum_1gm1gm);
  
  for(int a = 0; a<Ns; a++)
  {
    double drdtau = gamma_dot_0/mm/g_np1*pow(fabs(tau_np1[a]/g_np1), 1.0/mm - 1.0);
    double drdg   = -drdtau*tau_np1[a]/g_np1;

    double R2_a = fabs(gamma_dots[a])/tau_np1[a];
    
    compute_P_alpha(m,a,&Pa);
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
  Matrix_inv(U, UI);

  Matrix_AxB(*dMdu,1.0,0.0,UI,0,V,0);   
  Matrix_Vec2Mat(*dMdu,3,3);
  
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

int plasticity_model_slip_system(Constitutive_model *m)
{
  int err = 0;
  int N_SYS = 12; // depended on slip system

  int j_max = 3;  
  Matrix_redim((m->vars).Fs[TENSOR_P_sys], N_SYS*j_max*j_max, 1);
  double *P_sys = (m->vars).Fs[TENSOR_P_sys].m_pdata;
  
  Vec_v((m->vars).state_vars[0], VAR_Ns+1)  = N_SYS;

	for (int k = 0; k<N_SYS; k++)
	{
		for (int j = 0; j<j_max; j++)
		{
			for (int i = 0; i<j_max; i++)
				P_sys[Index_3D(k,j,i,j_max,j_max)] = 0.0;
		}
		P_sys_sp(k, P_sys);
	}
	
	return err;	
}