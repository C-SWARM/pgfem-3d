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

  sprintf((*info)->F_names[TENSOR_Fn],        "Fn");
  sprintf((*info)->F_names[TENSOR_pFn],       "pFp");
  sprintf((*info)->F_names[TENSOR_Fnp1],      "Fnp1");
  sprintf((*info)->F_names[TENSOR_pFnp1],     "pFpnp1");    
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
  double *P_sys = ((m->param)->Psys)->m_pdata;
  Matrix_redim(*Pa, 3,3);
  for(int a=0; a<3; a++)
  {
    for(int b=0;b<3;b++)
      Mat_v(*Pa,a+1,b+1) = P_sys[Index_3D(alpha,a,b,3,3)];
  }
  return 0;
}

int compute_C_D_alpha(Constitutive_model *m, Matrix(double) *aC, Matrix(double) *aD, 
                    Matrix(double) *eFn, Matrix(double) *eFnp1,Matrix(double) *M, Matrix(double) *Pa, Matrix(double) *S,
                    Matrix(double) *L, Matrix(double) *C)
{
  Matrix_init(*aC, 0.0);  
  Matrix_init(*aD, 0.0);
  
  Matrix(double) LC, AA, CAA, eFnp1AA, eFnp1AAMT, MI;
  Matrix_construct_redim(double,LC,       3,3);
  Matrix_construct_redim(double,AA,       3,3);
  Matrix_construct_redim(double,CAA,      3,3);
  Matrix_construct_redim(double,eFnp1AA,  3,3);
  Matrix_construct_redim(double,eFnp1AAMT,3,3);    
  Matrix_construct_redim(double,MI,       3,3);      

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

int compute_dMdu(Constitutive_model *m, Matrix(double) *dMdu, 
                 Matrix(double) *Grad_du, Matrix(double) *eFn, Matrix(double) *eFnp1, Matrix(double) *M,
                 Matrix(double) *S, Matrix(double) *L, double dt)                 
{
  // compute dMdu:U = -grad(du):B
  // Grad_du = Grad(du)

  double *state_var = (m->vars).state_vars[0].m_pdata;
  
  int N_SYS             = (m->param)->N_SYS;;
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
  Matrix_AxB(C,1.0,0.0,*eFnp1,1,*eFnp1,0);
  
  Matrix(double) U,UI,II,B,aCxPa,CxP,aDxPa,DxP;
  
  Matrix_construct_redim(double, U,     81,1); // 3x3x3x3 tensor
  Matrix_construct_redim(double, UI,    9,9);
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
  Matrix_construct_init(double, V, 3,3,0.0);
  Matrix_Tns2_dd_Tns4(V,*Grad_du,B);
  Matrix_AeqB(V, -1.0, V);
  
  Matrix_Tns4_mat_9x9(U);
  Matrix_inv(U, UI);
  
  dMdu->m_row = 1;
  dMdu->m_col = 9;
  Matrix_Mat2Vec(V);        
  Matrix_AxB(*dMdu,1.0,0.0,V,1,UI,0);   
  Matrix_Vec2Mat(*dMdu,3,3);*/
  
  Matrix(double) V;
  Matrix_construct_redim(double, V, 3,3);
  Matrix_Tns4_dd_Tns2(V,B,*Grad_du);
  Matrix_AeqB(V, -1.0, V);
  
  Matrix_Tns4_mat_9x9(U);
  Matrix_inv(U, UI);

  Matrix_Mat2Vec(V);
  Matrix_Mat2Vec(*dMdu);
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

int plasticity_model_slip_system(Matrix_double *P)
{
  int err = 0;
  int N_SYS = 12; // depended on slip system

  int j_max = 3;  
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

void elastic_stress(MaterialProperties *Props, double *F, double *S)
{
  Matrix(double) Se, Fe;
  Matrix_construct_redim(double,Se,3,3);
  Matrix_construct(double, Fe);   
  Matrix_init_w_array(Fe,3,3,F);  
  constitutive_model_update_elasticity(Props->cm, &Fe, 0.0, NULL, &Se, 0);

  for(int a=0;a<9;a++)
    S[a] = Se.m_pdata[a];
    
  Matrix_cleanup(Se);
  Matrix_cleanup(Fe);   
}

void elastic_tangent(MaterialProperties *Props, double *F, double *L)
{
  Matrix(double) Se, Fe, Le;
  Matrix_construct_redim(double,Se,3,3); 
  Matrix_construct(double, Fe);  
  Matrix_construct_redim(double,Le,81,1);     
  Matrix_init_w_array(Fe,3,3,F);  
  constitutive_model_update_elasticity(Props->cm, &Fe, 0.0, &Le, &Se, 1);

  for(int a=0;a<81;a++)
    L[a] = Le.m_pdata[a];
    
  Matrix_cleanup(Se);
  Matrix_cleanup(Le);  
  Matrix_cleanup(Fe);   
}

int plasticity_model_integration_ip(Matrix_double *pFnp1, Constitutive_model *m, Matrix_double *Fnp1, Matrix_double *Fe_n, double dt)
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
	
	int j_max = 3;

  Matrix(double) Fe_I,Fp_np1k;
	
  Matrix_construct_redim(double,Fe_I,3,3);
  Matrix_construct_redim(double,Fp_np1k,3,3);

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
  Matrix_construct_redim(double,S_n,3,3);       
  Matrix_construct_redim(double,pFnp1_I,3,3);       
  Matrix_construct_redim(double,eFnp1,3,3);           
	
	Matrix_inv(*pFnp1,pFnp1_I);
	Matrix_AxB(eFnp1,1.0,0.0,*Fnp1,0,pFnp1_I,0);
	
	elastic_stress(&Props, eFnp1.m_pdata, S_n.m_pdata);  

//  double gamma_dot = 0.0;    //////////////////////// compute g_np1
  for (int k = 0; k<N_SYS; k++)
	{
	  double tau_k = Tau_Rhs_sp(k, P_sys, eFnp1.m_pdata, S_n.m_pdata);
    Vec_v(Fs[TENSOR_tau],k+1) = tau_k;
    Vec_v(Fs[TENSOR_gamma_dot],k+1) = gamma_Rate_PL(&Params,g_n,tau_k);
//    gamma_dot += fabs(Vec_v(Fs[TENSOR_gamma_dot],k+1)); //////////////////////// compute g_np1
  }
  double g_Rhs = g_Rate_VK(&Params,&Struc,g_n, Fs[TENSOR_gamma_dot].m_pdata);
  state_var[VAR_g_np1] =  g_n + dt*g_Rhs;
  state_var[VAR_L_np1] =  L_np1;
  
  
  //////////////////////// compute g_np1
//  double gm_gms   = gamma_dot/Params.gamsdot;
//  double gs_np1 = 0.0;   
//  if(fabs(gm_gms)>1.0e-15)
//    gs_np1 = Params.gs0*pow(fabs(gm_gms),Params.w);
  
    
//  double gg = ((gs_np1-Params.g0)*g_n + dt*Params.G0*gs_np1*gamma_dot)/(gs_np1 - Params.g0 + dt*Params.G0*gamma_dot);
//  printf("%e %e %e\n", g_n + dt*g_Rhs, gg, g_n + fabs(dt*g_Rhs - gg));  
//  state_var[VAR_g_np1] = gg;
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
  
  state_var[VAR_gamma_dot_0] = 1.0;
  state_var[VAR_gamma_dot_s] = 50.0e+9;  
  state_var[VAR_m]           = 0.05;
  state_var[VAR_g0]          = 210.0;  
  state_var[VAR_G0]          = 200.0;    
  state_var[VAR_gs_0]        = 330.0;
  state_var[VAR_w]           = 0.005;
  
    
  state_var[VAR_g_n] = state_var[VAR_g0];
  state_var[VAR_g_np1] = state_var[VAR_g0];    
  state_var[VAR_L_n] = 0.0;  
  state_var[VAR_L_np1] = 0.0;
  
  int N_SYS = (m->param)->N_SYS;
  Matrix(double) *Fs = (m->vars).Fs;
  Matrix_redim(Fs[TENSOR_tau],N_SYS, 1);
  Matrix_redim(Fs[TENSOR_gamma_dot],N_SYS, 1);    
  return err;
}

int plasticity_model_test(const HOMMAT *hmat, Matrix(double) *L_in, int Print_results)
{  
  Constitutive_model m;
  Model_parameters p;
  
  constitutive_model_construct(&m);
  model_parameters_construct(&p);  
  model_parameters_initialize(&p, NULL, NULL, hmat, CRYSTAL_PLASTICITY);
  constitutive_model_initialize(&m, &p);
  
  p.N_SYS = plasticity_model_slip_system(p.Psys);
  plasticity_model_read_parameters(&m);  
  
  int err = 0;
  int N_SYS = p.N_SYS;
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
	Props.cm = &m;
	Props.compute_elastic_stress = elastic_stress;
	Props.compute_elastic_tangent = elastic_tangent;	

	/*--------MaterialParameters_Settings--------*/
	MaterialParameters Param;
	Param.Model_Type = PL_VK;
	Param.Parameters_Count = 7;
	Param.gam0dot = 1.0;
	Param.m_matl = 0.05;
	Param.gs0 = 330.0;
	Param.gamsdot = 50000000000.0;
	Param.w = 0.005;
	Param.G0 = 200.0;
	Param.g0 = 210.0;
	
	double T_Initial = 0.0;
	double T_Final = 1.0;
	double dt = 0.0001;
  double Load_History = 1.0;
	int Load_Type = UNIAXIAL_COMPRESSION;	
	
	int Num_Steps=(int)(ceil(T_Final/dt));  
  
  int j_max = 3;
  Matrix(double) *P_sys = p.Psys;
  
  Matrix(double) F_np1,F_n,Fe_n,Fe_I,Fp_n,Fp_np1,Fp_np1k;
  Matrix(double) F_I,S_n,T_n,L,tmp,R_tmp,Fp_I;
  Matrix(double) Tau_Array, gamma_RateArray, tmp_Eigen;
  
  Matrix_construct_redim(double,F_np1  ,j_max,j_max);
  Matrix_construct_redim(double,F_n    ,j_max,j_max);
  Matrix_construct_redim(double,Fe_n   ,j_max,j_max);
  Matrix_construct_redim(double,Fe_I   ,j_max,j_max);
  Matrix_construct_redim(double,Fp_n   ,j_max,j_max);
  Matrix_construct_redim(double,Fp_np1 ,j_max,j_max); 
  Matrix_construct_redim(double,Fp_np1k,j_max,j_max); 
  Matrix_construct_redim(double,F_I    ,j_max,j_max);
  Matrix_construct_redim(double,S_n    ,j_max,j_max);
  Matrix_construct_redim(double,T_n    ,j_max,j_max);
  Matrix_construct_redim(double,L      ,j_max,j_max);
  Matrix_construct_redim(double,tmp    ,j_max,j_max);
  Matrix_construct_redim(double,R_tmp  ,j_max,j_max);
  Matrix_construct_redim(double,Fp_I   ,j_max,j_max);
  
  Matrix_construct_redim(double,Tau_Array      , N_SYS,1);
  Matrix_construct_redim(double,gamma_RateArray, N_SYS,1);
  Matrix_construct_redim(double,tmp_Eigen,j_max,1);

  double g0=Param.g0;
  double g_n = g0;
	double g_np1 = 0.0;
	  
	switch(Load_Type)
	{
    case UNIAXIAL_TENSION:
		  // Tension
  		Matrix_eye(L,3);
  		Mat_v(L,1,1) = -(0.5)*Load_History;
  		Mat_v(L,2,2) = -(0.5)*Load_History;
  		Mat_v(L,3,3) = +(1.0)*Load_History;
      break;
    case SIMPLE_SHEAR: 
      // shear
      Matrix_init(L, 0.0);
      Mat_v(L, 1,2) = Load_History;
      break;
    case PLAIN_STRAIN_COMPRESSION:     
		  // plain_strain_compression
  		Matrix_eye(L,3);
  		Mat_v(L,1,1) = (1.0)*Load_History;
  		Mat_v(L,2,2) = (0.0)*Load_History;
  		Mat_v(L,3,3) = -(1.0)*Load_History;
      break;      
    case UNIAXIAL_COMPRESSION:
      // compressions use same setting, contiue             
    case CYCLIC_LOADING:
      // compressions use same setting, contiue             
    case STRESS_RELAXATION:
  		Matrix_eye(L,3);
  		Mat_v(L,1,1) = (0.5)*Load_History;
  		Mat_v(L,2,2) = (0.5)*Load_History;
  		Mat_v(L,3,3) = -(1.0)*Load_History;
      break;      
    default:
      break;
	}	
	
	if(L_in != NULL)
	  Matrix_AeqB(L, 1.0, *L_in);
	  
	Matrix_eye(F_n,3);
	Matrix_init(F_np1, 0.0);
	F_Implicit_sp(dt, F_n.m_pdata, L.m_pdata, F_np1.m_pdata);
	Matrix_print(F_np1);

  Matrix_AeqB(Fe_n,1.0,F_n);
	Matrix_eye(Fp_n,3);
	Matrix_inv(Fe_n, Fe_I);

	Matrix_AxB(Fp_np1k,1.0,0.0,Fe_I,0,F_np1,0);

	double L_np1k = 0.001;
	double L_np1;
	double L_n = 0.0;

	Matrix_init(S_n, 0.0);
	elastic_stress(&Props,Fe_n.m_pdata, S_n.m_pdata);

	double g_np1k2 = 0.0;
	double g_RESIDUAL = 0.0;
	double g_np1k1 = g_n;
	double det_Fe;

  Matrix_det(Fe_n, det_Fe);
  Matrix(double) FeS;
  Matrix_construct_redim(double,FeS,j_max,j_max);
	Matrix_init(T_n, 0.0);
	Matrix_AxB(FeS,1.0,0.0,Fe_n,0,S_n,0);
	Matrix_AxB(T_n,1.0/det_Fe,0.0,FeS,0,Fe_n,1);

  if(Print_results)
  	printf("time \t Effective Stress \n");
	

  double t = T_Initial;

	double norm_T, norm_S;

  Matrix_ddot(T_n,T_n,norm_T);
  Matrix_ddot(S_n,S_n,norm_S);
    
  double s_eff=sqrt(3.0/2.0)*sqrt(norm_T);
	double PK2_eff = sqrt(3.0/2.0)*sqrt(norm_S);

  if(Print_results)
	  printf("%e \t %e %e\n",t,s_eff,PK2_eff);
 
	for (int kk=1; kk<Num_Steps+1; kk++)
	{
		Matrix_init(Fp_np1, 0.0);
		L_np1 = 0.0;

		Staggered_NewtonRapson_Testing_sp(&Props,&Param,&Struc,&Solver,
				                       t,dt,P_sys->m_pdata,g_n,Fp_n.m_pdata,F_np1.m_pdata,
				                       L_np1k,&L_np1,Fp_np1k.m_pdata,Fp_np1.m_pdata);

		Matrix_init(Fp_n, 0.0);
		Matrix_AeqB(Fp_n,1.0,Fp_np1);

		L_n = 0.0;
		L_n = L_np1;

		Matrix_init(Fp_I, 0.0);
		Matrix_inv(Fp_n, Fp_I);

		Matrix_init(Fe_n, 1.0);
		Matrix_AxB(Fe_n,1.0,0.0,F_np1,0,Fp_I,0);

		Matrix_init(Fe_I, 0.0);
		Matrix_inv(Fe_n, Fe_I);

		if(Load_Type==STRESS_RELAXATION && (t > 0.5) )
		{
			Mat_v(L,1,1) = Mat_v(L,2,2) = Mat_v(L,3,3) = 0.0;
		}
		if(Load_Type==CYCLIC_LOADING )
		{
			//plug_for_cyclic_loading:
			if ((t > 0.5) && (t <= 0.9))
			{
				double new_rate = -1.0;
				Mat_v(L,1,1)=-(0.5)*new_rate;
				Mat_v(L,2,2)=-(0.5)*new_rate;
				Mat_v(L,3,3)=+(1.0)*new_rate;
			}
			if (t>0.9 && (t <= 1.6))
			{
				double new_rate = 1.0;
				Mat_v(L,1,1)=-(0.5)*new_rate;
			  Mat_v(L,2,2)=-(0.5)*new_rate;
				Mat_v(L,3,3)=+(1.0)*new_rate;
			}
		}

		t=((double) (kk))*dt;

		//compute::F_np1
		Matrix_init(F_n, 0.0);
		Matrix_AeqB(F_n,1.0,F_np1);
		Matrix_init(F_np1, 0.0);
     	if( Load_Type==UNIAXIAL_COMPRESSION || 
     	    Load_Type==STRESS_RELAXATION || 
     	    Load_Type==CYCLIC_LOADING)
		{
			//CMP:compression
			Matrix_eye(L,3);
			Mat_v(L,1,1) = +(0.5)*Load_History;
			Mat_v(L,2,2) = +(0.5)*Load_History;
			Mat_v(L,3,3) = -(1.0)*Load_History;
		}
		
	  if(L_in != NULL)
	    Matrix_AeqB(L, 1.0, *L_in);
		
		F_Implicit_sp(dt, F_n.m_pdata, L.m_pdata, F_np1.m_pdata);		

		Matrix_init(Fp_np1k, 0.0);
		Matrix_AxB(Fp_np1k,1.0,0.0,Fe_I,0,F_np1,0);

		L_np1k = 0.0;
		L_np1k = L_n;

		Matrix_init(S_n,0.0);
		elastic_stress(&Props,Fe_n.m_pdata,S_n.m_pdata);

		Matrix_init(T_n,0.0);
    Matrix_det(Fe_n, det_Fe);
	  Matrix_AxB(FeS,1.0,0.0,Fe_n,0,S_n,0);
	  Matrix_AxB(T_n,1.0/det_Fe,0.0,FeS,0,Fe_n,1);

    Matrix_ddot(T_n,T_n,norm_T);
    Matrix_ddot(S_n,S_n,norm_S);
    
		s_eff=sqrt(3.0/2.0)*sqrt(norm_T);
		PK2_eff = sqrt(3.0/2.0)*sqrt(norm_S);

    if(Print_results)
  	  printf("%e \t %e %e\n",t,s_eff,PK2_eff);
  
		Matrix_init(Tau_Array, 0.0);
		Matrix_init(gamma_RateArray, 0.0);
		double gmdot_tmp = 0.0;
		for (int k = 0; k<N_SYS; k++)
		{
		  double tau_k = Tau_Rhs_sp(k, P_sys->m_pdata, Fe_n.m_pdata, S_n.m_pdata);
			Vec_v(Tau_Array,k+1) = tau_k;
			Vec_v(gamma_RateArray,k+1) = gamma_Rate_PL(&Param,g_n,tau_k);
			gmdot_tmp += fabs(Vec_v(gamma_RateArray,k+1));
		}

		double g_Rhs = g_Rate_VK(&Param,&Struc,g_n, gamma_RateArray.m_pdata);
		g_np1 = g_n + dt*(g_Rhs);
		g_Rhs = 0.0;

		//update_hardening
		g_n = g_np1;
		g_np1 = 0.0;
	}
	
  Matrix_cleanup(F_np1  );
  Matrix_cleanup(F_n    );
  Matrix_cleanup(Fe_n   );
  Matrix_cleanup(Fe_I   );
  Matrix_cleanup(Fp_n   );
  Matrix_cleanup(Fp_np1 ); 
  Matrix_cleanup(Fp_np1k); 
  Matrix_cleanup(F_I    );
  Matrix_cleanup(S_n    );
  Matrix_cleanup(T_n    );
  Matrix_cleanup(L      );
  Matrix_cleanup(tmp    );
  Matrix_cleanup(R_tmp  );
  Matrix_cleanup(Fp_I   );	

  Matrix_cleanup(Tau_Array);
  Matrix_cleanup(gamma_RateArray);
  Matrix_cleanup(tmp_Eigen);  

	
  Matrix_cleanup(FeS);

  constitutive_model_destroy(&m);
  model_parameters_destroy(&p);  	  
  return err;
}

