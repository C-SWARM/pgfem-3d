/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "constitutive_model.h"

#include "plasticity_model_none.h"
#include "plasticity_model.h"

#include "material.h"
#include "hommat.h"
#include "matgeom.h"
#include "PGFEM_io.h"
#include "data_structure_c.h"
#include "elem3d.h"

#include "CM.h"

Define_Matrix(double);

int constitutive_model_construct(Constitutive_model *m)
{
  int err = 0;
  m->param = NULL;
  err += state_variables_build(&(m->vars));
  return err;
}

int constitutive_model_initialize(Constitutive_model *m,
                                  const Model_parameters *p)
{
  int err = 0;
  if (p == NULL){
    err++;
  } else {
    m->param = p;
    Model_var_info *info = NULL;
    m->param->get_var_info(&info);
    err += state_variables_initialize(&(m->vars),info->n_Fs,info->n_vars);
    err += model_var_info_destroy(&info);
  }
  return err;
}

int constitutive_model_destroy(Constitutive_model *m)
{
  int err = 0;
  /* drop pointer to Model_parameters object (deallocated elsewhere) */
  m->param = NULL;
  err += state_variables_destroy(&(m->vars));
  return err;
}

int model_var_info_destroy(Model_var_info **info)
{
  int err = 0;
  Model_var_info *t_info = *info;
  /* invalidate pointer */
  *info = NULL;

  /* deallocate internal memory */
  for (size_t i=0, e=t_info->n_Fs; i < e; i++){
    if(t_info->F_names[i])
      free(t_info->F_names[i]);
  }
  if(t_info->F_names)
    free(t_info->F_names);

  for (size_t i=0, e=t_info->n_vars; i < e; i++){
    if(t_info->var_names[i])
      free(t_info->var_names[i]);
  }
  if(t_info->var_names)
    free(t_info->var_names);

  /* destroy memory for structure */
  free(t_info);

  return err;
}

int model_parameters_construct(Model_parameters *p)
{
  int err = 0;
  /* poison all values */
  p->p_mat = NULL;
  p->p_mgeom = NULL;
  p->p_hmat = NULL;
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_dudj = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->get_var_info = NULL;
  p->type = -1;
  return err;
}

int model_parameters_initialize(Model_parameters *p,
                          const MATERIAL *p_mat,
                          const MATGEOM_1 *p_mgeom,
                          const HOMMAT *p_hmat,
                          const size_t type)
{
  int err = 0;
  p->p_mat = p_mat;
  p->p_mgeom = p_mgeom;
  p->p_hmat = p_hmat;
  p->type = type;
  switch(type) {
  case HYPER_ELASTICITY:
    err += plasticity_model_none_initialize(p);
    break;
  case CRYSTAL_PLASTICITY:
    err += plasticity_model_initialize(p);
    break;
  case BPA_PLASTICITY:
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",type);
    err++;
    break;
  }
  return err;
}

int model_parameters_destroy(Model_parameters *p)
{
  int err = 0;
  /* drop pointer to material (material free'd elsewhere) */
  p->p_mat = NULL;
  p->p_mgeom = NULL;
  p->p_hmat = NULL;

  /* drop function pointers */
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_dudj = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->get_var_info = NULL;

  /* reset counters/flags */
  p->type = -1;

  return err;
}

int constitutive_model_update_elasticity(Constitutive_model *m, Matrix_double *Fe, double dt,
                               Matrix_double *L, Matrix_double *S, int compute_stiffness)
{
  int err = 0;
  void *ctx;
  double J;
  Matrix(double) C, CI;    
  Matrix_construct_redim(double,C,3,3); 
  Matrix_construct_redim(double,CI,3,3);     
  
  Matrix_AxB(C, 1.0, 0.0, *Fe, 1, *Fe, 0);
  Matrix_inv(C,CI);        
  Matrix_det(*Fe, J);

  switch(m->param->type)
  {
    case HYPER_ELASTICITY:
      err += plasticity_model_none_ctx_build(&ctx, C.m_pdata, &J);
      break;
    case CRYSTAL_PLASTICITY:
      
      err += plasticity_model_ctx_build(&ctx, C.m_pdata, &J);
      
      break;
    case BPA_PLASTICITY:
    default:
      PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",m->param->type);
      err++;
      break;
  }

  // compute stress
  double dudj = 0.0;
  double d2udj2 = 0.0;
  double nu = (m->param->p_hmat)->nu;
  double G = (m->param->p_hmat)->G;
  double kappa = ((2.0*G *(1.0+nu))/(3.0*(1.0-2.0*nu)));    
  err += m->param->compute_dev_stress(m, ctx, S);
  err += m->param->compute_dudj(m,ctx,&dudj);    
  Matrix_AplusB(*S, kappa*J*dudj,CI,1.0,*S);
  //compute stiffness
  if(compute_stiffness)
  {  
    Matrix(double) CIoxCI, CICI, SoxS;    
    Matrix_construct_redim(double,CIoxCI,81,1);
    Matrix_construct_redim(double,CICI,81,1);             
    Matrix_construct_redim(double,SoxS,81,1);
    
    err += m->param->compute_dev_tangent(m, ctx, L);              
    err += m->param->compute_d2udj2(m,ctx,&d2udj2);
  
    for(int I=1; I<=3; I++)
    {
      for(int JJ=1; JJ<=3; JJ++)
      {
        for(int P=1; P<=3; P++)
        {
          for(int Q=1; Q<=3; Q++)
          {
            Tns4_v(CIoxCI,I,JJ,P,Q) = Mat_v(CI,I,JJ)*Mat_v(CI,P,Q);
            Tns4_v(SoxS,I,JJ,P,Q) = Mat_v(CI,I,JJ)*Mat_v(CI,P,Q);            
            Tns4_v(CICI,I,JJ,P,Q) = Mat_v(CI,I,P)*Mat_v(CI,Q,JJ);
          }
        }
      }
    }
  
    double H = 0.0; /* use stored damage evolution parameter */
    for(int I=1; I<=81; I++)
    {
      Vec_v(*L, I) += kappa*(J*dudj + J*J*d2udj2)*Vec_v(CIoxCI, I)
                   - 2.0*kappa*J*dudj*Vec_v(CICI, I) - H*Vec_v(SoxS, I);
    }    
    Matrix_cleanup(CIoxCI);
    Matrix_cleanup(CICI);
    Matrix_cleanup(SoxS);             
  }

  switch(m->param->type)
  {
    case HYPER_ELASTICITY:
      err += plasticity_model_none_ctx_destroy(&ctx);
      break;
    case CRYSTAL_PLASTICITY:
      err += plasticity_model_ctx_destroy(&ctx);
      break;
    case BPA_PLASTICITY:
    default:
      PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",m->param->type);
      err++;
      break;
  }
  
  Matrix_cleanup(C);
  Matrix_cleanup(CI); 
  return err; 
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

int integration_ip(Matrix_double *pFnp1, Constitutive_model *m, Matrix_double *Fe_n, double dt)
{
  const HOMMAT *hmat = m->param->p_hmat;
  double *state_var = (m->vars).state_vars[0].m_pdata;
  Matrix(double) *Fs = (m->vars).Fs;
  
  int err = 0;
  int N_SYS = (int) state_var[VAR_Ns];
  double *P_sys = Fs[TENSOR_P_sys].m_pdata;
  double *Fp    = Fs[TENSOR_Fp].m_pdata;
  double *F     = Fs[TENSOR_F].m_pdata;
  
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
  state_var[VAR_gamma_dot_0] = 1.0;
  state_var[VAR_gamma_dot_s] = 50.0e+9;  
  state_var[VAR_m]           = 0.05;
  state_var[VAR_g0]          = 0.21;  
  state_var[VAR_G0]          = 0.2;    
  state_var[VAR_gs_0]        = 0.33;
  state_var[VAR_w]           = 0.005;  
  
//  Matrix_print((m->vars).state_vars[0]);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
	
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
  
  double L_np1k = state_var[VAR_L_np1];
  double L_np1 = 0.0;
  Matrix_inv(*Fe_n, Fe_I);
  Matrix_AxB(Fp_np1k,1.0,0.0,Fe_I,0,Fs[TENSOR_F],0);
	
	Staggered_NewtonRapson_Testing_sp(&Props,&Params,&Struc,&Solver,0.0,dt,
									   P_sys,g_n,Fp,F,L_np1k,&L_np1,Fp_np1k.m_pdata,pFnp1->m_pdata);

  Matrix(double) S_n;
  Matrix_construct_redim(double,S_n,3,3);       
	
	elastic_stress(&Props, Fe_n->m_pdata, S_n.m_pdata);  
  for (int k = 0; k<N_SYS; k++)
	{
	  double tau_k = Tau_Rhs_sp(k, P_sys, Fe_n->m_pdata, S_n.m_pdata);
    Vec_v(Fs[TENSOR_tau],k+1) = tau_k;
    Vec_v(Fs[TENSOR_gamma_dot],k+1) = gamma_Rate_PL(&Params,g_n,tau_k);
  }

  double g_Rhs = g_Rate_VK(&Params,&Struc,g_n, Fs[TENSOR_gamma_dot].m_pdata);
  state_var[VAR_g_np1] =  g_n + dt*g_Rhs;
  state_var[VAR_L_np1] =  L_np1;

  Matrix_cleanup(Fe_I);
  Matrix_cleanup(Fp_np1k);
  Matrix_cleanup(S_n);
  
////////////////////////////////////////////////////////////////////////////  
  double dd;
  Matrix_det(*pFnp1, dd);
//  Matrix_print(*pFnp1);  
  printf("g_n+1: %e, dd: %e, L_n+1: %e\n", state_var[VAR_g_np1], dd, state_var[VAR_L_np1]);
/////////////////////////////////////////////////////////////////////////////
  
	return err;  
  
}

int constitutive_model_update_plasticity(Matrix_double *pFnp1,
                                         Matrix_double *eFn,
                                         Matrix_double *pFn,
                                         Constitutive_model *m, double dt)
{
  int err = 0;  

  switch(m->param->type)
  {
    case HYPER_ELASTICITY:
      Matrix_eye(*pFnp1,3);
      return err;
    
    case CRYSTAL_PLASTICITY:
    {  
      Matrix_AeqB(m->vars.Fs[1],1.0,*pFn);

/////////////////////////////////////////////////////////////////////  
     dt = 1.0;
      Matrix_eye(*pFnp1,3);
/////////////////////////////////////////////////////////////////////        
                                                 
      integration_ip(pFnp1, m, eFn, dt);
/////////////////////////////////////////////////////////////////////      
    //  Matrix_print(*pFnp1);
      Matrix_eye(*pFnp1,3);
/////////////////////////////////////////////////////////////////////
      return err;
    }  
    case BPA_PLASTICITY:
      return err;
    
    default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",m->param->type);
    err++;
    break;
  }      
  return err;  
}

int constitutive_model_update_dMdu(Constitutive_model *m, Matrix_double *dMdu, Matrix_double *Fe, Matrix_double *S, Matrix_double *L, Matrix_double *Grad_du, double dt)
{
  int err = 0;
  switch(m->param->type) 
  {
  case HYPER_ELASTICITY:
    Matrix_init(*dMdu, 0.0);
    return err;
  case CRYSTAL_PLASTICITY:
  {
    void *ctx;
    double J;
    Matrix(double) C;
    Matrix_construct_redim(double,C,3,3); 
    Matrix_AxB(C, 1.0, 0.0, *Fe, 1, *Fe, 0);
    Matrix_det(*Fe, J);
    err += plasticity_model_ctx_build(&ctx, C.m_pdata, &J);
    err += compute_dMdu(m,dMdu,Grad_du,Fe,S,L,dt);
    err += plasticity_model_ctx_destroy(&ctx);
    Matrix_cleanup(C);
    /////////////////////////////////////////////
    Matrix_init(*dMdu, 0.0);    
    ////////////////////////////////////////////
    break;    
  }  
  case BPA_PLASTICITY:
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",m->param->type);
    err++;
    break;
  }
  return err; 
}


int build_model_parameters_list(Model_parameters **param_list,
                                const int n_mat,
                                const MATGEOM_1 *p_mgeom,
                                const HOMMAT *hmat_list)
{
  int err = 0;
  if (n_mat <= 0) return 1;
  (*param_list) = malloc(n_mat*sizeof(**param_list));

  /* for now set all model type to HYPER_ELASTIC. See issue #22 */
  int type = HYPER_ELASTICITY;
  for (int i = 0; i < n_mat; i++) {
    err += model_parameters_construct(& ((*param_list)[i]) );
    err += model_parameters_initialize(& ((*param_list)[i]),
                                       NULL,
                                       p_mgeom,
                                       hmat_list + i,
                                       type);
  }

  return err;
}

int destroy_model_parameters_list(const int n_mat,
                                  Model_parameters *param_list)
{
  int err = 0;
  if (param_list == NULL) return 0;

  for (int i = 0; i < n_mat; i++) {
    err += model_parameters_destroy(param_list + i);
  }
  free(param_list);
  return 0;
}

int init_all_constitutive_model(EPS *eps,
                                const int ne,
                                const ELEMENT *elem,
                                const Model_parameters *param_list)
{
  int err = 0;
  if (ne <= 0) return 1;

  for (int i = 0; i < ne; i++) {
    /* aliases */
    EPS *p_eps = eps + i;
    const ELEMENT *p_el = elem + i;
    const Model_parameters *p_param = param_list + (p_el->mat[2]);

    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++) {
      err += constitutive_model_initialize((p_eps->model) + j, p_param);
    }
  }
  return err;
}

int constitutive_model_test(const HOMMAT *hmat)
{
 
  Constitutive_model m;
  Model_parameters p;
  
  constitutive_model_construct(&m);
  model_parameters_construct(&p);  
  model_parameters_initialize(&p, NULL, NULL, hmat, CRYSTAL_PLASTICITY);
  constitutive_model_initialize(&m, &p);
    
  int err = 0;
  int N_SYS = 12;
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

  printf("%e\n", Props.Lame_I);
/*	
	Props.Lame_I =75600.0;
	Props.Lame_II =26100.0;
	Props.Modulus_Elastic = (Props.Lame_II*((3.0*Props.Lame_I+2.0*Props.Lame_II)/(Props.Lame_I+Props.Lame_II)));
	Props.Modulus_Shear = Props.Lame_II;
	Props.Poissons_Ratio= ((Props.Modulus_Elastic/(2.0*Props.Lame_II))-1.0);
	Props.Modulus_Bulk = (Props.Modulus_Elastic/(3.0-6.0*Props.Poissons_Ratio));*/
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
	double dt = 0.001;
	double Load_History = 1.0;
	int Load_Type = UNIAXIAL_COMPRESSION;	
	
	int Num_Steps=(int)(ceil(T_Final/dt));  
  
  int j_max = 3;
  Matrix(double) P_sys;  
  Matrix_construct_redim(double, P_sys, N_SYS*j_max*j_max, 1);

	for (int k = 0; k<N_SYS; k++)
	{
		for (int j = 0; j<j_max; j++)
		{
			for (int i = 0; i<j_max; i++)
				P_sys.m_pdata[Index_3D(k,j,i,j_max,j_max)] = 0.0;
		}
		P_sys_sp(k, P_sys.m_pdata);
	}

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

	Matrix_eye(F_n,3);
	Matrix_init(F_np1, 0.0);
	F_Implicit_sp(dt, F_n.m_pdata, L.m_pdata, F_np1.m_pdata);

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


	printf("time \t Effective Stress \n");
	double norm_T;
	Matrix_ddot(T_n,T_n,norm_T);
	double s_eff=sqrt(3.0/2.0)*sqrt(norm_T);

  double t = T_Initial;
	printf("%e \t %e\n",t,s_eff);
	
	for (int kk=1; kk<Num_Steps+1; kk++)
	{
		Matrix_init(Fp_np1, 0.0);
		L_np1 = 0.0;

		Staggered_NewtonRapson_Testing_sp(&Props,&Param,&Struc,&Solver,
				                       t,dt,P_sys.m_pdata,g_n,Fp_n.m_pdata,F_np1.m_pdata,
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
		s_eff=sqrt(3.0/2.0)*sqrt(norm_T);
    printf("%e \t %e\n",t,s_eff);
  
		Matrix_init(Tau_Array, 0.0);
		Matrix_init(gamma_RateArray, 0.0);
		double gmdot_tmp = 0.0;
		for (int k = 0; k<N_SYS; k++)
		{
		  double tau_k = Tau_Rhs_sp(k, P_sys.m_pdata, Fe_n.m_pdata, S_n.m_pdata);
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

  printf("integration test is successfull\n");
	
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
  Matrix_cleanup(P_sys);

  constitutive_model_destroy(&m);
  model_parameters_destroy(&p);  	  
  return err;
}

