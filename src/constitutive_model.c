/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "constitutive_model.h"

#include "material.h"
#include "hommat.h"
#include "matgeom.h"
#include "PGFEM_io.h"
#include "data_structure_c.h"

/* model headers */
#include "plasticity_model_none.h"
#include "plasticity_model.h"

#include "CM.h"
#include "physics_helper.h"

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
  Matrix_det(m->vars.Fs[0], J);

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
      for(int J=1; J<=3; J++)
      {
        for(int P=1; P<=3; P++)
        {
          for(int Q=1; Q<=3; Q++)
          {
            Tns4_v(CIoxCI,I,J,P,Q) = Mat_v(CI,I,J)*Mat_v(CI,P,Q);
            Tns4_v(SoxS,I,J,P,Q) = Mat_v(CI,I,J)*Mat_v(CI,P,Q);            
            Tns4_v(CICI,I,J,P,Q) = Mat_v(CI,I,P)*Mat_v(CI,Q,J);
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
  Matrix_construct_redim(double,Le,81,1);     
  Matrix_init_w_array(Fe,3,3,F);  
  constitutive_model_update_elasticity(Props->cm, &Fe, 0.0, &Le, &Se, 1);

  for(int a=0;a<81;a++)
    L[a] = Le.m_pdata[a];
    
  Matrix_cleanup(Se);
  Matrix_cleanup(Le);  
  Matrix_cleanup(Fe);   
}

int integration_ip(Matrix_double *pFnp1, Constitutive_model *m, Matrix_double *Fe_n, Matrix_double *S_n, double dt)
{
  int err = 0;
  int N_SYS=Vec_v((m->vars).state_vars[0], VAR_Ns+1);;
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
	Props.Lame_I =1.2/*Lambda*/;
	Props.Lame_II =1.2/*Shear*/;
	Props.Modulus_Elastic = (Props.Lame_II*((3.0*Props.Lame_I+2.0*Props.Lame_II)/(Props.Lame_I+Props.Lame_II)));
	Props.Modulus_Shear = Props.Lame_II;
	Props.Poissons_Ratio= ((Props.Modulus_Elastic/(2.0*Props.Lame_II))-1.0);
	Props.Modulus_Bulk = (Props.Modulus_Elastic/(3.0-6.0*Props.Poissons_Ratio));
	Props.use_hyperelastic = 0;
	Props.cm = m;
	Props.compute_elastic_stress = elastic_stress;
	Props.compute_elastic_tangent = elastic_tangent;

	/*--------MaterialParameters_Settings--------*/
	MaterialParameters Params;
	Params.Model_Type = PL_VK;
	Params.Parameters_Count = 7;
	Params.gam0dot = 1.0;
	Params.m_matl = 0.05;
	Params.gs0 = 0.330;
	Params.gamsdot = 50000000000.0;
	Params.w = 0.005;
	Params.G0 = 0.2;
	Params.g0 = 0.21;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
  Vec_v((m->vars).state_vars[0], VAR_g_n+1) = 210.0;
  Vec_v((m->vars).state_vars[0], VAR_L_n+1) = 0.001;
  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
	
	int j_max = 3;

	/*needed param info*/
		
	double *P_nsysR = (double *) malloc(sizeof(double)*N_SYS*j_max*j_max);

	for (int k = 0; k<N_SYS; k++)
	{
		for (int j = 0; j<j_max; j++)
		{
			for (int i = 0; i<j_max; i++)
				P_nsysR[Index_3D(k,j,i,j_max,j_max)] = 0.0;
		}
		P_sys_sp(k, P_nsysR);
	}
	
  Matrix(double) Fe_I,Fp_np1k;
	
  Matrix_construct_redim(double,Fe_I,3,3);
  Matrix_construct_redim(double,Fp_np1k,3,3);
  
  double g_n = Vec_v((m->vars).state_vars[0], VAR_g_n+1);
  // F_np1 = m->vars.Fs[TENSOR_F]: F
  // F_pn = m->vars.Fs[TENSOR_Fp]: Fp
  
  double L_np1k = Vec_v((m->vars).state_vars[0], VAR_L_n+1);
  double L_np1 = 0.0;
  Matrix_inv(*Fe_n, Fe_I);
  Matrix_AxB(Fp_np1k,1.0,0.0,Fe_I,0,m->vars.Fs[0],0);
  
	Staggered_NewtonRapson_Testing_sp(&Props,&Params,&Struc,&Solver,0.0,dt,
									   P_nsysR,g_n,m->vars.Fs[TENSOR_Fp].m_pdata,m->vars.Fs[TENSOR_F].m_pdata
									   ,L_np1k,&L_np1,Fp_np1k.m_pdata,pFnp1->m_pdata);

  Matrix(double) gamma_RateArray,Tau_Array;	
  Matrix_construct_redim(double,Tau_Array,N_SYS,1);
  Matrix_construct_redim(double,gamma_RateArray,N_SYS,1);
  
//  Matrix_print(*pFnp1);  
//  printf("L_n+1: %e\n", L_np1);

  for (int k = 0; k<N_SYS; k++)
	{
    Vec_v(Tau_Array,k+1) = Tau_Rhs_sp(k, P_nsysR, Fe_n->m_pdata, S_n->m_pdata);
    Vec_v(gamma_RateArray,k+1) = gamma_Rate_PL(&Params,g_n,Vec_v(Tau_Array,k+1));
  }

  double g_Rhs = g_Rate_VK(&Params,&Struc,g_n, gamma_RateArray.m_pdata);
  Vec_v((m->vars).state_vars[0], VAR_g_n+1) =  g_n + dt*g_Rhs;
  Vec_v((m->vars).state_vars[0], VAR_L_n+1) =  L_np1;

  Matrix_cleanup(Fe_I);
  Matrix_cleanup(Fp_np1k);
  
  Matrix_cleanup(gamma_RateArray);
  Matrix_cleanup(Tau_Array);
	return err;  
  
}

int constitutive_model_update_plasticity(Matrix_double *pFnp1,Matrix_double *eFn,Matrix_double *pFn,Constitutive_model *m, double dt)
{
  int err = 0;  

  switch(m->param->type)
  {
    case HYPER_ELASTICITY:
      Matrix_eye(*pFnp1,3);
      return err;
    
    case CRYSTAL_PLASTICITY:
    {  
      Matrix(double) S;
      Matrix_construct_redim(double,S,3,3);       

      Matrix_AeqB(m->vars.Fs[1],1.0,*pFn);

/////////////////////////////////////////////////////////////////////  
     dt = 1.0;
      Matrix_eye(*pFnp1,3);
/////////////////////////////////////////////////////////////////////        
      constitutive_model_update_elasticity(m,eFn,dt,NULL,&S,0);                                           
      integration_ip(pFnp1, m, eFn, &S, dt);
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
    Matrix_AxB(C, 1.0, 0.0, m->vars.Fs[0], 1, m->vars.Fs[0], 0);
    Matrix_det(m->vars.Fs[0], J);
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
