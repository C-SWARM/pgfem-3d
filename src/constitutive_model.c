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
#include "femlib.h"

#include "CM.h"

#ifndef _Matrix_double
Define_Matrix(double);
#define _Matrix_double 1
#endif

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

int model_var_info_print(FILE *f,
                         const Model_var_info * info)
{
  int err = 0;
  fprintf(f,"F names: ");
  for(int i = 0, e = info->n_Fs; i < e; i++) fprintf(f,"%s ",info->F_names[i]);
  fprintf(f,"\nVar names: ");
  for(int i = 0, e = info->n_vars; i < e; i++) fprintf(f,"%s ",info->var_names[i]);
  fprintf(f,"\n");
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
  p->Psys = NULL;
  p->N_SYS = 0;
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
  p->Psys = NULL;
  p->N_SYS = 0;
  switch(type) {
  case HYPER_ELASTICITY:
    err += plasticity_model_none_initialize(p);
    break;
  case CRYSTAL_PLASTICITY:
  {
    p->Psys = malloc(sizeof(*(p->Psys)));
    Matrix_construct(double,*(p->Psys));
    err += plasticity_model_initialize(p);
    break;
  }
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
  p->N_SYS = 0;
  
  /* drop function pointers */
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_dudj = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->get_var_info = NULL;
  if(p->Psys){
    Matrix_cleanup(*(p->Psys));
    free(p->Psys);
  }
  p->Psys = NULL;
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
            Tns4_v(SoxS,I,JJ,P,Q) = Mat_v(*S,I,JJ)*Mat_v(*S,P,Q);            
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

int constitutive_model_update_plasticity(Matrix_double *pFnp1,
                                         Matrix_double *Fnp1,
                                         Matrix_double *eFn,
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
      plasticity_model_integration_ip(pFnp1, m, Fnp1, eFn, dt);
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

int constitutive_model_update_dMdu(Constitutive_model *m, Matrix_double *dMdu, Matrix_double *eFn, Matrix_double *eFnp1, Matrix_double *M, 
                                   Matrix_double *S, Matrix_double *L, Matrix_double *Grad_du, double dt)
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
    Matrix_AxB(C, 1.0, 0.0, *eFnp1, 1, *eFnp1, 0);
    Matrix_det(*eFnp1, J);
    err += plasticity_model_ctx_build(&ctx, C.m_pdata, &J);
    err += compute_dMdu(m,dMdu,Grad_du,eFn,eFnp1,M,S,L,dt);
    err += plasticity_model_ctx_destroy(&ctx);
    Matrix_cleanup(C);
        
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
////////////////////////////////////////////////////////////////////
  type = CRYSTAL_PLASTICITY; // for the test
////////////////////////////////////////////////////////////////////  
  
  for (int i = 0; i < n_mat; i++) {
    err += model_parameters_construct(&((*param_list)[i]) );
    err += model_parameters_initialize(&((*param_list)[i]),
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

int read_constitutive_model_parameters(EPS *eps,
                                const int ne,
                                const ELEMENT *elem,
                                const int n_mat,
                                Model_parameters *param_list)
{ 
  int err = 0;
  int type = HYPER_ELASTICITY;
////////////////////////////////////////////////////////////////////
  type = CRYSTAL_PLASTICITY; // for the test
////////////////////////////////////////////////////////////////////   
  switch(type)
  {
    case HYPER_ELASTICITY:
      return err;
    
    case CRYSTAL_PLASTICITY:
    { 
      for(int a=0; a<n_mat; a++)
      {
        Matrix(double) *P = param_list[a].Psys;
        param_list[a].N_SYS = plasticity_model_slip_system(P);
      } 
      
      for(int a=0; a<ne; a++)
      {
        long n_ip = 0;
        int_point(elem[a].toe,&n_ip);
        for(int ip=0; ip<n_ip; ip++)
        {
          Constitutive_model *m = &(eps[a].model[ip]);
          plasticity_model_read_parameters(m); 
        }
      }
      return err;
    }  
    case BPA_PLASTICITY:
      return err;
    
    default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",type);
    err++;
    break;
  }      
  return err;   
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

int constitutive_model_update_time_steps(EPS *eps, const int ne, const ELEMENT *elem)
{
  int err = 0;
  if (ne <= 0) return 1;

  for (int i = 0; i < ne; i++) 
  {
    const ELEMENT *p_el = elem + i;
    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++)     
    {
      Constitutive_model *m = &(eps[i].model[j]);
      Matrix(double) *Fs = (m->vars).Fs;
      double *state_var = (m->vars).state_vars[0].m_pdata;
      Matrix_AeqB(Fs[TENSOR_Fn], 1.0,Fs[TENSOR_Fnp1]);
      Matrix_AeqB(Fs[TENSOR_pFn],1.0,Fs[TENSOR_pFnp1]);
      state_var[VAR_g_n] = state_var[VAR_g_np1];
      state_var[VAR_L_n] = state_var[VAR_L_np1];
    }
  }
  return err;  
}

int constitutive_model_update_time_steps_test(ELEMENT *elem, NODE *node, HOMMAT *hommat, EPS *eps, 
                                        const int ne, const int nn, const int ndofn,
                                        double* r, double dt)
{
  int nsd = 3;
  int total_Lagrangian = 0;
  int err = 0;
  if (ne <= 0) return 1;

  Matrix(double) Fn, pFn, pFnI, eFn, Fnp1, pFnp1, Fr;
  Matrix_construct_redim(double,Fn,3,3);
  Matrix_construct_redim(double,pFn,3,3);
  Matrix_construct_redim(double,pFnI,3,3);
  Matrix_construct_redim(double,eFn,3,3);    
  Matrix_construct_redim(double,Fnp1,3,3); 
  Matrix_construct_redim(double,pFnp1,3,3);      
  Matrix_construct_redim(double,Fr,3,3);    
    
  for (int e = 0; e < ne; e++) 
  {
    
    int intg_order = 1;    
    if (PGFEM3D_DEV_TEST)
      intg_order = 0;

    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe, e, elem, node, intg_order,total_Lagrangian);
    int nne = fe.nne;
    
    Matrix(double) u;  
    Matrix_construct_init(double,u,nne*nsd,1,0.0);
                        
    for(int a = 0; a<nne; a++)
    {      
      int nid = Vec_v(fe.node_id, a+1);
      for(int b=0; b<nsd; b++)
      {
        Vec_v(u, a*nsd+b+1) = r[nid*ndofn + b];
      }
    }    
    
    for (int ip = 1; ip <=fe.nint; ip++)     
    {
      FEMLIB_elem_basis_V(&fe, ip);
      FEMLIB_update_shape_tensor(&fe);
      FEMLIB_update_deformation_gradient(&fe,nsd,u.m_pdata,Fr);
            
      Constitutive_model *m = &(eps[e].model[ip-1]);
      Matrix(double) *Fs = (m->vars).Fs;
      double *state_var = (m->vars).state_vars[0].m_pdata;
      
   
      Matrix_AeqB(Fn,1.0,Fs[TENSOR_Fn]);
      Matrix_AeqB(pFn,1.0,Fs[TENSOR_pFn]);      
        
      Matrix_inv(pFn, pFnI);
      Matrix_AxB(eFn,1.0,0.0,Fn,0,pFnI,0); 
    
      // --> update plasticity part
      if(total_Lagrangian)
      {    
        Matrix_AeqB(Fnp1,1.0,Fr);  // Fn+1 
      }
      else
      {
        Matrix_AxB(Fnp1,1.0,0.0,Fr,0,Fn,0);  // Fn+1    
      }   
    
      constitutive_model_update_plasticity(&pFnp1,&Fnp1,&eFn,m,dt);     


      Matrix_AeqB(Fs[TENSOR_Fn], 1.0,Fs[TENSOR_Fnp1]);
      Matrix_AeqB(Fs[TENSOR_pFn],1.0,Fs[TENSOR_pFnp1]);
      state_var[VAR_g_n] = state_var[VAR_g_np1];
      state_var[VAR_L_n] = state_var[VAR_L_np1];
    }
  }
  
  
  /*********************/
  /* Coordinate update */
  /*********************/
   for(int n = 0;n<nn; n++)
   {
    for(int a=0;a<nsd;a++)
    {
      int II = node[n].id[a];
      if (II != 0)
      {
	      if (a == 0) node[n].x1 = node[n].x1_fd + r[n*ndofn + a];
	      else if (a == 1) node[n].x2 = node[n].x2_fd + r[n*ndofn + a];
	      else if (a == 2) node[n].x3 = node[n].x3_fd + r[n*ndofn + a];
      }
    }
  }/* end n < nn */  
  
  Matrix_cleanup(Fn);
  Matrix_cleanup(pFn);
  Matrix_cleanup(pFnI);
  Matrix_cleanup(eFn);    
  Matrix_cleanup(Fnp1);
  Matrix_cleanup(pFnp1);  
  Matrix_cleanup(Fr);   
  return err;  
}

int constitutive_model_test(const HOMMAT *hmat, Matrix_double *L_in, int Print_results)
{
  int err = plasticity_model_test(hmat, L_in, Print_results);
  return err;
}
