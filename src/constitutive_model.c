/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */


#include "CM.h"
#include "constitutive_model.h"
#include "plasticity_model_none.h"
#include "plasticity_model.h"
#include "PGFEM_io.h"
#include "data_structure_c.h"
#include "elem3d.h"
#include "femlib.h"
#include "index_macros.h"

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
  p->get_pF = NULL;
  p->get_pFn = NULL;
  p->get_eF = NULL;
  p->get_eFn = NULL;
  p->destroy_ctx = NULL;
  p->compute_dMdu = NULL;
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
  p->get_Fn = NULL;
  p->get_pF = NULL;
  p->get_pFn = NULL;
  p->get_eF = NULL;
  p->get_eFn = NULL;
  p->destroy_ctx = NULL;
  p->compute_dMdu = NULL;
  if(p->Psys){
    Matrix_cleanup(*(p->Psys));
    free(p->Psys);
  }
  p->Psys = NULL;
  /* reset counters/flags */
  p->type = -1;

  return err;
}

int constitutive_model_update_elasticity(const Constitutive_model *m,
                                         const Matrix_double *F,
                                         const double dt,
                                         Matrix_double *L,
                                         Matrix_double *S,
                                         const int compute_stiffness)
{
  int err = 0;
  void *ctx;
  const Model_parameters *func = m->param;
  double J;
  Matrix(double) Fe,C, CI;    
  Matrix_construct_redim(double,C,3,3); 
  Matrix_construct_redim(double,CI,3,3);     
  Matrix_construct_redim(double, Fe, 3,3);

  /* this is TEMPORARY until expand hyperelastic to maintain an
     elastic deformation gradient */
  switch (m->param->type){
  case HYPER_ELASTICITY:
    Matrix_AeqB(Fe,1.0,*F);
    break;
  default:
    func->get_eF(m,&Fe);
    break;
  }

  Matrix_AxB(C, 1.0, 0.0, Fe, 1, Fe, 0);
  Matrix_inv(C,CI);        
  Matrix_det(Fe, J);

  switch(m->param->type)
  {
    case HYPER_ELASTICITY:
      err += plasticity_model_none_ctx_build(&ctx, C.m_pdata, &J);
      break;
    case CRYSTAL_PLASTICITY:
      
      err += plasticity_model_ctx_build(&ctx, F->m_pdata, dt);
      
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
  err += func->compute_dev_stress(m, ctx, S);
  err += func->compute_dudj(m,ctx,&dudj);    
  Matrix_AplusB(*S, kappa*J*dudj,CI,1.0,*S);
  //compute stiffness
  if(compute_stiffness)
  {  
    Matrix(double) CIoxCI, CICI, SoxS;    
    Matrix_construct_redim(double,CIoxCI,81,1);
    Matrix_construct_redim(double,CICI,81,1);             
    Matrix_construct_redim(double,SoxS,81,1);
    
    err += func->compute_dev_tangent(m, ctx, L);              
    err += func->compute_d2udj2(m,ctx,&d2udj2);
  
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

  func->destroy_ctx(&ctx);  
  Matrix_cleanup(C);
  Matrix_cleanup(CI); 
  return err; 
}

int constitutive_model_update_plasticity(Matrix_double *pFnp1,
                                         const Matrix_double *Fnp1,
                                         const Matrix_double *eFn,
                                         Constitutive_model *m,
                                         const double dt)
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

int build_model_parameters_list(Model_parameters **param_list,
                                const int n_mat,
                                const MATGEOM_1 *p_mgeom,
                                const HOMMAT *hmat_list, const int type)
{
  int err = 0;
  if(type<0)
    return err;  
  if (n_mat <= 0) return 1;
  
  (*param_list) = malloc(n_mat*sizeof(**param_list));

  /* for now set all model type to HYPER_ELASTIC. See issue #22 */
  
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
                                Model_parameters *param_list,
                                const int type)
{ 
  int err = 0;
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
  if(param_list==NULL)
    return err;

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

int constitutive_model_update_time_steps(EPS *eps,
                                         const int ne,
                                         const ELEMENT *elem)
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
      m->param->update_state_vars(m);
    }
  }
  return err;  
}

int constitutive_model_update_time_steps_test(ELEMENT *elem,
                                              NODE *node,
                                              HOMMAT *hommat,
                                              EPS *eps,
                                              const int ne,
                                              const int nn,
                                              const int ndofn,
                                              double* r,
                                              double dt)
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
      m->param->update_state_vars(m);
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
               if (a == 0)      node[n].x1 += r[n*ndofn + a];
               else if (a == 1) node[n].x2 += r[n*ndofn + a];
               else if (a == 2) node[n].x3 += r[n*ndofn + a];
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

int stiffness_el_hyper_elasticity(double *lk,
        const int ii,
        const int ndofn,
        const int nne,
        const int nsd,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        MATGEOM matgeom,
        const long *nod,
        const NODE *node,
        double dt,
        SIG *sig,
        EPS *eps,
        const SUPP sup,
        double *r_e)
{
  int err = 0;
  int total_Lagrangian = 1;
  int compute_stiffness = 1;
  
  double *u;
  u = (double *) malloc(sizeof(double)*nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  Matrix(double) F,L,S,ST_ab,ST_wg,AA,BB,CC,sAA,sBB,sCC,LsBB;

  Matrix_construct_redim(double, F   ,3,3);
  Matrix_construct_redim(double, L   ,81,1);  
  Matrix_construct_redim(double, S   ,3,3);    
  Matrix_construct_redim(double, AA  ,3,3);
  Matrix_construct_redim(double, BB  ,3,3);
  Matrix_construct_redim(double, CC  ,3,3);  
  Matrix_construct_redim(double, sAA ,3,3);
  Matrix_construct_redim(double, sBB ,3,3);  
  Matrix_construct_redim(double, sCC ,3,3);    
  Matrix_construct_redim(double, LsBB,3,3);
  
  Matrix_construct(double,ST_ab);
  Matrix_construct(double,ST_wg);
  
  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);   
  
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,F);
    
    Constitutive_model *m = &(eps[ii].model[ip-1]);
    
    Matrix_init(L,0.0);
    Matrix_init(S,0.0);    
    
    constitutive_model_update_elasticity(m,&F,dt,&L,&S,compute_stiffness);
    
    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];

        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                    3,3,3,1.0,F.m_pdata,3,ptrST_ab,3,0.0,AA.m_pdata,3);
        Matrix_symmetric(AA,sAA);
        
        for(int w=0; w<nne; w++)
        {
          for(int g=0; g<nsd; g++)
          { 
            const double* const ptrST_wg = &(fe.ST)[idx_4_gen(w,g,0,0,
                                                      nne,nsd,nsd,nsd)];             
            for(int _i=0; _i<nsd; _i++)
            {
                for(int _j=0; _j<nsd; _j++)
                {
                    BB.m_pdata[idx_2(_i,_j)] = CC.m_pdata[idx_2(_i,_j)] = 0.0;
                    for(int _k=0; _k<nsd; _k++)
                    {
                        /* BB = F' ST_wg */
                        BB.m_pdata[idx_2(_i,_j)] += F.m_pdata[idx_2(_k,_i)]*ptrST_wg[idx_2(_k,_j)];
                        
                        /* CC = ST_wg' ST_ab */
                        CC.m_pdata[idx_2(_i,_j)] += ptrST_wg[idx_2(_k,_i)]*ptrST_ab[idx_2(_k,_j)];
                    }
                }
            }           
            
            Matrix_symmetric(BB,sBB);
            Matrix_symmetric(CC,sCC);
            
            Matrix_init(LsBB,0.0);
            for(int _i=0; _i<9; _i++)
            {
              if(fabs(sBB.m_pdata[_i]) < 1.0e-15) continue;
              for(int _j=0; _j<9; _j++)
                LsBB.m_pdata[_j] += L.m_pdata[9*_j+_i]*sBB.m_pdata[_i];
            }            
            
            //Matrix_Tns4_dd_Tns2(LsBB,L,sBB);

            const int lk_idx = idx_K(a,b,w,g,nne,nsd);  
            for(int n=0; n<9; n++)
            {                 
              lk[lk_idx] += fe.detJxW*(sAA.m_pdata[n]*LsBB.m_pdata[n]
                                       + S.m_pdata[n]*sCC.m_pdata[n]);
            }
          }
        }
      }
    }
  }
  free(u);
  
  Matrix_cleanup(F);
  Matrix_cleanup(L);  
  Matrix_cleanup(S);    
  Matrix_cleanup(ST_ab);
  Matrix_cleanup(ST_wg);
  Matrix_cleanup(AA);
  Matrix_cleanup(BB);
  Matrix_cleanup(CC);  
  Matrix_cleanup(sAA);
  Matrix_cleanup(sBB);  
  Matrix_cleanup(sCC);    
  Matrix_cleanup(LsBB);
    
  FEMLIB_destruct(&fe);
  
  return err;
}

int residuals_el_hyper_elasticity(double *f,
        const int ii,
        const int ndofn,
        const int nne,
        const int nsd,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        MATGEOM matgeom,
        const long *nod,
        const NODE *node,
        double dt,
        SIG *sig,
        EPS *eps,
        const SUPP sup,
        double *r_e)
{
  int err = 0;
  int total_Lagrangian = 1;
  int compute_stiffness = 0;  
  double *u;
  u = (double *) malloc(sizeof(double)*nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  Matrix(double) F,S,ST_ab,AA,sAA;

  Matrix_construct_redim(double, F   ,3,3);
  Matrix_construct_redim(double, S   ,3,3);    
  Matrix_construct_redim(double, AA  ,3,3);
  Matrix_construct_redim(double, sAA ,3,3);
  
  Matrix_construct(double,ST_ab);
  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);      

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,F);
    
    Constitutive_model *m = &(eps[ii].model[ip-1]);
    Matrix_init(S,0.0);    
    
    constitutive_model_update_elasticity(m,&F,dt,NULL,&S,compute_stiffness);
              
    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];


        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                    3,3,3,1.0,F.m_pdata,3,ptrST_ab,3,0.0,AA.m_pdata,3);
//        Matrix_init_w_array(ST_ab,3,3,ptrST_ab);
//        Matrix_AxB(AA,1.0,0.0,F,1,ST_ab,0); 
        Matrix_symmetric(AA,sAA);

//        Matrix_init_w_array(ST_ab,3,3,ptrST_ab);
//        Matrix_AxB(AA,1.0,0.0,F,1,ST_ab,0); 
//        Matrix_symmetric(AA,sAA);
        
        int fe_id = a*ndofn + b;              
                
        for(int n=0; n<9; n++)
          f[fe_id] += fe.detJxW*sAA.m_pdata[n]*S.m_pdata[n];
      }
    }
  }
  
  free(u);
  Matrix_cleanup(F);
  Matrix_cleanup(S);
  Matrix_cleanup(AA);
  Matrix_cleanup(sAA);
  Matrix_cleanup(ST_ab);  
    
  FEMLIB_destruct(&fe);
  return err;
}        

int stiffness_el_crystal_plasticity(double *lk,
                                    const int ii,
                                    const int ndofn,
                                    const int nne,
                                    const int nsd,
                                    const ELEMENT *elem,
                                    const HOMMAT *hommat,
                                    MATGEOM matgeom,
                                    const long *nod,
                                    const NODE *node,
                                    double dt,
                                    SIG *sig,
                                    EPS *eps,
                                    const SUPP sup,
                                    double *r_e)
{
  static const int total_Lagrangian = 0;
  
  int err = 0;

  double *u = malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = malloc(sizeof(*dMdu_all) * 9 * nne * nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  Matrix(double) Fn, Fr, Fnp1;
  Matrix(double) pFn;
  Matrix(double) FreFn,eFnp1,pFnp1,pFnp1_I,L,dMdu,S;  
  Matrix(double) pFnI, eFn, M, eFnM;  
  Matrix(double) ST_ab, ST_wg, AA, BB, CC;
  Matrix(double) sAA, sBB, sCC;
  Matrix(double) MTeFnT_sAA, MTeFnT_sAA_eFn,MTeFnT_sAA_eFnM,FrTFr,MTeFnT_FrTFr,MTeFnT_FrTFreFn,MTeFnT_FrTFreFndMdu,dCdu,MTeFnT_sBB;
  Matrix(double) L_dCdu,MTeFnT_sCC,MTeFnT_sCC_eFnM,MTeFnT_sAA_eFndMdu,sMTeFnT_sAA_eFndMdu;
      
  Matrix_construct_redim(double,Fn ,3,3);
  Matrix_construct_redim(double,Fr ,3,3);
  Matrix_construct_redim(double,Fnp1 ,3,3);      

  Matrix_construct_redim(double,pFn ,3,3);
  
  Matrix_construct_redim(double,L ,81,1);  
  Matrix_construct_redim(double,dMdu ,3,3);
  Matrix_construct_redim(double,FreFn,3,3);      
  Matrix_construct_redim(double,eFnp1,3,3);
  Matrix_construct_redim(double,pFnp1,3,3);
  Matrix_construct_redim(double,pFnp1_I,3,3);    
  Matrix_construct_redim(double,S    ,3,3);     
  Matrix_construct_redim(double,pFnI,3,3); 
  Matrix_construct_redim(double, eFn,3,3); 
  Matrix_construct_redim(double,   M,3,3);
  Matrix_construct_redim(double,eFnM,3,3); 
  Matrix_construct(double,ST_ab);
  Matrix_construct(double,ST_wg);
  Matrix_construct_redim(double,AA,3,3);
  Matrix_construct_redim(double,BB,3,3);
  Matrix_construct_redim(double,CC,3,3);  
  Matrix_construct_redim(double,sAA, 3,3);
  Matrix_construct_redim(double,sBB, 3,3);  
  Matrix_construct_redim(double,sCC, 3,3);    
  Matrix_construct_redim(double,MTeFnT_sAA        ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFn     ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFnM    ,3,3);
  Matrix_construct_redim(double,FrTFr            ,3,3);
  Matrix_construct_redim(double,MTeFnT_FrTFr      ,3,3);
  Matrix_construct_redim(double,MTeFnT_FrTFreFn    ,3,3);
  Matrix_construct_redim(double,MTeFnT_FrTFreFndMdu,3,3);
  Matrix_construct_redim(double,dCdu             ,3,3);
  Matrix_construct_redim(double,MTeFnT_sBB        ,3,3);
  Matrix_construct_redim(double,L_dCdu           ,3,3);
  Matrix_construct_redim(double,MTeFnT_sCC        ,3,3);
  Matrix_construct_redim(double,MTeFnT_sCC_eFnM    ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFndMdu ,3,3);
  Matrix_construct_redim(double,sMTeFnT_sAA_eFndMdu,3,3);  
  
  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);
  int compute_stiffness = 1;      

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,Fr);
    
    Constitutive_model *m = &(eps[ii].model[ip-1]);
    //Matrix(double) *Fs = (m->vars).Fs;

    /* get a shortened pointer for simplified CM function calls */
    const Model_parameters *func = m->param;
    err += func->get_Fn(m,&Fn);
    err += func->get_pFn(m,&pFn);
    err += func->get_eFn(m,&eFn);
    
    // --> update plasticity part
    if(total_Lagrangian)
    { 
      Matrix(double) FnI;
      Matrix_construct_redim(double, FnI,3,3);
      Matrix_inv(Fn,FnI);
      Matrix_AeqB(Fnp1,1.0,Fr);  // Fn+1 
      Matrix_AxB(Fr,1.0,0.0,Fnp1,0,FnI,0);  // Fn+1          
      Matrix_cleanup(FnI);         
    }
    else
    {
      Matrix_AxB(Fnp1,1.0,0.0,Fr,0,Fn,0);  // Fn+1    
    }   
    Matrix_AxB(FrTFr,1.0,0.0,Fr,1,Fr,0); 

    /* need to have called the integration algorithm. This should be
       done OUTSIDE of the stiffness/tangent functions */
    constitutive_model_update_plasticity(&pFnp1,&Fnp1,&eFn,m,dt);
    err += func->get_pF(m,&pFnp1);

    Matrix_inv(pFnp1, pFnp1_I);
    Matrix_AxB(M,1.0,0.0,pFn,0,pFnp1_I,0);
    // <-- update plasticity part 

    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(S,0.0);    
    
    Matrix_AxB(eFnp1,1.0,0.0,Fnp1,0,pFnp1_I,0);
    constitutive_model_update_elasticity(m,&eFnp1,dt,&L,&S,compute_stiffness);
    // <-- update elasticity part

    // --> start computing tagent
    Matrix_AxB(FreFn,1.0,0.0,Fr,0,eFn,0);
    Matrix_AxB(eFnM,1.0,0.0,eFn,0,M,0);
    
    double Jn; Matrix_det(Fn, Jn);

    void *tmp_ctx = NULL;
    err += plasticity_model_ctx_build(&tmp_ctx, Fnp1.m_pdata, dt);
    err += func->compute_dMdu(m, tmp_ctx, fe.ST, nne, ndofn, dMdu_all);
    err += func->destroy_ctx(&tmp_ctx);

    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];
        Matrix_init_w_array(ST_ab,3,3,ptrST_ab);
        Matrix_AxB(AA,1.0,0.0,Fr,1,ST_ab,0); 
        Matrix_symmetric(AA,sAA);
        
        Matrix_AxB(MTeFnT_sAA,1.0,0.0,eFnM,1,sAA,0);
        Matrix_AxB(MTeFnT_sAA_eFn,1.0,0.0,MTeFnT_sAA,0,eFn,0);
        Matrix_AxB(MTeFnT_sAA_eFnM,1.0,0.0,MTeFnT_sAA_eFn,0,M,0);  

        for(int w=0; w<nne; w++)
        {
          for(int g=0; g<nsd; g++)
          { 
            const double* const ptrST_wg = &(fe.ST)[idx_4_gen(w,g,0,0,
                                                      nne,nsd,nsd,nsd)]; 
            const double* const ptr_dMdu_wg = &(dMdu_all[idx_4_gen(w,g,0,0,
                                                                   nne,nsd,nsd,nsd)]);

            Matrix_init_w_array(ST_wg,3,3,ptrST_wg); 
            Matrix_init_w_array(dMdu,3,3,ptr_dMdu_wg);

            Matrix_AxB(BB,1.0,0.0,Fr,1,ST_wg,0); 
            Matrix_symmetric(BB,sBB);
            Matrix_AxB(CC, 1.0,0.0,ST_ab,1,ST_wg,0);
            Matrix_symmetric(CC,sCC);
            
            // compute dCdu
            Matrix_AxB(MTeFnT_FrTFr,1.0,0.0,eFnM,1,FrTFr,0);
            Matrix_AxB(MTeFnT_FrTFreFn,1.0,0.0,MTeFnT_FrTFr,0,eFn,0);
            Matrix_AxB(MTeFnT_FrTFreFndMdu,1.0,0.0,MTeFnT_FrTFreFn,0,dMdu,0);            
            Matrix_symmetric(MTeFnT_FrTFreFndMdu,dCdu);
            
            Matrix_AxB(MTeFnT_sBB,1.0,0.0,eFnM,1,sBB,0);
            Matrix_AxB(dCdu,1.0,1.0,MTeFnT_sBB,0,eFnM,0);
            
            // compute MTeFnT_sAA_eFnM:L:dCdu
            Matrix_Tns4_dd_Tns2(L_dCdu,L,dCdu);
            double MTeFnT_sAA_eFnM_L_dCdu = 0.0;
            Matrix_ddot(MTeFnT_sAA_eFnM,L_dCdu,MTeFnT_sAA_eFnM_L_dCdu);
            
            // compute MTeFnT_sCC_eFnM
            Matrix_AxB(MTeFnT_sCC,1.0,0.0,eFnM,1,sCC,0);
            Matrix_AxB(MTeFnT_sCC_eFnM,1.0,0.0,MTeFnT_sCC,0,eFnM,0);
            
            // compute MTeFnT_sCC_eFnM:S
            double MTeFnT_sCC_eFnM_S = 0.0;
            Matrix_ddot(MTeFnT_sCC_eFnM,S,MTeFnT_sCC_eFnM_S);
            
            // compute MTeFnT_sAA_eFndMdu
            Matrix_AxB(MTeFnT_sAA_eFndMdu,1.0,0.0,MTeFnT_sAA_eFn,0,dMdu,0);    
            Matrix_symmetric(MTeFnT_sAA_eFndMdu, sMTeFnT_sAA_eFndMdu);        

            // compute MTeFnT_sAA_eFndMdu:S
            double sMTeFnT_sAA_eFndMdu_S = 0.0;            
            Matrix_ddot(sMTeFnT_sAA_eFndMdu,S,sMTeFnT_sAA_eFndMdu_S);
            
            const int lk_idx = idx_K(a,b,w,g,nne,nsd);  
                      
            lk[lk_idx] += 1.0/Jn*fe.detJxW*(MTeFnT_sAA_eFnM_L_dCdu + 2.0*sMTeFnT_sAA_eFndMdu_S + MTeFnT_sCC_eFnM_S);            
          }
        }
      }
    }
  }
  free(u);
  
  Matrix_cleanup(Fn);
  Matrix_cleanup(Fr);
  Matrix_cleanup(Fnp1);
  Matrix_cleanup(pFn);      
  Matrix_cleanup(L);  
  Matrix_cleanup(dMdu);
  Matrix_cleanup(FreFn);      
  Matrix_cleanup(eFnp1);
  Matrix_cleanup(pFnp1);
  Matrix_cleanup(pFnp1_I);
  Matrix_cleanup(S);    
  Matrix_cleanup(pFnI); 
  Matrix_cleanup(eFn); 
  Matrix_cleanup(M);
  Matrix_cleanup(eFnM); 
  Matrix_cleanup(ST_ab);
  Matrix_cleanup(ST_wg);
  Matrix_cleanup(AA);
  Matrix_cleanup(BB);
  Matrix_cleanup(CC);  
  Matrix_cleanup(sAA);
  Matrix_cleanup(sBB);  
  Matrix_cleanup(sCC);    
  Matrix_cleanup(MTeFnT_sAA);
  Matrix_cleanup(MTeFnT_sAA_eFn);
  Matrix_cleanup(MTeFnT_sAA_eFnM);
  Matrix_cleanup(FrTFr);
  Matrix_cleanup(MTeFnT_FrTFr);
  Matrix_cleanup(MTeFnT_FrTFreFn);
  Matrix_cleanup(MTeFnT_FrTFreFndMdu);
  Matrix_cleanup(dCdu);
  Matrix_cleanup(MTeFnT_sBB);
  Matrix_cleanup(L_dCdu);
  Matrix_cleanup(MTeFnT_sCC);
  Matrix_cleanup(MTeFnT_sCC_eFnM);
  Matrix_cleanup(MTeFnT_sAA_eFndMdu);
  Matrix_cleanup(sMTeFnT_sAA_eFndMdu);
    
  FEMLIB_destruct(&fe);
  free(dMdu_all);
 
  return err;
}        

int residuals_el_crystal_plasticity(double *f,
        const int ii,
        const int ndofn,
        const int nne,
        const int nsd,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        MATGEOM matgeom,
        const long *nod,
        const NODE *node,
        double dt,
        SIG *sig,
        EPS *eps,
        const SUPP sup,
        double *r_e)
{
  int total_Lagrangian = 0;
  
  int err = 0;
    
  double *u;
  u = (double *) malloc(sizeof(double)*nne*nsd);
  
  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }

  Matrix(double) Fn, Fr, Fnp1, pFn;
  Matrix_construct_redim(double,Fn ,3,3);
  Matrix_construct_redim(double,Fr ,3,3);
  Matrix_construct_redim(double,Fnp1 ,3,3); 
  
  Matrix_construct_redim(double,pFn ,3,3);     
       
  Matrix(double) eFnp1,pFnp1,pFnp1_I,L,S;  
  Matrix_construct_redim(double,eFnp1,3,3);  
  Matrix_construct_redim(double,pFnp1,3,3);
  Matrix_construct_redim(double,pFnp1_I,3,3);      
  Matrix_construct_redim(double,L ,81,1);  
  Matrix_construct_redim(double,S    ,3,3);
  
  Matrix(double) pFnI, eFn, M, eFnM;  
  Matrix_construct_redim(double,pFnI,3,3); 
  Matrix_construct_redim(double, eFn,3,3); 
  Matrix_construct_redim(double,   M,3,3);
  Matrix_construct_redim(double,eFnM,3,3);
            
  Matrix(double) ST_ab, AA, sAA;
  Matrix_construct(double,ST_ab);
  Matrix_construct_redim(double,AA,3,3);
  Matrix_construct_redim(double,sAA, 3,3);

  Matrix(double) MTeFnT_sAA, MTeFnT_sAA_eFnM;
  Matrix_construct_redim(double,MTeFnT_sAA        ,3,3);
  Matrix_construct_redim(double,MTeFnT_sAA_eFnM    ,3,3);

  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);      
  int compute_stiffness = 0;
   
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,Fr);
    
    Constitutive_model *m = &(eps[ii].model[ip-1]);
    Matrix(double) *Fs = (m->vars).Fs;    


//    if(updated_Lagrangian)
//    {
      Matrix_AeqB(Fn,1.0,Fs[TENSOR_Fn]);
      Matrix_AeqB(pFn,1.0,Fs[TENSOR_pFn]);      
//    }   
//    else
//    {  
//      Matrix_eye(Fn,3);
//      Matrix_eye(pFn,3);
//    } 
    

    Matrix_inv(pFn, pFnI);
    Matrix_AxB(eFn,1.0,0.0,Fn,0,pFnI,0); 
   
    // --> update plasticity part
    if(total_Lagrangian)
    {
      Matrix(double) FnI;
      Matrix_construct_redim(double, FnI,3,3);
      Matrix_inv(Fn,FnI);
      Matrix_AeqB(Fnp1,1.0,Fr);  // Fn+1 
      Matrix_AxB(Fr,1.0,0.0,Fnp1,0,FnI,0);  // Fn+1          
      Matrix_cleanup(FnI);          
    }
    else
    {
      Matrix_AxB(Fnp1,1.0,0.0,Fr,0,Fn,0);  // Fn+1    
    }      
    constitutive_model_update_plasticity(&pFnp1,&Fnp1,&eFn,m,dt);
    Matrix_AxB(M,1.0,0.0,pFnI,0,pFnp1,0);    
    // <-- update plasticity part

    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(S,0.0);    
    
    Matrix_inv(pFnp1, pFnp1_I);
    Matrix_AxB(eFnp1,1.0,0.0,Fnp1,0,pFnp1_I,0);
    constitutive_model_update_elasticity(m,&eFnp1,dt,&L,&S,compute_stiffness);
    // <-- update elasticity part
            
    Matrix_AxB(eFnM,1.0,0.0,eFn,0,M,0);
    double Jn; Matrix_det(Fn, Jn);
    
    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];
        Matrix_init_w_array(ST_ab,3,3,ptrST_ab);
        Matrix_AxB(AA,1.0,0.0,Fr,1,ST_ab,0); 
        Matrix_symmetric(AA,sAA);

        Matrix_AxB(MTeFnT_sAA,1.0,0.0,eFnM,1,sAA,0);
        Matrix_AxB(MTeFnT_sAA_eFnM,1.0,0.0,MTeFnT_sAA,0,eFnM,0); 
        double MTeFnT_sAA_eFnM_S = 0.0; 
        Matrix_ddot(MTeFnT_sAA_eFnM,S,MTeFnT_sAA_eFnM_S);          
        
        int fe_id = a*ndofn + b;              
        f[fe_id] += 1.0/Jn*fe.detJxW*MTeFnT_sAA_eFnM_S;
      }
    }       
  }
  
  free(u);
  Matrix_cleanup(Fn);
  Matrix_cleanup(Fr);
  Matrix_cleanup(Fnp1);
  Matrix_cleanup(pFn);

  Matrix_cleanup(eFnp1);    
  Matrix_cleanup(pFnp1);
  Matrix_cleanup(pFnp1_I);        
  Matrix_cleanup(L);  
  Matrix_cleanup(S);

  Matrix_cleanup(pFnI);
  Matrix_cleanup(eFn);  
  Matrix_cleanup(M);
  Matrix_cleanup(eFnM);  
  
  Matrix_cleanup(ST_ab);
  Matrix_cleanup(AA);
  Matrix_cleanup(sAA);
  Matrix_cleanup(MTeFnT_sAA);
  Matrix_cleanup(MTeFnT_sAA_eFnM);
    
  FEMLIB_destruct(&fe);
  return err;
}
