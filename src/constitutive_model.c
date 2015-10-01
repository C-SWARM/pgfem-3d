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
#include "plasticity_model_BPA.h"

#include "material.h"
#include "hommat.h"
#include "matgeom.h"
#include "PGFEM_io.h"
#include "data_structure_c.h"
#include "elem3d.h"
#include "femlib.h"
#include "index_macros.h"

#ifndef _Matrix_double
Define_Matrix(double);
#define _Matrix_double 1
#endif

#define DIM 3
#define TENSOR_LEN 9

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
  p->get_Fn = NULL;
  p->get_pF = NULL;
  p->get_pFn = NULL;
  p->get_eF = NULL;
  p->get_eFn = NULL;
  p->get_hardening = NULL;
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
  case TESTING:
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
    err += plasticity_model_BPA_initialize(p);
    break;
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
  p->get_hardening = NULL;
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

int constitutive_model_update_elasticity_npa(const Constitutive_model *m,
                                         const Matrix_double *F,
                                         const double dt,
                                         Matrix_double *L,
                                         Matrix_double *S,
                                         const int compute_stiffness, double alpha)
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
  case CRYSTAL_PLASTICITY:
    if(alpha<0)
      func->get_eF(m,&Fe);
    else
      Matrix_AeqB(Fe,1.0,*F);
    break;
  default:    
    func->get_eF(m,&Fe);
    break;
  }

  Matrix_AxB(C, 1.0, 0.0, Fe, 1, Fe, 0);
  Matrix_inv(C,CI);        
  Matrix_det(Fe, J);

  switch(m->param->type) {
  case TESTING:
  case HYPER_ELASTICITY:
    err += plasticity_model_none_ctx_build(&ctx, F->m_pdata);
    break;
  case CRYSTAL_PLASTICITY:
    err += plasticity_model_ctx_build(&ctx, F->m_pdata, dt,alpha);
    break;
  case BPA_PLASTICITY:
    err += plasticity_model_BPA_ctx_build(&ctx, F->m_pdata, dt);
    break;
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
  Matrix_cleanup(Fe);
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
  double alpha = -1.0;   // if alpha < 0, no inertia
  
  err += constitutive_model_update_elasticity_npa(m,F,dt,L,S,compute_stiffness,alpha);
  
  return err;
  
}
int build_model_parameters_list(Model_parameters **param_list,
                                const int n_mat,
                                const MATGEOM_1 *p_mgeom,
                                const HOMMAT *hmat_list,
                                const int type) /* <-- see #22 */
{
  int err = 0;
  if(type<0)
    return err;  
  if (n_mat <= 0) return 1;
  
  (*param_list) = malloc(n_mat*sizeof(**param_list));
  
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
  if(type<0)
    return err;
  switch (type) {
  case TESTING:
  case HYPER_ELASTICITY: /* nothing required */ break;
  case CRYSTAL_PLASTICITY:
    plasticity_model_read_parameters(eps, ne, elem, n_mat, param_list);
    break;

  case BPA_PLASTICITY:
    {
      for (int a = 0; a < ne; a++) {
        long n_ip = 0;
        int_point(elem[a].toe, &n_ip);
        for (int ip = 0; ip < n_ip; ip++) {
          Constitutive_model *m = &(eps[a].model[ip]);
          err += plasticity_model_BPA_set_initial_values(m);
        }
      }
    }
    break;

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

static int constitutive_model_update_time_steps(EPS *eps,
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

int constitutive_model_update_time_steps_test(const ELEMENT *elem,
                                              NODE *node,
                                              EPS *eps,
                                              const int ne,
                                              const int nn,
                                              const int ndofn,
                                              const double* r,
                                              const double dt,
                                              const int total_Lagrangian)
{
  int nsd = 3;
  int err = 0;
  if (ne <= 0) return 1;

  err += constitutive_model_update_time_steps(eps,ne,elem);
  
  /*********************/
  /* Coordinate update */
  /*********************/
  if(total_Lagrangian) {
    for(int n = 0;n<nn; n++) {
      for(int a=0;a<nsd;a++) {
        int II = node[n].id[a];
        if (II != 0){
          if (a == 0)      node[n].x1 = node[n].x1_fd + r[n*ndofn + a];
          else if (a == 1) node[n].x2 = node[n].x2_fd + r[n*ndofn + a];
          else if (a == 2) node[n].x3 = node[n].x3_fd + r[n*ndofn + a];
        }
      }
    }/* end n < nn */
  } else {
    for(int n = 0;n<nn; n++) {
      for(int a=0;a<nsd;a++) {
        int II = node[n].id[a];
        if (II != 0){
          if (a == 0)      node[n].x1 += r[n*ndofn + a];
          else if (a == 1) node[n].x2 += r[n*ndofn + a];
          else if (a == 2) node[n].x3 += r[n*ndofn + a];
        }
      }
    }/* end n < nn */
  }
  
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
                                  const long *nod,
                                  const NODE *node,
                                  const double dt,
                                  EPS *eps,
                                  const SUPP sup,
                                  const double *r_e)
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
                                  const long *nod,
                                  const NODE *node,
                                  const double dt,
                                  EPS *eps,
                                  const SUPP sup,
                                  const double *r_e)
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
                                    const long *nod,
                                    const NODE *node,
                                    const double dt,
                                    EPS *eps,
                                    const SUPP sup,
                                    const double *r_e,
                                    const int total_Lagrangian)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia

  double *u = malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = malloc(sizeof(*dMdu_all) * 9 * nne * nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  enum {Fn,Fr,Fnp1,pFn,dMdu,eFnp1,pFnp1,pFnp1_I,S,
    eFn,M,eFnM,ST_ab,ST_wg,AA,BB,CC,sAA,sBB,sCC,
    MTeFnT_sAA_eFn,MTeFnT_sAA_eFnM,FrTFr,
    MTeFnT_FrTFreFn,MTeFnT_FrTFreFndMdu,dCdu,
    L_dCdu,MTeFnT_sCC_eFnM,MTeFnT_sAA_eFndMdu,
    sMTeFnT_sAA_eFndMdu,eFnMT,Fend};
  
  Matrix(double) L;  
  Matrix_construct_redim(double,L ,81,1);  

  /* list of second-order tensors */
  Matrix(double) *F2 = malloc(Fend*sizeof(Matrix(double)));
  
  for (int a = 0; a < Fend; a++) {
    Matrix_construct_redim(double, F2[a],3,3);
  }

  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);
  int compute_stiffness = 1;      

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,F2[Fr]);
    
    Constitutive_model *m = &(eps[ii].model[ip-1]);

    /* get a shortened pointer for simplified CM function calls */
    const Model_parameters *func = m->param;
    err += func->get_Fn(m,&F2[Fn]);
    err += func->get_pFn(m,&F2[pFn]);
    err += func->get_eFn(m,&F2[eFn]);
    
    // --> update plasticity part
    if(total_Lagrangian)
    { 
      Matrix(double) FnI;
      Matrix_construct_redim(double, FnI,3,3);
      Matrix_inv(F2[Fn],FnI);
      Matrix_AeqB(F2[Fnp1],1.0,F2[Fr]);  /* set F2[Fnp1] */
      Matrix_AxB(F2[Fr],1.0,0.0,F2[Fnp1],0,FnI,0);  /* recompute F2[Fr] */         
      Matrix_cleanup(FnI);         
    }
    else
    {
      Matrix_AxB(F2[Fnp1],1.0,0.0,F2[Fr],0,F2[Fn],0);  // F2[Fn]+1    
    }   
    Matrix_AxB(F2[FrTFr],1.0,0.0,F2[Fr],1,F2[Fr],0); 

    /* need to have called the integration algorithm. This should be
       done OUTSIDE of the stiffness/tangent functions */
    void *ctx = NULL;
    switch (m->param->type){
    case TESTING:
      err += plasticity_model_none_ctx_build(&ctx,F2[Fnp1].m_pdata);
      break;
    case CRYSTAL_PLASTICITY:
      err += plasticity_model_ctx_build(&ctx,F2[Fnp1].m_pdata,dt,alpha);
      break;
    case BPA_PLASTICITY:
      err += plasticity_model_BPA_ctx_build(&ctx,F2[Fnp1].m_pdata,dt);
      break;      
    default: assert(0); break;
    }
    err += func->compute_dMdu(m, ctx, fe.ST, nne, ndofn, dMdu_all);
    err += func->destroy_ctx(&ctx);
    err += func->get_pF(m,&F2[pFnp1]);
    err += func->get_eF(m,&F2[eFnp1]);

    Matrix_inv(F2[pFnp1], F2[pFnp1_I]);
    Matrix_AxB(F2[M],1.0,0.0,F2[pFn],0,F2[pFnp1_I],0);
    // <-- update plasticity part 

    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(F2[S],0.0);    
    
    constitutive_model_update_elasticity(m,&F2[eFnp1],dt,&L,&F2[S],compute_stiffness);
    // <-- update elasticity part

    // --> start computing tagent
    Matrix_AxB(F2[eFnM],1.0,0.0,F2[eFn],0,F2[M],0);
    Matrix_AeqBT(F2[eFnMT],1.0,F2[eFnM]);
        
    double Jn; Matrix_det(F2[Fn], Jn);

    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];
        Matrix_init_w_array(F2[ST_ab],3,3,ptrST_ab);
        Matrix_AxB(F2[AA],1.0,0.0,F2[Fr],1,F2[ST_ab],0); 
        Matrix_symmetric(F2[AA],F2[sAA]);
        
        Matrix_Tns2_AxBxC(F2[MTeFnT_sAA_eFn],1.0,0.0,F2[eFnMT],F2[sAA],F2[eFn]);        
        Matrix_AxB(F2[MTeFnT_sAA_eFnM],1.0,0.0,F2[MTeFnT_sAA_eFn],0,F2[M],0);  

        for(int w=0; w<nne; w++)
        {
          for(int g=0; g<nsd; g++)
          { 
            const double* const ptrST_wg = &(fe.ST)[idx_4_gen(w,g,0,0,
                                                      nne,nsd,nsd,nsd)]; 
            const double* const ptr_dMdu_wg = &(dMdu_all[idx_4_gen(w,g,0,0,
                                                                   nne,nsd,nsd,nsd)]);

            Matrix_init_w_array(F2[ST_wg],3,3,ptrST_wg); 
            Matrix_init_w_array(F2[dMdu],3,3,ptr_dMdu_wg);

            Matrix_AxB(F2[BB],1.0,0.0,F2[Fr],1,F2[ST_wg],0); 
            Matrix_symmetric(F2[BB],F2[sBB]);
            Matrix_AxB(F2[CC], 1.0,0.0,F2[ST_ab],1,F2[ST_wg],0);
            Matrix_symmetric(F2[CC],F2[sCC]);
            
            // compute F2[dCdu]
            Matrix_Tns2_AxBxC(F2[MTeFnT_FrTFreFn],1.0,0.0,F2[eFnMT],F2[FrTFr],F2[eFn]);            
            Matrix_AxB(F2[MTeFnT_FrTFreFndMdu],1.0,0.0,F2[MTeFnT_FrTFreFn],0,F2[dMdu],0);                                    
            Matrix_symmetric(F2[MTeFnT_FrTFreFndMdu],F2[dCdu]);
                        
            Matrix_Tns2_AxBxC(F2[dCdu],1.0,1.0,F2[eFnMT],F2[sBB],F2[eFnM]);            
            
            // compute F2[MTeFnT_sAA_eFnM]:L:F2[dCdu]
            Matrix_Tns4_dd_Tns2(F2[L_dCdu],L,F2[dCdu]);
            double MTeFnT_sAA_eFnM_L_dCdu = 0.0;
            Matrix_ddot(F2[MTeFnT_sAA_eFnM],F2[L_dCdu],MTeFnT_sAA_eFnM_L_dCdu);
            
            // compute F2[MTeFnT_sCC_eFnM]
            Matrix_Tns2_AxBxC(F2[MTeFnT_sCC_eFnM],1.0,0.0,F2[eFnMT],F2[sCC],F2[eFnM]);
                        
            // compute F2[MTeFnT_sCC_eFnM]:F2[S]
            double MTeFnT_sCC_eFnM_S = 0.0;
            Matrix_ddot(F2[MTeFnT_sCC_eFnM],F2[S],MTeFnT_sCC_eFnM_S);
            
            // compute F2[MTeFnT_sAA_eFndMdu]
            Matrix_AxB(F2[MTeFnT_sAA_eFndMdu],1.0,0.0,F2[MTeFnT_sAA_eFn],0,F2[dMdu],0);    
            Matrix_symmetric(F2[MTeFnT_sAA_eFndMdu], F2[sMTeFnT_sAA_eFndMdu]);        

            // compute F2[MTeFnT_sAA_eFndMdu]:F2[S]
            double sMTeFnT_sAA_eFndMdu_S = 0.0;            
            Matrix_ddot(F2[sMTeFnT_sAA_eFndMdu],F2[S],sMTeFnT_sAA_eFndMdu_S);
            
            const int lk_idx = idx_K(a,b,w,g,nne,nsd);  
                      
            lk[lk_idx] += 1.0/Jn*fe.detJxW*(MTeFnT_sAA_eFnM_L_dCdu + 2.0*sMTeFnT_sAA_eFndMdu_S + MTeFnT_sCC_eFnM_S);            
          }
        }
      }
    }
  }
  free(u);
  
  Matrix_cleanup(L); 

  /* destroy second-order tensors */
  for (int a = 0; a < Fend; a++) {
    Matrix_cleanup(F2[a]);
  }
  free(F2);

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
                                    const long *nod,
                                    const NODE *node,
                                    const double dt,
                                    EPS *eps,
                                    const SUPP sup,
                                    const double *r_e,
                                    const int total_Lagrangian)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia
    
  double *u = (double *) malloc(sizeof(double)*nne*nsd);
  
  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  enum {Fn,Fr,Fnp1,pFn,eFnp1,pFnp1,pFnp1_I,
    L,S,pFnI,eFn,M,eFnM,ST_ab,AA,sAA,MTeFnT_sAA_eFnM,eFnMT,Fend}; 
  
  /* second-order tensors */
  Matrix(double) *F2 = malloc(Fend*sizeof(Matrix(double)));
  for (int a = 0; a < Fend; a++) {
    Matrix_construct_redim(double, F2[a],3,3);
  }

  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);      
  int compute_stiffness = 0;
   
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,F2[Fr]);

    Constitutive_model *m = &(eps[ii].model[ip-1]);

    m->param->get_Fn(m,&F2[Fn]);
    m->param->get_pFn(m,&F2[pFn]);    

    Matrix_inv(F2[pFn], F2[pFnI]);
    m->param->get_eFn(m,&F2[eFn]);
   
    // --> update plasticity part
    if(total_Lagrangian)
    {
      Matrix(double) FnI;
      Matrix_construct_redim(double, FnI,3,3);
      Matrix_inv(F2[Fn],FnI);
      Matrix_AeqB(F2[Fnp1],1.0,F2[Fr]);  /* set F2[Fnp1] */
      Matrix_AxB(F2[Fr],1.0,0.0,F2[Fnp1],0,FnI,0);  /* recompute F2[Fr] */         
      Matrix_cleanup(FnI);          
    }
    else
    {
      Matrix_AxB(F2[Fnp1],1.0,0.0,F2[Fr],0,F2[Fn],0);  /* compute F2[Fnp1] */    
    }      

    void *ctx = NULL;
    switch (m->param->type){
    case TESTING:
      err += plasticity_model_none_ctx_build(&ctx,F2[Fnp1].m_pdata);
      break;
    case CRYSTAL_PLASTICITY:
      err += plasticity_model_ctx_build(&ctx,F2[Fnp1].m_pdata,dt,alpha);
      break;
    case BPA_PLASTICITY:
      err += plasticity_model_BPA_ctx_build(&ctx,F2[Fnp1].m_pdata,dt);
      break;      
    default: assert(0); break;
    }
    err += m->param->integration_algorithm(m,ctx);
    err += m->param->destroy_ctx(&ctx);
    err += m->param->get_pF(m,&F2[pFnp1]);

    Matrix_AxB(F2[M],1.0,0.0,F2[pFnI],0,F2[pFnp1],0);    
    // <-- update plasticity part

    // --> update elasticity part
    Matrix_init(F2[S],0.0);    
    
    Matrix_inv(F2[pFnp1], F2[pFnp1_I]);
    Matrix_AxB(F2[eFnp1],1.0,0.0,F2[Fnp1],0,F2[pFnp1_I],0);
    constitutive_model_update_elasticity(m,&F2[eFnp1],dt,NULL,&F2[S],compute_stiffness);
    // <-- update elasticity part
            
    Matrix_AxB(F2[eFnM],1.0,0.0,F2[eFn],0,F2[M],0);
    Matrix_AeqBT(F2[eFnMT],1.0,F2[eFnM]);
    double Jn; Matrix_det(F2[Fn], Jn);
    
    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];
        Matrix_init_w_array(F2[ST_ab],3,3,ptrST_ab);
        Matrix_AxB(F2[AA],1.0,0.0,F2[Fr],1,F2[ST_ab],0); 
        Matrix_symmetric(F2[AA],F2[sAA]);

        Matrix_Tns2_AxBxC(F2[MTeFnT_sAA_eFnM],1.0,0.0,F2[eFnMT],F2[sAA],F2[eFnM]);
        double MTeFnT_sAA_eFnM_S = 0.0; 
        Matrix_ddot(F2[MTeFnT_sAA_eFnM],F2[S],MTeFnT_sAA_eFnM_S);          
        
        int fe_id = a*ndofn + b;              
        f[fe_id] += 1.0/Jn*fe.detJxW*MTeFnT_sAA_eFnM_S;
      }
    }       
  }
  
  free(u);
  
  /* destroyu list of second-order tenosors */
  for(int a = 0; a < Fend; a++) {
    Matrix_cleanup(F2[a]);
  }
  free(F2);

  FEMLIB_destruct(&fe);
  return err;
}

int constitutive_model_update_output_variables(SIG *sig,
                                               EPS *eps,
                                               const int ne,
                                               const double dt)
{
  int err = 0;

  static const double eye[TENSOR_LEN] = {[0] = 1.0, [4] = 1.0, [8] = 1.0};
  /* *** ASSUME LINEAR ELEMENTS -- 1 INTEGRATION POINT *** */
  const int ip = 0;

  /* deformation gradient */
  Matrix_double F, eF, pF, S;
  Matrix_construct_redim(double, F, DIM, DIM);
  Matrix_construct_redim(double, eF, DIM, DIM);
  Matrix_construct_redim(double, pF, DIM, DIM);
  Matrix_construct_redim(double, S, DIM, DIM);

  for (int i = 0; i < ne; i++) {
    /* *** ASSUME LINEAR ELEMENTS -- 1 INTEGRATION POINT *** */
    const Constitutive_model *m = eps[i].model;
    const Model_parameters *func = m->param;

    err += func->get_Fn(m, &F);
    err += func->get_eFn(m, &eF);
    err += func->get_pFn(m, &pF);
    err += constitutive_model_update_elasticity(m,&F,dt,NULL,&S,0);

    /* get aliases to Matrix data for simpler access */
    const double *Sd = S.m_pdata;
    const double *eFd = eF.m_pdata;
    const double *pFd = pF.m_pdata;
    const double eJ = det3x3(eFd);

    /* store symmetric part of S (PK2) */
    sig[i].il[ip].o[0] = Sd[idx_2(0,0)]; /* XX */
    sig[i].il[ip].o[1] = Sd[idx_2(1,1)]; /* YY */
    sig[i].il[ip].o[2] = Sd[idx_2(2,2)]; /* ZZ */
    sig[i].il[ip].o[3] = Sd[idx_2(1,2)]; /* YZ */
    sig[i].il[ip].o[4] = Sd[idx_2(0,2)]; /* XZ */
    sig[i].il[ip].o[5] = Sd[idx_2(0,1)]; /* XY */

    /* store elastic deformation */
    memcpy(eps[i].il[ip].F, eFd, TENSOR_LEN * sizeof(*eFd));

    /* store the hardening parameter */
    err += func->get_hardening(m, &eps[i].dam[ip].wn);

    /* compute/store the plastic stretch */
    eps[i].dam[ip].Xn = sqrt( cblas_ddot(TENSOR_LEN, pFd, 1, pFd, 1) / 3.0 );

    /* Compute the Cauchy Stress sigma = 1/eJ eF S eF' */
    double sigma[TENSOR_LEN] = {};
    double temp[TENSOR_LEN] = {};
    double temp_I[TENSOR_LEN] = {};
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                DIM,DIM,DIM, 1.0 / eJ, eFd, DIM, Sd, DIM,
                0.0, temp,DIM);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                DIM, DIM, DIM, 1.0, temp, DIM, eFd, DIM,
                0.0, sigma, DIM);

    /* store symmetric part */
    sig[i].el.o[0] = sigma[idx_2(0,0)]; /* XX */
    sig[i].el.o[1] = sigma[idx_2(1,1)]; /* YY */
    sig[i].el.o[2] = sigma[idx_2(2,2)]; /* ZZ */
    sig[i].el.o[3] = sigma[idx_2(1,2)]; /* YZ */
    sig[i].el.o[4] = sigma[idx_2(0,2)]; /* XZ */
    sig[i].el.o[5] = sigma[idx_2(0,1)]; /* XY */

    /* Compute the logarithmic strain e = 1/2(I - inv(FF'))*/
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                DIM, DIM, DIM, 1.0, eFd, DIM, eFd, DIM,
                0.0, temp, DIM);
    err += inv3x3(temp, temp_I);
    /* e <-- temp is the Euler strain */
    for(int i = 0; i < TENSOR_LEN; i++){
      temp[i] = 0.5*(eye[i]-temp_I[i]);
    }

    /* store symmetric part (also Eng. strain) */
    eps[i].el.o[0] = temp[idx_2(0,0)];
    eps[i].el.o[1] = temp[idx_2(1,1)];
    eps[i].el.o[2] = temp[idx_2(2,2)];

    eps[i].el.o[3] = 2. * temp[idx_2(1,2)];
    eps[i].el.o[4] = 2. * temp[idx_2(0,2)];
    eps[i].el.o[5] = 2. * temp[idx_2(0,1)];

  }

  Matrix_cleanup(F);
  Matrix_cleanup(eF);
  Matrix_cleanup(pF);
  Matrix_cleanup(S);

  /* Compute equivalent stress and strain */
  Mises (ne, sig, eps, 0);

  return err;
}



int stiffness_el_crystal_plasticity_w_inertia(double *lk,
                                    const int ii,
                                    const int ndofn,
                                    const int nne,
                                    const int nsd,
                                    const ELEMENT *elem,
                                    const long *nod,
                                    const NODE *node,
                                    const double dt,
                                    EPS *eps,
                                    const SUPP sup,
                                    const double *r_e,
                                    double alpha)
{
  int err = 0;
  int total_Lagrangian = 1;
  double *u = malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = malloc(sizeof(*dMdu_all) * 9 * nne * nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  enum {Fn,Fr,Fnp1,pFn,dMdu,eFnp1,pFnp1,pFnp1_I,S,
    eFn,M,Mn,Mnp1,eFnM,ST_ab,ST_wg,AA,BB,CC,sAA,sBB,sCC,
    MTeFnT_sAA_eFn,MTeFnT_sAA_eFnM,FrTFr,
    MTeFnT_FrTFreFn,MTeFnT_FrTFreFndMdu,dCdu,
    L_dCdu,MTeFnT_sCC_eFnM,MTeFnT_sAA_eFndMdu,
    sMTeFnT_sAA_eFndMdu,eFnMT,Fend};
  
  Matrix(double) L;  
  Matrix_construct_redim(double,L ,81,1);  

  /* list of second-order tensors */
  Matrix(double) *F2 = malloc(Fend*sizeof(Matrix(double)));
  
  for (int a = 0; a < Fend; a++) {
    Matrix_construct_redim(double, F2[a],3,3);
  }

  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);
  int compute_stiffness = 1;      

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,F2[Fnp1]);
    
    Constitutive_model *m = &(eps[ii].model[ip-1]);

    err += m->param->get_pF(m,&F2[pFnp1]);
    err += m->param->get_pFn(m,&F2[pFn]);
    err += m->param->get_Fn(m,&F2[Fn]);
        
    Matrix_inv(F2[pFnp1], F2[Mnp1]);
    Matrix_inv(F2[pFn],   F2[Mn]);
    
    mid_point_rule(F2[M].m_pdata, F2[Mn].m_pdata, F2[Mnp1].m_pdata, alpha, nsd*nsd);
    mid_point_rule(F2[Fr].m_pdata, F2[Fn].m_pdata, F2[Fnp1].m_pdata, alpha, nsd*nsd);    
    
    Matrix_AxB(F2[FrTFr],1.0,0.0,F2[Fr],1,F2[Fr],0);  

    /* need to have called the integration algorithm. This should be
       done OUTSIDE of the stiffness/tangent functions */
    void *ctx = NULL;
    switch (m->param->type){
    case TESTING:
      err += plasticity_model_none_ctx_build(&ctx,F2[Fnp1].m_pdata);
      break;
    case CRYSTAL_PLASTICITY:
      err += plasticity_model_ctx_build(&ctx,F2[Fnp1].m_pdata,dt,alpha);
      break;
    case BPA_PLASTICITY:
      err += plasticity_model_BPA_ctx_build(&ctx,F2[Fnp1].m_pdata,dt);
      break;      
    default: assert(0); break;
    }
        
    err += m->param->compute_dMdu(m, ctx, fe.ST, nne, ndofn, dMdu_all);
    err += m->param->destroy_ctx(&ctx);

    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(F2[S],0.0);    
    
    constitutive_model_update_elasticity_npa(m,&F2[Fr],dt,&L,&F2[S],compute_stiffness, alpha);
    // <-- update elasticity part

    // --> start computing tagent
    //Matrix_AxB(F2[eFnM],1.0,0.0,F2[eFn],0,F2[M],0);
    Matrix_AeqBT(F2[eFnM],1.0,F2[M]); // eFn = 1 for total_Lagrangian
    Matrix_AeqBT(F2[eFnMT],1.0,F2[eFnM]);
        
    double Jn; // Matrix_det(F2[Fn], Jn);
    Jn = 1.0;

    for(int a=0; a<nne; a++)
    {
      for(int b=0; b<nsd; b++)
      {
        const double* const ptrST_ab = &(fe.ST)[idx_4_gen(a,b,0,0,
                                                nne,nsd,nsd,nsd)];
        Matrix_init_w_array(F2[ST_ab],3,3,ptrST_ab);
        Matrix_AxB(F2[AA],1.0,0.0,F2[Fr],1,F2[ST_ab],0);
        Matrix_symmetric(F2[AA],F2[sAA]);

//       Matrix_Tns2_AxBxC(F2[MTeFnT_sAA_eFn],1.0,0.0,F2[eFnMT],F2[sAA],F2[eFn]);
//       Matrix_AxB(F2[MTeFnT_sAA_eFnM],1.0,0.0,F2[MTeFnT_sAA_eFn],0,F2[M],0);
        Matrix_Tns2_AxBxC(F2[MTeFnT_sAA_eFnM],1.0,0.0,F2[eFnMT],F2[sAA],F2[M]);

        for(int w=0; w<nne; w++)
        {
          for(int g=0; g<nsd; g++)
          {
            const double* const ptrST_wg = &(fe.ST)[idx_4_gen(w,g,0,0,
                                                      nne,nsd,nsd,nsd)];
            const double* const ptr_dMdu_wg = &(dMdu_all[idx_4_gen(w,g,0,0,
                                                                   nne,nsd,nsd,nsd)]);

            Matrix_init_w_array(F2[ST_wg],3,3,ptrST_wg);
            Matrix_init_w_array(F2[dMdu],3,3,ptr_dMdu_wg);

            Matrix_AxB(F2[BB],1.0,0.0,F2[Fr],1,F2[ST_wg],0);
            Matrix_symmetric(F2[BB],F2[sBB]);
            Matrix_AxB(F2[CC], 1.0,0.0,F2[ST_ab],1,F2[ST_wg],0);
            Matrix_symmetric(F2[CC],F2[sCC]);

            // compute F2[dCdu]
//            Matrix_Tns2_AxBxC(F2[MTeFnT_FrTFreFn],1.0,0.0,F2[eFnMT],F2[FrTFr],F2[eFn]);
            Matrix_Tns2_AxBxC(F2[MTeFnT_FrTFreFndMdu],1.0,0.0,F2[eFnMT],F2[FrTFr],F2[dMdu]);
//            Matrix_AxB(F2[MTeFnT_FrTFreFndMdu],1.0,0.0,F2[MTeFnT_FrTFreFn],0,F2[dMdu],0);
            Matrix_symmetric(F2[MTeFnT_FrTFreFndMdu],F2[dCdu]);
                        
            Matrix_Tns2_AxBxC(F2[dCdu],1.0,1.0,F2[eFnMT],F2[sBB],F2[eFnM]);
            
            // compute F2[MTeFnT_sAA_eFnM]:L:F2[dCdu]
            Matrix_Tns4_dd_Tns2(F2[L_dCdu],L,F2[dCdu]);
            double MTeFnT_sAA_eFnM_L_dCdu = 0.0;
            Matrix_ddot(F2[MTeFnT_sAA_eFnM],F2[L_dCdu],MTeFnT_sAA_eFnM_L_dCdu);
            
            // compute F2[MTeFnT_sCC_eFnM]
            Matrix_Tns2_AxBxC(F2[MTeFnT_sCC_eFnM],1.0,0.0,F2[eFnMT],F2[sCC],F2[eFnM]);
                        
            // compute F2[MTeFnT_sCC_eFnM]:F2[S]
            double MTeFnT_sCC_eFnM_S = 0.0;
            Matrix_ddot(F2[MTeFnT_sCC_eFnM],F2[S],MTeFnT_sCC_eFnM_S);
            
            // compute F2[MTeFnT_sAA_eFndMdu]
            Matrix_AxB(F2[MTeFnT_sAA_eFndMdu],1.0,0.0,F2[MTeFnT_sAA_eFn],0,F2[dMdu],0);
            Matrix_symmetric(F2[MTeFnT_sAA_eFndMdu], F2[sMTeFnT_sAA_eFndMdu]);

            // compute F2[MTeFnT_sAA_eFndMdu]:F2[S]
            double sMTeFnT_sAA_eFndMdu_S = 0.0;
            Matrix_ddot(F2[sMTeFnT_sAA_eFndMdu],F2[S],sMTeFnT_sAA_eFndMdu_S);
            
            const int lk_idx = idx_K(a,b,w,g,nne,nsd);

            lk[lk_idx] += 1.0/Jn*fe.detJxW*(MTeFnT_sAA_eFnM_L_dCdu + 2.0*sMTeFnT_sAA_eFndMdu_S + MTeFnT_sCC_eFnM_S);
          }
        }
      }
    }
  }
  free(u);
  
  Matrix_cleanup(L);

  /* destroy second-order tensors */
  for (int a = 0; a < Fend; a++) {
    Matrix_cleanup(F2[a]);
  }
  free(F2);

  FEMLIB_destruct(&fe);
  free(dMdu_all);
 
  return err;
}

int residuals_el_crystal_plasticity_n_plus_alpha(double *f,
                                    const Constitutive_model *m,
                                    const int ii,
                                    const int ndofn,                                    
                                    const int nne,
                                    const int nsd,
                                    const ELEMENT *elem,
                                    const long *nod,
                                    const NODE *node,
                                    const double dt,
                                    Matrix(double) *Mnp1,
                                    Matrix(double) *Mn,
                                    Matrix(double) *Fnp1,
                                    Matrix(double) *Fn,
                                    double alpha,
                                    double dt_alpha_1_minus_alpha,
                                    FEMLIB *fe)
{
  // Total Lagrangian based
  int err = 0;   
  enum {M,F,S,ST_ab,AA,sAA,MT,MT_sAA_M,Fend}; 
  
  // second-order tensors
  Matrix(double) *F2 = malloc(Fend*sizeof(Matrix(double)));
  for (int a = 0; a < Fend; a++) {
    Matrix_construct_redim(double, F2[a],nsd,nsd);
  }

  int compute_stiffness = 0;
   
  mid_point_rule(F2[M].m_pdata, Mn->m_pdata, Mnp1->m_pdata, alpha, nsd*nsd);
  mid_point_rule(F2[F].m_pdata, Fn->m_pdata, Fnp1->m_pdata, alpha, nsd*nsd);
    
  Matrix_AeqBT(F2[MT],1.0,F2[M]);  
  Matrix_init(F2[S],0.0);    
  constitutive_model_update_elasticity_npa(m,&F2[F],dt,NULL,&F2[S],compute_stiffness, alpha);
  
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &(fe->ST)[idx_4_gen(a,b,0,0,
                                              nne,nsd,nsd,nsd)];
      Matrix_init_w_array(F2[ST_ab],3,3,ptrST_ab);
      Matrix_AxB(F2[AA],1.0,0.0,F2[F],1,F2[ST_ab],0); 
      Matrix_symmetric(F2[AA],F2[sAA]);

      Matrix_Tns2_AxBxC(F2[MT_sAA_M],1.0,0.0,F2[MT],F2[sAA],F2[M]);
      double MT_sAA_M_S = 0.0; 
      Matrix_ddot(F2[MT_sAA_M],F2[S],MT_sAA_M_S);          
      
      int fe_id = a*ndofn + b;              
      f[fe_id] += dt_alpha_1_minus_alpha*fe->detJxW*MT_sAA_M_S;
    }
  }         
  
  // destroyu list of second-order tenosors
  for(int a = 0; a < Fend; a++) {
    Matrix_cleanup(F2[a]);
  }
  free(F2);
  return err;
}

int residuals_el_crystal_plasticity_w_inertia(double *f,
                                    const int ii,
                                    const int ndofn,
                                    const int nne,
                                    const int nsd,
                                    const ELEMENT *elem,
                                    const long *nod,
                                    const NODE *node,
                                    const double dt,
                                    EPS *eps,
                                    const SUPP sup,
                                    const double *r_e,
                                    const double alpha)
{
  int err = 0;
  int total_Lagrangian = 1;
    
  double *u       = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_npa   = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_nm1pa = (double *) malloc(sizeof(double)*nne*nsd);    
    
  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  enum {Fnp1,Fn,Fnm1,pFnp1,pFn,pFnm1,Mnp1,Mn,Mnm1,Fend}; 
  
  // second-order tensors
  Matrix(double) *F2 = malloc(Fend*sizeof(Matrix(double)));
  for (int a = 0; a < Fend; a++) {
    Matrix_construct_redim(double, F2[a],3,3);
  }

  FEMLIB fe;
  FEMLIB_initialization_by_elem(&fe, ii, elem, node, 0,total_Lagrangian);      
  int compute_stiffness = 0;
  
  memset(f_npa, 0, sizeof(double)*nne*ndofn);   
  memset(f_nm1pa, 0, sizeof(double)*nne*ndofn);   
    
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);
    FEMLIB_update_shape_tensor(&fe);  
    FEMLIB_update_deformation_gradient(&fe,ndofn,u,F2[Fnp1]);

    Constitutive_model *m = &(eps[ii].model[ip-1]);

    void *ctx = NULL;
    switch (m->param->type){
    case TESTING:
      err += plasticity_model_none_ctx_build(&ctx,F2[Fnp1].m_pdata);
      break;
    case CRYSTAL_PLASTICITY:
      err += plasticity_model_ctx_build(&ctx,F2[Fnp1].m_pdata,dt,alpha);
      break;
    case BPA_PLASTICITY:
      err += plasticity_model_BPA_ctx_build(&ctx,F2[Fnp1].m_pdata,dt);
      break;      
    default: assert(0); break;
    }
    err += m->param->integration_algorithm(m,ctx);
    err += m->param->destroy_ctx(&ctx);
    err += m->param->get_pF(m,&F2[pFnp1]);
    err += m->param->get_pFn(m,&F2[pFn]);
    err += m->param->get_pFnm1(m,&F2[pFnm1]); 

    err += m->param->get_Fn(m,  &F2[Fn]);
    err += m->param->get_Fnm1(m,&F2[Fnm1]);
            
    Matrix_inv(F2[pFnp1], F2[Mnp1]);
    Matrix_inv(F2[pFn],   F2[Mn]);
    Matrix_inv(F2[pFnm1], F2[Mnm1]);
    
    double dt_1_minus_alpha = -dt*(1.0-alpha);
    err += residuals_el_crystal_plasticity_n_plus_alpha(f_npa,m,ii,ndofn,nne,nsd,
                                elem,nod,node,dt,
                                &F2[Mnp1],&F2[Mn],&F2[Fnp1],&F2[Fn],
                                alpha, dt_1_minus_alpha,&fe);
                                
    double dt_alpha = -dt*alpha;

    err += residuals_el_crystal_plasticity_n_plus_alpha(f_nm1pa,m,ii,ndofn,nne,nsd,
                                elem,nod,node,dt,
                                &F2[Mn],&F2[Mnm1],&F2[Fn],&F2[Fnm1],
                                alpha, dt_alpha,&fe);
  }
  
  for(int a=0; a<nne*nsd; a++)
    f[a] += f_npa[a] + f_nm1pa[a];
  
  
  free(u);
  free(f_npa);
  free(f_nm1pa);
  
  // destroyu list of second-order tenosors
  for(int a = 0; a < Fend; a++) {
    Matrix_cleanup(F2[a]);
  }
  free(F2);

  FEMLIB_destruct(&fe);
  return err;
}