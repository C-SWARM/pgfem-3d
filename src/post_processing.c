#include "post_processing.h"
#include "enumerations.h"
#include "utils.h"
#include "allocation.h"
#include "PGFem3D_to_VTK.hpp"

#include "constitutive_model.h"
#include "plasticity_model.h"

int read_from_VTK(const PGFem3D_opt *opts, int myrank, int step, double *u)
{
  int err = 0;
  char filename[1024];
  sprintf(filename,"%s/VTK/STEP_%.5d/%s_%d_%d.vtu",opts->opath,step,opts->ofname,myrank, step);   
  err += read_VTK_file(filename, u);      
  return err;
}

void post_processing_compute_stress_disp_ip(FEMLIB *fe, int e, Matrix(double) S, HOMMAT *hommat, ELEMENT *elem, 
                          Matrix(double) F, double Pn)
{
  Matrix(double) C,CI,devS;
  Matrix_construct_init(double,C,3,3,0.0);
  Matrix_construct_init(double,CI,3,3,0.0);
  Matrix_construct_init(double,devS,3,3,0.0);

  Matrix_AxB(C,1.0,0.0,F,1,F,0);
  double J;
  Matrix_det(F, J);
  Matrix_inv(C, CI);
  
  int mat = elem[e].mat[2];
  double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu)); 

  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
  Stress(C.m_pdata,&hommat[mat],devS.m_pdata);

  dUdJFuncPtr DUDJ = getDUdJFunc(1,&hommat[mat]);
  double dUdJ = 0.0;
  DUDJ(J,&hommat[mat],&dUdJ);
  double kappaJdUdJ = kappa*J*dUdJ;
  
//  for(int a = 0; a<9; a++)
//    S.m_pdata[a] = devS.m_pdata[a] + kappaJdUdJ*CI.m_pdata[a];
  Matrix_AplusB(S,1.0,devS,kappaJdUdJ,CI);
            
  Matrix_cleanup(C);
  Matrix_cleanup(CI);
  Matrix_cleanup(devS);        
}

void post_processing_compute_stress_3f_ip(FEMLIB *fe, int e, Matrix(double) S, HOMMAT *hommat, ELEMENT *elem, 
                          Matrix(double) F, double Pn)
{
  Matrix(double) C,CI,devS;
  Matrix_construct_init(double,C,3,3,0.0);
  Matrix_construct_init(double,CI,3,3,0.0);
  Matrix_construct_init(double,devS,3,3,0.0);
  
/*  for(int a=0; a<9; a++)
  {
    if(F.m_pdata[a]<1.0e-6)
      F.m_pdata[a] = 0.0;
  }*/    
  
  Matrix_AxB(C,1.0,0.0,F,1,F,0);
  double J;
  Matrix_det(F, J);
  Matrix_inv(C, CI);
  int mat = elem[e].mat[2];
  double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));      			      
  
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
  dUdJFuncPtr UP = getDUdJFunc(1, &hommat[mat]);
  Stress(C.m_pdata,&hommat[mat],devS.m_pdata);
  
  double Up = 0.0;
  UP(J,&hommat[mat],&Up);
//  printf("%e, %e\n", Pn, kappa*Up);
  double JPn = J*Pn;
  
//  for(int a = 0; a<9; a++)
//    S.m_pdata[a] = devS.m_pdata[a] + JPn*CI.m_pdata[a]; 
  Matrix_AplusB(S,1.0,devS,JPn,CI);
          
  Matrix_cleanup(C);
  Matrix_cleanup(CI);
  Matrix_cleanup(devS);        
}

void post_processing_compute_stress_plasticity_ip(FEMLIB *fe, int e, int ip, Matrix(double) *S, HOMMAT *hommat, ELEMENT *elem, EPS *eps,
                          Matrix(double) *Fr, double Pn)
{
  int updated_Lagrangian = 0;
  int compute_stiffness = 0;
  
  Constitutive_model *m = &(eps[e].model[ip-1]);
  Matrix(double) *Fs = (m->vars).Fs;
  Matrix(double) Fn, pFnp1_I, Fnp1, eFnp1;
  Matrix_construct_init(double,Fn,3,3,0.0);
  Matrix_construct_init(double,Fnp1,3,3,0.0);  
  Matrix_construct_init(double,pFnp1_I,3,3,0.0);
  Matrix_construct_init(double,eFnp1,3,3,0.0);  
  if(updated_Lagrangian)
    Matrix_AeqB(Fn,1.0,Fs[TENSOR_Fn]);
  else
    Matrix_eye(Fn,3);
  Matrix_AxB(Fnp1,1.0,0.0,*Fr,0,Fn,0);
  Matrix_inv(Fs[TENSOR_pFnp1], pFnp1_I);  
  Matrix_AxB(eFnp1,1.0,0.0,Fnp1,0,pFnp1_I,0); 
  constitutive_model_update_elasticity(m,&eFnp1,0.0,NULL,S,compute_stiffness);
  double Jn; Matrix_det(Fn, Jn); 
  
  Matrix_cleanup(Fn);
  Matrix_cleanup(Fnp1);  
  Matrix_cleanup(pFnp1_I);
  Matrix_cleanup(eFnp1);   
}                          
void post_processing_compute_stress(double *GS, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, int analysis)
{
  int nsd = 3;
  Matrix(double) F,S,LS;
  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,S,3,3,0.0);
  Matrix_construct_init(double,LS,3,3,0.0);
  double LV = 0.0;
  double GV = 0.0;
  
  for(int e = 0; e<ne; e++)
  {    
    int intg_order = 1;    
    if (PGFEM3D_DEV_TEST)
      intg_order = 0;

    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe, e, elem, node, intg_order,1);
    int nne = fe.nne;
    
    Matrix(double) Np, u, P;  
    Matrix_construct_init(double,Np,npres,1,0.0);
    Matrix_construct_init(double,u,nne*nsd,1,0.0);
    Matrix_construct_init(double,P,npres,1,0.0);        
                        
    for(int a = 0; a<nne; a++)
    {      
      int nid = Vec_v(fe.node_id, a+1);
      for(int b=0; b<nsd; b++)
      {
        Vec_v(u, a*nsd+b+1) = r[nid*ndofn + b];
      }
    }
    
    if(analysis==TF)
    {     	
      if(npres==1)
        Vec_v(P, 1) = eps[e].d_T[0];
      else
      {
        for(int a = 0; a<nne; a++)
        {
          int nid = Vec_v(fe.node_id, a+1);
          Vec_v(P, a+1) = r[nid*ndofn + 3];
        }
      }          
    }
    
    for(int ip = 1; ip<=fe.nint; ip++)
    {
      FEMLIB_elem_basis_V(&fe, ip);  
      FEMLIB_update_shape_tensor(&fe);
      FEMLIB_update_deformation_gradient(&fe,nsd,u.m_pdata,F);
      double Pn = 0.0;          
      if(analysis==TF)
      {  
        FEMLIB_elem_shape_function(&fe,ip,npres, Np);

        for(int a=1; a<=npres; a++)
          Pn += Vec_v(Np,a)*Vec_v(P,a);

        post_processing_compute_stress_3f_ip(&fe,e,S,hommat,elem,F,Pn);
        
        double J;
        Matrix_det(F, J);
//        printf("%e, %e\n", J, eps[e].T[0]);
      }
      else
      {
        if (PGFEM3D_DEV_TEST)
          post_processing_compute_stress_plasticity_ip(&fe,e,ip,&S,hommat,elem,eps,&F,Pn);
        else
          post_processing_compute_stress_disp_ip(&fe,e,S,hommat,elem,F,Pn);
          
      }  
      
      LV += fe.detJxW;

      for(int a=0; a<9; a++)
        LS.m_pdata[a] += S.m_pdata[a]*fe.detJxW;
	  }
    FEMLIB_destruct(&fe);
    Matrix_cleanup(Np);
    Matrix_cleanup(u);
    Matrix_cleanup(P);    

  }
  
  Matrix_cleanup(F);
  Matrix_cleanup(S);
    
  MPI_Allreduce(LS.m_pdata,GS,9,MPI_DOUBLE,MPI_SUM,mpi_comm);
  MPI_Allreduce(&LV,&GV,1,MPI_DOUBLE,MPI_SUM,mpi_comm);  
  Matrix_cleanup(LS);
  
  printf("Total V = %e\n", GV);
  for(int a=0; a<9; a++)
    GS[a] = GS[a]/GV;
};

