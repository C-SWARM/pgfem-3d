#include "post_processing.h"
#include "enumerations.h"
#include "utils.h"
#include "allocation.h"
#include "PGFem3D_to_VTK.hpp"

#include "constitutive_model.h"

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

void post_processing_compute_stress4CM(FEMLIB *fe, int e, int ip, Matrix(double) *S, double *Jnp1, 
                                                  HOMMAT *hommat, ELEMENT *elem, EPS *eps)
{
  int total_Lagrangian = 0;
  int compute_stiffness = 0;
  
  Constitutive_model *m = &(eps[e].model[ip-1]);
  Matrix(double) Fnp1, eFnp1;
  Matrix_construct_redim(double,eFnp1,3,3);
  Matrix_construct_redim(double,Fnp1,3,3);
  /* after update (i.e., converged step) the *Fn = *Fnp1 */
  m->param->get_Fn(m,&Fnp1);
  m->param->get_eFn(m,&eFnp1);

  constitutive_model_update_elasticity(m,&eFnp1,0.0,NULL,S,compute_stiffness);  
  Matrix_det(Fnp1, *Jnp1);
  
  Matrix_cleanup(Fnp1);
  Matrix_cleanup(eFnp1);   
}                          
void post_processing_compute_stress(double *GS, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts)
                    
{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM)
    intg_order = 0;
    
  if(opts->analysis_type==CM && opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;
  
  int nsd = 3;
  Matrix(double) F,S,LS;
  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,S,3,3,0.0);
  Matrix_construct_init(double,LS,3,3,0.0);  

  double LV = 0.0;
  double GV = 0.0;
  
  for(int e = 0; e<ne; e++)
  {        
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe, e, elem, node, intg_order,total_Lagrangian);
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

    if(opts->analysis_type==TF)
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
      double Jnp1 = 1.0;
      switch(opts->analysis_type)
      {
        case DISP:
          post_processing_compute_stress_disp_ip(&fe,e,S,hommat,elem,F,Pn);
          break;
        case TF:
        {
          FEMLIB_elem_shape_function(&fe,ip,npres, Np);
  
          for(int a=1; a<=npres; a++)
            Pn += Vec_v(Np,a)*Vec_v(P,a);
  
          post_processing_compute_stress_3f_ip(&fe,e,S,hommat,elem,F,Pn);
          break;
        }           
        case CM:
        {
          post_processing_compute_stress4CM(&fe,e,ip,&S,&Jnp1,hommat,elem,eps);
          break;
        }      
        default:
          break;
      }
      
      LV += fe.detJxW/Jnp1;

      for(int a=0; a<9; a++)
        LS.m_pdata[a] += S.m_pdata[a]*fe.detJxW/Jnp1;
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
  
  for(int a=0; a<9; a++)
    GS[a] = GS[a]/GV;    
}

void post_processing_deformation_gradient(double *GF, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts)
                    
{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM)
    intg_order = 0;
    
  if(opts->analysis_type==CM && opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;    
  
  int nsd = 3;
  Matrix(double) F,LF;
  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,LF,3,3,0.0);  

  double LV = 0.0;
  double GV = 0.0;
  
  for(int e = 0; e<ne; e++)
  {        
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
    
    for(int ip = 1; ip<=fe.nint; ip++)
    {      
      FEMLIB_elem_basis_V(&fe, ip);  
      FEMLIB_update_shape_tensor(&fe);
      FEMLIB_update_deformation_gradient(&fe,nsd,u.m_pdata,F);
      
      if(opts->analysis_type==CM)
      { 
        Constitutive_model *m = &(eps[e].model[ip-1]);
        double Jnp1 = 1.0;
        Matrix(double) Fnp1;
        Matrix_construct_redim(double,Fnp1,3,3);
        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_Fn(m,&Fnp1);
        if(!total_Lagrangian)
          Matrix_det(Fnp1, Jnp1);
        
        LV += fe.detJxW/Jnp1;
        for(int a=0; a<9; a++)
          LF.m_pdata[a] += Fnp1.m_pdata[a]*fe.detJxW/Jnp1;

        Matrix_cleanup(Fnp1);
      }
      else
      {        
        LV += fe.detJxW;
        for(int a=0; a<9; a++)
          LF.m_pdata[a] += F.m_pdata[a]*fe.detJxW;
      }
    }
    FEMLIB_destruct(&fe);
    Matrix_cleanup(u);
  }
      
  MPI_Allreduce(LF.m_pdata,GF,9,MPI_DOUBLE,MPI_SUM,mpi_comm);
  MPI_Allreduce(&LV,&GV,1,MPI_DOUBLE,MPI_SUM,mpi_comm);    
  Matrix_cleanup(F);
  Matrix_cleanup(LF);
  
  for(int a=0; a<9; a++)
    GF[a] = GF[a]/GV;
}

void post_processing_deformation_gradient_elastic_part(double *GF, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts)
                    
{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM)
    intg_order = 0;
    
  if(opts->analysis_type==CM && opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;    
  
  int nsd = 3;
  Matrix(double) F,LF;
  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,LF,3,3,0.0);  

  double LV = 0.0;
  double GV = 0.0;
  
  for(int e = 0; e<ne; e++)
  {        
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
    
    for(int ip = 1; ip<=fe.nint; ip++)
    {      
      FEMLIB_elem_basis_V(&fe, ip);  
      FEMLIB_update_shape_tensor(&fe);
      FEMLIB_update_deformation_gradient(&fe,nsd,u.m_pdata,F);
      
      if(opts->analysis_type==CM)
      { 
        Constitutive_model *m = &(eps[e].model[ip-1]);
        double Jnp1 = 1.0;
        Matrix(double) Fnp1, eFnp1;
        Matrix_construct_redim(double,Fnp1,3,3);
        Matrix_construct_redim(double,eFnp1,3,3);        
        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_Fn(m,&Fnp1);
        m->param->get_eFn(m,&eFnp1);
        if(!total_Lagrangian)        
          Matrix_det(Fnp1, Jnp1);
        
        LV += fe.detJxW/Jnp1;
        for(int a=0; a<9; a++)
          LF.m_pdata[a] += eFnp1.m_pdata[a]*fe.detJxW/Jnp1;

        Matrix_cleanup(Fnp1);
        Matrix_cleanup(eFnp1);        
      }
      else
      {        
        LV += fe.detJxW;
        for(int a=0; a<9; a++)
          LF.m_pdata[a] += F.m_pdata[a]*fe.detJxW;
      }
    }
    FEMLIB_destruct(&fe);
    Matrix_cleanup(u);
  }
      
  MPI_Allreduce(LF.m_pdata,GF,9,MPI_DOUBLE,MPI_SUM,mpi_comm);
  MPI_Allreduce(&LV,&GV,1,MPI_DOUBLE,MPI_SUM,mpi_comm);    
  Matrix_cleanup(F);
  Matrix_cleanup(LF);
  
  for(int a=0; a<9; a++)
    GF[a] = GF[a]/GV;
}


void post_processing_plastic_hardness(double *G_gn, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts)
                    
{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM)
    intg_order = 0;
    
  if(opts->analysis_type==CM && opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;    
  
  int nsd = 3;
  double L_gn = 0.0;

  double LV = 0.0;
  double GV = 0.0;
  
  for(int e = 0; e<ne; e++)
  {        
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe, e, elem, node, intg_order,total_Lagrangian);
    int nne = fe.nne;
        
    for(int ip = 1; ip<=fe.nint; ip++)
    {      
      FEMLIB_elem_basis_V(&fe, ip);  
      FEMLIB_update_shape_tensor(&fe);
      
      Constitutive_model *m = &(eps[e].model[ip-1]);        
      double g_n = 0.0;      
      m->param->get_hardening(m,&g_n);      
      
      
      double Jnp1 = 0.0;
      Matrix(double) Fnp1;
      Matrix_construct_redim(double,Fnp1,3,3);
      /* after update (i.e., converged step) the *Fn = *Fnp1 */
      m->param->get_Fn(m,&Fnp1);
      Matrix_det(Fnp1, Jnp1);
      
      LV += fe.detJxW/Jnp1;
      L_gn += g_n*fe.detJxW/Jnp1;

      Matrix_cleanup(Fnp1);
    }
    FEMLIB_destruct(&fe);
  }
      
  MPI_Allreduce(&L_gn,G_gn,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  MPI_Allreduce(&LV,&GV,1,MPI_DOUBLE,MPI_SUM,mpi_comm);    
  
  (*G_gn) = (*G_gn)/GV;
}

void post_processing_potential_energy(double *GE, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts)
                    
{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM)
    intg_order = 0;
    
  if(opts->analysis_type==CM && opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0; 
  
  int nsd = 3;
  Matrix(double) F,C;
  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,C,3,3,0.0);  

  double LE = 0.0;

  double LV = 0.0;
  double GV = 0.0;
   
  for(int e = 0; e<ne; e++)
  { 
    int mat = elem[e].mat[2];
    double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));    

    devPotentialFuncPtr dev_func = getDevPotentialFunc(1,&hommat[mat]);
    UFuncPtr              U_func =            getUFunc(1,&hommat[mat]);

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
    
    for(int ip = 1; ip<=fe.nint; ip++)
    {      
      FEMLIB_elem_basis_V(&fe, ip);  
      FEMLIB_update_shape_tensor(&fe);
      FEMLIB_update_deformation_gradient(&fe,nsd,u.m_pdata,F);
      
      if(opts->analysis_type==CM)
      { 
        Constitutive_model *m = &(eps[e].model[ip-1]);
        double Jnp1 = 1.0;
        Matrix(double) Fnp1, eFnp1;
        Matrix_construct_redim(double,Fnp1,3,3);
        Matrix_construct_redim(double,eFnp1,3,3);
               
        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_Fn(m,&Fnp1);
        m->param->get_eFn(m,&eFnp1);
        if(!total_Lagrangian)        
          Matrix_det(Fnp1, Jnp1);
        
        Matrix_AxB(C,1.0,0.0,eFnp1,1,eFnp1,0);        
        double W = 0.0;
        dev_func(C.m_pdata,&hommat[mat],&W);          
        double U = 0.0;
        double J = 0.0;
        Matrix_det(eFnp1, J);
        U_func(J,&hommat[mat],&U);

        LE += (W+kappa*U)*fe.detJxW/Jnp1;
        
        Matrix_cleanup(Fnp1);
        Matrix_cleanup(eFnp1);
                
      }
      else
      {        
        Matrix_AxB(C,1.0,0.0,F,1,F,0);        
        double W = 0.0;
        dev_func(C.m_pdata,&hommat[mat],&W);          
        double U = 0.0;
        double J = 0.0;
        Matrix_det(F, J);
        U_func(J,&hommat[mat],&U);

        LE += (W+kappa*U)*fe.detJxW;
      }
    }
    FEMLIB_destruct(&fe);
    Matrix_cleanup(u);    
  }
      
  MPI_Allreduce(&LE,GE,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  Matrix_cleanup(F);
  Matrix_cleanup(C);
}

void post_processing_deformed_volume(double *GV, ELEMENT *elem, long ne, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts)
                    
{
  int total_Lagrangian = 1;
  int intg_order = 1;
  
  if(opts->analysis_type==CM && opts->cm==CRYSTAL_PLASTICITY)
  {  
    total_Lagrangian = PLASTICITY_TOTAL_LAGRANGIAN;
    intg_order = 0;
  }     
  
  int nsd = 3;

  double LV = 0.0;
  *GV = 0.0;
  
  for(int e = 0; e<ne; e++)
  {        
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
    
    for(int ip = 1; ip<=fe.nint; ip++)
    {      
      FEMLIB_elem_basis_V(&fe, ip);  
      FEMLIB_update_shape_tensor(&fe);
      
      if(opts->analysis_type==CM && opts->cm==CRYSTAL_PLASTICITY)
      { 
        Constitutive_model *m = &(eps[e].model[ip-1]);
        double Jnp1 = 1.0;
        Matrix(double) Fnp1;
        Matrix_construct_redim(double,Fnp1,3,3);
        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_Fn(m,&Fnp1);
        if(!total_Lagrangian)
          Matrix_det(Fnp1, Jnp1);
        
        LV += fe.detJxW/Jnp1;
        Matrix_cleanup(Fnp1);
      }
      else
      {        
        LV += fe.detJxW;
      }
    }
    FEMLIB_destruct(&fe);
    Matrix_cleanup(u);
  }
      
  MPI_Allreduce(&LV,GV,1,MPI_DOUBLE,MPI_SUM,mpi_comm);    
  
}

