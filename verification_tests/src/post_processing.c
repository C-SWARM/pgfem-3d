#include "post_processing.h"
#include "enumerations.h"
#include "utils.h"
#include "allocation.h"


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

void post_processing_compute_stress(double *GS, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, int analysis)
{
  int nsd;
  Matrix(double) F,S,LS;
  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,S,3,3,0.0);
  Matrix_construct_init(double,LS,3,3,0.0);
  double LV = 0.0;
  double GV = 0.0;
  
  for(int e = 0; e<ne; e++)
  {
    
    int nne = elem[e].toe;
    long *nod = aloc1l (nne);
    elemnodes(e,nne,nod,elem);
        
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe, e, elem, node);
    Matrix(double) u;  
    Matrix_construct_init(double,u,nne*ndofn,1,0.0);
                        
    for(int a = 0; a<nne; a++)
    {      
      int nid = Vec_v(fe.node_id, a+1);
      for(int b=0; b<ndofn; b++)
        Vec_v(u, a*ndofn+b+1) = r[nid*ndofn + b];
    }
    
    if(analysis!=DISP)
    {  
      printf("Only displacement based analysis is now supported\n");
      return;
    }
    
    for(int ip = 1; ip<=fe.nint; ip++)
    {
      FEMLIB_elem_basis_V(&fe, ip);  
      FEMLIB_update_shape_tensor(&fe);
      FEMLIB_update_deformation_gradient(&fe,ndofn,u.m_pdata,F);

      double Pn = 0.0;          
      post_processing_compute_stress_disp_ip(&fe,e,S,hommat,elem,F,Pn);
        
      
      LV += fe.detJxW;
                                                          
      for(int a=0; a<9; a++)
        LS.m_pdata[a] += S.m_pdata[a]*fe.detJxW;
	  }

    FEMLIB_destruct(&fe);
    Matrix_cleanup(u);

  }
  
  Matrix_cleanup(F);
  Matrix_cleanup(S);
    
  MPI_Allreduce(LS.m_pdata,GS,9,MPI_DOUBLE,MPI_SUM,mpi_comm);
  MPI_Allreduce(&LV,&GV,1,MPI_DOUBLE,MPI_SUM,mpi_comm);  
  Matrix_cleanup(LS);
  
  for(int a=0; a<9; a++)
    GS[a] = GS[a]/GV;
};