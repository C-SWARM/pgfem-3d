#ifndef __MMS_L2NORM_h__
#define __MMS_L2NORM_h__

#include<math.h>
#include "enumerations.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "elem3d.h"

#include "MMS.h"
#include "femlib.h"

void compute_L2_error(double *GL2_err, ELEMENT *elem, long ne, NODE *node, double* r, double *Ph, double *Vh, double t, MPI_Comm mpi_comm, const PGFem3D_opt *opts, const HOMMAT *hommat)
{
  int myrank = 0;
  MPI_Comm_rank (mpi_comm,&myrank);
  int ndofn = 3;

  double L2_err[3];
   L2_err[0] =  L2_err[1] =  L2_err[2] = 0.0;
  GL2_err[0] = GL2_err[1] = GL2_err[2]= 0.0; 
  
  for(long e = 0; e<ne; e++)
  {  
    
    const int mat = elem[e].mat[2];    
    long nne = elem[e].toe;
    long *nod = aloc1l(nne);
    elemnodes(e,nne,nod,elem);
    
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe, e, elem, node, 3, 1);
                 
    for(int ip = 1; ip<=fe.nint; ip++)
    {
      FEMLIB_elem_basis_V(&fe, ip); 
      double u[3], uh[3], du[3];
      
      double x_ip = Vec_v(fe.x_ip,1);
      double y_ip = Vec_v(fe.x_ip,2);
      double z_ip = Vec_v(fe.x_ip,3);
                
      MMS_displacement(u, t, x_ip, y_ip, z_ip);      
      
      uh[0] = uh[1] = uh[2] = 0.0;                  
      for(long a = 1; a<=nne; a++)
      {        
        long nid = nod[a-1];     
        
        uh[0] += Vec_v(fe.N,a)*r[nid*ndofn + 0];
        uh[1] += Vec_v(fe.N,a)*r[nid*ndofn + 1];
        uh[2] += Vec_v(fe.N,a)*r[nid*ndofn + 2];        
      }
      
      double P = 0;
      double V = 0;
      if(opts->analysis_type==TF)
      {
        MMS_pressure_volume(&P, &V, &hommat[mat], t, x_ip, y_ip, z_ip);
      }         
      
      du[0] = u[0] - uh[0];
      du[1] = u[1] - uh[1];
      du[2] = u[2] - uh[2];
                          
      L2_err[0] += (du[0]*du[0] + du[1]*du[1] + du[2]*du[2])*fe.detJxW;
      L2_err[1] += (Ph[e]-P)*(Ph[e]-P)*fe.detJxW;      
      L2_err[2] += (Vh[e]-V)*(Vh[e]-V)*fe.detJxW;            
    }
    
    FEMLIB_destruct(&fe);
    dealoc1(nod);
  }
  
  MPI_Allreduce(L2_err,GL2_err,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
};  

#endif