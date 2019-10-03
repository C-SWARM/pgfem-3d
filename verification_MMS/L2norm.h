#ifndef PGFEM3D_MMS_L2NORM_H
#define PGFEM3D_MMS_L2NORM_H

#include "allocation.h"
#include "elem3d.h"
#include "enumerations.h"
#include "femlib.h"
#include "MMS.h"
#include "utils.h"
#include "pgfem3d/Communication.hpp"
#include <cmath>

using namespace multiscale::net;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double AnalyticalTemperature(const double t, const double x, const double y, const double z){
  return exp(-t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
}

void compute_L2_Thermal(double *GL2_err, 
                        Element *elem,
                        long ne,
                        Node *node,
                        double* r,
                        double t, 
                        const pgfem3d::CommunicationStructure *com)
{
  double L2_err = 0.0;
  *GL2_err = 0.0;
  
  for(long e = 0; e<ne; e++)
  {    
    FEMLIB fe(e,elem,node, 1, 1);

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);

      double x_ip = fe.x_ip(0);
      double y_ip = fe.x_ip(1);
      double z_ip = fe.x_ip(2);

      double T = AnalyticalTemperature(t, x_ip, y_ip, z_ip);
      double Th = 0.0;

      for(long a = 0; a<fe.nne; a++)
        Th += fe.N(a)*r[fe.node_id(a)];

      double dT = T - Th;

      L2_err += dT*dT*fe.detJxW;
    }
  }

  com->net->allreduce(&L2_err,GL2_err,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
}

void compute_L2_Mechanical_MMS(double *GL2_err, Element *elem, long ne, Node *node, double* r,
		      double *Ph, double *Vh, double t,
		      const pgfem3d::CommunicationStructure *com,
		      const PGFem3D_opt *opts, const HOMMAT *hommat,
		      const bool is4cm)
{
  int ndofn = 3;

  double L2_err[3];
  L2_err[0] =  L2_err[1] =  L2_err[2] = 0.0;
  GL2_err[0] = GL2_err[1] = GL2_err[2]= 0.0;
  
  for(long e = 0; e<ne; e++)
  {

    const int mat = elem[e].mat[2];

    MATERIAL_ELASTICITY mat_e;

    set_properties_using_E_and_nu(&mat_e,hommat[mat].E,hommat[mat].nu);
    mat_e.m01   = hommat[mat].m01;
    mat_e.m10   = hommat[mat].m10;
    mat_e.G     = hommat[mat].G;
    mat_e.kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));
    mat_e.devPotFlag = hommat[mat].devPotFlag;
    mat_e.volPotFlag = hommat[mat].volPotFlag;
    
    HyperElasticity elast;
    elast.construct_elasticity(&mat_e, true);
    
    long nne = elem[e].toe;
    long *nod = aloc1l(nne);
    elemnodes(e,nne,nod,elem);

    FEMLIB fe(e,elem,node, 1, 1);

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      double u[3], uh[3], du[3];

      double x_ip = fe.x_ip(0);
      double y_ip = fe.x_ip(1);
      double z_ip = fe.x_ip(2);

      MMS_displacement(u, t, x_ip, y_ip, z_ip, is4cm);

      uh[0] = uh[1] = uh[2] = 0.0;
      for(long a = 0; a<nne; a++)
      {
        long nid = nod[a];

        uh[0] += fe.N(a)*r[nid*ndofn + 0];
        uh[1] += fe.N(a)*r[nid*ndofn + 1];
        uh[2] += fe.N(a)*r[nid*ndofn + 2];
      }

      double P = 0;
      double V = 0;
      if(opts->analysis_type==TF || opts->analysis_type==CM3F)
      {
        MMS_pressure_volume(&P, &V, &elast, t, x_ip, y_ip, z_ip, is4cm);
      }

      du[0] = u[0] - uh[0];
      du[1] = u[1] - uh[1];
      du[2] = u[2] - uh[2];

      L2_err[0] += (du[0]*du[0] + du[1]*du[1] + du[2]*du[2])*fe.detJxW;
      L2_err[1] += (Ph[e]-P)*(Ph[e]-P)*fe.detJxW;
      L2_err[2] += (Vh[e]-V)*(Vh[e]-V)*fe.detJxW;
    }

    dealoc1(nod);
  }

  com->net->allreduce(L2_err,GL2_err,3,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
}

#endif // #define PGFEM3D_MMS_L2NORM_H
