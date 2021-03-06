#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "post_processing.h"
#include "PGFem3D_to_VTK.hpp"
#include "allocation.h"
#include "constitutive_model.h"
#include "enumerations.h"
#include "utils.h"
#include <ttl/ttl.h>
#include "PGFem3D_data_structure.h"

using namespace pgfem3d;
using namespace multiscale::net;

//ttl declarations
namespace {
using namespace ttl;
constexpr Index<'I'> I;
constexpr Index<'J'> J;
constexpr Index<'K'> K;
constexpr Index<'L'> L;

template<class T1, class T2>
inline int inv(T1 &A, T2 &AI) {
  int err = inv3x3(A.data, AI.data);
  return err;
}
}

int read_from_VTK(const PGFem3D_opt *opts, int myrank, int step, double *u)
{
  int err = 0;
  char filename[1024];
  snprintf(filename, sizeof(filename), "%s/VTK/STEP_%.5d/%s_%d_%d.vtu",
           opts->opath, step,opts->ofname,myrank, step);
  err += read_VTK_file(filename, u);
  return err;
}

void post_processing_compute_stress_disp_ip(FEMLIB *fe, int e, double *S_in, HOMMAT *hommat, Element *elem,
                                            double *F_in, double Pn)
{
  Tensor<2,3,double*> F(F_in), S(S_in);
  Tensor<2,3,double> C,CI,devS;

  C = F(K,I)*F(K,J);
  double detF = det(F);
  inv(C,CI);

  int mat = elem[e].mat[2];
  double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
  Stress(C.data,&hommat[mat],devS.data);

  dUdJFuncPtr DUDJ = getDUdJFunc(1,&hommat[mat]);
  double dUdJ = 0.0;
  DUDJ(detF,&hommat[mat],&dUdJ);

  //  for(int a = 0; a<9; a++)
  //    S.m_pdata[a] = devS.m_pdata[a] + kappaJdUdJ*CI.m_pdata[a];
  S(I,J) = devS(I,J) + kappa*detF*dUdJ*CI(I,J);
}

void post_processing_compute_stress_3f_ip(FEMLIB *fe, int e, double *S_in, HOMMAT *hommat, Element *elem,
                                          double *F_in, double theta)
{
  Tensor<2,3,double*> F(F_in), S(S_in);
  Tensor<2,3,double> C,CI,devS;

  C = F(K,I)*F(K,J);
  inv(C,CI);

  int mat = elem[e].mat[2];
  double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));

  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
  dUdJFuncPtr UP = getDUdJFunc(1, &hommat[mat]);
  Stress(C.data,&hommat[mat],devS.data);

  double Up = 0.0;
  UP(theta,&hommat[mat],&Up);
  double JPn = kappa*theta*Up;

  S(I,J) = devS(I,J) + JPn*CI(I,J);
}

void post_processing_compute_stress4CM(FEMLIB *fe, int e, int ip, double *S, double *Jnp1,
                                       HOMMAT *hommat, Element *elem, EPS *eps)
{
  bool compute_stiffness = false;

  Constitutive_model *m = &(eps[e].model[ip]);
  Tensor<2,3,double> Fnp1, eFnp1;
  m->param->get_F(m,Fnp1.data,1);
  m->param->get_eF(m,eFnp1.data,1);

  constitutive_model_default_update_elasticity(m,eFnp1.data,NULL,S,compute_stiffness);
  *Jnp1 = det(Fnp1);
}
void post_processing_compute_stress(double *GS, Element *elem, HOMMAT *hommat, long ne, int npres, Node *node, EPS *eps,
                                    double* r, double *Vnp1, int ndofn,
				    const CommunicationStructure *com, const PGFem3D_opt *opts)

{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM || opts->analysis_type==CM3F)
    intg_order = 0;

  if((opts->analysis_type==CM || opts->analysis_type==CM3F)&& opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  int nsd = 3;
  Tensor<2,3,double> F, S = {}, LS = {};

  double LV = 0.0;
  double GV = 0.0;

  for(int e = 0; e<ne; e++)
  {
    FEMLIB fe(e, elem, node, intg_order,total_Lagrangian);
    int nne = fe.nne;

    Matrix<double> Np, u;
    Np.initialization(npres,1,0.0);
    u.initialization(nne*nsd,1,0.0);

    for(int a = 0; a<nne; a++)
    {
      int nid = fe.node_id(a);
      for(int b=0; b<nsd; b++)
      {
        u(a*nsd+b) = r[nid*ndofn + b];
      }
    }

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();
      fe.update_deformation_gradient(nsd,u.m_pdata,F.data);
      double Pn = 0.0;
      double Jnp1 = 1.0;
      switch(opts->analysis_type)
      {
       case DISP:
        post_processing_compute_stress_disp_ip(&fe,e,S.data,hommat,elem,F.data,Pn);
        break;
       case TF:
       {
         double theta = Vnp1[e];
         fe.elem_shape_function(ip,npres, Np.m_pdata);
         post_processing_compute_stress_3f_ip(&fe,e,S.data,hommat,elem,F.data,theta);
         break;
       }
       case CM:   // intended to flow
       case CM3F:
         {
           post_processing_compute_stress4CM(&fe,e,ip,S.data,&Jnp1,hommat,elem,eps);
           break;
         }
       default:
        break;
      }

      LV += fe.detJxW/Jnp1;

      for(int a=0; a<9; a++)
        LS.data[a] += S.data[a]*fe.detJxW/Jnp1;
    }
  }
  com->net->allreduce(LS.data,GS,9,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  com->net->allreduce(&LV,&GV,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);

  for(int a=0; a<9; a++)
    GS[a] = GS[a]/GV;
}

void post_processing_deformation_gradient(double *GF, Element *elem, HOMMAT *hommat, long ne, int npres, Node *node, EPS *eps,
                                          double* r, int ndofn, const CommunicationStructure *com, const PGFem3D_opt *opts)

{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM || opts->analysis_type==CM3F)
    intg_order = 0;

  if((opts->analysis_type==CM || opts->analysis_type==CM3F) && opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  int nsd = 3;
  Tensor<2,3,double> F,LF = {};

  double LV = 0.0;
  double GV = 0.0;

  for(int e = 0; e<ne; e++)
  {
    FEMLIB fe(e, elem, node, intg_order,total_Lagrangian);
    int nne = fe.nne;

    Matrix<double> u(nne,nsd);

    for(int a = 0; a<nne; a++)
    {
      int nid = fe.node_id(a);
      for(int b=0; b<nsd; b++)
      {
        u(a,b) = r[nid*ndofn+b];
      }
    }

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();
      fe.update_deformation_gradient(nsd,u.m_pdata,F.data);

      if(opts->analysis_type==CM || opts->analysis_type==CM3F)
      {
        Constitutive_model *m = &(eps[e].model[ip]);
        double Jnp1 = 1.0;
        Tensor<2,3,double> Fnp1 = {};
        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_F(m,Fnp1.data,1);
        if(!total_Lagrangian)
          Jnp1 = det(Fnp1);

        LV += fe.detJxW/Jnp1;
        for(int a=0; a<9; a++)
          LF.data[a] += Fnp1.data[a]*fe.detJxW/Jnp1;
      }
      else
      {
        LV += fe.detJxW;
        for(int a=0; a<9; a++)
          LF.data[a] += F.data[a]*fe.detJxW;
      }
    }
  }

  com->net->allreduce(LF.data,GF,9,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  com->net->allreduce(&LV,&GV,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);

  for(int a=0; a<9; a++)
    GF[a] = GF[a]/GV;
}

void post_processing_deformation_gradient_elastic_part(double *GF, Element *elem, HOMMAT *hommat, long ne, int npres, Node *node, EPS *eps,
                                                       double* r, int ndofn, const CommunicationStructure *com, const PGFem3D_opt *opts)

{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM || opts->analysis_type==CM3F)
    intg_order = 0;

  if((opts->analysis_type==CM || opts->analysis_type==CM3F)&& opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  int nsd = 3;
  Tensor<2,3,double> F,LF = {};

  double LV = 0.0;
  double GV = 0.0;

  for(int e = 0; e<ne; e++)
  {
    FEMLIB fe(e, elem, node, intg_order,total_Lagrangian);
    int nne = fe.nne;

    Matrix<double> u(nne,nsd);

    for(int a = 0; a<nne; a++)
    {
      int nid = fe.node_id(a);
      for(int b=0; b<nsd; b++)
      {
        u(a,b) = r[nid*ndofn+b];
      }
    }

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();
      fe.update_deformation_gradient(nsd,u.m_pdata,F.data);

      if(opts->analysis_type==CM || opts->analysis_type==CM3F)
      {
        Constitutive_model *m = &(eps[e].model[ip]);
        double Jnp1 = 1.0;
        Tensor<2,3,double> Fnp1 = {}, eFnp1 = {};

        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_F(m,Fnp1.data,1);
        m->param->get_eF(m,eFnp1.data,1);
        if(!total_Lagrangian)
          Jnp1 = det(Fnp1);

        LV += fe.detJxW/Jnp1;
        for(int a=0; a<9; a++)
          LF.data[a] += eFnp1.data[a]*fe.detJxW/Jnp1;
      }
      else
      {
        LV += fe.detJxW;
        for(int a=0; a<9; a++)
          LF.data[a] += F.data[a]*fe.detJxW;
      }
    }
  }

  com->net->allreduce(LF.data,GF,9,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  com->net->allreduce(&LV,&GV,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);

  for(int a=0; a<9; a++)
    GF[a] = GF[a]/GV;
}


void post_processing_plastic_hardness(double *G_gn, Element *elem, HOMMAT *hommat, long ne, int npres, Node *node, EPS *eps,
                                      double* r, int ndofn, const CommunicationStructure *com, const PGFem3D_opt *opts)

{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM || opts->analysis_type==CM3F)
    intg_order = 0;

  if((opts->analysis_type==CM || opts->analysis_type==CM3F) && opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  double L_gn = 0.0;

  double LV = 0.0;
  double GV = 0.0;

  for(int e = 0; e<ne; e++)
  {
    FEMLIB fe(e, elem, node, intg_order,total_Lagrangian);

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();

      Constitutive_model *m = &(eps[e].model[ip]);
      double g_n = 0.0;
      m->param->get_hardening(m,&g_n,2);

      Tensor<2,3,double> Fnp1;
      /* after update (i.e., converged step) the *Fn = *Fnp1 */
      m->param->get_F(m,Fnp1.data,1);
      double Jnp1 = det(Fnp1);

      LV += fe.detJxW/Jnp1;
      L_gn += g_n*fe.detJxW/Jnp1;
    }
  }

  com->net->allreduce(&L_gn,G_gn,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  com->net->allreduce(&LV,&GV,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);

  (*G_gn) = (*G_gn)/GV;
}

void post_processing_potential_energy(double *GE, Element *elem, HOMMAT *hommat, long ne, int npres, Node *node, EPS *eps,
                                      double* r, int ndofn, const CommunicationStructure *com, const PGFem3D_opt *opts)

{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM || opts->analysis_type==CM3F)
    intg_order = 0;

  if((opts->analysis_type==CM || opts->analysis_type==CM3F)&& opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  int nsd = 3;
  Tensor<2,3,double> F,C;

  double LE = 0.0;

  for(int e = 0; e<ne; e++)
  {
    int mat = elem[e].mat[2];
    double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));

    devPotentialFuncPtr dev_func = getDevPotentialFunc(1,&hommat[mat]);
    UFuncPtr              U_func =            getUFunc(1,&hommat[mat]);

    FEMLIB fe(e, elem, node, intg_order,total_Lagrangian);
    int nne = fe.nne;

    Matrix<double> u(nne,nsd);

    for(int a = 0; a<nne; a++)
    {
      int nid = fe.node_id(a);
      for(int b=0; b<nsd; b++)
      {
        u(a,b) = r[nid*ndofn+b];
      }
    }

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();
      fe.update_deformation_gradient(nsd,u.m_pdata,F.data);

      if(opts->analysis_type==CM || opts->analysis_type==CM3F)
      {
        Constitutive_model *m = &(eps[e].model[ip]);
        double Jnp1 = 1.0;
        Tensor<2,3,double> Fnp1, eFnp1;

        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_F(m,Fnp1.data,1);
        m->param->get_eF(m,eFnp1.data,1);
        if(!total_Lagrangian)
          Jnp1 = det(Fnp1);

        C = eFnp1(K,I)*eFnp1(K,L);
        double W = 0.0;
        dev_func(C.data,&hommat[mat],&W);
        double U = 0.0;
        double detF = det(eFnp1);
        U_func(detF,&hommat[mat],&U);

        LE += (W+kappa*U)*fe.detJxW/Jnp1;
      }
      else
      {
        C = F(K,I)*F(K,J);
        double W = 0.0;
        dev_func(C.data,&hommat[mat],&W);
        double U = 0.0;
        double detF = det(F);
        U_func(detF,&hommat[mat],&U);

        LE += (W+kappa*U)*fe.detJxW;
      }
    }
  }

  com->net->allreduce(&LE,GE,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
}

void post_processing_deformed_volume(double *GV, Element *elem, long ne, Node *node, EPS *eps,
                                     double* r, int ndofn, const CommunicationStructure *com,
				     const PGFem3D_opt *opts)

{
  int total_Lagrangian = 1;
  int intg_order = 1;

  if(opts->analysis_type==CM || opts->analysis_type==CM3F)
    intg_order = 0;

  if((opts->analysis_type==CM || opts->analysis_type==CM3F)&& opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  int nsd = 3;

  double LV = 0.0;
  *GV = 0.0;

  for(int e = 0; e<ne; e++)
  {
    FEMLIB fe(e, elem, node, intg_order,total_Lagrangian);
    int nne = fe.nne;

    Matrix<double> u(nne,nsd);

    for(int a = 0; a<nne; a++)
    {
      int nid = fe.node_id(a);
      for(int b=0; b<nsd; b++)
      {
        u(a,b) = r[nid*ndofn+b];
      }
    }

    for(int ip = 0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();

      if(opts->analysis_type==CM && opts->cm==CRYSTAL_PLASTICITY)
      {
        Constitutive_model *m = &(eps[e].model[ip]);
        double Jnp1 = 1.0;
        Tensor<2,3,double> Fnp1;

        /* after update (i.e., converged step) the *Fn = *Fnp1 */
        m->param->get_F(m,Fnp1.data,1);
        if(!total_Lagrangian)
          Jnp1 = det(Fnp1);

        LV += fe.detJxW/Jnp1;
      }
      else
      {
        LV += fe.detJxW;
      }
    }
  }

  com->net->allreduce(&LV,GV,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
}

/// Compute and print maximum element pressure.
/// Each process compute maximum pressure for each material and gather 
/// all maximum pressure and print maximum pressure with coodinates where
/// maximum pressure measured.
///
/// \param[in] grid     an object containing all mesh info
/// \param[in] mat      a material object
/// \param[in] fv       object for field variables
/// \param[in] com      communication object
void post_processing_max_pressure(const Grid &grid,
                                  const MaterialProperty &mat,
                                  const CommunicationStructure *com,
                                  const FieldVariables &fv)
{
  int myrank = com->rank;
  
  Matrix<double> LP(com->nproc, mat.nmat, 0.0);
  Matrix<double> GP(com->nproc, mat.nmat, 0.0);
  Matrix<int> eid(mat.nmat, 1, 0);
  
  for(int ia=0; ia<grid.ne; ia++)
  {
    int mat_id = (grid.element[ia]).mat[0];
    double P = (fv.sig[ia].el.o[0] + fv.sig[ia].el.o[1] + fv.sig[ia].el.o[2])/3.0;
    
    if(fabs(LP(myrank, mat_id)) < fabs(P))
    {  
        LP(myrank, mat_id) = P;
        eid(mat_id) = ia;
    }
  }
  
  com->net->allreduce(LP.m_pdata,GP.m_pdata,com->nproc*mat.nmat,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  
  Matrix<bool> have_max(mat.nmat, 1, true);
  for(int ia=0; ia<com->nproc; ia++)
  {
    if(ia==myrank)
      continue;
    
    for(int ib=0; ib<mat.nmat; ib++)
    {    
      if(fabs(GP(ia, ib)) > fabs(LP(myrank, ib)))
        have_max(ib) = false;
    }
  }
  
  for(int ia=0; ia<mat.nmat; ia++)
  {  
    if(have_max(ia))
    {
      FEMLIB fe(eid(ia), grid.element, grid.node, 0,1);
      double x[3];
      x[0] = x[1] = x[2] = 0.0;
    
      for(int ib=0; ib<fe.nne; ib++)
      {
        for(int ic=0; ic<3; ic++)
          x[ic] += fe.node_coord(ib, ic);
      }
    
      x[0] /= fe.nne;
      x[1] /= fe.nne;
      x[2] /= fe.nne;
                          
      PGFEM_printf("max. pressure for (mat=%d): %e at process = %d, element = %d (xyz = %e, %e, %e)\n", 
                    ia, LP(myrank, ia), myrank, eid(ia), x[0], x[1], x[2]);
    }
  }
}

/// compute surface area of face on which 
/// Neumann boundary (NB) conditions are applied
///
/// \param[in] grid an object containing all mesh info
/// \param[in] com  communication object
void post_processing_compute_NBE_area(Grid &grid,
                                      const pgfem3d::CommunicationStructure *com){

  if(grid.NBE.m_row*grid.NBE.m_col == 0)
    return;
          
  for(int mp_id = 0; mp_id<grid.NBE.m_row*grid.NBE.m_col; ++mp_id){

    double lA = 0.0;
    double gA = 0.0;
              
    for(int nbe_id = 0; nbe_id < grid.NBE(mp_id).nbe_no; ++nbe_id)
    {
      const int eid = grid.NBE(mp_id).element_ids(nbe_id, 0);
      FEMLIB fe(eid, grid.element, grid.node, 1, 1);

      // number of boundary element to be integrated
      const int bnd_elem_no = grid.NBE(mp_id).element_ids(nbe_id, 1);      
      
      for(int iA = 0; iA<bnd_elem_no; ++iA){      
        const int face_id = grid.NBE(mp_id).bnd_elements(nbe_id)(iA, 0);      
        FemLibBoundary fes(&fe, face_id, fe.intg_order);
      
        // apply quadrature rule 
        for(int ip=0; ip<fes.nint; ++ip){
          fes.elem_basis_S(ip);
          lA += fes.detJxW;
        }
      }
    }
    
  
    com->net->allreduce(&lA,&gA,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
    
    if(com->rank==0)
      std::cout << "Surface area (Neumann boundary (NB) conditions for physics " << mp_id << ") = " << gA << endl;
  }
}
                                    
