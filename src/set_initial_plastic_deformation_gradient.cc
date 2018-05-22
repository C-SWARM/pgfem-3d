#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "femlib.h"
#include "constitutive_model.h"
#include "set_initial_plastic_deformation_gradient.h"
#include "PGFem3D_data_structure.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "PLoc_Sparse.h"
#include "utils.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

/// Fill the B matrix, size of B have to be created before call this function.
///
/// \param[out] &B  filled B matrix
/// \param[in]  &dN derivative of shape function
/// \param[in]  a   index of node (local id of an element)
void get_B(Matrix<double> &B,
           Matrix<double> &dN,
           const int a)
{
  B.set_values(0.0);
  B(1, 1) = dN(a, 1);
  B(2, 2) = dN(a, 2);
  B(3, 3) = dN(a, 3);
  B(4, 2) = dN(a, 3);
  B(4, 3) = dN(a, 2);
  B(5, 1) = dN(a, 3);
  B(5, 3) = dN(a, 1);
  B(6, 1) = dN(a, 2);
  B(6, 2) = dN(a, 1);
}

/// Fill the transposed B matrix, size of B have to be created before call this function.
///
/// \param[out] &B  filled B matrix
/// \param[in]  &dN derivative of shape function
/// \param[in]  a   index of node (local id of an element)
void get_BT(Matrix<double> &B,
            Matrix<double> &dN,
            const int a)
{
  B.set_values(0.0);
  B(1, 1) = dN(a, 1);
  B(2, 2) = dN(a, 2);
  B(3, 3) = dN(a, 3);
  B(2, 4) = dN(a, 3);
  B(3, 4) = dN(a, 2);
  B(1, 5) = dN(a, 3);
  B(3, 5) = dN(a, 1);
  B(1, 6) = dN(a, 2);
  B(2, 6) = dN(a, 1);
}

/// Compute element level left and right hand side matrix and vector
/// for updating geometry by given plastic deformation gradient
///
/// \param[in]  *fe         container of finite element resources
/// \param[out] **Lk        computed stiffness matrix
/// \param[in]  *grid       mesh object
/// \param[in]  *fv         object for field variables
/// \param[in]  *sol        object for solution scheme
/// \param[in]  *load       object for loading
/// \param[in]  *com        communication object
/// \param[in]  *Ddof       shifted degree of freedom (global end dof id in each domain)
/// \param[in]  myrank      current process rank
/// \param[in]  interior    if 1, local element
///                            0, element on the communication boundary
/// \param[in]  *opts       structure PGFem3D option
/// \param[in]  mp_id       mutiphysics id
/// \param[in]  do_assemble if yes, assmeble computed element stiffness matrix into
///                                 sparce solver matrix or compute element stiffness only
/// \param[out] *f          computed right hand side vector
/// \param[out] *f_disp     computed right hand side vector due to prescribed values
///                         [Kii Kio]<ui>   <bi>
///                         [Koi Koo]<uo> = <bo>
///                         [Kii][ui] = <bi> - [Kio]<uo>, f_disp = [Kio]<uo>
/// \param[in] &D           6x6 stiffness matrix
/// \param[in] *pF_in       given plastic deformation gradient
/// \return non-zero on internal error
int compute_LHS_RHS_pK_elem(FEMLIB *fe,
                            double **Lk,
                            const Grid *grid,
                            const FieldVariables *fv,
                            const Solver *sol,
                            const LoadingSteps *load,
                            const CommunicationStructure *com,
                            int *Ddof,
                            const int myrank,
                            const int interior,
                            const PGFem3D_opt *opts,
                            const int mp_id,
                            const int do_assemble,
                            double *f,
                            double *f_disp,
                            const Matrix<double> &D,
                            double *pF_in)
{
  int err = 0;
  SUPP sup = load->sups[mp_id];

  long *nod = fe->node_id.m_pdata;
  int ndofn  = fv->ndofn;
  int nne    = fe->nne;
  int nsd    = fe->nsd;
  long ndofe = nne*ndofn;

  int nint   = fe->nint;

  Matrix<long> cnL(ndofe, 1);
  Matrix<long> cnG(ndofe, 1);

  get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,grid->node,cnL.m_pdata,mp_id);
  get_dof_ids_on_elem_nodes(1,nne,ndofn,nod,grid->node,cnG.m_pdata,mp_id);

  Matrix<double> lk(ndofe,ndofe,0.0);
  Matrix<double> fi(ndofe, 1, 0.0);
  Matrix<double> pF;
  Matrix<double> du(nsd, nne, 0.0);
  Matrix<double> u(ndofe, 1, 0.0);
  Matrix<double> x(nsd, nne);
  x.trans(fe->node_coord);

  pF.use_reference(nsd,nsd,pF_in);
  du.prod(pF,x);
  du.sub(x);

  Matrix<double> BaT(3,6,0.0), Bb(6,3,0.0), BaTD(3,6,0.0), BaTDBb(3,3);
  for(int ip = 1; ip<=nint; ip++)
  {
    // Udate basis functions at the integration points.
    fe->elem_basis_V(ip);

    for(int na = 1; na<=nne; na++)
    {
      get_BT(BaT, fe->dN, na);
      BaTD.prod(BaT, D);
      for(int nb = 1; nb<=nne; nb++)
      {
        get_B(Bb, fe->dN, nb);
        BaTDBb.prod(BaTD, Bb);

        for(int ia = 1; ia<=nsd; ia++)
        {
          int p = nsd*(na-1)+ia;
          for(int ja = 1; ja<=nsd; ja++)
          {
            int q = nsd*(nb-1)+ja;
            lk(p,q) += BaTDBb(ia, ja)*fe->detJxW;
          }
        }
      }
    }
  }


  for(int ia=1; ia<=nne; ia++)
  {
    for(int ib=1; ib<=nsd; ib++)
      u((ia-1)*nsd+ib) = du(ib, ia);
  }

  fi.prod(lk, u);

  for(int ia = 0; ia<nne; ia++)
  {
    for(int ib = 0; ib<nsd; ib++)
    {
      int II = grid->node[nod[ia]].id_map[mp_id].id[ib] - 1;
      if (II < 0) continue;
      f[II] += fi.m_pdata[ia*nsd + ib];
    }
  }

  // Assemble
  if(do_assemble)
  {
    PLoc_Sparse(Lk,lk.m_pdata,
                com->Ai,
                com->Ap,
                cnL.m_pdata,cnG.m_pdata,ndofe,Ddof,
                com->GDof,
                myrank,
                com->nproc,
                com->spc,
                interior,
                sol->system,
                opts->analysis_type);
  }


  bool compute_load4pBCs = false;
  for(int ia=0; ia<ndofe; ia++){
    if(cnL.m_pdata[ia] <= -1)
    {
      compute_load4pBCs = true;
      break;
    }
  }

  if(compute_load4pBCs)
  {
    Matrix<double> disp((nne)*ndofn,1,0.0), f_loc((nne)*ndofn,1,0.0);

    // get the bc increment
    for(int ia=0; ia<nne; ia++)
    {
      for(int ib=0; ib<ndofn; ib++)
      {
        int id = ia*ndofn + ib;
        if(cnL.m_pdata[id] <= -1)
          disp.m_pdata[id] = sup->defl[abs(cnL.m_pdata[id])-1];
        else
          disp.m_pdata[id] = 0.0;
      }
    }
    f_loc.prod(lk,disp);

    // element -> localization
    for(int ia=0; ia<nne; ia++)
    {
      for(int ib=0; ib<ndofn; ib++)
      {
        int id_e = ia*ndofn + ib;
        int id_l = grid->node[nod[ia]].id_map[mp_id].id[ib]-1;
        if (id_l < 0)  continue;
        f_disp[id_l] += f_loc.m_pdata[id_e];
      }
    }
  }

  return err;
}

/// Compute left and right hand side matrix and vector
/// for updating geometry by given plastic deformation gradient. Left hand side matrix
/// will be stored soler object (sparce matrix).
///
/// \param[in]       *grid    mesh object
/// \param[in]       *fv      object for field variables
/// \param[in]       *mat     a material object
/// \param[in, out]  *sol     object for solution scheme
/// \param[in]       *load    object for loading
/// \param[in]       *com     communication object
/// \param[in]       *opts    structure PGFem3D option
/// \param[in]       *mp      mutiphysics object
/// \param[in]       mp_id    mutiphysics id
/// \param[in]       myrank   current process rank
/// \param[out]      *f       computed right hand side vector
/// \param[in]       *pF      given plastic deformation gradient
/// \return non-zero on internal error
int compute_LHS_and_RHS_of_pK(const Grid *grid,
                              const FieldVariables *fv,
                              const MaterialProperty *mat,
                              Solver *sol,
                              const LoadingSteps *load,
                              const CommunicationStructure *com,
                              const PGFem3D_opt *opts,
                              const Multiphysics& mp,
                              const int mp_id,
                              const int myrank,
                              double *f)
{
  int err = 0;
  int total_Lagrangian  = 1;
  int intg_order        = 0;
  int do_assemble       = 1; // udate stiffness matrix

  double I[9] = {1.0,0.0,0.0,
                 0.0,1.0,0.0,
                 0.0,0.0,1.0};
                 
  Matrix<double> f_disp(fv->ndofd, 1, 0.0);

  // arbitrary values for elastic deformation, will not effect on results.
  double E = 100.0;
  double nu = 0.3;
  Matrix<double> D(6,6,0.0);

  double lambda = E*nu/(1.0+nu)/(1.0-2.0*nu);
  double mu     = E/2.0/(1.0+nu);
  D(1,1) = D(2,2) = D(3,3) = lambda + 2.0*mu;
  D(1,2) = D(1,3) = D(2,1) = D(2,3) = D(3,1) = D(3,2) = lambda;
  D(4,4) = D(5,5) = D(6,6) = mu;

  double **Lk,**recieve;
  com->spc->post_stiffmat(&Lk,&recieve);
  
  Matrix<int> Ddof(com->nproc,1);

  Ddof.m_pdata[0] = com->DomDof[0];
  for (int ia=1; ia<com->nproc; ia++)
    Ddof.m_pdata[ia] = Ddof.m_pdata[ia-1] + com->DomDof[ia];

  for(int eid=0; eid<com->nbndel; eid++)
  {
    // construct finite element library
    // it provide element wise integration info
    // such as basis function, weights, ...
    
    double pFI[9];
    memcpy(pFI,I,sizeof(double)*9);

    const int mat_id = grid->element[com->bndel[eid]].mat[2];
    double *pF = mat->hommat[mat_id].param->pF;
    if(pF!=NULL)
      inv3x3(pF, pFI);
    
    FEMLIB fe(com->bndel[eid],grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 0;
    err += compute_LHS_RHS_pK_elem(&fe,Lk,grid,fv,sol,load,com,Ddof.m_pdata,
                                   myrank,interior,opts,
                                   mp_id,do_assemble,f,f_disp.m_pdata,D,pFI);

    if(err != 0)
      break;
  }
  com->spc->send_stiffmat();
  
  int skip = 0;
  int idx  = 0;
  
  for(int eid=0; eid<grid->ne; eid++)
  {
    int is_it_in = is_element_interior(eid,&idx,&skip,com->nbndel,com->bndel,myrank);

    if(is_it_in==-1)
    {
      err = 1;
      break;
    }

    if(is_it_in==0)
      continue;
      
      
    double pFI[9];
    memcpy(pFI,I,sizeof(double)*9);

    const int mat_id = grid->element[eid].mat[2];
    double *pF = mat->hommat[mat_id].param->pF;
    if(pF!=NULL)
      inv3x3(pF, pFI);

    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 1;

    err += compute_LHS_RHS_pK_elem(&fe,Lk,grid,fv,sol,load,com,Ddof.m_pdata,
                                   myrank,interior,opts,
                                   mp_id,do_assemble,f,f_disp.m_pdata,D,pFI);

    if(err != 0)
      break;
  }

  com->spc->assemble_nonlocal_stiffmat(sol->system);
  com->spc->finalize_stiffmat();
  sol->system->assemble(); 
  
  for(int ia=0; ia<fv->ndofd; ia++)
    f[ia] -= f_disp.m_pdata[ia];

  return err;
}

/// Compute maximum boundary values due to given deformation gradient.
///
/// \param[out] *bcv_in  computed maximum bc values
/// \param[in]  *grid    mesh object
/// \param[in]  *mat     a material object
/// \param[in]  *fv      object for field variables
/// \param[in]  *load    object for loading
/// \param[in]  *com     communication object
/// \param[in]  mp_id    mutiphysics id
void compute_maximum_BC_values(double *bcv_in,
                               const Grid *grid,
                               const MaterialProperty *mat,
                               const FieldVariables *fv,
                               const LoadingSteps *load,
			       const CommunicationStructure *com,
                               const int mp_id)
{
  SUPP sup = load->sups[mp_id];
  int npd = sup->npd;

  Matrix<double> bcv, max_disp(npd, 1, 1.0e-15), Max_disp(npd, 1), Min_disp(npd, 1);

  bcv.use_reference(npd, 1, bcv_in);

  const int ne = sup->nde;
  const long *el_id = sup->lepd;

  for(int iA=0; iA<ne; iA++)
  {
    // construct finite element library
    // it provide element wise integration info
    // such as basis function, weights, ...
    
    long eid = el_id[iA];
    const int mat_id = grid->element[eid].mat[2];
    double *pF = mat->hommat[mat_id].param->pF;
    if(pF==NULL)
      continue;
      
    double pFI[9];
    inv3x3(pF, pFI);    
    
    FEMLIB fe(eid,grid->element,grid->node,1,1);

    long *nod  = fe.node_id.m_pdata;
    int ndofn  = fv->ndofn;
    int nne    = fe.nne;
    int nsd    = fe.nsd;
    long ndofe = nne*ndofn;
  
    Matrix<long> cnL(ndofe, 1);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,grid->node,cnL.m_pdata,mp_id);
    
    for(int ia=0; ia<nne; ia++){  
      for(int ib=0; ib<nsd; ib++){
        int id = cnL.m_pdata[ia*ndofn + ib];
        if( id >= 0)
          continue;
          
        double x[3], X[3], u[3];
        x[0] = fe.node_coord(ia+1,1); 
        x[1] = fe.node_coord(ia+1,2);
        x[2] = fe.node_coord(ia+1,3);

        X[0] = pFI[0]*x[0] + pFI[1]*x[1] + pFI[2]*x[2];
        X[1] = pFI[3]*x[0] + pFI[4]*x[1] + pFI[5]*x[2];
        X[2] = pFI[6]*x[0] + pFI[7]*x[1] + pFI[8]*x[2];
        
        u[0] = X[0] - x[0];
        u[1] = X[1] - x[1];
        u[2] = X[2] - x[2];
        
        if(fabs(max_disp(abs(id)))<fabs(u[ib]))
          max_disp(abs(id)) = u[ib];
      }
    }
  }

  com->net->allreduce(max_disp.m_pdata,Max_disp.m_pdata,npd,NET_DT_DOUBLE,NET_OP_MAX,com->comm);
  com->net->allreduce(max_disp.m_pdata,Min_disp.m_pdata,npd,NET_DT_DOUBLE,NET_OP_MIN,com->comm);
    
  for(int ib=1; ib<=npd; ib++)
  {
    bcv(ib) = Max_disp(ib);
    if(fabs(Max_disp(ib))<fabs(Min_disp(ib)))
      bcv(ib) = Min_disp(ib);
  }
}

/// Compute initial variance due to given deformation gradient and update geometry to
/// set initial plastic deformation gradient. In setting initial deformation matching
/// plastic deformation gradient, zero elastic deformation will be imposed.
///
/// \param[in, out]  *grid    mesh object
/// \param[in]       *fv      object for field variables
/// \param[in]       *mat     a material object
/// \param[in, out]  *sol     object for solution scheme
/// \param[in, out]  *load    object for loading
/// \param[in]       *com     communication object
/// \param[in]       *opts    structure PGFem3D option
/// \param[in]       *mp      mutiphysics object
/// \param[in]       mp_id    mutiphysics id
/// \param[in]       *pF      given plastic deformation gradient
/// \param[out]      *bcv     prescribed boundary values due to initial plastic deformation
void update_geometry_for_inital_pF(Grid *grid,
                                   const FieldVariables *fv,
                                   const MaterialProperty *mat,
                                   Solver *sol,
                                   LoadingSteps *load,
                                   const CommunicationStructure *com,
                                   const PGFem3D_opt *opts,
                                   const Multiphysics& mp,
                                   const int mp_id,
                                   const double *bcv)
{
  int myrank = com->rank;
  int npd = load->sups[mp_id]->npd;
  for(int ia=0; ia<npd; ia++)
    (load->sups[mp_id])->defl[ia] = bcv[ia];

  Matrix<double> lF(fv->ndofd, 1, 0.0);
  Matrix<double> GF(com->DomDof[myrank], 1, 0.0);
  Matrix<double>  U(com->DomDof[myrank], 1, 0.0);
  SOLVER_INFO s_info;

  sol->system->zero();

  compute_LHS_and_RHS_of_pK(grid,fv,mat,sol,load,com,
			    opts, mp, mp_id, myrank, lF.m_pdata);

  LToG(lF.m_pdata,GF.m_pdata,fv->ndofd,com);

  double hypre_time = sol->system->solveSystem(opts, GF.m_pdata, U.m_pdata, 0, 0,
                                               com->DomDof, &s_info);
  
  if(myrank == 0)
    printf("time = %e: solver: nor = %8.8e || iter = %d\n",hypre_time,s_info.res_norm,s_info.n_iter);

  GToL(U.m_pdata,fv->u_np1,fv->ndofd,com);

  sol->system->zero();

  for(int ic=0; ic<grid->nn; ic++)
  {
    for(int im=0; im<3; im++)
    {
      long id = grid->node[ic].id_map[mp_id].id[im];
      double du = 0.0;

      if(id==0)
        continue;

      if(id<0)
        du = (load->sups[mp_id])->defl[abs(id)-1];
      else
        du = fv->u_np1[id-1];

      if(im==0)
        grid->node[ic].x1 = grid->node[ic].x1_fd = grid->node[ic].x1_fd + du;
      if(im==1)
        grid->node[ic].x2 = grid->node[ic].x2_fd = grid->node[ic].x2_fd + du;
      if(im==2)
        grid->node[ic].x3 = grid->node[ic].x3_fd = grid->node[ic].x3_fd + du;

      if(opts->restart < 0)
      {
        fv->u_n[ic*3+im] = fv->u_nm1[ic*3+im] = -du;

        if(id>0)
          fv->u_np1[id-1] = -du;
      }
    }
  }

  for(int ib=0; ib<npd; ib++)
    (load->sups[mp_id])->defl[ib] = -(load->sups[mp_id])->defl[ib];


}

/// Update initial deformation gradients accordingly after updating geometry.
///
/// \param[in] *grid       mesh object
/// \param[in] *fv         object for field variables
/// \param[in] PGFem3D_opt *opts
void update_element_deformation_gradient(const Grid *grid,
                                         FieldVariables *fv,
                                         const PGFem3D_opt *opts)
{
  int total_Lagrangian = 1;
  int intg_order = 0;

  for(int eid=0; eid<grid->ne; eid++)
  {
    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);

    Matrix<double> d(fe.nne*fe.nsd, 1,0.0);
    for(int ic=0; ic<fe.nne; ic++)
    {
      long nid = fe.node_id(ic+1);
      for(int id=0; id<fe.nsd;id++)
        d.m_pdata[ic*fe.nsd+id] = fv->u_n[nid*fv->ndofn + id];
    }
    
    double v = 0.0;
    double theta = 0.0;
    double one_over_pJ = 0.0;
    Matrix<double> pJ(fe.nint, 1, 0.0);
    for(int ip=0; ip<fe.nint; ip++)
    {
      double Fnp1[9];
      fe.elem_basis_V(ip+1);
      fe.update_shape_tensor();
      fe.update_deformation_gradient(fe.nsd,d.m_pdata,Fnp1);
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      double *tempF = m->param->pF;
      m->param->p_hmat->param->pF = Fnp1;
      m->param->set_init_vals(m);
      m->param->p_hmat->param->pF = tempF;
      
      double J = det3x3(Fnp1);
      v += fe.detJxW;
      theta += J*fe.detJxW;
      one_over_pJ += 1.0/J*fe.detJxW;
      pJ(ip+1) = J;
    }
    theta       /= v;
    one_over_pJ /= v;
    
    if(opts->analysis_type == CM3F)
    {      
      double eJ = theta*one_over_pJ;
      const Model_parameters *mp = fv->eps[eid].model[0].param;
      double P = mp->compute_dudj(fv->eps[eid].model + 0, eJ, -1, 0.5);
      fv->tf.V_np1(eid+1, 1) = fv->tf.V_n(eid+1, 1) = fv->tf.V_nm1(eid+1, 1) = theta;
      fv->tf.P_np1(eid+1, 1) = fv->tf.P_n(eid+1, 1) = fv->tf.P_nm1(eid+1, 1) = P;
    }      
  }
}

/// Compute boundary values and displacement due to given initial plastic deformation gradient(inverse of pF), and
/// set initial values such that when total deformation gradient is computed, initial condition is imposed to have
/// conditions as F = pF, and eF = 1.
///
/// \param[in, out]  *grid    mesh objectn
/// \param[in]       *fv      object for field variables
/// \param[in]       *mat     a material object
/// \param[in]       *sol     object for solution scheme
/// \param[in, out]  *load    object for loading
/// \param[in]       *com     communication object
/// \param[in]       *opts    structure PGFem3D option
/// \param[in]       *mp      mutiphysics object
/// \param[in]       mp_id    mutiphysics id
/// \param[in]       myrank   current process rank
/// \param[in]       *pF      given plastic deformation gradient
/// \return non-zero on internal error
int set_initial_plastic_deformation_gradient_(Grid *grid,
                                             FieldVariables *fv,
                                             MaterialProperty *mat,                                             
                                             pgfem3d::Solver *sol,
                                             LoadingSteps *load,
                                             const CommunicationStructure *com,
                                             const PGFem3D_opt *opts,
                                             const Multiphysics& mp,
                                             const int mp_id)
{
  int err = 0;
  int myrank = com->rank;
  int G_do_intialization = 0;
  int L_do_intialization = 0;
  
  for(int ia=0; ia<mat->nhommat; ia++){
    if(NULL != mat->hommat[ia].param->pF){
      L_do_intialization = 1;
      break;
    }
  }
  
  com->net->allreduce(&L_do_intialization, &G_do_intialization, 1, NET_DT_INT,
		      NET_OP_SUM,com->comm);
  if(G_do_intialization==0)
    return 0;

  int npd = load->sups[mp_id]->npd;
  Matrix<double> bcv(npd, 1, 0.0);
  compute_maximum_BC_values(bcv.m_pdata, grid, mat, fv, load, com, mp_id);

  // print setting 
  if(myrank==0)
  {
    printf("maximum displacement to be set for updating initial geometry: ");
    for(int ia=1; ia<=npd; ia++)
      printf("(%d)=%e ", ia, bcv(ia));

    printf("\n");

    for(int ia=0; ia<mat->nhommat; ia++)
    {
      double *pF = mat->hommat[ia].param->pF;
      if(NULL == pF)
        continue;

      printf("set initial plastic deformation for homat[%d]:\npF=[", ia);
      printf("%e %e %e\n%e %e %e\n%e %e %e]\n", pF[0], pF[1], pF[2]
                                              , pF[3], pF[4], pF[5]
                                              , pF[6], pF[7], pF[8]);
    }
  }

  // Going back to pFI configuration to have initial pF when BC is applied.
  // So, use pFI
  update_geometry_for_inital_pF(grid, fv, mat, sol, load, com,
                                opts, mp, mp_id, bcv.m_pdata);
  
  update_element_deformation_gradient(grid, fv, opts);

  return err;
}

/// Compute boundary values and displacement due to given initial plastic deformation gradient(inverse of pF), and
/// set initial values such that when total deformation gradient is computed, initial condition is imposed to have
/// conditions as F = pF, and eF = 1.
///
/// \param[in, out]  *grid    mesh object
/// \param[in]       *fv      object for field variables
/// \param[in]       *mat     a material object
/// \param[in]       *sol     object for solution scheme
/// \param[in, out]  *load    object for loading
/// \param[in]       *com     communication object
/// \param[in]       mpi_comm MPI_COMM_WORLD
/// \param[in]       *opts    structure PGFem3D option
/// \param[in]       *mp      mutiphysics object
/// \param[in]       mp_id    mutiphysics id
/// \return non-zero on internal error
int set_initial_plastic_deformation_gradient(Grid *grid,
                                             FieldVariables *fv,
                                             MaterialProperty *mat,                                             
                                             pgfem3d::Solver *sol,
                                             LoadingSteps *load,
                                             const CommunicationStructure *com,
                                             const PGFem3D_opt *opts,
                                             const Multiphysics& mp,
                                             const int mp_id){
  int total_Lagrangian = 1;
  int intg_order = 0;

  for(int eid=0; eid<grid->ne; eid++)
  {
    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);

    for(int ip=0; ip<fe.nint; ip++)
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->set_init_vals(m);
    }
    if(opts->analysis_type == CM3F)
    {      
      double eJ = 1.0;
      const Model_parameters *mp = fv->eps[eid].model[0].param;
      double P = mp->compute_dudj(fv->eps[eid].model + 0, eJ, -1, 0.5);
      fv->tf.P_np1(eid+1, 1) = fv->tf.P_n(eid+1, 1) = fv->tf.P_nm1(eid+1, 1) = P;
    }      
  }
  return 0;
}
