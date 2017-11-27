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

using pgfem3d::Solver;


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
                com->comm,
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
/// \param[in, out]  *sol     object for solution scheme
/// \param[in]       *load    object for loading
/// \param[in]       *com     communication object
/// \param[in]       mpi_comm MPI_COMM_WORLD
/// \param[in]       *opts    structure PGFem3D option
/// \param[in]       *mp      mutiphysics object
/// \param[in]       mp_id    mutiphysics id
/// \param[in]       myrank   current process rank
/// \param[out]      *f       computed right hand side vector
/// \param[in]       *pF      given plastic deformation gradient
/// \return non-zero on internal error
int compute_LHS_and_RHS_of_pK(const Grid *grid,
                              const FieldVariables *fv,
                              Solver *sol,
                              const LoadingSteps *load,
                              const CommunicationStructure *com,
                              const MPI_Comm mpi_comm,
                              const PGFem3D_opt *opts,
                              const Multiphysics *mp,
                              const int mp_id,
                              const int myrank,
                              double *f,
                              double *pF)
{
  int err = 0;
  int total_Lagrangian  = 1;
  int intg_order        = 0;
  int do_assemble       = 1; // udate stiffness matrix

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
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;

  err += init_and_post_stiffmat_comm(&Lk,&recieve,&req_r,&sta_r,
                                     mpi_comm,com->comm);

  Matrix<int> Ddof(com->nproc,1);

  Ddof.m_pdata[0] = com->DomDof[0];
  for (int ia=1; ia<com->nproc; ia++)
    Ddof.m_pdata[ia] = Ddof.m_pdata[ia-1] + com->DomDof[ia];

  for(int eid=0; eid<com->nbndel; eid++)
  {
    // construct finite element library
    // it provide element wise integration info
    // such as basis function, weights, ...
    FEMLIB fe(com->bndel[eid],grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 0;
    err += compute_LHS_RHS_pK_elem(&fe,Lk,grid,fv,sol,load,com,Ddof.m_pdata,
                                   myrank,interior,opts,
                                   mp_id,do_assemble,f,f_disp.m_pdata,D,pF);

    if(err != 0)
      break;
  }
  err += send_stiffmat_comm(&sta_s,&req_s,Lk,mpi_comm,com->comm);

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

    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 1;

    err += compute_LHS_RHS_pK_elem(&fe,Lk,grid,fv,sol,load,com,Ddof.m_pdata,
                                   myrank,interior,opts,
                                   mp_id,do_assemble,f,f_disp.m_pdata,D,pF);

    if(err != 0)
      break;
  }

  err += assemble_nonlocal_stiffmat(com->comm,sta_r,req_r,sol->system,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,com->comm);
  sol->system->assemble(); 
  
  for(int ia=0; ia<fv->ndofd; ia++)
    f[ia] -= f_disp.m_pdata[ia]; 

  // stiffnes build is completed
  // deallocate memory
  for(int ia=0; ia<com->nproc; ia++)
  {
    free(recieve[ia]);
    free(Lk[ia]);
  }
  free(recieve);
  free(Lk);
  free(sta_s);
  free(sta_r);
  free(req_s);
  free(req_r);

  return err;
}

/// Compute maximum boundary values due to given deformation gradient.
///
/// \param[out] *bcv_in  computed maximum bc values
/// \param[in]  *pF      given deformation gradient
/// \param[in]  *grid    mesh object
/// \param[in]  *fv      object for field variables
/// \param[in]  *load    object for loading
/// \param[in]  mpi_comm MPI_COMM_WORLD
/// \param[in]  myrank   current process rank
/// \param[in]  mp_id    mutiphysics id
void compute_maximum_BC_values(double *bcv_in,
                               const double *pF,
                               const Grid *grid,
                               const FieldVariables *fv,
                               const LoadingSteps *load,
                               const MPI_Comm mpi_comm,
                               const int mp_id,
                               const int myrank)
{
  int npd = load->sups[mp_id]->npd;
  int ndn = load->sups[mp_id]->ndn;
  Matrix<double> bcv, max_disp(npd, 1, 1.0e-15), Max_disp(npd, 1), Min_disp(npd, 1);
  
  bcv.use_reference(npd, 1, bcv_in);
  
  for(int ib=0; ib<ndn; ib++)
  {
    int nid = load->sups[mp_id]->lnpd[ib];
    double x[3], X[3], u[3];
    x[0] = grid->node[nid].x1_fd;
    x[1] = grid->node[nid].x2_fd;
    x[2] = grid->node[nid].x3_fd;
    
    X[0] = pF[0]*x[0] + pF[1]*x[1] + pF[2]*x[2];
    X[1] = pF[3]*x[0] + pF[4]*x[1] + pF[5]*x[2];
    X[2] = pF[6]*x[0] + pF[7]*x[1] + pF[8]*x[2];
    
    u[0] = X[0] - x[0];
    u[1] = X[1] - x[1];
    u[2] = X[2] - x[2];
  
    for(int ic=0; ic<grid->nsd; ic++)
    {
      long id = grid->node[nid].id_map[mp_id].id[ic];
      if(id<0)
      {  
        if(fabs(max_disp(abs(id)))<fabs(u[ic]))
          max_disp(abs(id)) = u[ic];
      }
    }
  }
  
  MPI_Allreduce(max_disp.m_pdata,Max_disp.m_pdata,npd,MPI_DOUBLE,MPI_MAX,mpi_comm);
  MPI_Allreduce(max_disp.m_pdata,Min_disp.m_pdata,npd,MPI_DOUBLE,MPI_MIN,mpi_comm);
    
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
/// \param[in, out]  *sol     object for solution scheme
/// \param[in, out]  *load    object for loading
/// \param[in]       *com     communication object
/// \param[in]       mpi_comm MPI_COMM_WORLD
/// \param[in]       *opts    structure PGFem3D option
/// \param[in]       *mp      mutiphysics object
/// \param[in]       mp_id    mutiphysics id
/// \param[in]       myrank   current process rank
/// \param[in]       *pF      given plastic deformation gradient
/// \param[out]      *bcv     prescribed boundary values due to initial plastic deformation
void update_geometry_for_inital_pF(Grid *grid,
                                   const FieldVariables *fv,
                                   Solver *sol,
                                   LoadingSteps *load,
                                   const CommunicationStructure *com,
                                   const MPI_Comm mpi_comm,
                                   const PGFem3D_opt *opts,
                                   const Multiphysics *mp,
                                   const int mp_id,
                                   const int myrank,
                                   double *pF,
                                   const double *bcv)
{
  int npd = load->sups[mp_id]->npd;
  for(int ia=0; ia<npd; ia++)
    (load->sups[mp_id])->defl[ia] = bcv[ia];
    
  Matrix<double> lF(fv->ndofd, 1, 0.0);
  Matrix<double> GF(fv->ndofd, 1, 0.0);          
  Matrix<double>  U(fv->ndofd, 1, 0.0);
  SOLVER_INFO s_info;

  sol->system->zero();
        
  compute_LHS_and_RHS_of_pK(grid,fv,sol,load,com,
                            mpi_comm, opts, mp, mp_id, myrank, lF.m_pdata, pF);

  LToG(lF.m_pdata,GF.m_pdata,myrank,com->nproc,
       fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);

  double hypre_time = sol->system->solveSystem(opts, GF.m_pdata, U.m_pdata, 0, 0,
                                               com->DomDof, &s_info);

  if(myrank == 0)
    printf("time = %e: solver: nor = %8.8e || iter = %d\n",hypre_time,s_info.res_norm,s_info.n_iter);


  GToL(U.m_pdata,fv->u_np1,myrank,com->nproc,
       fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);

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
/// \param[in]       *grid    mesh object
/// \param[in]       *fv      object for field variables
void update_element_deformation_gradient(const Grid *grid,
                                         FieldVariables *fv)
{
  int total_Lagrangian = 1; 
  int intg_order = 0;
  
  for(int eid=0; eid<grid->ne; eid++)
  {
    FEMLIB fe;
    fe.initialization(eid,grid->element,grid->node,intg_order,total_Lagrangian);
                
    Matrix<double> d(fe.nne*fe.nsd, 1,0.0);
    for(int ic=0; ic<fe.nne; ic++)
    { 
      long nid = fe.node_id(ic+1);
      for(int id=0; id<fe.nsd;id++)
        d.m_pdata[ic*fe.nsd+id] = fv->u_n[nid*fv->ndofn + id];
    }            
  
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
    }
  }
}

/// Compute boundary values and displacement due to given initial plastic deformation gradient(inverse of pF), and
/// set initial values such that when total deformation gradient is computed, initial condition is imposed to have
/// conditions as F = pF, and eF = 1.
///
/// \param[in, out]  *grid    mesh object
/// \param[in]       *fv      object for field variables
/// \param[in]       *sol     object for solution scheme
/// \param[in, out]  *load    object for loading
/// \param[in]       *com     communication object
/// \param[in]       mpi_comm MPI_COMM_WORLD
/// \param[in]       *opts    structure PGFem3D option
/// \param[in]       *mp      mutiphysics object
/// \param[in]       mp_id    mutiphysics id
/// \param[in]       myrank   current process rank
/// \param[in]       *pF      given plastic deformation gradient
/// \return non-zero on internal error
int set_initial_plastic_deformation_gradient(Grid *grid,
                                             FieldVariables *fv,
                                             pgfem3d::Solver *sol,
                                             LoadingSteps *load,
                                             const CommunicationStructure *com,
                                             const MPI_Comm mpi_comm,
                                             const PGFem3D_opt *opts,
                                             const Multiphysics *mp,
                                             const int mp_id,
                                             const int myrank,
                                             double *pF)
{
  int err = 0;
  
  double pFI[9];
  inv3x3(pF, pFI);  

  int npd = load->sups[mp_id]->npd;
  Matrix<double> bcv(npd, 1, 0.0);
  compute_maximum_BC_values(bcv.m_pdata, pFI, grid, fv, load, mpi_comm, mp_id, myrank);
  
  // print setting 
  if(myrank==0)
  {
    printf("maximum displacement to be set for updating initial geometry: ");
    for(int ia=1; ia<=npd; ia++)
      printf("(%d)=%e ", ia, bcv(ia));

    printf("\n");    
    printf("set initial plastic deformation:\npF=[");
    printf("%e %e %e\n%e %e %e\n%e %e %e]\n", pF[0], pF[1], pF[2]
                                            , pF[3], pF[4], pF[5]
                                            , pF[6], pF[7], pF[8]);
  }

  // Going back to pFI configuration to have initial pF when BC is applied.
  // So, use pFI
  update_geometry_for_inital_pF(grid, fv, sol, load, com, 
                                mpi_comm, opts, mp, mp_id, myrank, pFI, bcv.m_pdata);
                                
  update_element_deformation_gradient(grid, fv);  
  
  return err;
}                              