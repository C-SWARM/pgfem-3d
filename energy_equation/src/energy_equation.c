/// Define energy equation: function for computing stiffness matrix and residual vector
/// 
/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#include "energy_equation.h"
#include "femlib.h"
#include "PGFem3D_data_structure.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "PLoc_Sparse.h"
#include <math.h>
#include "constitutive_model.h"
#include "utils.h"
#include "material_properties.h" // <= constitutive model material properties
#include "hyperelasticity.h"     // <= constitutive model elasticity

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

int compute_deF_over_dhF(Matrix(double) *dF,
                         Matrix(double) *F,
                         Matrix(double) *pFI, 
                         Matrix(double) *hFI)
{
  int err = 0;  
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          for(int M=1; M<=DIM_3; M++)
          {
            for(int N=1; N<=DIM_3; N++)
              Tns4_v(*dF,I,J,K,L) = -Mat_v(*F,I,M)*Mat_v(*hFI,M,K)*Mat_v(*hFI,L,N)*Mat_v(*pFI,N,J);
          }
        }
      }
    }
  }
  return err;
}                         

int compute_xB(Matrix(double) *xB,
               Matrix(double) *eF,
               Matrix(double) *dFdF,
               Matrix(double) *S,
               Matrix(double) *L4)
{
  int err = 0;
  Matrix(double) eB;
  Matrix_construct_redim(double, eB,DIM_3x3x3x3,1);

  for(int I=1; I<=DIM_3; I++)
  {
    for(int K=1; K<=DIM_3; K++)
    {
      for(int J=1; J<=DIM_3; J++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          for(int M=1; M<=DIM_3; M++)
            Tns4_v(eB,I,K,J,L) = Mat_v(*S,L,K)*(I==J) + Mat_v(*eF,I,J)*Tns4_v(*L4,J,K,L,M)*Mat_v(*eF,J,M);
        }
      }
    }
  }  
    
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          for(int M=1; M<=DIM_3; M++)
          {
            for(int N=1; N<=DIM_3; N++)
            {
              for(int Q=1; Q<=DIM_3; Q++)
              {
                for(int R=1; R<=DIM_3; R++)                
                  Tns4_v(*xB,I,J,K,L) = Tns4_v(*dFdF,I,J,M,N)*Tns4_v(eB,N,N,Q,R)*Tns4_v(*dFdF,Q,R,K,L);
              }
            }
          }
        }
      }
    }
  }
  Matrix_cleanup(eB);
  return err;
} 
/// compute heat generation due to plastic flow
///
/// \param[in] mat MATERIAL_PROPERTY object
/// \param[in] fv_m mechanical FIELD_VARIABLES object
/// \param[in] dT temperature change
/// \param[in] dt time step size
/// \param[in] eid element id
/// \param[in] ip integration point id
/// \param[in] mat_id material id
/// \return heat rate
double compute_mechanical_heat_gen_plastic(const MATERIAL_PROPERTY *mat,
                                           const FIELD_VARIABLES *fv_m,
                                           const double dT,
                                           const double dt,
                                           const int eid,
                                           const int ip,
                                           const int mat_id)
{
  int err = 0;
  double Q_m = 0.0;
  int compute_stiffness = 0;

  Matrix(double) F, eF, pF, pFI, pFn, pFdot, S, eP, pP;
  Matrix_construct_init(double, F,    DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, eF,    DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, pF,    DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, pFI,   DIM_3,DIM_3,0.0);  
  Matrix_construct_init(double, pFn,   DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, pFdot, DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, S,     DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, eP,    DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, pP,    DIM_3,DIM_3,0.0);    
  
  Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
  const Model_parameters *func = m->param;

  err += func->get_Fn(   m, &F);  
  err += func->get_eFn(  m, &eF);
  err += func->get_pFn(  m, &pF);
  err += func->get_pFnm1(m, &pFn);
  
  Matrix_inv(pF, pFI);

  for(int ia=0; ia<9; ia++)
    pFdot.m_pdata[ia] = (pF.m_pdata[ia] - pFn.m_pdata[ia])/dt;   
    
  ELASTICITY *elast = (m->param)->cm_elast;

  double *tempS = elast->S; // temporal pointer to update *L, and *S using elast
  double *tempL = elast->L; 
  elast->S = S.m_pdata;
  elast->L = NULL;       
  
  elast->update_elasticity(elast,eF.m_pdata,compute_stiffness);
  Matrix_AxB(eP,1.0,0.0,eF,0,S,0);
  for(int ik = 1; ik<=DIM_3; ik++)
  {
    for(int N = 1; N<=DIM_3; N++)
    {
      for(int ia = 1; ia<=DIM_3; ia++)
      {
        for(int ja = 1; ja<=DIM_3; ja++)
        {
          for(int M = 1; M<=DIM_3; M++)
            Mat_v(pP,ik,N) += -Mat_v(eP,ia,ja)*Mat_v(F,ia,M)*Mat_v(pFI,M,ik)*Mat_v(pFI,N,ja);
        }
      }
    }
  }  

  Matrix_ddot(pP, pFdot, Q_m);
  
  elast->S = tempS;
  elast->L = tempL;        

  Matrix_cleanup(F);
  Matrix_cleanup(eF);  
  Matrix_cleanup(pF);
  Matrix_cleanup(pFI);  
  Matrix_cleanup(pFn);
  Matrix_cleanup(pFdot);
  Matrix_cleanup(S);
  Matrix_cleanup(eP);
  Matrix_cleanup(pP);

  return -Q_m;
}

/// compute heat generation due to elastic flow
///
/// \param[in] mat MATERIAL_PROPERTY object
/// \param[in] fv_m mechanical FIELD_VARIABLES object
/// \param[in] dT temperature change
/// \param[in] dt time step size
/// \param[in] eid element id
/// \param[in] ip integration point id
/// \param[in] mat_id material id
/// \return heat rate
double compute_mechanical_heat_gen_elastic(const MATERIAL_PROPERTY *mat,
                                           const FIELD_VARIABLES *fv_m,
                                           const double T,
                                           const double dT,
                                           const double dt,
                                           const int eid,
                                           const int ip,
                                           const int mat_id)
{
  int err = 0;
  double Q_m = 0.0;
  int compute_stiffness = 1;
  
  double Tdot = dT/dt;

  Matrix(double) L, deF_dhF, hB;
  Matrix_construct_init(double, L,       DIM_3x3x3x3, 1, 0.0);
  Matrix_construct_init(double, deF_dhF, DIM_3x3x3x3, 1, 0.0);
  Matrix_construct_init(double, hB,      DIM_3x3x3x3, 1, 0.0);    
      
  enum {F,eF,pF,pFI,hF,hFI,hFp,hFp_d,hFpp,S,eP,hP,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3 ,0.0);
  }  

  // get mechanical part of deformation gradient
  Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
  const Model_parameters *func = m->param;
  ELASTICITY *elast = (m->param)->cm_elast;
  err += func->get_Fn(   m, F2+F);  
  err += func->get_eFn(  m, F2+eF);
  err += func->get_pFn(  m, F2+pF);

  // compute thermal part of deformation gradient    
  double ax = mat->mater[mat_id].ax;
  double ay = mat->mater[mat_id].ay;
  double az = mat->mater[mat_id].az;  

  Mat_v(F2[hF],  1,1) = 1.0 + ax*dT;  
  Mat_v(F2[hF],  2,2) = 1.0 + ay*dT;
  Mat_v(F2[hF],  3,3) = 1.0 + az*dT;
  
  Mat_v(F2[hFp], 1,1) = ax;
  Mat_v(F2[hFp], 2,2) = ay;
  Mat_v(F2[hFp], 3,3) = az;
  
  Mat_v(F2[hFpp], 1,1) = Mat_v(F2[hFpp], 2,2) = Mat_v(F2[hFpp], 3,3) = 0.0;

  Matrix_AeqB(F2[hFp_d], Tdot, F2[hFpp]);
  Matrix_inv(F2[hF], F2[hFI]);
  Matrix_inv(F2[pF], F2[pFI]);  


  double *tempS = elast->S; // temporal pointer to update *L, and *S using elast
  double *tempL = elast->L; 
  elast->S = F2[S].m_pdata;
  elast->L = L.m_pdata;
  
  elast->update_elasticity(elast,F2[eF].m_pdata,compute_stiffness);
  Matrix_AxB(F2[eP],1.0,0.0,F2[eF],0,F2[S],0);  
  
  err += compute_deF_over_dhF(&deF_dhF,F2+F,F2+pFI,F2+hFI);  
  err += compute_xB(&hB,F2+eF,&deF_dhF,F2+S,&L);

  
  for(int I = 1; I<=DIM_3; I++)
  {
    for(int J = 1; J<=DIM_3; J++)
    {
      Q_m += Mat_v(F2[hP],I,J)*Mat_v(F2[hFpp],I,J);
      for(int K = 1; K<=DIM_3; K++)
      {
        for(int L = 1; L<=DIM_3; L++)
          Q_m += Mat_v(F2[hFp],I,J)*Tns4_v(hB,I,J,K,L)*Mat_v(F2[hFp],I,J);
      }
    }
  } 

  Q_m *= T*Tdot;
  
  elast->S = tempS;
  elast->L = tempL;        
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);   
  free(F2);
  
  Matrix_cleanup(L);
  Matrix_cleanup(deF_dhF);
  Matrix_cleanup(hB);  

  return Q_m;
}

/// compute temperature on node in a element
///
/// \param[in] cn local DOF ids for field variables
/// \param[in] ndofe number of degree of freedom for element
/// \param[in] T nodal temperature for all node
/// \param[in] dT temperature increment
/// \param[in] elem ELEMENT object
/// \param[in] node NODE object
/// \param[in] sup struct for BC's
/// \param[in] T_e computed nodal temperature for current element
/// \param[in] T0 reference temperature
/// \return non-zero on internal error
int get_temperature_elem(const long *cn,
                         const long ndofe,
                         const double *T,
                         const double *dT,
                         const ELEMENT *elem,
                         const NODE *node,
                         const SUPP sup,
                         double *T_e,
                         double T0)
{
  int err = 0;
  for(int i=0; i< ndofe; i++){
    const int id = cn[i];
    const int aid = abs(id) - 1;

    if (id == 0){
      T_e[i] = T0;
    } else if (id > 0){
      T_e[i] = T[aid] + dT[aid];
    } else {
      T_e[i] = T0 + sup->defl[aid] + sup->defl_d[aid];
    }
  }
  return err;
}

/// assemble residuals for heat conduction problem
/// 
/// element-wize residuals are merged into global residual vector
///
/// \param[in] fe container of finite element resources
/// \param[in] fi local residual(element-wize)
/// \param[in] grid a mesh object
/// \param[in,out] fv field variable object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int energy_equation_residuals_assemble(FEMLIB *fe,
                                       double *fi,
                                       GRID *grid,
                                       FIELD_VARIABLES *fv,
                                       const int mp_id)
{
  int err = 0;
  long *nod = fe->node_id.m_pdata;
   
  for(int ia = 0; ia<fe->nne; ia++) 
  {
    for(int ib = 0; ib<fv->ndofn; ib++)
    {
      int II = grid->node[nod[ia]].id_map[mp_id].id[ib] - 1;
      if (II < 0) continue;
      fv->f_u[II] += fi[ia*(fv->ndofn) + ib];
    }
  }
  return err;
}

/// compute residuals for heat conduction problem at the element level
///
/// Actual computation for constructing residual vector takes place here.
/// Memory of the residual vector should be allocated before this function calls.
///
/// \param[in] fe container of finite element resources
/// \param[out] fi_in local residual vector (element level) 
/// \param[in] du temperature increment
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] FV array of field variable object 
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] dt_in time step size
/// \return non-zero on internal error
int energy_equation_compute_residuals_elem(FEMLIB *fe,
                                           double *fi_in,
                                           double *du,
                                           GRID *grid,
                                           MATERIAL_PROPERTY *mat,
                                           FIELD_VARIABLES *fv,
                                           LOADING_STEPS *load,
                                           int mp_id,
                                           double dt_in)
                                           
{
  int err = 0;
  
  int eid = fe->curt_elem_id;  
  const int hmat_id = (grid->element[eid]).mat[2];  
  double rho_0 = (mat->hommat[hmat_id]).density;
  const int mat_id = (mat->hommat[hmat_id]).mat_id;
  
  MATERIAL_THERMAL *thermal = (mat->thermal) + mat_id;  
  Matrix(double) k;
  k.m_pdata = thermal->k;
  k.m_row = k.m_col = DIM_3;
  
  double cp = thermal->cp;
  double Q = 0.0;
  
  double dt = 1.0; // for the quasi steady state
  if(rho_0>0)
    dt = dt_in;  
  
  SUPP sup = load->sups[mp_id];

  long *nod = fe->node_id.m_pdata; // use only address, no need deallication
  int ndofn = fv->ndofn;  
  long ndofe = (fe->nne)*ndofn;
  long *cnL = aloc1l(ndofe);
  
  Matrix(double) fi;
  fi.m_pdata = fi_in;
  fi.m_row = ndofe;  
  
  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cnL,mp_id); 
  
  Matrix(double) q, Tnp1, Tn;  
  Matrix_construct_redim(double, q,grid->nsd,1);
  Matrix_construct_init(double, Tnp1,fe->nne,1,0.0); 
  Matrix_construct_init(double, Tn,  fe->nne,1,0.0);
  
  //compute nodal value
  get_temperature_elem(cnL,ndofe,fv->u_np1,du,grid->element,grid->node,sup,Tnp1.m_pdata,fv->u0);
  
  int myrank = 0;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpi_comm,&myrank);

  for(int ia=0; ia<fe->nne; ia++)
    Vec_v(Tn,   ia+1) = fv->u_n[nod[ia]];
  
  for(int ip = 1; ip<=fe->nint; ip++)
  {
    // Udate basis functions at the integration points.
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);
    
    double Temp = 0.0;
    double dT   = 0.0;
    Matrix_init( q,0.0);
            
    // compute varialbes at the integration point
    for(int ia=1; ia<=fe->nne; ia++)
    { 
      // k = [nsd, nsd], dN = [nne, nsd], Tnp1  = [nne, 1],
      // q = [nsd, 1] = k*dN'*T
     
      for(int ib=1; ib<=grid->nsd; ib++)
      {
        for(int ic=1; ic<grid->nsd; ic++)        
          Vec_v(q,ib) += Mat_v(k, ib, ic)*Mat_v(fe->dN,ia,ic)*Vec_v(Tnp1, ia);
      }
        
      Temp += Vec_v(fe->N,ia)*Vec_v(Tnp1, ia);
      dT   += Vec_v(fe->N,ia)*(Vec_v(Tnp1, ia)-Vec_v(Tn, ia));
    }
    
    if(fv->n_coupled > 0)
    {  
      FIELD_VARIABLES *fv_m = fv->fvs[0];
      switch(fv->coupled_physics_ids[0])
      {
        case MULTIPHYSICS_MECHANICAL:
        {
          Q += compute_mechanical_heat_gen_elastic(mat,fv_m,Temp,dT,dt,eid,ip,mat_id);
          Q += compute_mechanical_heat_gen_plastic(mat,fv_m,dT,dt,eid,ip,mat_id);
          break;
        }  
        
        default:
          Q += 0.0;
      }    
    }
      
    // R = rho_0*cp*dT + dt*grad.q - dt*Q = 0;
    for(int ia=1; ia<=fe->nne; ia++)
    {      
      Vec_v(fi,ia) += rho_0*cp*Vec_v(fe->N,ia)*(dT - dt*Q)*(fe->detJxW);
      for(int ib=1; ib<=grid->nsd; ib++)
        Vec_v(fi,ia) += dt*Mat_v(fe->dN,ia,ib)*Vec_v(q,ib)*(fe->detJxW);
    }  
  }
  
  free(cnL);
  Matrix_cleanup(q);
  Matrix_cleanup(Tnp1);
  Matrix_cleanup(Tn);
  
  return err;
}

/// compute residuals for heat conduction problem
///
/// Every element will be visited to compute residuals and local residuals are
/// assembled to global residual vector. The actual computation for constructing
/// residuals takes place in energy_equation_compute_residuals_elem function.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv field variable object 
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] use_updated if use_updated=1, compute residuals updated temperature
///                           use_updated=0, compute residuals using temporal temperature
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_residuals(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      LOADING_STEPS *load,
                                      const int mp_id,
                                      int use_updated,
                                      double dt)
{
  int err = 0;
  int total_Lagrangian = 0;
  int intg_order = 0;
  
  double *du;
  if(use_updated)
    du = fv->f;
  else
    du = fv->d_u;  
    
  for(int eid=0; eid<grid->ne; eid++)
  {
    // Construct finite element library.
    // It provide element wise integration info 
    // such as basis function, weights, ...
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);
    
    // do volume integration at an element, 
    // fe needs to be updated by integration points 
    Matrix(double) fi;
    Matrix_construct_init(double, fi, fe.nne, 1, 0.0);
    err += energy_equation_compute_residuals_elem(&fe,fi.m_pdata,du,grid,mat,fv,load,mp_id,dt);
    err += energy_equation_residuals_assemble(&fe,fi.m_pdata,grid,fv,mp_id);
    
    Matrix_cleanup(fi);
    FEMLIB_destruct(&fe);
  }
  return err;
}

/// compute stiffness for heat conduction problem at the element level
///
/// Actual computation for constructing stiffness matrix takes place here.
/// Local memory for the element level stiffness matrix is created and 
/// merged into memory for assembling and communications.
///
/// \param[in] fe container of finite element resources
/// \param[out] LK local stiffness matrix for assmebling and communication 
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] fv field variable object 
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com object for communications
/// \param[in] Ddof number of degree of freedoms accumulated through processes
/// \param[in] myrank current process rank
/// \param[in] interior indentifier to distinguish element in interior or
///             communication boundary.
/// \param[in] mp_id mutiphysics id
/// \param[in] do_assemble if yes, assmeble computed element stiffness matrix into 
///            sparce solver matrix or compute element stiffness only
/// \param[in] compute_load4pBCs, if yes, compute external flux due to Dirichlet BCs
///                               if no, no compute external flux due to Dirichlet BCs
/// \param[in] dt_in time step size
/// \return non-zero on internal error
int energy_equation_compute_stiffness_elem(FEMLIB *fe,
                                           double **Lk,                                           
                                           GRID *grid,
                                           MATERIAL_PROPERTY *mat,
                                           FIELD_VARIABLES *fv,
                                           SOLVER_OPTIONS *sol,
                                           LOADING_STEPS *load,
                                           COMMUNICATION_STRUCTURE *com,
                                           int *Ddof,
                                           int myrank,
                                           int interior,
                                           const PGFem3D_opt *opts,
                                           int mp_id,
                                           int do_assemble,
                                           int compute_load4pBCs,
                                           double dt_in)
                                           
{  
  int err = 0;
  int eid = fe->curt_elem_id;  
  const int hmat_id = (grid->element[eid]).mat[2];  
  double rho_0 = (mat->hommat[hmat_id]).density;
  const int mat_id = (mat->hommat[hmat_id]).mat_id;
  
  MATERIAL_THERMAL *thermal = (mat->thermal) + mat_id;  
  Matrix(double) k;
  k.m_pdata = thermal->k;
  k.m_row = k.m_col = DIM_3;
  
  double cp = thermal->cp;
  double Q = 0.0;
  
  double dt = 1.0; // for the quasi steady state
  if(rho_0>0)
    dt = dt_in;    
  
  long *nod = fe->node_id.m_pdata;
  int ndofn = fv->ndofn;  
  long ndofe = (fe->nne)*ndofn;
  long *cnL = aloc1l(ndofe);
  long *cnG = aloc1l(ndofe);
  
  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cnL,mp_id);
  get_dof_ids_on_elem_nodes(1,fe->nne,ndofn,nod,grid->node,cnG,mp_id);
  
  Matrix(double) lk;
  Matrix_construct_init(double,lk,ndofe,ndofe,0.0);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    // Udate basis functions at the integration points.
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);
    for(int ia=1; ia<=fe->nne; ia++)
    {
      for(int ib=1; ib<=fe->nne; ib++)
      {
        Mat_v(lk,ia,ib) += rho_0*cp*Vec_v(fe->N,ia)*Vec_v(fe->N,ib)*(fe->detJxW);
        for(int im = 1; im<=grid->nsd; im++)
        {
          for(int in = 1; in<=grid->nsd; in++)
            Mat_v(lk,ia,ib) += dt*Mat_v(fe->dN,ia,in)*Mat_v(k, in, im)*Mat_v(fe->dN,ib,im)*(fe->detJxW);
        }    
      }
    }  
  } 

  // Assemble
  if(do_assemble)
  {  
    PLoc_Sparse(Lk,lk.m_pdata,
                com->Ai,
                com->Ap,
                cnL,cnG,ndofe,Ddof,
                com->GDof,
                myrank,
                com->nproc,
                com->comm,
                interior,
                sol->PGFEM_hypre,
                opts->analysis_type);                               
  }
  
  if(compute_load4pBCs)
  {
    int ndofn = 1;
    int k = 0;
    int jj = 0;
    Matrix(double) u, f_loc;    
    Matrix_construct_init(double,u    ,(fe->nne)*ndofn,1,0.0);
    Matrix_construct_init(double,f_loc,(fe->nne)*ndofn,1,0.0);
    
    // get the bc increment
    for(int ia=0; ia<fe->nne; ia++)
    {
      for(int ib=0; ib<ndofn; ib++)
      {
        int id = ia*ndofn + ib;
        if(cnL[id] <= -1)
          u.m_pdata[id] = load->sups[mp_id]->defl_d[abs(cnL[id])-1];
        else
          u.m_pdata[id] = 0.0;
      }
    }
    Matrix_AxB(f_loc,1.0,0.0,lk,0,u,0);
    
    // element -> localization
    for(int ia=0; ia<fe->nne; ia++)
    {
      for(int ib=0; ib<ndofn; ib++)
      {
        int id_e = ia*ndofn + ib;
        int id_l = grid->node[nod[ia]].id_map[mp_id].id[ib]-1;
        if (id_l < 0)  continue;
          fv->f_defl[id_l] += f_loc.m_pdata[id_e];
      }
    }
    
    Matrix_cleanup(u);
    Matrix_cleanup(f_loc);        
  }
  
  Matrix_cleanup(lk);
    
  free(cnL);
  free(cnG);
  return err;
}

/// compute stiffness for heat conduction problem
///
/// In building global stiffness matrix, elements on communcation boundary are
/// first computed and interior elements are visited later. As soon as stiffness is 
/// computed on the communcation boundary, communication is established and the interior stiffness
/// is commputed. In this way, the communcation and the computation are overlaid.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] fv field variable object 
/// \param[in] sol object for solution scheme
/// \param[in] com object for communications
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_stiffness(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      COMMUNICATION_STRUCTURE *com,
                                      MPI_Comm mpi_comm,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt)
{
  int err = 0;
  int total_Lagrangian  = 0;
  int intg_order        = 0;
  int do_assemble       = 1; // udate stiffness matrix
  int compute_load4pBCs = 0; // if 1, compute load due to Dirichlet BCs. 
                             // if 0, update only stiffness

  double **Lk,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;  

  err += init_and_post_stiffmat_comm(&Lk,&recieve,&req_r,&sta_r,
                                     mpi_comm,com->comm);  

  Matrix(int) Ddof;
  Matrix_construct_redim(int, Ddof,com->nproc,1);
  
  Ddof.m_pdata[0] = com->DomDof[0];
  for (int ia=1; ia<com->nproc; ia++)
    Ddof.m_pdata[ia] = Ddof.m_pdata[ia-1] + com->DomDof[ia];
  
  for(int eid=0; eid<com->nbndel; eid++)
  {
    // construct finite element library
    // it provide element wise integration info 
    // such as basis function, weights, ... 
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,com->bndel[eid],grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 0;
    err += energy_equation_compute_stiffness_elem(&fe,Lk,grid,mat,fv,sol,NULL,com,Ddof.m_pdata,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt);
    
    FEMLIB_destruct(&fe);
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
    
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);
    
    // do volume integration at an element
    int interior = 1;
    
    err += energy_equation_compute_stiffness_elem(&fe,Lk,grid,mat,fv,sol,NULL,com,Ddof.m_pdata,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt);
      
    FEMLIB_destruct(&fe);
    if(err != 0)
      break;
  }        

  err += assemble_nonlocal_stiffmat(com->comm,sta_r,req_r,sol->PGFEM_hypre,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,com->comm);
    
  // stiffnes build is completed
  // deallocate memory
  for(int ia=0; ia<com->nproc; ia++)
  {
    free (recieve[ia]);
    free(Lk[ia]);
  }
  free (recieve);
  free (Lk);
  free (sta_s);
  free (sta_r);
  free (req_s);
  free (req_r);  
  
  Matrix_cleanup(Ddof);
  return err;
}

/// compute flux due to Dirichlet BCs
///
/// Compute flux vector for prescribed BCs(Dirichlet)
/// This compute load, f, as below:
/// [Kii Kio]<ui>   <bi>
/// [Koi Koo]<uo> = <bo>
/// [Kii][ui] = <bi> - [Kio]<uo>
/// where f = [Kio]<uo>, uo is Drichlet BCs  
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] fv field variable object 
/// \param[in] sol object for solution scheme
/// \param[in] com object for communications
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_load4pBCs(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      LOADING_STEPS *load,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt)
{
  int err = 0;
  int total_Lagrangian = 0;
  int intg_order       = 0;
  int interior         = 1;
  int do_assemble      = 0;
  int compute_load4pBCs= 1;  

  for(int ia=0; ia<load->sups[mp_id]->nde; ia++)
  {
    int eid = load->sups[mp_id]->lepd[ia];
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    err += energy_equation_compute_stiffness_elem(&fe,NULL,grid,mat,fv,sol,load,NULL,NULL,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt);    
    FEMLIB_destruct(&fe);
    if(err != 0)
      break;      
  }
    
  return err;
}

