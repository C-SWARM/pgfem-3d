#include "energy_equation.h"
#include "femlib.h"
#include "PGFem3D_data_structure.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "PLoc_Sparse.h"
#include <math.h>

int is_element_interior(int eid, int *idx, int *skip, COMMUNICATION_STRUCTURE *com,
                        int myrank)
{ 
  int is_it_in = 1;
  if(com->nbndel > 0) // most of time it is ture
  {
    if(*idx < com->nbndel-1)
    {
      if(eid == 0 && *idx == 0 && com->bndel[*idx] == 0)
      {
        (*idx)++;
        (*skip)++;
        is_it_in = 0;
      } 
      else if(eid == com->bndel[*idx])
      {
        (*idx)++;
        (*skip)++;
        is_it_in = 0;
      } 
      else if (*idx == 0 && eid < com->bndel[*idx])
        is_it_in = 1;
      else if (*idx > 0 && com->bndel[*idx-1] < eid && eid < com->bndel[*idx])
        is_it_in = 1;
      else 
      {
        is_it_in = -1;
        PGFEM_printf("[%d]ERROR: problem in determining if element %ld"
                     " is on interior.\n", myrank, eid);
      }
    } 
    else if(eid == com->bndel[com->nbndel-1])
      is_it_in = 0;
  }
  
  return is_it_in;    
} 

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

int energy_equation_compute_residuals_elem(FEMLIB *fe,
                                           double *fi_in,
                                           double *du,
                                           GRID *grid,
                                           MATERIAL_PROPERTY *mat,
                                           FIELD_VARIABLES *fv,
                                           LOADING_STEPS *load,
                                           int mp_id)
                                           
{
  int err = 0;
  double Q = 0.0;
  double k = 1.0;
  
  SUPP sup = load->sups[mp_id];

  long *nod = fe->node_id.m_pdata; // use only address, no need deallication
  int ndofn = fv->ndofn;  
  long ndofe = (fe->nne)*ndofn;
  long *cnL = aloc1l(ndofe);
  
  Matrix(double) fi;
  fi.m_pdata = fi_in;
  fi.m_row = ndofe;  
  
  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cnL,mp_id); 
  
  Matrix(double) q, dT, T;  
  Matrix_construct_redim(double, q,grid->nsd,1);
  Matrix_construct_redim(double,dT,grid->nsd,1);
  Matrix_construct_init(double, T,fe->nne,1,0.0); 
  
  //compute nodal value
  def_elem_total(cnL,ndofe,fv->u_np1,du,grid->element,grid->node,sup,T.m_pdata);  
  
  for(int ip = 1; ip<=fe->nint; ip++)
  {
    // Udate basis functions at the integration points.
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);
    
    double Temp = 0.0;
    Matrix_init( q,0.0);
    Matrix_init(dT,0.0);
            
    // compute varialbes at the integration point
    for(int ia=1; ia<=fe->nne; ia++)
    {       
      for(int ib=1; ib<=grid->nsd; ib++)
        Vec_v(q,ib) += k*Mat_v(fe->dN,ia,ib)*Vec_v(T, ia);
        
      Temp += Vec_v(fe->N,ia)*Vec_v(T, ia);
    }
    
    for(int ia=1; ia<=fe->nne; ia++)
    {
      Vec_v(fi,ia) += Vec_v(fe->N,ia)*Q*(fe->detJxW);
      for(int ib=1; ib<=grid->nsd; ib++)
        Vec_v(fi,ia) += Mat_v(fe->dN,ia,ib)*Vec_v(q,ib)*(fe->detJxW);
    }  
  }
  
  free(cnL);
  Matrix_cleanup(q);
  Matrix_cleanup(dT);
  Matrix_cleanup(T);    
  
  return err;
}

int energy_equation_compute_residuals(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      LOADING_STEPS *load,
                                      const int mp_id,
                                      int use_updated)
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
    err += energy_equation_compute_residuals_elem(&fe,fi.m_pdata,du,grid,mat,fv,load,mp_id);
    err += energy_equation_residuals_assemble(&fe,fi.m_pdata,grid,fv,mp_id);
    
    Matrix_cleanup(fi);
    FEMLIB_destruct(&fe);
  }
  return err;
}

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
                                           int compute_load4pBCs)
                                           
{
  int err = 0;
  double k = 1.0; //temporal conductivity
  
  const int mat_id = grid->element[fe->curt_elem_id].mat[2];
  double rho = mat->hommat[mat_id].density;  
  
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
        Mat_v(lk,ia,ib) += Vec_v(fe->N,ia)*Vec_v(fe->N,ib)*(fe->detJxW);
        for(int im = 1; im<=grid->nsd; im++)
          Mat_v(lk,ia,ib) += k*Mat_v(fe->dN,ia,im)*Mat_v(fe->dN,ib,im)*(fe->detJxW);
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

///
///
///
/// \param[in] print_level print level for a summary of the entire function call
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int energy_equation_compute_stiffness(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      COMMUNICATION_STRUCTURE *com,
                                      MPI_Comm mpi_comm,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id)
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
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 0;
    err += energy_equation_compute_stiffness_elem(&fe,Lk,grid,mat,fv,sol,NULL,com,Ddof.m_pdata,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs);
    
    FEMLIB_destruct(&fe);
    if(err != 0)
      break;      
  }
  err += send_stiffmat_comm(&sta_s,&req_s,Lk,mpi_comm,com->comm);
    
  int skip = 0;
  int idx  = 0;

  for(int eid=0; eid<grid->ne; eid++)
  {
    int is_it_in = is_element_interior(eid,&idx,&skip,com,myrank);
    
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
                                                  do_assemble,compute_load4pBCs);
      
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

int energy_equation_compute_load4pBCs(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      LOADING_STEPS *load,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id)
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
                                                  do_assemble,compute_load4pBCs);    
    FEMLIB_destruct(&fe);
    if(err != 0)
      break;      
  }
    
  return err;
}

