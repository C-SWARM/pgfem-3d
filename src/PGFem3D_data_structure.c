#include "PGFem3D_data_structure.h"

#include "node.h"
#include "element.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "sig.h"
#include "eps.h"
#include "matgeom.h"
#include "hommat.h"
#include "hypre_global.h"
#include "pgfem_comm.h"
#include "comm_hints.h"

/// initialize time stepping variable
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
/// 
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int time_stepping_initialization(PGFem3D_TIME_STEPPING *ts)
{
  int err = 0;
  ts->nt     = 0;
  ts->tim    = 0;
  ts->times  = NULL;
  ts->dt_n   = 0.0;
  ts->dt_np1 = 0.0;
  ts->print  = NULL;
  return err;
}


/// destruct time stepping variable
/// free memory spaces for member arrays and structs
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int destruct_time_stepping(PGFem3D_TIME_STEPPING *ts)
{
  int err = 0;
  if(NULL != ts->times) free(ts->times);
  if(NULL != ts->print) free(ts->print);
  err += time_stepping_initialization(ts);  
  return err;
}


/// initialize mesh object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs
///                  except nsd = 3 (number of spatial dimension)
///
/// \param[in, out] grid an object containing all mesh data
/// \return non-zero on internal error
int grid_initialization(GRID *grid)
{
  int err = 0;
  grid->Gnn     = 0;
  grid->Gne     = 0;
  grid->Gnbndel = 0;
  grid->Gn_be   = 0;
  grid->ne      = 0;
  grid->nn      = 0;
  grid->nsd     = 3;
  grid->n_be    = 0;
  grid->nce     = 0;
  grid->Gnce    = 0;
  grid->node    = NULL;
  grid->element = NULL;
  grid->b_elems = NULL;
  grid->coel    = NULL;
  return err;
}


/// destruct of mesh
/// free memory spaces for member arrays and structs
///
/// \param[in, out] grid an object containing all mesh data
/// \return non-zero on internal error
int destruct_grid(GRID *grid, 
                  const PGFem3D_opt *opts)
{
  int err = 0;
  destroy_bounding_elements(grid->n_be,grid->b_elems);
  destroy_elem(grid->element,grid->ne);
  destroy_node(grid->nn,grid->node);
  
  if(opts->cohesive == 1)
    destroy_coel(grid->coel,grid->nce);
    
  err += grid_initialization(grid);
  return err;
}


/// initialize field variables
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
/// 
/// \param[in, out] fv an object containing all field variables
/// \return non-zero on internal error
int field_varialbe_initialization(FIELD_VARIABLES *fv)
{
  int err = 0;
  fv->ndofn  = 0;
  fv->ndofd  = 0;
  fv->npres  = 0;
  fv->n_concentrations = 0;
  fv->u_np1  = NULL;
  fv->u_n    = NULL;
  fv->u_nm1  = NULL;
  fv->d_u    = NULL;
  fv->dd_u   = NULL;
  fv->f      = NULL;
  fv->R      = NULL;
  fv->f_defl = NULL;
  fv->RR     = NULL;
  fv->f_u    = NULL;
  fv->RRn    = NULL;
  fv->pores  = 0.0;
  fv->BS_x   = NULL;
  fv->BS_f   = NULL;
  fv->BS_f_u = NULL;
  fv->BS_RR  = NULL;
  fv->NORM   = 0.0;
  fv->sig    = NULL;
  fv->eps    = NULL;
  return err;
}


/// construct field variables
/// create memory space for member arrays and structs
/// prior to run this function, ndofn, ndofd, and npres must be assigned properly
/// 
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \param[in] com an object for communication
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int construct_field_varialbe(FIELD_VARIABLES *fv, 
                             GRID *grid,
                             COMMUNICATION_STRUCTURE *com,
                             const PGFem3D_opt *opts,
                             int myrank)
{
  int err = 0;
  long DomDof_myrank = com->DomDof[myrank];
  fv->u_np1  = aloc1(fv->ndofd);
  fv->u_n    = aloc1(grid->nn*fv->ndofn);
  fv->u_nm1  = aloc1(grid->nn*fv->ndofn);
  fv->d_u    = aloc1(fv->ndofd);
  fv->dd_u   = aloc1(fv->ndofd);
  fv->f      = aloc1(fv->ndofd);
  fv->R      = aloc1(fv->ndofd);
  fv->f_defl = aloc1(fv->ndofd);
  fv->RR     = aloc1(fv->ndofd);
  fv->f_u    = aloc1(fv->ndofd);
  fv->RRn    = aloc1(fv->ndofd);
  fv->pores  = 0.0;
  fv->BS_x   = aloc1(DomDof_myrank);
  fv->BS_f   = aloc1(DomDof_myrank);
  fv->BS_f_u = aloc1(DomDof_myrank);
  fv->BS_RR  = aloc1(DomDof_myrank);
  fv->NORM   = 0.0;
  fv->sig = build_sig_il(grid->ne,opts->analysis_type,grid->element);
  fv->eps = build_eps_il(grid->ne,grid->element,opts->analysis_type);
  if (opts->smoothing == 0)
    fv->sig_n = build_sig_el(grid->nn);

  return err;
}


/// destruct field variables
/// free memory spaces for member arrays and structs
/// 
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int destruct_field_varialbe(FIELD_VARIABLES *fv, 
                            GRID *grid,
                            const PGFem3D_opt *opts)
{
  int err = 0;
  if(NULL != fv->u_np1)  free(fv->u_np1);
  if(NULL != fv->u_n)    free(fv->u_n);
  if(NULL != fv->u_nm1)  free(fv->u_nm1);
  if(NULL != fv->d_u)    free(fv->d_u);
  if(NULL != fv->dd_u)   free(fv->dd_u);
  if(NULL != fv->f)      free(fv->f);
  if(NULL != fv->R)      free(fv->R);
  if(NULL != fv->f_defl) free(fv->f_defl);
  if(NULL != fv->RR)     free(fv->RR);
  if(NULL != fv->f_u)    free(fv->f_u);
  if(NULL != fv->RRn)    free(fv->RRn);
  if(NULL != fv->BS_x)   free(fv->BS_x);
  if(NULL != fv->BS_f)   free(fv->BS_f);
  if(NULL != fv->BS_f_u) free(fv->BS_f_u);
  if(NULL != fv->BS_RR)  free(fv->BS_RR);
    
  destroy_eps_il(fv->eps,grid->element,grid->ne,opts->analysis_type);
  destroy_sig_il(fv->sig,grid->element,grid->ne,opts->analysis_type);
  if(opts->smoothing == 0)
    destroy_sig_el(fv->sig_n, grid->nn);
    
  err += field_varialbe_initialization(fv);  
  return err;
}


/// initialize material properties
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] mat an object containing all material parameters
/// \return non-zero on internal error
int material_initialization(MATERIAL_PROPERTY *mat)
{
  int err = 0;
  mat->mater    = NULL;
  mat->hommat   = NULL;
  mat->matgeom  = NULL;
  mat->co_props = NULL;      
  mat->nhommat    = 0;
  mat->nmat       = 0;
  mat->n_orient   = 0;
  mat->n_co_props = 0;
  return err;
}


/// destruct material properties
/// free memory spaces for member arrays and structs
///
/// \param[in, out] mat an object containing all material parameters
/// \return non-zero on internal error
int destruct_material(MATERIAL_PROPERTY *mat,
                      const PGFem3D_opt *opts)
{
  int err = 0;
  if(NULL != mat->mater) free(mat->mater);
  destroy_matgeom(mat->matgeom,mat->n_orient);
  destroy_hommat(mat->hommat,mat->nhommat);

  if(opts->cohesive == 1)
    destroy_cohesive_props(mat->n_co_props,mat->co_props);

  err += material_initialization(mat);
  return err;
}


/// initialize iterative solver object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] sol an object containing data for linear solver
/// \return non-zero on internal error
int solution_scheme_initialization(SOLVER_OPTIONS *sol)
{
  int err            = 0;
  sol->n_step        = 0;
  sol->nor_min       = 0.0;
  sol->iter_max      = 0;
  sol->alpha         = 0.5; /// midpoint rule alpha default is 2nd order
  sol->microscale    = NULL;
  sol->PGFEM_hypre   = NULL;
  sol->FNR           = 0;
  sol->gama          = 0.0;
  sol->err           = 0.0;
  sol->iter_max_sol  = 0;
  sol->computer_zero = 0.0;
  return err;
}

/// initialize loading steps object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] load an object containing boundary increments
/// \return non-zero on internal error
int loading_steps_initialization(LOADING_STEPS *load)
{
  int err = 0;
  load->sup         = NULL;
  load->sup_defl    = NULL;
  load->nln         = 0;
  load->nle_s       = 0;
  load->nle_v       = 0;
  load->znod        = NULL;
  load->zele_s      = NULL;
  load->zele_v      = NULL;
  load->tim_load    = NULL;
  load->solver_file = NULL;
  return err;
}


/// destruct loading steps object
/// free memory spaces for member arrays and structs
///
/// \param[in, out] load an object containing boundary increments
/// \return non-zero on internal error
int destruct_loading_steps(LOADING_STEPS *load)
{
  int err = 0;  
  if(NULL != load->sup_defl) free(load->sup_defl);
  if(NULL != load->tim_load) free(load->tim_load);
  if(NULL != load->solver_file) fclose(load->solver_file);
          
  destroy_zatnode(load->znod,   load->nln);
  destroy_zatelem(load->zele_s, load->nle_s);
  destroy_zatelem(load->zele_v, load->nle_v);

  destroy_supp(load->sup);
  
  err += loading_steps_initialization(load);
  return err;
}


/// initialize communication structures
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int communication_structure_initialization(COMMUNICATION_STRUCTURE *com)
{
  int err = 0;
  com->nproc  = 0;
  com->Ap     = NULL;
	com->Ai     = NULL;
  com->DomDof = NULL;
  com->nbndel = 0;
  com->bndel  = NULL;
  com->comm   = NULL;
  com->GDof   = 0;
  com->NBN    = 0;
  com->hints  = NULL;
  return err;
}


/// destruct communication structures
/// free memory spaces for member arrays and structs
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int destruct_communication_structure(COMMUNICATION_STRUCTURE *com)
{
  int err = 0;
  if(NULL != com->Ap) free(com->Ap);
	if(NULL != com->Ai) free(com->Ai);
  if(NULL != com->DomDof) free(com->DomDof);
  if(NULL != com->bndel) free(com->bndel);
  destroy_commun(com->comm ,com->nproc);  
  Comm_hints_destroy(com->hints);
  
  err += communication_structure_initialization(com);
  return err;
}
