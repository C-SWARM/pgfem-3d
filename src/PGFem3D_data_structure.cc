#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFem3D_data_structure.h"
#include "allocation.h"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "constitutive_model.h"
#include "elem3d.h"
#include "element.h"
#include "eps.h"
#include "femlib.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "sig.h"
#include "utils.h"
#include "vtk_output.h"
#include <cassert>

#include <sstream>
#include <fstream>

using namespace pgfem3d;

/// initialize time stepping variable
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int time_stepping_initialization(TimeStepping *ts)
{
  int err = 0;
  ts->nt     = 0;
  ts->tim    = 0;
  ts->times  = NULL;
  ts->dt_n   = 0.0;
  ts->dt_np1 = 0.0;
  ts->print  = NULL;
  ts->tns    = NULL;
  return err;
}


/// destruct time stepping variable
/// free memory spaces for member arrays and structs
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int destruct_time_stepping(TimeStepping *ts)
{
  int err = 0;
  if(NULL != ts->times) PGFEM_free(ts->times);
  if(NULL != ts->print) PGFEM_free(ts->print);
  if(NULL != ts->tns)   PGFEM_free(ts->tns);
  err += time_stepping_initialization(ts);
  return err;
}


/// initialize mesh object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs
///                  except nsd = 3 (number of spatial dimension)
///
/// \param[in, out] grid an object containing all mesh data
/// \return non-zero on internal error
int grid_initialization(Grid *grid)
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
/// \param[in] mp multiphysics object
/// \return non-zero on internal error
int destruct_grid(Grid *grid,
                  const PGFem3D_opt *opts,
                  const Multiphysics& mp)
{
  int err = 0;
  destroy_bounding_elements(grid->n_be,grid->b_elems);
  destroy_elem(grid->element,grid->ne);
  destroy_node_multi_physics(grid->nn,grid->node, mp.physicsno);

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
int field_varialbe_initialization(FieldVariables *fv)
{
  int err = 0;
  fv->apply_initial_velocity = false;
  fv->u0     = 0.0;
  fv->Gndof  = 0;
  fv->ndofn  = 0;
  fv->ndofd  = 0;
  fv->npres  = 0;
  fv->nVol   = 0;
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
  fv->sig_n  = NULL;
  fv->n_coupled = 0;
  fv->coupled_physics_ids = NULL;
  fv->fvs    = NULL;
  fv->temporal = NULL;
  fv->statv_list = NULL;
  fv->subdivision_factor_n   = 0.0;
  fv->subdivision_factor_np1 = 1.0;
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
/// \param[in] mp mutiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id physics id
/// \return non-zero on internal error
int construct_field_varialbe(FieldVariables *fv,
                             Grid *grid,
                             CommunicationStructure *com,
                             const PGFem3D_opt *opts,
                             const Multiphysics& mp,
                             int myrank,
                             int mp_id)
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
  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
  {
    if(opts->analysis_type == CM || opts->analysis_type == CM3F)
    {
      const Element *elem = grid->element;
      int n_state_varialbles = 0;
      for(long eid=0; eid<grid->ne; eid++)
      {
        long nne = elem[eid].toe;
        long nint = 0;
        int_point(nne,&nint);
        n_state_varialbles += nint;
      }

      fv->statv_list = new State_variables[n_state_varialbles];
    }

    fv->sig = build_sig_il(grid->ne,opts->analysis_type,grid->element);
    fv->eps = build_eps_il(grid->ne,grid->element,opts->analysis_type,&(fv->statv_list));
    if (opts->smoothing == 0)
      fv->sig_n = build_sig_el(grid->nn);
  }
  else
    fv->eps = build_eps_il(grid->ne,grid->element,-1,NULL);

  return err;
}


/// destruct field variables
/// free memory spaces for member arrays and structs
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id physics id
/// \return non-zero on internal error
int destruct_field_varialbe(FieldVariables *fv,
                            Grid *grid,
                            const PGFem3D_opt *opts,
                            const Multiphysics& mp,
                            int mp_id)
{
  int err = 0;
  if(NULL != fv->u_np1)  PGFEM_free(fv->u_np1);
  if(NULL != fv->u_n)    PGFEM_free(fv->u_n);
  if(NULL != fv->u_nm1)  PGFEM_free(fv->u_nm1);
  if(NULL != fv->d_u)    PGFEM_free(fv->d_u);
  if(NULL != fv->dd_u)   PGFEM_free(fv->dd_u);
  if(NULL != fv->f)      PGFEM_free(fv->f);
  if(NULL != fv->R)      PGFEM_free(fv->R);
  if(NULL != fv->f_defl) PGFEM_free(fv->f_defl);
  if(NULL != fv->RR)     PGFEM_free(fv->RR);
  if(NULL != fv->f_u)    PGFEM_free(fv->f_u);
  if(NULL != fv->RRn)    PGFEM_free(fv->RRn);
  if(NULL != fv->BS_x)   PGFEM_free(fv->BS_x);
  if(NULL != fv->BS_f)   PGFEM_free(fv->BS_f);
  if(NULL != fv->BS_f_u) PGFEM_free(fv->BS_f_u);
  if(NULL != fv->BS_RR)  PGFEM_free(fv->BS_RR);
  if(NULL != fv->coupled_physics_ids) free(fv->coupled_physics_ids);
  if(NULL != fv->fvs)    free(fv->fvs);

  destroy_eps_il(fv->eps,grid->element,grid->ne,opts->analysis_type);
  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
  {
    destroy_sig_il(fv->sig,grid->element,grid->ne,opts->analysis_type);
    if(opts->smoothing == 0)
      destroy_sig_el(fv->sig_n, grid->nn);
  }

  if(NULL != fv->statv_list)
    delete[] fv->statv_list;

  err += field_varialbe_initialization(fv);
  return err;
}


/// initialize field variables thermal part
///
/// \param[in, out] fv an object containing all field variables for thermal
/// \return non-zero on internal error
int thermal_field_varialbe_initialization(FieldVariablesThermal *fv)
{
  int err = 0;
  fv->Gndof = 0;
  fv->ndofn = 1;
  fv->ndofd = 0;
  fv->T_np1 = NULL;
  fv->T_n   = NULL;
  fv->T_nm1 = NULL;
  fv->dT    = NULL;
  fv->dd_T  = NULL;
  fv->NORM  = 0;
  return err;
}

/// prepare temporal varialbes for staggering Newton Raphson iterations
///
/// Before call this function, physics coupling should be defined in
/// fv->n_coupled and fv->coupled_physics_ids
///
/// \param[in, out] fv an object containing all field variables for thermal
/// \param[in] grid an object containing all mesh data
/// \param[in] is_for_Mechanical if yes, prepare constitutive models
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int prepare_temporal_field_varialbes(FieldVariables *fv,
                                     Grid *grid,
                                     int is_for_Mechanical,
                                     const PGFem3D_opt *opts)
{
  int err =0;
  fv->temporal = new FieldVariablesTemporal{};
  fv->temporal->u_n   = aloc1(grid->nn*fv->ndofn);
  fv->temporal->u_nm1 = aloc1(grid->nn*fv->ndofn);
  if(is_for_Mechanical)
  {
    if(opts->analysis_type==CM3F)
      fv->temporal->tf.construct(grid->ne, fv->nVol, fv->npres, true);

    const Element *elem = grid->element;

    int n_state_varialbles = 0;
    for(int eid=0; eid<grid->ne; eid++)
    {
      int nne = elem[eid].toe;
      long nint = 0;
      int_point(nne,&nint);
      n_state_varialbles += nint;
    }
    fv->temporal->element_variable_no = n_state_varialbles;
    fv->temporal->var = new State_variables[n_state_varialbles];

    for(int eid=0; eid<grid->ne; eid++)
    {
      int nne = elem[eid].toe;
      long nint = 0;
      int_point(nne,&nint);

      for(int ip=0; ip<nint; ip++)
      {
        Constitutive_model *m = &(fv->eps[eid].model[ip]);
        Model_var_info info;
        m->param->get_var_info(info);
        err += fv->temporal->var[m->model_id].initialization(info.n_Fs,
                                                             info.n_vars,
                                                             info.n_flags);
      }
    }
  }
  return err;
}

/// destory temporal varialbes for staggering Newton Raphson iterations
///
/// should be called before destroying fv
///
/// \param[in, out] fv an object containing all field variables for thermal
/// \param[in] is_for_Mechanical if yes, prepare constitutive models
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int destory_temporal_field_varialbes(FieldVariables *fv,
                                     int is_for_Mechanical,
                                     const PGFem3D_opt *opts)
{
  int err =0;
  PGFEM_free(fv->temporal->u_n);
  PGFEM_free(fv->temporal->u_nm1);
  fv->temporal->u_n = NULL;
  fv->temporal->u_nm1 = NULL;

  if(is_for_Mechanical)
  {
    delete[] fv->temporal->var;

    fv->temporal->var = NULL;
    fv->temporal->element_variable_no = 0;
  }
  delete fv->temporal;
  fv->temporal = NULL;
  return err;
}

/// initialize material properties
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] mat an object containing all material parameters
/// \return non-zero on internal error
int material_initialization(MaterialProperty *mat)
{
  int err = 0;
  mat->density  = NULL;
  mat->mater    = NULL;
  mat->thermal  = NULL;
  mat->hommat   = NULL;
  mat->matgeom  = NULL;
  mat->nhommat    = 0;
  mat->nmat       = 0;
  mat->n_orient   = 0;
  mat->n_co_props = 0;
  mat->co_props = NULL;
  return err;
}


/// destruct material properties
/// free memory spaces for member arrays and structs
///
/// \param[in, out] mat an object containing all material parameters
/// \return non-zero on internal error
int destruct_material(MaterialProperty *mat, const PGFem3D_opt *opts)
{
  int err = 0;
  if(NULL != mat->density) free(mat->density);
  if(NULL != mat->mater)   PGFEM_free(mat->mater);
  delete [] mat->thermal;
  destroy_matgeom(mat->matgeom,mat->n_orient);
  destroy_hommat(mat->hommat,mat->nhommat);

  if(opts->cohesive == 1)
    destroy_cohesive_props(mat->n_co_props,mat->co_props);

  err += material_initialization(mat);
  return err;
}

/// initialize loading steps object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] load an object containing boundary increments
/// \return non-zero on internal error
int loading_steps_initialization(LoadingSteps *load)
{
  int err = 0;
  load->sups        = NULL;
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

/// construct loading steps object
/// free memory spaces for member arrays and structs
///
/// \param[in, out] load an object containing boundary increments
/// \param[in] mp multiphysics object
/// \return non-zero on internal error
int construct_loading_steps(LoadingSteps *load, const Multiphysics& mp)
{
  int err = 0;

  load->sups        = (SUPP *)   malloc(sizeof(SUPP)*mp.physicsno);
  load->sup_defl    = (double **) malloc(sizeof(double *)*mp.physicsno);
  load->tim_load    = (long **)   malloc(sizeof(long *)*mp.physicsno);
  load->solver_file = (FILE **)   malloc(sizeof(FILE *)*mp.physicsno);

  return err;
}

/// destruct loading steps object
/// free memory spaces for member arrays and structs
///
/// \param[in, out] load an object containing boundary increments
/// \param[in] mp multiphysics object
/// \return non-zero on internal error
int destruct_loading_steps(LoadingSteps *load, const Multiphysics& mp)
{
  int err = 0;

  destroy_zatnode(load->znod,   load->nln);
  destroy_zatelem(load->zele_s, load->nle_s);
  destroy_zatelem(load->zele_v, load->nle_v);

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    destroy_supp(load->sups[ia]);
    if(NULL != load->tim_load[ia]) PGFEM_free(load->tim_load[ia]);
    if(NULL != load->sup_defl[ia]) PGFEM_free(load->sup_defl[ia]);
    if(NULL != load->solver_file[ia]) fclose(load->solver_file[ia]);
  }
  free(load->sups);
  free(load->tim_load);
  free(load->sup_defl);
  free(load->solver_file);


  err += loading_steps_initialization(load);
  return err;
}


/// initialize communication structures
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int communication_structure_initialization(CommunicationStructure *com)
{
  int err = 0;
  com->nproc  = 0;
  com->rank   = 0;
  com->Ap     = NULL;
    com->Ai     = NULL;
  com->DomDof = NULL;
  com->nbndel = 0;
  com->bndel  = NULL;
  com->GDof   = 0;
  com->NBN    = 0;
  com->comm   = 0;
  com->net    = NULL;
  com->boot   = NULL;
  com->spc    = NULL;
  com->hints  = NULL;
  return err;
}


/// destruct communication structures
/// free memory spaces for member arrays and structs
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int destruct_communication_structure(CommunicationStructure *com)
{
  int err = 0;
  if(NULL != com->Ap) PGFEM_free(com->Ap);
  if(NULL != com->Ai) PGFEM_free(com->Ai);
  if(NULL != com->DomDof) PGFEM_free(com->DomDof);
  if(NULL != com->bndel)  PGFEM_free(com->bndel);

  err += communication_structure_initialization(com);
  return err;
}

/// initialize multiphysics object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs)
///
/// \param[in, out] mp an object for multiphysics stepping
/// \return non-zero on internal error
int multiphysics_initialization(Multiphysics& mp)
{
  int err = 0;
  mp.physicsno   = 0;
  mp.physicsname = NULL;
  mp.physics_ids = NULL;
  mp.ndim        = NULL;
  mp.write_no    = NULL;
  mp.total_write_no = 0;
  return err = 0;
}


/// construct multiphysics object
/// create memory space for member arrays and structs
///
/// \param[in, out] mp an object for multiphysics stepping
/// \param[in] physicsno number of physics
/// \return non-zero on internal error
int construct_multiphysics(Multiphysics& mp,
                           int physicsno)
{
  int err = 0;
  mp.physicsno   = physicsno;
  mp.physicsname = (char **) malloc(sizeof(char *)*physicsno);
  mp.physics_ids = (int*) malloc(sizeof(int)*physicsno);
  mp.ndim        = (int*) malloc(sizeof(int)*physicsno);
  mp.write_no    = (int*) malloc(sizeof(int)*physicsno);
  mp.write_ids.resize(physicsno);
  mp.coupled_ids.resize(physicsno);

  for(int ia=0; ia<physicsno; ia++)
  {
    mp.physicsname[ia] = (char *) malloc(sizeof(char)*1024);
    mp.physics_ids[ia] = 0;
    mp.ndim[ia]        = 0;
    mp.write_no[ia]    = 0;
  }
  mp.total_write_no  = 0;
  return err = 0;
}

/// set a physics
/// set physics id, number of degree freedom, and name
///
/// \param[in, out] mp an object for multiphysics
/// \param[in] obj_id id to access each physics
/// \param[in] mp_id multiphysics id
/// \param[in] n_dof number of degree freedom of the physics
/// \param[in] name physics name
/// \return non-zero on internal error
static int set_a_physics(const Multiphysics& mp,
                  int obj_id,
                  int mp_id,
                  int n_dof,
                  const char *name)
{
  int err = 0;
  mp.physics_ids[obj_id] = mp_id;
  mp.ndim[obj_id]        = n_dof;
  sprintf(mp.physicsname[obj_id], "%s", name);
  return err = 0;
}

/// destruct multiphysics object
/// free memory spaces for member arrays and structs
///
/// \param[in, out] mp an object for multiphysics stepping
/// \return non-zero on internal error
int destruct_multiphysics(Multiphysics& mp)
{
  int err = 0;
  if(NULL != mp.physicsname)
  {
    for(int ia=0; ia<mp.physicsno; ia++)
      if(NULL != mp.physicsname[ia]) free(mp.physicsname[ia]);

    free(mp.physicsname);
  }

  if(NULL != mp.physics_ids) free(mp.physics_ids);
  if(NULL != mp.ndim)        free(mp.ndim);
  if(NULL != mp.write_no)    free(mp.write_no);
  multiphysics_initialization(mp);
  return err = 0;
}

/// read and construct multiphysics
/// if multiphysics.in is not provied, default(mechanical) will be set
/// the file format is as below:
///
/// multiphysics.in
/// # <= comments
/// # Number of physic
/// 2
/// ######################################################
/// # Physics setting for thermal part
/// # [Physics id] = 0: Mechanical
/// #                1: Thermal
/// #                2: Chemistry (Not available)
/// #
/// # [Physics id] [Physics name] [degree of freedom]
/// 1 Thermal 1
/// #
/// # Print result
/// # [number of varialbs to be printed]
/// 2
/// # [data id] = 0: temperature
/// #             1: heat flux
/// #             2: heat generation
/// # curretly temperature and heat flux are supported
/// 0 1
/// # end of thermal
/// ######################################################
/// # Physics setting for mechanical part
/// 0 Mechanical 3
/// #
/// # Print result
/// 3
/// # [data id] = 0: Displacement
/// #             1: CauchyStress
/// #             2: EulerStrain
/// #             3: EffectiveStrain
/// #             4: EffectiveStress
/// #             5: CellProperty
/// #             6: Damage
/// #             7: Chi for damage model
/// #             8: F(total deformation gradient)
/// #             9: P(1st Piola Kirchhoff stress)
/// #             10: W(strain energy density)
/// 0 1 5
///
/// \param[in, out] mp an object for multiphysics stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_multiphysics_settings(Multiphysics& mp,
                               const PGFem3D_opt *opts,
                               int myrank)
{
  int err = 0;
  int physicsno = 0;

  char filename[1024];
  sprintf(filename,"%s/multiphysics.in",opts->ipath);
  FILE *in = NULL;
  in = fopen(filename,"r");

  if(in==NULL) // check file is readable
  {
    if(myrank==0)
    {
      PGFEM_printf("no [%s/multiphysics.in] is provided\n", opts->ipath);
      PGFEM_printf("Use default setting (Mechanical only).\n");
    }
  }
  else
  {
    err += scan_for_valid_line(in);
    CHECK_SCANF(in, "%d", &physicsno);
    if(physicsno>0)
    {
      err += construct_multiphysics(mp, physicsno);

      int physics_id, ndof;
      int n_couple;

      char name[1024];
      int cnt_pmr = 0;
      for(int ia=0; ia<physicsno; ia++)
      {
        // read physics id, physics name and number of degree of freedons on node
        err += scan_for_valid_line(in);
        CHECK_SCANF(in, "%d%s%d", &physics_id,name,&ndof);
        err += set_a_physics(mp, ia,physics_id,ndof, name);

        // read ids for coupling
        err += scan_for_valid_line(in);
        CHECK_SCANF(in, "%d", &n_couple);

        if(n_couple<0)
          n_couple = 0;

        mp.coupled_ids[ia].resize(n_couple + 1);
        mp.coupled_ids[ia][0] = n_couple;
        for(int ib = 0; ib<n_couple; ib++)
          CHECK_SCANF(in, "%d", &(mp.coupled_ids[ia][ib+1]));

        // read ids for writing results
        err += scan_for_valid_line(in);
        CHECK_SCANF(in, "%d", mp.write_no+ia);


        if(mp.write_no[ia]>0)
        {
          // read from file for writing results
          mp.write_ids[ia].resize(mp.write_no[ia]);
          err += scan_for_valid_line(in);
          cnt_pmr += mp.write_no[ia];
          for(int ib=0; ib<mp.write_no[ia]; ib++)
          CHECK_SCANF(in, "%d", &(mp.write_ids[ia][ib]));
        }
        if(mp.write_no[ia]==-1)
        {
          switch(physics_id)
          {
            case MULTIPHYSICS_MECHANICAL:
              mp.write_no[ia] = MECHANICAL_Var_NO;
              break;
            case MULTIPHYSICS_THERMAL:
              mp.write_no[ia] = Thermal_Var_NO;
              break;
            case MULTIPHYSICS_CHEMICAL:
              mp.write_no[ia] = CHEMICAL_Var_NO;
              break;
            default:
              mp.write_no[ia] = MECHANICAL_Var_NO;
          }

          // set default : all outputs
          mp.write_ids[ia].resize(mp.write_no[ia]);
          err += scan_for_valid_line(in);
          cnt_pmr += mp.write_no[ia];
          for(int ib=0; ib<mp.write_no[ia]; ib++)
            mp.write_ids[ia][ib] = ib;
        }
      }
      mp.total_write_no  = cnt_pmr;
    }
    fclose(in); // close file
  }

  if(physicsno<=0)
  {
    err += construct_multiphysics(mp, 1);
    err += set_a_physics(mp, 0, MULTIPHYSICS_MECHANICAL, 3, "Mechanical");

    mp.coupled_ids[0].resize(1);
    mp.coupled_ids[0][0] = 0;

    mp.write_no[0] = MECHANICAL_Var_NO;
    mp.write_ids[0].resize(mp.write_no[0]);
    for(int ib=0; ib<mp.write_no[0]; ib++)
      mp.write_ids[0][ib] = ib;
    mp.total_write_no  = MECHANICAL_Var_NO;
  }
  // print multiphysics setting
  if(myrank==0)
  {
    PGFEM_printf("Total number of physics: %d\n", mp.physicsno);
    for(int ia=0; ia<mp.physicsno; ia++)
    {
      PGFEM_printf("%d. physics name \t\t= %s\n", ia, mp.physicsname[ia]);
      PGFEM_printf("   # of unknown on node \t= %d\n", mp.ndim[ia]);
      PGFEM_printf("   # of physics to be coupled \t= %d", mp.coupled_ids[ia][0]);
      PGFEM_printf(", ids = ");
      for(int ib=0; ib<mp.coupled_ids[ia][0]; ib++)
        PGFEM_printf("%d ", mp.coupled_ids[ia][ib+1]);

      PGFEM_printf("\n");

      PGFEM_printf("   # of output variables \t= %d", mp.write_no[ia]);
      PGFEM_printf(", ids = ");

      assert((int)mp.write_ids[ia].size() >= mp.write_no[ia] && "attempted mp.write_ids[ia] out of bounds access");
      for(int ib=0; ib<mp.write_no[ia]; ib++){
        PGFEM_printf("%d ", mp.write_ids[ia][ib]);
      }

      PGFEM_printf("\n\n");
    }
  }
  return err;
}

/// construct element object in NBE contributing for NBC by features
///
/// \param[in] elem   element objct, used to identify boundary elements
/// \param[in] elemno total number of elements
/// \param[in] nsd    number of spatial dimensions
void
NeumannBoundaryElement::construct_list_of_boundary_elements(Element *elem,
                                                            const int elemno,
                                                            const int nsd){
  this->nbe_no = 0;
  this->element_check.initialization(elemno, 2, 0);
  
  // count number of element with NBE
  for(int eid=0; eid<elemno; ++eid){
    FEMLIB fe;
    // number of boundary element in element[eid]
    int bnd_elem_no = fe.number_of_boundary_elements(elem[eid].toe, nsd);
        
    // check this element has features for NB
    for(int ia=0; ia<this->feature_no; ++ia){
      for(int ib = 0; ib<bnd_elem_no; ++ib){
        if(this->bnd_feature(ia) == elem[eid].bnd_type[ib] &&
           this->bnd_feature_id(ia) == elem[eid].bnd_id[ib]){
          this->element_check(eid, 0) = 1;
          break;
        }
      }
      if(this->element_check(eid, 0) == 1){
        ++(this->nbe_no); 
        break;
      }
    }
  }
  
  // set size of lists 
  this->element_ids.initialization(this->nbe_no, 2, 0);
  this->bnd_elements.initialization(this->nbe_no, 1);
  // re-visit element to set list of NBE
  int nbe_cnt = 0;
  for(int eid=0; eid<elemno; ++eid){
    if(this->element_check(eid, 0)==0)
      continue;
      
    FEMLIB fe;
    // number of boundary element in element[eid]
    int bnd_elem_no = fe.number_of_boundary_elements(elem[eid].toe, nsd);
    
    int bnd_elems[MAX_BND_ELEMENT_NO] = {};
    int feature_idx[MAX_BND_ELEMENT_NO] = {};
    
    int cnt_bnd_elem = 0;
    for(int ia=0; ia<this->feature_no; ++ia){
      for(int ib = 0; ib<bnd_elem_no; ++ib){
        if(this->bnd_feature(ia) == elem[eid].bnd_type[ib] &&
           this->bnd_feature_id(ia) == elem[eid].bnd_id[ib]){
          bnd_elems[ib]   = 1;
          feature_idx[ib] = ia;
          ++cnt_bnd_elem;
        }
      }
    }
    if(cnt_bnd_elem>0){ // if this element has NB elements
                        // set size and feature ids
      this->element_check(eid, 1)   = nbe_cnt;
      this->element_ids(nbe_cnt, 0) = eid;
      this->element_ids(nbe_cnt, 1) = cnt_bnd_elem;      
      this->bnd_elements(nbe_cnt).initialization(cnt_bnd_elem, 2, 0);
      for(int ia=0, ib=0; ia<bnd_elem_no; ++ia){
        if(bnd_elems[ia] == 1){
          this->bnd_elements(nbe_cnt)(ib,   0) = ia;
          this->bnd_elements(nbe_cnt)(ib++, 1) = feature_idx[ia];
        }
      }
      ++nbe_cnt; // increase number of boundary elements for next count
    }
  }      
}


/// print all feature values
void 
NeumannBoundaryElement::print_feature(void){
  features.print("features", "%d");
  for(int ia=0, nf = feature_no; ia<nf; ia++){
    std::cout << "Load type(" << load_type(ia) << "): ";
    for(int ib=0; ib<load_type(ia); ++ib)
      std::cout << load(ia, ib) << " ";

    std::cout << endl;
  }
}

/// print all boundary elements values
void 
NeumannBoundaryElement::print_bnd_elements(void){
  for(int ia=0; ia<nbe_no; ++ia){
    for(int ib=0; ib<element_ids(ia, 1); ++ib){
      std::cout << "eid : " << element_ids(ia, 0);
      std::cout << " face no : " << element_ids(ia, 1);
      std::cout << " face id : " << bnd_elements(ia)(ib, 0);
      std::cout << " load id : " << bnd_elements(ia)(ib, 1) << endl;
    }
  }
}

/// read load (math expression) from in_dir/load_name and save
/// at NBE.load(feature_id, load_id)
///
/// \param[in] in_dir     NBE directory path
/// \param[in] load_fn    file name of math expression for NBC load
/// \param[in] feature_id feature id of this load
/// \param[in] load_id    load id of a feature
///                       load_id = 0 to 2: for condition (load type = 3) : [tx, ty, tz]
///                               = 0 only  for condition (load type = 1) : t.N
void
NeumannBoundaryElement::read_loads(const std::string &in_dir,
                                   const std::string &load_fn, 
                                   const int feature_id,
                                   const int load_id){
  std::stringstream fn;
  fn << in_dir << "/" << load_fn;
  std::ifstream ifs;
  ifs.open(fn.str());
  if(ifs.is_open())
    std::getline(ifs, load(feature_id, load_id));
  else{
    std::cout << "Cannot open [" << fn.str() << "]" << endl;
    abort();
  }
}

/// compute load values from string_function_of_time using string in load
///
/// \param[in] tnp1 time at n+1
/// \param[in] tn   time at n
/// \param[in] tnm1 time at n-1
void
NeumannBoundaryElement::compute_load_values(double tnp1,
                                            double tn,
                                            double tnm1){

  //cannot guarantee calling function with negative time is valid
  if(tn<0)
    tn = 0.0;
    
  if(tnm1<0)
    tnm1 = 0.0;
    
  for(int ia=0; ia<feature_no; ++ia){
    for(int ib=0; ib<load_type(ia); ++ib){
      load_values_np1(ia, ib) = string_function_of_time(load(ia, ib), tnp1);
      load_values_n(  ia, ib) = string_function_of_time(load(ia, ib), tn);
      load_values_nm1(ia, ib) = string_function_of_time(load(ia, ib), tnm1);
    }
  }
}