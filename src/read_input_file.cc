/// Functions are defined for reading input files
///
/// Authors:
///   Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
///   Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
///   Sangmin Lee, University of Notre Dame <slee43 [at] nd.edu>
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "read_input_file.h"
#include "Arc_length.h"
#include "PGFem3D_data_structure.h"
#include "allocation.h"
#include "constitutive_model.h"
#include "enumerations.h"
#include "gen_path.h"
#include "in.h"
#include "load.h"
#include "read_cryst_plast.h"
#include "restart.h"
#include "utils.h"
#include <cstring>
#include "three_field_element.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

/// read mechanical part of material properties
///
/// \param[in] fp file pointer for reading mechanical part of material properties
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \return non-zero on interal error
int read_material_for_Mechanical(FILE *fp,
                                 MaterialProperty *mat,
                                 const PGFem3D_opt *opts)
{
  int err = 0;
  for(int ia=0; ia<mat->nmat; ia++)
  {
    scan_for_valid_line(fp);
    if(read_material(fp,ia,mat->mater,opts->legacy))
      PGFEM_Abort();
  }
  return err;
}

/// read mechanical part of material properties
///
/// \param[in] fp file pointer for reading mechanical part of material properties
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \return non-zero on interal error
int read_material_for_Thermal(FILE *fp,
                              MaterialProperty *mat,
                              const PGFem3D_opt *opts)
{
  int err = 0;
  int param_in = 10;

  MaterialThermal *thermal = new MaterialThermal[mat->nmat];

  for(int ia=0; ia<mat->nmat; ia++)
  {
    double cp;
    double k[9];
    double FHS_MW = 1.0;

    int match = 0;
    scan_for_valid_line(fp);
    match += fscanf(fp, "%lf", &cp);

    scan_for_valid_line(fp);
    match += fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", k+0, k+1, k+2
                    , k+3, k+4, k+5
                    , k+6, k+7, k+8);

    if(match != param_in)
      PGFEM_Abort();

    thermal[ia].cp = cp;
    for(int ib=0; ib<9; ib++)
      thermal[ia].k[ib] = k[ib];

    // read optional value (fraction of heat source from mechanical work)
    scan_for_valid_line(fp);
    int read_no = fscanf(fp, "%lf", &FHS_MW);

    if(read_no == 1)
      thermal[ia].FHS_MW = FHS_MW;
    else
      thermal[ia].FHS_MW = 1.0;
  }

  mat->thermal = thermal;
  return err;
}

/// read material properties for multiphysics problem
///
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \param[in] mp multiphysics object
/// \return non-zero on interal error
int read_multiphysics_material_properties(MaterialProperty *mat,
                                          const PGFem3D_opt *opts,
                                          const Multiphysics& mp,
					  int myrank)
{
  int err = 0;
  char dirname[1024], fn[2048];
  sprintf(dirname,"%s/Material",opts->ipath);

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    sprintf(fn,"%s/%s.mat",dirname,mp.physicsname[ia]);

    FILE *fp = NULL;
    fp = fopen(fn, "r");

    if(fp==NULL)
    {
      if(myrank==0)
        PGFEM_printf("No [%s] exists.\n", fn);

      continue;
    }

    switch(mp.physics_ids[ia])
    {
     case MULTIPHYSICS_MECHANICAL:
      err += read_material_for_Mechanical(fp,mat,opts);
      break;
     case MULTIPHYSICS_THERMAL:
      err += read_material_for_Thermal(fp,mat,opts);
      break;
     case MULTIPHYSICS_CHEMICAL:
      break;
     default:
      break;
    }
    fclose(fp);
  }

  // read and set general material properties e.g. density

  mat->density = (double *) malloc(sizeof(double)*(mat->nmat));
  double *d = mat->density;
  for(int ia=0; ia<mat->nmat; ia++)
    d[ia] = 0.0;

  sprintf(fn,"%s/material.mat",dirname);
  FILE *fp = NULL;
  fp = fopen(fn, "r");

  if(fp == NULL)
  {
    if(myrank==0)
      PGFEM_printf("No [%s] exists. \nDensity is set to zero.\n", fn);

    return err;
  }

  int match = 0;
  for(int ia=0; ia<mat->nmat; ia++)
  {
    scan_for_valid_line(fp);
    match += fscanf(fp, "%lf", d+ia);
  }

  fclose(fp);

  // check number of densities that is read.
  if(match != mat->nmat)
  {
    if(myrank==0)
      PGFEM_printf("Material density is not read as many as number of materials.\n");

    PGFEM_Abort();
  }
  else
  {
    for(int ia = 0; ia<mat->nhommat; ia++)
    {
      (mat->hommat[ia]).density = d[(mat->hommat[ia]).mat_id];
      if(myrank==0)
        PGFEM_printf("Density(%d), %e\n", ia, (mat->hommat[ia]).density);
    }
  }
  return err;
}

/// count number of ranges that are seperated by comma
///
/// \param[in] str a string containing ranges
/// \return a integer number (counted number of ranges)
int count_number_of_ranges(char str[])
{
  int charno = strlen(str);
  int rangeno = 1;

  for(int a=0; a<=charno; a++)
  {
    if(str[a]==',')
      rangeno++;
  }

  return rangeno;
}

/// get list of ranges
/// \param[in] str_in a string containing ranges
///      str_in format:
///      ---------------------------------------------------------------
///      srt_in = "0:0.001:0.1"
///         or
///      str_in = "0:0.001:0.1,0.105:0.005:0.2"
///      ---------------------------------------------------------------
///      range format: start_value:step_size:end_value
///      ex) 1:1:10 => 1,2,3,4,5,6,7,8,9,10
///      ex) 0:0.2:1 => 0,0.2,0.4,0.6,0.8,1
/// \param[out] ranges double array of ranges read from the input string
/// \return non-zero on internal error
int interpret_ranges(double *ranges, char str_in[])
{
  int err = 0;
  int rangeno = count_number_of_ranges(str_in);
  int charno = strlen(str_in);

  char str[1024], *pch;
  memcpy(str, str_in, charno+1);

  pch = strtok (str,",");
  for(int a=0; a<rangeno; a++)
  {
    sscanf(pch, "%lf:%lf:%lf", ranges+a*3+0, ranges+a*3+1, ranges+a*3+2);
    pch = strtok(NULL, ",");
  }

  return err;
}

int read_input_file(const PGFem3D_opt *opts,
		    const CommunicationStructure *com,
                    long *nn,
                    long *Gnn,
                    long *ndofn,
                    long *ne,
                    long *lin_maxit,
                    double *lin_err,
                    double *lim_zero,
                    long *nmat,
                    long *n_concentrations,
                    long *n_orient,
                    Node **node,
                    Element **elem,
                    Material **material,
                    MATGEOM *matgeom,
                    SUPP *sup,
                    long *nln,
                    ZATNODE **znod,
                    long *nel_s,
                    ZATELEM **zelem_s,
                    long *nel_v,
                    ZATELEM **zelem_v,
                    const int *fv_ndofn,
                    const int physicsno,
                    const int *ndim,
                    char **physicsnames)
{
  int err = 0;
  int myrank = com->rank;

  /* compute filename and open file */
  char *filename = PGFEM_calloc(char, 500);
  sprintf(filename,"%s/%s%d.in",opts->ipath,opts->ifname,myrank);
  FILE *in = PGFEM_fopen(filename,"r");

  /* read header lines */
  CHECK_SCANF (in,"%ld %ld %ld",nn,ndofn,ne);
  CHECK_SCANF (in,"%ld %lf %lf",lin_maxit,lin_err,lim_zero);
  CHECK_SCANF (in,"%ld %ld %ld",nmat,n_concentrations,n_orient);

  (*node) = build_node_multi_physics(*nn,fv_ndofn,physicsno);
  (*elem) = build_elem(in,*ne,opts->analysis_type);
  (*material) = PGFEM_calloc(Material, *nmat);
  (*matgeom) = build_matgeom(*n_concentrations,*n_orient);

  *Gnn = read_nodes(in,*nn,*node,opts->legacy,com);
  /* NOTE: Supports assume only ndim supported dofs per node! */

  char BC[1024];
  sprintf(BC,"%s/BC",opts->ipath);

  if(is_directory_exist(BC))
  {
    if(myrank==0)
      PGFEM_printf("BC exists skip BC from filebase_*.in instead read boundary conditions from BC\n");

    // skip reading support from lagacy inputs
    int nbc; // temporal, don't need here
    int n[4];
    CHECK_SCANF(in, "%d", &nbc);
    for(int ia=0; ia<nbc; ia++)
      CHECK_SCANF(in, "%d %d %d %d", n+0,n+1,n+2,n+3);

    CHECK_SCANF(in, "%d", &nbc);
    double v;

    for(int ia=0; ia<nbc; ia++)
      CHECK_SCANF(in, "%lf", &v);

    // read boundary conditions if BC diretory exists
    char fn_bc[2048];
    char fn_bcv[2048];

    for(int ia=0; ia<physicsno; ia++)
    {
      sprintf(fn_bc,"%s/%s_%d.bc",BC,physicsnames[ia],myrank);
      sprintf(fn_bcv,"%s/%s.bcv",BC,physicsnames[ia]);

      FILE *fp = NULL;
      fp = fopen(fn_bc, "r");
      sup[ia] = read_Dirichlet_BCs(fp,*nn,ndim[ia],*node,ia);
      if(fp!=NULL) fclose(fp);

      fp = NULL;
      fp = fopen(fn_bcv, "r");
      if(fp==NULL)
      {
        PGFEM_printf("ERROR: Cannot open %s file. Exit.\n", fn_bcv);
        PGFEM_Abort();
      }
      err += read_Dirichlet_BCs_values(fp,*nn,ndim[ia],*node,sup[ia],ia);
      fclose(fp);
    }
  }
  else
  {
    for(int ia=0; ia<physicsno; ia++)
      sup[ia] = read_supports(in,*nn,ndim[ia],*node, ia);
  }

  read_elem(in,*ne,*elem,*sup,opts->legacy);
  for(int i=0, e=*nmat; i<e; i++){
    if ( read_material(in,i,*material,opts->legacy) ){
      PGFEM_Abort();
    }
  }

  if (feof(in)) {
    PGFEM_printerr("ERROR: prematurely reached EOF in %s(%s)\n",
                   __func__,__FILE__);
    PGFEM_Abort();
  }

  if ( override_material_properties(*nmat,opts,*material) ) {
    PGFEM_Abort();
  }

  read_matgeom(in,*n_concentrations,*n_orient,*matgeom);

  /* NOTE: Node/Element loading assumes forces only in ndim
     directions */
  /* node */
  CHECK_SCANF(in,"%ld",nln);
  *znod = build_zatnode (*ndofn,*nln);
  read_nodal_load (in,*nln,*ndofn,*znod);
  /* surface */
  CHECK_SCANF (in,"%ld",nel_s);
  *zelem_s = build_zatelem (*ndofn,*nel_s);
  read_elem_surface_load (in,*nel_s,*ndofn,*elem,*zelem_s);
  /* volume */
  CHECK_SCANF (in,"%ld",nel_v);
  *zelem_v = build_zatelem (*ndofn,*nel_v);

  /* check the ferror bit */
  if(ferror(in)) err++;

  /* free local memory and close file */
  free(filename);
  fclose(in);
  return err;
}

/// This function read mesh info, boundary conditions, and material properties
/// from main input files (*.in). While reading inputs, node, element, material,
/// and support objects are constructed. This function still utilizes
/// lagacy function (read_input_file). Later, this function and read_input_file need to
/// be combined.
///
/// \param[out] grid a mesh object
/// \param[out] mat a material object
/// \param[out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[in] mp multiphysics object
/// \param[in] com handle for communication
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int read_mesh_file(Grid *grid,
                   MaterialProperty *mat,
                   FieldVariables *FV,
                   Solver *SOL,
                   LoadingSteps *load,
                   const Multiphysics& mp,
		   const CommunicationStructure *com,
                   const PGFem3D_opt *opts)
{
  long ndofn;
  int myrank = com->rank;

  int *fv_ndofn = (int *) malloc(mp.physicsno*sizeof(int));

  for(int iA=0; iA<mp.physicsno; iA++)
    fv_ndofn[iA] = FV[iA].ndofn;

  int err = read_input_file(opts,
                            com,
                            &(grid->nn),
                            &(grid->Gnn),
                            &ndofn,
                            &(grid->ne),
                            &(SOL[0].iter_max_sol),
                            &(SOL[0].err),
                            &(SOL[0].computer_zero),
                            &(mat->nmat),
                            &(FV[0].n_concentrations),
                            &(mat->n_orient),
                            &(grid->node),
                            &(grid->element),
                            &(mat->mater),
                            &(mat->matgeom),
                            load->sups,
                            &(load->nln),
                            &(load->znod),
                            &(load->nle_s),
                            &(load->zele_s),
                            &(load->nle_v),
                            &(load->zele_v),
                            fv_ndofn,
                            mp.physicsno,
                            mp.ndim,
                            mp.physicsname);
  free(fv_ndofn);

  // read multiphysics material properties
  err += read_multiphysics_material_properties(mat,opts,mp,myrank);

  // update numerical solution scheme parameters
  FV[0].NORM = SOL[0].computer_zero;
  for(int iA=0; iA<mp.physicsno; iA++)
  {
    SOL[iA].iter_max_sol  = SOL[0].iter_max_sol;
    SOL[iA].err           = SOL[0].err;
    SOL[iA].computer_zero = SOL[0].computer_zero;
    FV[iA].n_concentrations = FV[0].n_concentrations;
    FV[iA].NORM             = FV[0].NORM;
  }
  // need to update number of elements that have prescribed BCs (supported)
  for (long ia=0;ia<grid->ne;ia++)
  {
    for(int iA = 1; iA<mp.physicsno; iA++) // iA = 0 is alreaded accounted in read_elem in read_input_file
    {
      const long *nod = grid->element[ia].nod;
      const int nne = grid->element[ia].toe;
      int is_it_supp = 0;
      for(long ja=0; ja<load->sups[iA]->ndn; ja++)
      {
        for (int ka=0; ka<nne; ka++)
        {
          if(load->sups[iA]->lnpd[ja] == nod[ka])
          {
            is_it_supp = 1;
            break;
          }
        }
        if(is_it_supp)
          break;
      }
      if(is_it_supp)
        (load->sups[iA]->nde)++;
    }
  }
  return err;
}

int read_time_steps(FILE *fp, TimeStepping *ts)
{
  int err = 0;

  // read number of computational times
  CHECK_SCANF(fp,"%ld",&(ts->nt));

  // read times
  ts->times = aloc1(ts->nt+1);
  for (long ia=0; ia<ts->nt+1; ia++)
    CHECK_SCANF(fp,"%lf",(ts->times)+ia);

  long n_p = 0;
  // read times for output
  CHECK_SCANF(fp,"%ld",&n_p);

  //Times for printing
  ts->print = times_print(fp,ts->nt,n_p);

  return err;
}

/// Read solver file for time stepping.
/// If command line includes override solver file option, all solver files will be overrided.
/// At the end of this function, file pointer is stored in LoadingSteps
/// in order to read load increments as time elapses.
/// Detailed slover file format can be found at the following link:
/// https://wiki-cswarm.crc.nd.edu/foswiki/pub/Main/CodeDevelopment/PGFem3DQuickStarts/generate_input_file.pdf
///
/// \param[out] time_steps object for time stepping
/// \param[out] mat a material object
/// \param[out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] crpl object for lagcy crystal plasticity
/// \param[in] mp multiphysics object
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_solver_file(TimeStepping *ts,
                     MaterialProperty *mat,
                     FieldVariables *FV,
                     Solver *SOL,
                     LoadingSteps *load,
                     CRPL *crpl,
                     const Multiphysics& mp,
                     const PGFem3D_opt *opts,
                     int myrank)
{
  int err = 0;
  // READ SOLVER FILE
  // override the default solver file with one specified
  // at commandline
  char filename[2048];
  char in_dat[1024];
  FILE *fp = NULL;
  if(opts->override_solver_file)
  {
    if(myrank == 0)
      PGFEM_printf("Overriding the default solver file with:\n%s\n", opts->solver_file);

    fp = fopen(opts->solver_file,"r");
  }
  else
  {
    // use the default file/filename
    sprintf(in_dat,"%s/%s",opts->ipath,opts->ifname);
    sprintf(filename,"%s%d.in.st",in_dat,myrank);
    fp = fopen(filename,"r");
    if(fp==NULL)
    {
      sprintf(filename,"%s%d.in.st",in_dat,0);
      fp = fopen(filename,"r");
    }
  }

  scan_for_valid_line(fp);

  CHECK_SCANF (fp,"%lf %ld %ld %ld",&(SOL[0].nor_min),&(SOL[0].iter_max),&(FV[0].npres),&(SOL[0].FNR));
  if(SOL[0].FNR == 2 || SOL[0].FNR == 3)
  {
    SOL[0].arc = (ARC_LENGTH_VARIABLES *) malloc(sizeof(ARC_LENGTH_VARIABLES));
    err += arc_length_variable_initialization(SOL[0].arc);
    CHECK_SCANF (fp,"%lf %lf",&(SOL[0].arc->dAL0),&(SOL[0].arc->dALMAX));
  }

  if(SOL[0].FNR == 4)
  {
    int physicsno = 0;
    CHECK_SCANF(fp, "%d %d", &physicsno, &(SOL[0].max_NR_staggering));

    if(physicsno != mp.physicsno)
    {
      if(myrank==0)
        PGFEM_printf("ERROR: Number of physics for setting parameters for the solver is not correct. Abort\n");

      PGFEM_Abort();
    }

    scan_for_valid_line(fp);

    for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
    {
      SOL[mp_id].FNR = 1;
      if(mp_id>0)
      {
        SOL[mp_id].max_NR_staggering = SOL[0].max_NR_staggering;
        SOL[mp_id].nor_min           = SOL[0].nor_min;
        SOL[mp_id].iter_max          = SOL[0].iter_max;
        FV[mp_id].npres              = FV[0].npres;
      }

      CHECK_SCANF (fp,"%d %d",&(SOL[mp_id].max_subdivision), &(SOL[mp_id].set_initial_residual));
      if(SOL[mp_id].set_initial_residual)
        CHECK_SCANF (fp,"%lf",&(SOL[mp_id].du));
    }
  }

  scan_for_valid_line(fp);

  // CRYSTAL PLASTICITY
  if(opts->analysis_type == FS_CRPL) {
    crpl = PGFEM_calloc (CRPL, mat->nmat);
    read_cryst_plast(fp,mat->nmat,crpl,opts->plc);
  }

  err += read_time_steps(fp, ts);

  // loading history exists in load directory
  char load_path[1024];
  char load_fn[2048];
  sprintf(load_path,"%s/load",opts->ipath);

  ts->tns = aloc1(mp.physicsno);
  int is_load_exist = 0;
  for(int ia=0; ia<mp.physicsno; ia++)
  {
    ts->tns[ia] = ts->times[0]; // set t(n) for individual physics

    sprintf(load_fn,"%s/%s.load",load_path,mp.physicsname[ia]);
    load->solver_file[ia] = NULL;
    load->solver_file[ia] = fopen(load_fn, "r"); // Load increments are needed to be read
    // while time is elapsing.
    // This file point needs to be freed end of the simulation
    // by calling destruction of the LoadingSteps

    if(load->solver_file[ia]==NULL)
      continue;

    is_load_exist = 1;
    long nlod_tim = 0;
    CHECK_SCANF(load->solver_file[ia],"%ld",&nlod_tim);

    // read times dependent load
    load->tim_load[ia] = NULL;
    load->tim_load[ia] = compute_times_load(load->solver_file[ia],ts->nt,nlod_tim);
  }

  // if no load directory exists,
  // loading history will be read and saved only in the 1st physics
  // from the lagacy solver file.
  if(is_load_exist==0)
  {
    long nlod_tim = 0;
    CHECK_SCANF (fp,"%ld",&nlod_tim);
    /* read times dependent load */
    load->tim_load[0] = compute_times_load(fp,ts->nt,nlod_tim);
    load->solver_file[0] = fp; // load increments are still need to be read
                               // while time is elapsing
                               // this file point needs to be freed end of the simulation
                               // by calling destruction of the LoadingSteps
  }

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    // if loading is not defined, set default (all zeros)
    if(load->solver_file[ia] == NULL)
    {
      if(ts->nt == 0)
        load->tim_load[ia] = aloc1l(1);
      else
        load->tim_load[ia] = aloc1l(ts->nt);
    }
  }
  return err;
}

/// Read initial conditions from lagcy format.
///
/// If no restart, this function reads initial conditions from *.initial files.
/// The file format can be found at the following link:
/// https://gitlab-cswarm.crc.nd.edu/pgfem_3d/pgfem_3d/wikis/how-to-set-initial-values
/// If restart is set from the commend line, *.inital files are used only for reading
/// mid point rule and material densities. Initial conditions are set by reading restart files.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv object for field variables
/// \param[out] sol object for solution scheme
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[out] tnm1 if restart, read time step info from the previous run
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_initial_values_lagcy(Grid *grid,
                              MaterialProperty *mat,
                              FieldVariables *fv,
                              Solver *sol,
                              LoadingSteps *load,
                              TimeStepping *ts,
                              PGFem3D_opt *opts,
                              const Multiphysics& mp,
                              double *tnm1,
                              int myrank)
{
  int err = 0;
  int mp_id = 0;

  char filename[1024];
  char line[1024];
  double dt = ts->times[1] - ts->times[0];;

  sprintf(filename,"%s/%s%d.initial",opts->ipath,opts->ifname,0);

  // restart option from command line is -1
  // check restart form initial file.
  if(opts->restart < 0)
  {
    FILE *fp_0 = fopen(filename,"r");

    if(fp_0 != NULL)
    {
      while(fgets(line, 1024, fp_0)!=NULL)
      {
        if(line[0]=='#')
          continue;

        sscanf(line, "%d", &(opts->restart));
        break;
      }
      fclose(fp_0);
    }
  }

  // check restart and read values
  if(opts->restart >= 0)
    err += read_restart(grid,fv,ts,load,opts,mp,tnm1,myrank);

  sprintf(filename,"%s/%s%d.initial",opts->ipath,opts->ifname,myrank);
  FILE *fp = fopen(filename,"r");

  int read_initial_file = 1;
  if(fp == NULL)
  {
    read_initial_file = 0;
    if(myrank==0)
      PGFEM_printf("Fail to open file [%s]. Quasi steady state\n", filename);
  }
  else
  {
    if((opts->analysis_type == CM || opts->analysis_type == CM3F) && opts->cm == UPDATED_LAGRANGIAN)
    {
      opts->cm = TOTAL_LAGRANGIAN;
      if(myrank==0)
      {
        PGFEM_printf("Updated Lagrangian is currently unavailable with inertia.\n");
        PGFEM_printf("Forced to Total Lagrangian (-cm = %d)\n", TOTAL_LAGRANGIAN);
      }
    }
  }

  if(read_initial_file)
  {
    if(myrank==0)
    {
      while(fgets(line, 1024, fp)!=NULL)
      {
        if(line[0]=='#')
          continue;

        double temp;
        sscanf(line, "%lf", &temp);
        break;
      }
    }

    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;

      sscanf(line, "%lf", &(sol[mp_id].alpha));
      break;
    }

    // read material density
    double *rho = (double *) malloc(sizeof(double)*mat->nmat);
    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;
      for(int a=0; a<mat->nmat; a++)
      {
        sscanf(line, "%lf", rho+a);
        if(a<mat->nmat-1)
          if (fgets(line, 1024, fp) == nullptr) {
            abort();
          }
      }
      break;
    }

    for(int ia = 0; ia<mat->nhommat; ia++)
    {
      (mat->hommat[ia]).density = rho[(mat->hommat[ia]).mat_id];
      if(myrank==0)
        PGFEM_printf("Density(%d), %e\n", ia, rho[(mat->hommat[ia]).mat_id]);
    }

    free(rho);

    if(opts->restart<0)
    {
      while(fgets(line, 1024, fp)!=NULL)
      {
        if(line[0]=='#')
          continue;

        long nid;
        double u[3], v[3];
        sscanf(line, "%ld %lf %lf %lf %lf %lf %lf", &nid, u+0, u+1, u+2, v+0, v+1, v+2);

        fv[mp_id].u_n[nid*3+0] = u[0];
        fv[mp_id].u_n[nid*3+1] = u[1];
        fv[mp_id].u_n[nid*3+2] = u[2];
        fv[mp_id].u_nm1[nid*3+0] = u[0]-dt*v[0];
        fv[mp_id].u_nm1[nid*3+1] = u[1]-dt*v[1];
        fv[mp_id].u_nm1[nid*3+2] = u[2]-dt*v[2];
      }
    }
    fclose(fp);
  }

  for(long idx_a = 0; idx_a<grid->nn; idx_a++)
  {
    for(long idx_b = 0; idx_b<fv[mp_id].ndofn; idx_b++)
    {
      long id = grid->node[idx_a].id_map[mp_id].id[idx_b];
      if(id>0)
        fv[mp_id].u_np1[id-1] = fv[mp_id].u_n[idx_a*fv[mp_id].ndofn + idx_b];
    }
  }


  return err;
}

/// Read initial conditions for mechanical problem
///
/// # can be used to add comments. The initial condition file should
/// provide Material density as many as number of materials.
/// and disp. in x, y, z velocity x, y, z at t(n=0) followed by node id.
/// e.g
/// # reference temperature
/// 1000.0
/// # initial temperature
/// 0 1.0 1.0 1.0 100.0 0.0 0.0
/// 2 1.0 1.0 1.0 100.0 0.0 0.0
/// :  :
///
/// \param[in] fp file pointer for reading IC for thermal
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in, out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int read_initial_for_Mechanical(FILE *fp,
                                Grid *grid,
                                MaterialProperty *mat,
                                FieldVariables *fv,
                                Solver *sol,
                                LoadingSteps *load,
                                TimeStepping *ts,
                                PGFem3D_opt *opts,
                                const Multiphysics& mp,
                                int myrank,
                                int mp_id)
{
  int err = 0;
  char line[1024];
  double dt = ts->times[1] - ts->times[0];;

  if((opts->analysis_type == CM || opts->analysis_type == CM3F) && opts->cm == UPDATED_LAGRANGIAN)
  {
    opts->cm = TOTAL_LAGRANGIAN;
    if(myrank==0)
    {
      PGFEM_printf("Updated Lagrangian is currently unavailable with inertia.\n");
      PGFEM_printf("Forced to Total Lagrangian (-cm = %d)\n", TOTAL_LAGRANGIAN);
    }
  }

  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;

    sscanf(line, "%lf", &(sol->alpha));
    break;
  }

  // read material density
  double *rho = (double *) malloc(sizeof(double)*mat->nmat);
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;
    for(int a=0; a<mat->nmat; a++)
    {
      sscanf(line, "%lf", rho+a);
      if(a<mat->nmat-1) {
        if (fgets(line, 1024, fp) == nullptr) {
          abort();
        }
      }
    }
    break;
  }

  for(int ia = 0; ia<mat->nhommat; ia++)
  {
    (mat->hommat[ia]).density = rho[(mat->hommat[ia]).mat_id];
    if(myrank==0)
      PGFEM_printf("Density(%d), %e\n", ia, rho[(mat->hommat[ia]).mat_id]);
  }

  free(rho);

  if(opts->restart < 0)
  {
    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;

      long nid;
      double u[3], v[3];
      sscanf(line, "%ld %lf %lf %lf %lf %lf %lf", &nid, u+0, u+1, u+2, v+0, v+1, v+2);

      fv->u_n[nid*3+0] = u[0];
      fv->u_n[nid*3+1] = u[1];
      fv->u_n[nid*3+2] = u[2];
      fv[mp_id].u_nm1[nid*3+0] = u[0]-dt*v[0];
      fv[mp_id].u_nm1[nid*3+1] = u[1]-dt*v[1];
      fv[mp_id].u_nm1[nid*3+2] = u[2]-dt*v[2];
      if(fabs(v[0]) > sol->computer_zero ||
         fabs(v[1]) > sol->computer_zero ||
         fabs(v[2]) > sol->computer_zero)
         fv->apply_initial_velocity = true;
    }
    if(opts->analysis_type == TF)
      compute_3f_initial_conditions(grid, mat, fv);

    if(opts->analysis_type == CM3F || opts->analysis_type == CM3F)
      compute_cm_initial_conditions(grid, mat, fv, load, mp, mp_id, opts->analysis_type);
  }

  for(long idx_a = 0; idx_a<grid->nn; idx_a++)
  {
    for(long idx_b = 0; idx_b<fv->ndofn; idx_b++)
    {
      long id = grid->node[idx_a].id_map[mp_id].id[idx_b];
      if(id>0)
        fv->u_np1[id-1] = fv->u_n[idx_a*fv->ndofn + idx_b];
    }
  }

  return err;
}


/// Read initial conditions for thermal problem
///
/// Initial condition includes reference temperature(T0, default = 300)
/// and temperature at t(n=0) followed by node id. # can be used to add comments.
/// e.g
/// # reference temperature
/// 300
/// # initial temperature
/// 0 310
/// 1 310
/// :  :
///
///
/// \param[in] fp file pointer for reading IC for thermal
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in, out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int read_initial_for_Thermal(FILE *fp,
                             Grid *grid,
                             MaterialProperty *mat,
                             FieldVariables *fv,
                             Solver *sol,
                             TimeStepping *ts,
                             PGFem3D_opt *opts,
                             int myrank,
                             int mp_id)
{
  int err = 0;
  char line[1024];

  double T0 = 300.0;
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;

    sscanf(line, "%lf", &T0);
    fv->u0 = T0; // set reference temperature
    if(myrank==0)
      PGFEM_printf("Default initial temperature: %e\n", T0);

    break;
  }

  if(opts->restart < 0)
  {
    // set default
    for(int ia=0; ia<grid->nn; ia++)
    {
      fv->u_nm1[ia] = T0;
      fv->u_n[ia] = T0;
    }

    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;

      long nid;
      double u;
      sscanf(line, "%ld %lf", &nid, &u);

      fv->u_n[nid] = u;
      fv->u_nm1[nid] = u;
    }
  }

  for(int ia = 0; ia<grid->nn; ia++)
  {
    long id = grid->node[ia].id_map[mp_id].id[0];
    if(id>0)
      fv->u_np1[id-1] = fv->u_n[ia];
  }

  return err;
}

/// Read initial conditions from lagcy format.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[out] tnm1 if restart, read time step info from the previous run
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_initial_values_IC(Grid *grid,
                           MaterialProperty *mat,
                           FieldVariables *FV,
                           Solver *SOL,
                           LoadingSteps *load,
                           TimeStepping *ts,
                           PGFem3D_opt *opts,
                           const Multiphysics& mp,
                           double *tnm1,
                           int myrank)
{
  int err = 0;

  // check restart and read restart values
  if(opts->restart >= 0)
    err += read_restart(grid,FV,ts,load,opts,mp,tnm1,myrank);

  char IC[1024];
  sprintf(IC,"%s/IC",opts->ipath);

  char fn_0[2048], fn[2048];

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    sprintf(fn_0,"%s/%s_0.initial",IC,mp.physicsname[ia]);
    sprintf(fn  ,"%s/%s_%d.initial",IC,mp.physicsname[ia], myrank);

    FILE *fp = NULL;
    fp = fopen(fn, "r");

    if(fp==NULL)
    {
      fp = fopen(fn_0, "r");
      if(fp==NULL)
      {
        if(myrank==0)
          PGFEM_printf("No [%s] exists. Use default ICs.\n", fn_0);

        continue;
      }
    }
    switch(mp.physics_ids[ia])
    {
     case MULTIPHYSICS_MECHANICAL:
      err += read_initial_for_Mechanical(fp,grid,mat,FV+ia,SOL+ia,load,ts,opts,mp,myrank,ia);
      break;
     case MULTIPHYSICS_THERMAL:
      err += read_initial_for_Thermal(fp,grid,mat,FV+ia,SOL+ia,ts,opts,myrank,ia);
      break;
     case MULTIPHYSICS_CHEMICAL:
      break;
     default:
      break;
    }

    fclose(fp);
  }
  return err;
}



void read_simulation_methods(char *filenameMS, PGFem3D_opt *opts) {
    int Nt,Ne,t,e,temp;

  FILE *fpms;
    fpms = fopen(filenameMS,"r");
  if(fpms == NULL)
  {
      PGFEM_printf("Failed to open file [%s]. \n", filenameMS);
  }
  fscanf(fpms,"%d",&Ne);                                    //number of cohesive elements
  fscanf(fpms,"%d",&Nt);                                    //number of time steps
  opts->methods = (int*) PGFEM_calloc (int,Nt*Ne);
      for(t = 0; t < Nt; t++) {
        for(e = 0; e < Ne; e++) {
          fscanf(fpms,"%d",&temp);
          opts->methods[e + Ne*t] = temp;
        }
      }
    fclose(fpms);
}


/// Read initial conditions.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[out] tnm1 if restart, read time step info from the previous run
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_initial_values(Grid *grid,
                        MaterialProperty *mat,
                        FieldVariables *FV,
                        Solver *SOL,
                        LoadingSteps *load,
                        TimeStepping *ts,
                        PGFem3D_opt *opts,
                        const Multiphysics& mp,
                        double *tnm1,
                        int myrank)
{
  int err = 0;
  char IC[1024];
  sprintf(IC,"%s/IC",opts->ipath);

  if(is_directory_exist(IC))
  {
    if(myrank==0)
      PGFEM_printf("IC directory exists, read initial conditions from IC\n");
    err += read_initial_values_IC(grid,mat,FV,SOL,load,ts,opts,mp,tnm1,myrank);
  }
  else
  {
    if(myrank==0)
      PGFEM_printf("No IC directory exists, read inital conditions from *.initial\n");

    err += read_initial_values_lagcy(grid,mat,FV+0,SOL+0,load,ts,opts,mp,tnm1,myrank);
  }
  return err;
}

/// Read loads increments.
/// As time is elapsing, loads increments are read from solver file which
/// file pointer is saved in LoadingSteps. Prior to run this function,
/// read_initial_values function shold be called which open and the solver file pointer.
/// The file pointer will be freed when LoadingSteps object is destoryed.
/// The number of loads increments should be exact as read before in read_initial_values.
///
/// \param[in] grid a mesh object
/// \param[in] fv object for field variables
/// \param[out] load object for loading
/// \param[in] mp multiphysics object
/// \param[in] tim time step ID
/// \param[in] com handle for communication
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_and_apply_load_increments(Grid *grid,
                                   FieldVariables *fv,
                                   LoadingSteps *load,
                                   const Multiphysics& mp,
                                   long tim,
				   const CommunicationStructure *com)
{
  int err = 0;
  int myrank = com->rank;
  
  //  read nodal prescribed boundary values
  for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
  {
    if(load->solver_file[mp_id]==NULL)
      continue;

    if(load->tim_load[mp_id][tim] == 1 && tim == 0)
    {
      if (myrank == 0)
        PGFEM_printf ("Incorrect load input for Time = 0\n");

      PGFEM_Comm_code_abort(com, 0);
    }

    if(load->tim_load[mp_id][tim] == 1 && tim != 0)
    {
      for(long ia=0;ia<load->sups[mp_id]->npd;ia++)
      {
        CHECK_SCANF(load->solver_file[mp_id],"%lf",(load->sups[mp_id])->defl_d + ia);
        (load->sup_defl[mp_id])[ia] = load->sups[mp_id]->defl_d[ia];
      }
      if(mp.physics_ids[mp_id]==MULTIPHYSICS_MECHANICAL)
      {
        // read nodal load in the subdomain
        read_nodal_load(load->solver_file[mp_id],load->nln,grid->nsd,load->znod);
        // read elem surface load */
        read_elem_surface_load(load->solver_file[mp_id],load->nle_s,grid->nsd,grid->element,load->zele_s);
        //  Node - generation of the load vector
        load_vec_node(fv[mp_id].R,load->nln,grid->nsd,load->znod,grid->node,MULTIPHYSICS_MECHANICAL);
        //  ELEMENT - generation of the load vector
        load_vec_elem_sur(fv[mp_id].R,load->nle_s,grid->nsd,grid->element,load->zele_s);
      }
    }
  }

  return err;
}

/// Read read cohesive elements.
/// This function was part of the main function and extracted to
/// Modularized reading cohesive elements and reorganize the main function structures.
/// The function call read not only cohesive elements but also cohesive properties (*.in.co_props).
/// The functin dependency includes lagacy code (read_cohe_elem) to read *.in.co files, too.
///
/// \param[out] grid a mesh object
/// \param[out] mat a material object
/// \param[in] opts structure PGFem3D option
/// \param[in] ensight object
/// \param[in] com handle for communication
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_cohesive_elements(Grid *grid,
                           MaterialProperty *mat,
                           const PGFem3D_opt *opts,
                           Ensight *ensight,
                           const CommunicationStructure *com)
{
  int err = 0;
  int myrank = com->rank;
  char in_dat[1024], filename[2048];

  sprintf(in_dat,"%s/%s",opts->ipath,opts->ifname);

  FILE *fp;
  // read cohesive properties
  sprintf(filename,"%s%d.in.co_props",in_dat,myrank);
  fp = fopen(filename,"r");
  read_cohesive_properties(fp,&(mat->n_co_props),&(mat->co_props),com);
  fclose(fp);

  /* read coheisve elements */
  sprintf(filename,"%s%d.in.co",in_dat,myrank);
  fp = fopen(filename,"r");

  long ncom = 0;
  double **comat = NULL;

  /* temporary leftovers from old file format */
  CHECK_SCANF(fp,"%ld\n",&ncom);

  /* to silence warning message. need to pull this legacy bit of
   * code out completely. Cohesive porperties provided in separate
   * file. This leads to *very* small memory leak */
  if(ncom <= 0)
    comat = aloc2(1,4);
  else
    comat = aloc2(ncom,4);


  /* read the cohesive element info */
  grid->coel = read_cohe_elem (fp,ncom,grid->nsd,grid->nn,grid->node,&(grid->nce),
                               comat,ensight,opts->vis_format,
                               myrank,mat->co_props);
  if(ncom <= 0)
    dealoc2(comat,ncom);
  else
    dealoc2(comat,1);

  fclose (fp);

  /* Global number of cohesive elements */
  com->net->allreduce(&(grid->nce),&(grid->Gnce),1,NET_DT_LONG,NET_OP_SUM,com->comm);
  return err;
}
