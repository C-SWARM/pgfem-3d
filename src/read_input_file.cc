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
#include "three_field_element.h"
#include "utils.h"

#include <cstring>
#include <limits>
#include <sstream>
#include <fstream>

using namespace pgfem3d;
using namespace multiscale::net;

/// read mechanical part of material properties
///
/// \param[in] fp file pointer for reading mechanical part of material properties
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \return non-zero on interal error
static int read_material_for_Mechanical(FILE *fp, MaterialProperty *mat,
                                        const PGFem3D_opt *opts)
{
  for (long i = 0, e = mat->nmat; i < e; ++i) {
    scan_for_valid_line(fp);
    if (read_material(fp, mat->mater[i], opts->legacy)) {
      PGFEM_Abort();
    }
  }
  return 0;
}

/// read mechanical part of material properties
///
/// \param[in] fp file pointer for reading mechanical part of material properties
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \return non-zero on interal error
static int read_material_for_Thermal(FILE *fp, MaterialProperty *mat,
                                     const PGFem3D_opt *opts)
{
  auto thermal = mat->thermal = new MaterialThermal[mat->nmat];

  for (long i = 0, e = mat->nmat; i < e; ++i) {
    scan_for_valid_line(fp);
    CHECK_SCANF(fp, "%lf", &thermal[i].cp);

    scan_for_valid_line(fp);
    auto k = thermal[i].k;
    CHECK_SCANF(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &k[0], &k[1], &k[2], &k[3], &k[4], &k[5], &k[6], &k[7], &k[8]);

    // read optional value (fraction of heat source from mechanical work)
    scan_for_valid_line(fp);
    thermal[i].FHS_MW = 1.0;                    // default
    fscanf(fp, "%lf", &thermal[i].FHS_MW);      // conditionally overwrite
  }

  return 0;
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
  sprintf(dirname, "%s/Material", opts->ipath);

  for(int i = 0, e = mp.physicsno; i < e; ++i) {
    sprintf(fn, "%s/%s.mat", dirname, mp.physicsname[i]);

    if (FILE *fp = fopen(fn, "r")) {
      switch(mp.physics_ids[i]) {
       case MULTIPHYSICS_MECHANICAL:
        err += read_material_for_Mechanical(fp, mat, opts);
        break;
       case MULTIPHYSICS_THERMAL:
        err += read_material_for_Thermal(fp, mat, opts);
        break;
       case MULTIPHYSICS_CHEMICAL:
        break;
       default:
        break;
      }
      fclose(fp);
    }
    else if (myrank == 0) {
      PGFEM_printf("No [%s] exists.\n", fn);
    }
  }

  // read and set general material properties e.g. density
  auto* d = mat->density = PGFEM_malloc<double>(mat->nmat);
  for (long i = 0, e = mat->nmat; i < e; ++i) {
    d[i] = 0.0;
  }

  int match = 0;
  sprintf(fn, "%s/material.mat", dirname);
  if (FILE *fp = fopen(fn, "r")) {
    for (long i = 0, e = mat->nmat; i < e; ++i) {
      scan_for_valid_line(fp);
      match += fscanf(fp, "%lf", &d[i]);
    }
    fclose(fp);
  }
  else {
    if (myrank == 0) {
      PGFEM_printf("No [%s] exists. \nDensity is set to zero.\n", fn);
    }

    return err;
  }

  // check number of densities that is read.
  if (match != mat->nmat) {
    if (myrank == 0) {
      PGFEM_printf("Material density is not read as many as number of "
                   "materials.\n");
    }
    PGFEM_Abort();
  }

  for(long i = 0, e = mat->nhommat; i < e; ++i) {
    mat->hommat[i].density = d[mat->hommat[i].mat_id];
    if (myrank == 0) {
      PGFEM_printf("Density(%d), %e\n", i, mat->hommat[i].density);
    }
  }

  return err;
}


/// read Neumann boundary conditions if NBC directory exists.
/// As NBC is read, NBC object in grid will be created 
///
/// \param[in,out] grid   mesh object, NBC member will be updated
/// \param[in]     mp     multiphysics object
/// \param[in]     opts   PGFem3D options
/// \param[in]     myrank current process rank
void read_Neumann_boundary_conditions(Grid *grid,
                                      const Multiphysics &mp,
                                      const PGFem3D_opt *opts,
                                      const int myrank){
  std::stringstream ss_nbe_dir;
  ss_nbe_dir << opts->ipath << "/" << "NBC";
  if(is_directory_exist(ss_nbe_dir.str().c_str())){
    if(myrank==0)
      PGFEM_printf("NBC directory exists, read Neumann boundary conditions from NBC\n");

    for(int ia = 0, e = mp.physicsno; ia < e; ++ia) {
      std::stringstream fn;
      fn << ss_nbe_dir.str() << "/" << mp.physicsname[ia] << ".nbc";
      
      std::ifstream ifs;
      ifs.open(fn.str());

      if (ifs.is_open()) {
        // read number of features and allocate data size
        ifs >> grid->NBE(ia).feature_no;
        
        if(grid->NBE(ia).feature_no > 0){
          grid->NBE(ia).set_feature_size(grid->NBE(ia).feature_no, mp.ndim[ia]);

          for(int ib=0, nf = grid->NBE(ia).feature_no; ib<nf; ib++){
            ifs >> grid->NBE(ia).features(ib, 0); // features(ib, 0): T3D feature type
            ifs >> grid->NBE(ia).features(ib, 1); // features(ib, 1): T3D feature id
            ifs >> grid->NBE(ia).load_type(ib);   // type is same as number of load values

            for(int ic=0; ic<grid->NBE(ia).load_type(ib); ++ic){
              std::string load_func_fn; // this is filename of individual load 
              ifs >> load_func_fn;
              // math expression in load_func_file is updated and saved in grid->NBE.load
              grid->NBE(ia).read_loads(ss_nbe_dir.str(), load_func_fn, ib, ic);
            }
          }
          
          // construct list of elemnets contributing for NBC
          grid->NBE(ia).construct_list_of_boundary_elements(grid->element, grid->ne, grid->nsd);

        }
      }
      else if (myrank == 0) {
        std::cout << "No " << fn.str() << "] exists." << endl;
      }
    }
  }
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

  *node = build_node_multi_physics(*nn,fv_ndofn,physicsno);
  *elem = build_elem(in,*ne,opts->analysis_type);
  *material = PGFEM_calloc(Material, *nmat);
  *matgeom = build_matgeom(*n_concentrations,*n_orient);

  *Gnn = read_nodes(in,*nn,*node,opts->legacy,com);
  /* NOTE: Supports assume only ndim supported dofs per node! */

  char BC[1024];
  sprintf(BC, "%s/BC", opts->ipath);

  if (is_directory_exist(BC)) {
    if (myrank == 0) {
      PGFEM_printf("BC exists skip BC from filebase_*.in instead read boundary "
                   "conditions from BC\n");
    }

    // skip some input data along this path
    int nbc;
    int n[4];
    CHECK_SCANF(in, "%d", &nbc);
    for (int i = 0; i < nbc; ++i) {
      CHECK_SCANF(in, "%d %d %d %d", n+0, n+1, n+2, n+3);
    }

    CHECK_SCANF(in, "%d", &nbc);
    double v;
    for (int i = 0; i < nbc; ++i) {
      CHECK_SCANF(in, "%lf", &v);
    }

    // read boundary conditions if BC diretory exists
    for (int i = 0; i < physicsno; ++i) {
      char fn_bc[2048];
      char fn_bcv[2048];
      sprintf(fn_bc, "%s/%s_%d.bc", BC, physicsnames[i], myrank);
      sprintf(fn_bcv, "%s/%s.bcv", BC, physicsnames[i]);

      // sup[i] is allocated in read_Dirichlet_BCs such thtat
      // read_Dirichlet_BCs function must be called even though fp is NULL.
      FILE *fp = fopen(fn_bc, "r");
      sup[i] = read_Dirichlet_BCs(fp, *nn, ndim[i], *node, i);

      if(fp!=NULL){
        fclose(fp);
      }

      if (FILE *fp = fopen(fn_bcv, "r")) {
        err += read_Dirichlet_BCs_values(fp, *nn, ndim[i], *node, sup[i], i);
        fclose(fp);
      }
      else {
        PGFEM_printf("ERROR: Cannot open %s file. Exit.\n", fn_bcv);
        PGFEM_Abort();
      }
    }
  }
  else {
    for (int i = 0; i < physicsno; ++i) {
      sup[i] = read_supports(in, *nn, ndim[i], *node, i);
    }
  }

  read_elem(in, *ne, *elem, *sup, opts->legacy);
  for (long i = 0, e = *nmat; i < e; ++i) {
    if (read_material(in, (*material)[i], opts->legacy)) {
      PGFEM_Abort();
    }
  }

  if (feof(in)) {
    PGFEM_printerr("ERROR: prematurely reached EOF in %s(%s)\n",
                   __func__,__FILE__);
    PGFEM_Abort();
  }

  if (override_material_properties(*nmat, opts, *material)) {
    PGFEM_Abort();
  }

  read_matgeom(in, *n_concentrations, *n_orient, *matgeom);

  /* NOTE: Node/Element loading assumes forces only in ndim
     directions */
  CHECK_SCANF(in, "%ld", nln);
  assert(0 <= *nln);
  *znod = build_zatnode(*ndofn, *nln);
  read_nodal_load(in, *nln, *ndofn, *znod);
  /* surface */
  CHECK_SCANF (in,"%ld", nel_s);
  assert(0 <= *nel_s);
  *zelem_s = build_zatelem(*ndofn, *nel_s);
  read_elem_surface_load(in, *nel_s, *ndofn, *elem, *zelem_s);

  /* volume */
  CHECK_SCANF (in, "%ld", nel_v);
  assert(0 <= *nel_v);
  *zelem_v = build_zatelem(*ndofn, *nel_v);

  if (ferror(in)) err++;
  PGFEM_free(filename);
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
  // Copies the ndofn from each multiphysics structure into a continuous array.
  std::vector<int> fv_ndofn(mp.physicsno);
  for (int i = 0, e = mp.physicsno; i < e; ++i) {
    assert(0 < FV[i].ndofn and FV[i].ndofn < std::numeric_limits<int>::max());
    fv_ndofn[i] = FV[i].ndofn;
  }

  long ndofn = 0;                               // output from read_input_file
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
                            &fv_ndofn[0],
                            mp.physicsno,
                            mp.ndim,
                            mp.physicsname);

  // read multiphysics material properties
  int myrank = com->rank;
  err += read_multiphysics_material_properties(mat, opts, mp, myrank);
  
  read_Neumann_boundary_conditions(grid, mp, opts, myrank);

  // update numerical solution scheme parameters
  FV[0].NORM = SOL[0].computer_zero;
  for (int i = 0, e = mp.physicsno; i < e; ++i) {
    SOL[i].iter_max_sol    = SOL[0].iter_max_sol;
    SOL[i].err             = SOL[0].err;
    SOL[i].computer_zero   = SOL[0].computer_zero;
    FV[i].n_concentrations = FV[0].n_concentrations;
    FV[i].NORM             = FV[0].NORM;
  }

  // need to update number of elements that have prescribed BCs (supported)
  for (long i = 0, e = grid->ne; i < e; ++i) {
    auto nod = grid->element[i].nod;
    auto nne = grid->element[i].toe;

    // @todo[ld] This isn't really ideal. In reality, it would be nice to
    //           provide some functionality in the element and/or SUPP to do
    //           this operation, rather than embedding it as a lambda.
    auto is_supported = [&](const SUPP_1& s) {
      for (long i = 0, e = s.ndn; i < e; ++i) {
        for (long j = 0, e = nne; j < e; ++j) {
          if (s.lnpd[i] == nod[j]) {
            return 1;
          }
        }
      }
      return 0;
    };

    // j = 0 is alreaded accounted in read_elem in read_input_file
    for (int j = 1, e = mp.physicsno; j < e; ++j) {
      load->sups[j]->nde += is_supported(*load->sups[j]);
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

int read_solver_file_multiscale(MultiscaleCommon *c,
                                MULTISCALE_SOLUTION *s,
                                SOLVER_FILE *solver_file,
                                const PGFem3D_opt *opts,
                                const int myrank)
{
  int err    = 0;
  int n_step = 0;
  double pores = 0.0;
  double *sup_defl = NULL;
  CRPL *crpl = NULL;
  s->tim = 0;

  // initialize and define multiphysics
  Multiphysics mp;

  int id = MULTIPHYSICS_MECHANICAL;
  int ndim = c->ndofn;
  int MECHANICAL_Var_NO = 0;
  int write_no = MECHANICAL_Var_NO;

  mp.write_ids.resize(MECHANICAL_Var_NO);
  vector<int> coupled_ids;
  char *physicsname = (char *) malloc(sizeof(char)*1024);
  {
    for(int ia=0; ia<MECHANICAL_Var_NO; ia++)
      mp.write_ids[ia].push_back(ia); //sets ia as the fist element of each vector

    coupled_ids.push_back(0);
    sprintf(physicsname, "Mechanical");

    mp.physicsno      = 1;
    mp.physicsname    = &physicsname;
    mp.physics_ids    = &id;
    mp.ndim           = &ndim;
    mp.write_no       = &write_no;
    mp.coupled_ids.push_back(coupled_ids);
    mp.total_write_no = MECHANICAL_Var_NO;
  }

  // initialize and define mesh object
  Grid grid;
  grid_initialization(&grid);
  {
    grid.ne          = c->ne;
    grid.nn          = c->nn;
    grid.element     = c->elem;
    grid.b_elems     = NULL;
    grid.node        = c->node;
    grid.nce         = c->nce;
    grid.coel        = c->coel;
  }

  // initialize and define field variables
  FieldVariables fv;
  {
    field_varialbe_initialization(&fv);
    fv.ndofn  = c->ndofn;
    fv.ndofd  = c->ndofd;
    fv.npres  = c->npres;
    fv.sig    = s->sig_e;
    fv.eps    = s->eps;
    fv.u_np1  = s->r;
    fv.f      = s->f;
    fv.d_u    = s->d_r;
    fv.dd_u   = s->rr;
    fv.R      = s->R;
    fv.f_defl = s->f_defl;
    fv.RR     = s->RR;
    fv.f_u    = s->f_u;
    fv.RRn    = s->RRn;
    fv.pores  = pores;
    fv.BS_x   = s->BS_x;
    fv.BS_f   = s->BS_f;
    fv.BS_RR  = s->BS_RR;
    fv.BS_f_u = s->BS_f_u;
    fv.NORM   = s->NORM;
  }

  /// initialize and define iterative solver object
  Solver sol{};
  {
    if(solver_file==NULL)
    {
      sol.nor_min  = c->lin_err;
      sol.FNR      = 1; // full NR
      sol.iter_max = c->maxit_nl;
    }
    else
    {
      sol.nor_min  = solver_file->nonlin_tol;
      sol.FNR      = solver_file->nonlin_method;
      sol.iter_max = solver_file->max_nonlin_iter;
    }
    sol.n_step     = n_step;
    sol.system     = c->SOLVER;
    sol.err        = c->lin_err;
  }

  // initialize and define loading steps object
  LoadingSteps load;
  {
    loading_steps_initialization(&load);
    load.sups     = &(c->supports);
    load.sup_defl = &sup_defl;
    load.solver_file = &solver_file->file;
    load.tim_load    = &solver_file->load_steps;
  }

  // initialize and define material properties
  MaterialProperty mat;
  {
    material_initialization(&mat);
    mat.hommat  = c->hommat;
    mat.matgeom = c->matgeom;
  }

  /// initialize and define time stepping variable
  TimeStepping ts;
  {
    time_stepping_initialization(&ts);
    if(solver_file==NULL)
      ts.nt = 1;
    else
      ts.nt  = solver_file->n_step;

    ts.tim = s->tim;
    ts.times = s->times;
    ts.dt_n   = 0.0;
    ts.dt_np1 = s->dt;
    ts.print  = NULL;
    ts.tns    = NULL;
  }

  err += read_solver_file(&ts, &mat, &fv, &sol, &load, crpl, mp, opts, myrank);

  /* assign the returned values */
  c->maxit_nl = sol.iter_max;
  c->nce = grid.nce;

  solver_file->nonlin_tol           = sol.nor_min;
  solver_file->nonlin_method        = sol.FNR;
  solver_file->max_nonlin_iter      = sol.iter_max;
  solver_file->n_pressure_nodes     = fv.npres;
  solver_file->nonlin_method        = sol.FNR;
  solver_file->n_step               = ts.nt;
  solver_file->times                = ts.times;
  solver_file->print_steps          = ts.print;

  /* nonlinear solving method options  */
  switch(solver_file->nonlin_method){
   case NEWTON_METHOD:
    break;

   case ARC_LENGTH_METHOD:
   case AUX_ARC_LENGTH_METHOD:
    /* additional options */
    solver_file->nonlin_method_opts[0] = sol.arc->dAL0;
    //solver_file->nonlin_method_opts[1] = sol.arc->dALMAX;
    break;

   case MULTIPHYSICS_NEWTON_METHOD:
    break;

   case TAYLOR_MODEL:
    break;      
  }

  free(physicsname);

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
      for(long a=0; a<mat->nmat; a++)
      {
        sscanf(line, "%lf", rho+a);
        if(a<mat->nmat-1)
          if (fgets(line, 1024, fp) == nullptr) {
            abort();
          }
      }
      break;
    }

    for(long ia = 0; ia<mat->nhommat; ia++)
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
/// \param[in] read_from_0 if ture, read information from filebase_0.intial
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
                                int mp_id,
                                const bool read_from_0)
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
    for(long a=0; a<mat->nmat; a++)
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

  for(long ia = 0; ia<mat->nhommat; ia++)
  {
    (mat->hommat[ia]).density = rho[(mat->hommat[ia]).mat_id];
    if(myrank==0)
      PGFEM_printf("Density(%d), %e\n", ia, rho[(mat->hommat[ia]).mat_id]);
  }

  free(rho);

  if(!read_from_0 && opts->restart < 0){
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
/// \param[in] read_from_0 if ture, read information from filebase_0.intial
/// \return non-zero on internal error
int read_initial_for_Thermal(FILE *fp,
                             Grid *grid,
                             MaterialProperty *mat,
                             FieldVariables *fv,
                             Solver *sol,
                             TimeStepping *ts,
                             PGFem3D_opt *opts,
                             int myrank,
                             int mp_id,
                             const bool read_from_0)
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
    for(long ia=0; ia<grid->nn; ia++)
    {
      fv->u_nm1[ia] = T0;
      fv->u_n[ia] = T0;
    }

    if(!read_from_0){

      while(fgets(line, 1024, fp)!=NULL){
        if(line[0]=='#')
          continue;

        long nid;
        double u;
        sscanf(line, "%ld %lf", &nid, &u);

        fv->u_n[nid] = u;
        fv->u_nm1[nid] = u;
      }
    }
  }

  for(long ia = 0; ia<grid->nn; ia++)
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
    bool read_from_0 = false;

    if(fp==NULL)
    {
      fp = fopen(fn_0, "r");
      read_from_0 = true;
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
      err += read_initial_for_Mechanical(fp,grid,mat,FV+ia,SOL+ia,load,ts,opts,mp,myrank,ia,read_from_0);
      break;
     case MULTIPHYSICS_THERMAL:
      err += read_initial_for_Thermal(fp,grid,mat,FV+ia,SOL+ia,ts,opts,myrank,ia,read_from_0);
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
  CHECK_SCANF(fpms,"%d",&Ne);                   //number of cohesive elements
  CHECK_SCANF(fpms,"%d",&Nt);                   //number of time steps
  opts->methods = (int*) PGFEM_calloc (int,Nt*Ne);
      for(t = 0; t < Nt; t++) {
        for(e = 0; e < Ne; e++) {
          CHECK_SCANF(fpms,"%d",&temp);
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

