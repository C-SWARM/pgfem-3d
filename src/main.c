
/* HEADER */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "PFEM3d.h"
#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

/* Standard headers/libs */
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>

/*=== PFEM3d headers ===*/
#include "PGFEM_io.h"
#include "allocation.h"
#include "Arc_length.h"
#include "build_distribution.h"
#include "homogen.h"
#include "hypre_global.h"
#include "in.h"
#include "load.h"
#include "matice.h"
#include "matrix_printing.h"
#include "Newton_Raphson.h"
#include "out.h"
#include "Printing.h"
#include "print_dist.h"
#include "profiler.h"
#include "Psparse_ApAi.h"
#include "read_cryst_plast.h"
#include "renumber_ID.h"
#include "set_fini_def.h"
#include "skyline.h"
#include "utils.h"
#include "SetGlobalNodeNumbers.h"
#include "interface_macro.h"
#include "computeMacroF.h"
#include "computeMacroS.h"
#include "vtk_output.h"
#include "PGFem3D_options.h"
#include "gen_path.h"
#include "initialize_damage.h"
#include "bounding_element.h"
#include "bounding_element_utils.h"
#include "element.h"
#include "node.h"
#include "generate_dof_ids.h"
#include "applied_traction.h"
#include "read_input_file.h"

#include "three_field_element.h"
#include "constitutive_model.h"
#include "post_processing.h"
#include "restart.h"

#include "comm_hints.h"
#include "fd_residuals.h"
#include "dynamics.h"

static const int periodic = 0;

/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/
#define SAVE_RESTART_FILE 1

/// Print PGFem3D running options
/// This will print out mesh info, analysis options
///
/// \param[in] argc number of arguments passed through command line
/// \param[in] argv arguments passed through command line
/// \param[in] fv field variables (FIELD_VARIABLES object)
/// \param[in] grid mesh info (GRID object)
/// \param[in] com commuincation info (COMMUNICATION_STRUCTURE object)
/// \param[in] load info for loading steps (LOADING_STEPS object)
/// \param[in] gem flag for Generalized finite element method
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int print_PGFem3D_run_info(int argc,char *argv[],
                           GRID *grid,
                           COMMUNICATION_STRUCTURE *com,
                           LOADING_STEPS *load,
                           long gem,
                           PGFem3D_opt *opts)
{
  int err = 0;

  PrintTitleV1();
  for (int ia = 0; ia < argc; ia++)
    PGFEM_printf("%s ",argv[ia]);
    
  PGFEM_printf("\n\n");
    
  switch(opts->analysis_type)
  {
    case ELASTIC:
      PGFEM_printf ("ELASTIC ANALYSIS\n");
      break;
    case TP_ELASTO_PLASTIC:
      PGFEM_printf ("TWO PHASE COMPOSITE SYSTEM : ELASTO-PLASTIC ANALYSIS\n");
      break;
    case FS_CRPL:
      PGFEM_printf ("FINITE STRAIN CRYSTAL ELASTO-PLASTICITY\n");
      break;
    case FINITE_STRAIN:
      if (opts->cohesive == 0) {
        PGFEM_printf ("FINITE STRAIN ELASTICITY\n");
      } else {
        PGFEM_printf ("FINITE STRAIN ELASTICITY WITH COHESIVE FRACTURE\n");
      }
      break;
    case STABILIZED:
      if (opts->cohesive == 0 && gem == 0) {
        PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION : stb = %12.5e\n",
                opts->stab);
      } else if( opts->cohesive == 1) {
        PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION"
                " WITH COHESIVE FRACTURE : stb = %12.5e\n",
                opts->stab);
      } else if ( gem == 1) {
        PGFEM_printf ("GENERALIZED FINITE ELEMENT METHOD\n");
        PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION : stb = %12.5e\n",
                opts->stab);
      }
      break;
    case MINI:
      PGFEM_printf("FINITE STRAIN HYPERELASTICITY W/ MINI ELEMENT\n");
      break;
    case MINI_3F:
      PGFEM_printf("FINITE STRAIN HYPERELASTICITY W/ MINI 3 FIELD ELEMENT\n");
      break;
    case DISP:
      PGFEM_printf("FINITE STRAIN DAMAGE HYPERELASTICITY:\n"
              "TOTAL LAGRANGIAN DISPLACEMENT-BASED ELEMENT\n");
      break;
    case TF:
      PGFEM_printf("FINITE STRAIN TREE FIELDS HYPERELASTICITY:\n"
              "TOTAL LAGRANGIAN TREE FIELDS-BASED ELEMENT\n");
      break;
    case CM:
    {
      PGFEM_printf("USE CONSTITUTIVE MODEL INTERFACE: ");
      switch(opts->cm)
      {
        case UPDATED_LAGRANGIAN:
          PGFEM_printf("UPDATED LAGRANGIAN\n");
          break;
        case TOTAL_LAGRANGIAN:
          PGFEM_printf("TOTAL LAGRANGIAN\n");
          break;
        case MIXED_ANALYSIS_MODE:
          PGFEM_printf("MIXED ANALYSIS MODE\n");
          break;
        default:
          PGFEM_printf("UPDATED LAGRANGIAN\n");
          break;
      }
      break;
    }
    default:
      PGFEM_printerr("ERROR: unrecognized analysis type!\n");
      PGFEM_Abort();
      break;
  }
  
  if(opts->multi_scale){
    if((load->sups[MULTIPHYSICS_MECHANICAL])->npd>= 9){
      PGFEM_printf("*** BULK Multiscale Modelling ***\n");
    } else {
      PGFEM_printf("*** INTERFACE Multiscale Modelling ***\n");
    }
  }
  
  PGFEM_printf ("\n");
  PGFEM_printf ("SolverPackage: ");
  assert(opts->solverpackage == HYPRE);
  switch(opts->solver){
    case HYPRE_GMRES: PGFEM_printf ("HYPRE - GMRES\n"); break;
    case HYPRE_BCG_STAB: PGFEM_printf ("HYPRE - BiCGSTAB\n"); break;
    case HYPRE_AMG: PGFEM_printf ("HYPRE - BoomerAMG\n"); break;
    case HYPRE_FLEX: PGFEM_printf ("HYPRE - FlexGMRES\n"); break;
    case HYPRE_HYBRID: PGFEM_printf ("HYPRE - Hybrid (GMRES)\n"); break;
    default:
      PGFEM_printerr("Unrecognized solver package!\n");
      PGFEM_Abort();
      break;
  }
  
  PGFEM_printf("Preconditioner: ");
  switch(opts->precond){
    case PARA_SAILS: PGFEM_printf ("HYPRE - PARASAILS\n"); break;
    case PILUT: PGFEM_printf ("HYPRE - PILUT\n"); break;
    case EUCLID: PGFEM_printf ("HYPRE - EUCLID\n"); break;
    case BOOMER: PGFEM_printf ("HYPRE - BoomerAMG\n"); break;
    case NONE: PGFEM_printf ("PGFEM3D - NONE\n"); break;
    case DIAG_SCALE: PGFEM_printf ("PGFEM3D - DIAGONAL SCALE\n"); break;
    case JACOBI: PGFEM_printf ("PGFEM3D - JACOBI\n"); break;
  }
  PGFEM_printf ("\n");
  PGFEM_printf ("Number of total nodes                    : %ld\n", grid->Gnn);
  PGFEM_printf ("Number of nodes on domain interfaces     : %ld\n", com->NBN);
  PGFEM_printf ("Total number of elements                 : %ld\n", grid->Gne);
  PGFEM_printf ("Number of elems on the COMM interfaces   : %ld\n", grid->Gnbndel);
  PGFEM_printf ("Total number of bounding (surf) elems    : %d\n",  grid->Gn_be);

  return err;
}

/// Print PGFem3D run time info
/// This will print out how much time takes for the finite element analysis
///
/// \param[in] total_time total simulation time
/// \param[in] hypre_time time took for linear system solver
/// \param[in] usage detailed time info
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int print_PGFem3D_final(double total_time,
                        double hypre_time,
                        struct rusage *usage,
                        int myrank)
{
  int err = 0;
  getrusage(RUSAGE_SELF, usage);
  PGFEM_printf("\n");
  PGFEM_printf("Time of analysis on processor [%d] - "
          " System %ld.%ld, User %ld.%ld.\n",
          myrank,(usage->ru_stime).tv_sec, (usage->ru_stime).tv_usec,
          (usage->ru_utime).tv_sec,(usage->ru_utime).tv_usec);
  
  PGFEM_printf("Total HYPRE solve time on processor [%d] - %f\n",
          myrank,hypre_time);
  PGFEM_printf("Total time (no MPI_Init()) - %f\n",total_time);
  return err;
}

int single_scale_main(int argc,char *argv[])
{
  int err = 0; 
  int mp_id_M = MULTIPHYSICS_MECHANICAL;
  int mp_id_T = MULTIPHYSICS_THERMAL;
  int mp_id_C = MULTIPHYSICS_CHEMICAL;
  
  static const int ndim = 3;
  /* Create MPI communicator. Currently aliased to MPI_COMM_WORLD but
   * may change */
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  //----------------------------------------------------------------------
  // create and initialization of PGFem3D objects
  //----------------------------------------------------------------------
  //---->
  // Multiphysics
  MULTIPHYSICS mp;
  int physicsno = 2;
  err += multiphysics_initialization(&mp);
  err += construct_multiphysics(&mp, physicsno);
  
  err += set_a_physics(&mp, 0, MULTIPHYSICS_MECHANICAL, 3, "Mechanical");
  err += set_a_physics(&mp, 1, MULTIPHYSICS_THERMAL,    1, "Thermal");    
  
  FIELD_VARIABLES *fv          =         (FIELD_VARIABLES *) malloc(physicsno*sizeof(FIELD_VARIABLES));
  SOLVER_OPTIONS  *sol         =          (SOLVER_OPTIONS *) malloc(physicsno*sizeof(SOLVER_OPTIONS));
  COMMUNICATION_STRUCTURE *com = (COMMUNICATION_STRUCTURE *) malloc(physicsno*sizeof(COMMUNICATION_STRUCTURE));

  for(int ia = 0; ia<physicsno; ia++)
  {
    err += field_varialbe_initialization(fv+ia);
    fv[ia].ndofn = mp.ndim[ia];    
    err += solution_scheme_initialization(sol+ia);
    err += communication_structure_initialization(com+ia);      
  }
    
  GRID grid;
  MATERIAL_PROPERTY mat;
  PGFem3D_TIME_STEPPING time_steps;
  LOADING_STEPS load;  
  ARC_LENGTH_VARIABLES arc;  

  err += time_stepping_initialization(&time_steps);
  err += grid_initialization(&grid); // grid.nsd = 3 is the default
  err += material_initialization(&mat);
  err += loading_steps_initialization(&load);
  load.sups = (SUPP *) malloc(sizeof(SUPP)*physicsno);
  err += arc_length_variable_initialization(&arc);

  //<---------------------------------------------------------------------
    
  struct rusage usage;

  double VVolume = 0.0;

  Model_parameters *param_list = NULL;

  char filename[500],in_dat[500],out_dat[500];
  
  /* boundary elements */
  
  long gem = 0;
  
  /* CRYSTAL PLASTICITY */
  CRPL *crpl = NULL;
  
  /* MI */
  double GVolume = 0.0;
  int namelen = 0;

  char processor_name[MPI_MAX_PROCESSOR_NAME];
      
  /* Ensight */
  ENSIGHT ensight;
  
  /* original volume */
  double oVolume = 0.0;
  
  double hypre_time = 0.0;
    
  /* ***** Set up debug log ***** */
  FILE *debug_log = NULL;
  /* debug_log = fopen("debug.log","w"); */
  debug_log = stdout;
  /* debug_log = stderr; */
  
  int myrank = 0;
  
  /*=== END INITIALIZATION === */
  
  int flag_MPI_Init;
  MPI_Initialized(&flag_MPI_Init);
  if(!flag_MPI_Init)
  {
    MPI_Init (&argc,&argv);
  }
  
  MPI_Comm_rank (mpi_comm,&myrank);
  MPI_Comm_size (mpi_comm,&(com[0].nproc)); 
  for(int ia=1; ia<physicsno; ia++)    
    com[ia].nproc = com[0].nproc;
    
  MPI_Get_processor_name (processor_name,&namelen);
  PGFEM_initialize_io(NULL,NULL);
  
  if(myrank == 0){
    PGFEM_printf("=== SINGLE SCALE ANALYSIS ===\n\n");
  }
  
  double total_time = 0.0;
  /* MPI_Barrier(mpi_comm); */
  total_time -= MPI_Wtime();
  
  MPI_Errhandler_set(mpi_comm,MPI_ERRORS_ARE_FATAL);
  
  if(myrank == 0) {
    PGFEM_printf("\n\nInitializing PFEM3d\n\n");
  }
  
  
  /*=== Parse the command line for options ===*/
  PGFem3D_opt options;
  if (argc <= 2){
    if(myrank == 0){
      print_usage(stdout);
    }
    exit(0);
  }
  set_default_options(&options);
  re_parse_command_line(myrank,2,argc,argv,&options);
  if(myrank == 0){
    print_options(stdout,&options);
  }
  
  /* set up solver variables */
  if(options.solverpackage == HYPRE)
  {
    for(int ia=0; ia<physicsno; ia++)
    {
      initialize_PGFEM_HYPRE_solve_info(&(sol[ia].PGFEM_hypre));
      (sol[ia].PGFEM_hypre)->solver_type = options.solver;
      (sol[ia].PGFEM_hypre)->precond_type = options.precond;
    }
  } else {
    if(myrank == 0){
      PGFEM_printerr("ERROR: Only HYPRE solvers are supported.\n");
    }
    PGFEM_Comm_code_abort(mpi_comm,0);
  }
  
  /* for attaching debugger */
  while (options.debug);
  
  /* visualization */
  switch(options.vis_format){
    case VIS_ENSIGHT: case VIS_VTK:
      ensight = PGFEM_calloc (1,sizeof(ENSIGHT_1));
      break;
    default: ensight = NULL; break;
  }
  
  /* abort early if unrecognized analysis type */
  if(options.analysis_type < 0
          || options.analysis_type >= ANALYSIS_MAX){
    if(myrank == 0){
      PGFEM_printerr("ERROR: Unregognized analysis type given (%d)!"
              " Please provide an analysis type (see the help menu).\n",
              options.analysis_type);
    }
    PGFEM_Abort();
  }
  
  if(options.restart > -1){
    if(myrank == 0)
      PGFEM_printerr("Restart from step number :%d.\n", options.restart);
  }
  
  sprintf(in_dat,"%s/%s",options.ipath,options.ifname);
  sprintf(out_dat,"%s/%s",options.opath,options.ofname);
  
  if(make_path(options.opath,DIR_MODE) != 0){
    if(myrank == 0){
      PGFEM_printf("Could not create path (%s)!\n"
              "Please check input and try again.\n\n",options.opath);
      print_usage(stdout);
    }
    PGFEM_Comm_code_abort(mpi_comm,-1);
  }

  //----------------------------------------------------------------------
  // read main input files ( *.in)
  //----------------------------------------------------------------------
  //---->  
  {
    int in_err = 0;            
    in_err = read_mesh_file(&grid,&mat,fv+0,sol+0,&load,&mp,mpi_comm,&options);
            
    if(in_err){
      PGFEM_printerr("[%d]ERROR: incorrectly formatted input file!\n",
              myrank);
      PGFEM_Abort();
    }    
  }
  //<---------------------------------------------------------------------
    
  /*=== READ COMM HINTS ===*/
  {
    char *fn = Comm_hints_filename(options.ipath, options.ifname, myrank);
    com[0].hints = Comm_hints_construct();
    int ch_err = Comm_hints_read_filename(com[0].hints, fn);
    MPI_Allreduce(MPI_IN_PLACE, &ch_err, 1, MPI_INT, MPI_SUM, mpi_comm);
    if (ch_err) {
      Comm_hints_destroy(com[0].hints);
      com[0].hints = NULL;
      if (myrank == 0) {
        PGFEM_printerr("WARNING: One or more procs could not load communication hints.\n"
                "Proceeding using fallback functions.\n");
      }
    }
    free(fn);
  }
  
  // use commuincation hints build at 0 for other physics (>0) 
  // memory will be deallocated once by checking is it NULL
  for(int ia=1; ia<physicsno; ia++)
    com[ia].hints = com[0].hints;
  
  for(int ia=0; ia<physicsno; ia++)
  {
    if(mp.physics_ids[ia]!=mp_id_M)
      continue;
      
    /*=== OVERRIDE PRESCRIBED DISPLACEMENTS ===*/
    if(options.override_pre_disp){
      if(override_prescribed_displacements(load.sups[mp_id_M],&options) != 0){
        PGFEM_printerr("[%d]ERROR: an error was encountered when"
                " reading the displacement override file.\n"
                "Be sure that there are enough prescribed"
                " displacements in the file.\n",myrank);
        PGFEM_Abort();
      }
    }
  
    /*=== MULTISCALE INFORMATION ===*/
    if(options.multi_scale){
      (load.sups[mp_id_M])->multi_scale = options.multi_scale;
      int ms_err = read_interface_macro_normal_lc(options.ipath,load.sups[mp_id_M]);
      if(ms_err != 0){
        PGFEM_printerr("[%d] ERROR: could not read normal from file!\n"
                "Check that the file \"%s/normal.in\""
                " exists and try again.\n",
                myrank,options.ipath);
        PGFEM_Abort();
      }
    }
  }
  
  /*=== BOUNDING ELEMENTS ===*/
  /* NOTE: These might be ripped out... */
  /* ADDED/tested 12/18/2012 MM */
  {
    char bnd_file[500];
    sprintf(bnd_file,"%s%d.in.bnd",in_dat,myrank);
    read_bounding_elements_fname(bnd_file,3,&(grid.n_be),&(grid.b_elems),mpi_comm);
    bounding_element_set_local_ids(grid.n_be,grid.b_elems,grid.element);
    bounding_element_reverse_mapping(grid.n_be,grid.b_elems,grid.element);
  }
  
  /*==== ADDITIONAL SETUP ===*/
  
  /* list of elements with prescribed deflection */
  list_el_prescribed_def(load.sups[mp_id_M],grid.node,grid.element,grid.b_elems,grid.ne,grid.n_be,grid.nn);
  
  /* list of elements on the COMMUNICATION boundary */
  //build for 0
  com[0].nbndel = 0;
  com[0].bndel = list_boundary_el(grid.ne,grid.element,grid.nn,grid.node,myrank,&(com[0].nbndel));

  for(int ia=0; ia<physicsno; ia++)
  {
    
    if(mp.physics_ids[ia]==mp_id_M)
    {
      //material matrices (Mechanical part) of the phases
      Mat_3D_orthotropic (mat.nmat,mat.mater,options.analysis_type);
  
      long ***a = NULL;
      a = aloc3l (mat.nmat,mat.nmat,fv[mp_id_M].n_concentrations);
      mat.nhommat = list (a,grid.ne,mat.nmat,fv[mp_id_M].n_concentrations,grid.element);
  
      //alocation of the material matrices
      mat.hommat = build_hommat (mat.nhommat);
  
      //creates material matrices of the homogeneous medium : LOCAL
      // COORDINATE SYSTEM
      hom_matrices (a,grid.ne,mat.nmat,fv[mp_id_M].n_concentrations,grid.element,mat.mater,mat.matgeom,
              mat.hommat,mat.matgeom->SH,options.analysis_type);
  
      dealoc3l(a,mat.nmat,mat.nmat);
    }
    
    // use commuincation boundary info for other physic using one build on 0
    // memory will be deallocated once by checking is it NULL    
    if(ia>0)
    {  
      com[ia].bndel  = com[0].bndel;
      com[ia].nbndel = com[0].nbndel;
    }

    // Create Graph for communication
    // GrComm = CreateGraph(nproc,myrank,nn,node);
    com[ia].DomDof = aloc1l(com[ia].nproc);

  }

  long *DomNe = aloc1l(com[0].nproc);
  long *DomNn = aloc1l(com[0].nproc);
  DomNe[myrank] = grid.ne;
  DomNn[myrank] = grid.nn;
  
  // Read cohesive elements
  if(options.cohesive == 1)
    err += read_cohesive_elements(&grid,&mat, &options, ensight, mpi_comm, myrank);

  /* Gather number of element from all domains */
  MPI_Gather (&(grid.ne),1,MPI_LONG,DomNe,1,MPI_LONG,0,mpi_comm);
  
  /* Total number of boundary elements */
  MPI_Reduce(&(com[0].nbndel),&(grid.Gnbndel),1,MPI_LONG,MPI_SUM,0,mpi_comm);
  
  /* Total number of bounding elements */
  MPI_Reduce(&(grid.n_be),&(grid.Gn_be),1,MPI_INT,MPI_SUM,0,mpi_comm);
  
  /* Gather number of nodes from all domains */
  MPI_Gather (&(grid.nn),1,MPI_LONG,DomNn,1,MPI_LONG,0,mpi_comm);
  
  if (myrank == 0 && PFEM_DEBUG)
    PGFEM_printf(" Done.\nRedistributing information...");                          
    

  for(int ia=0; ia<physicsno; ia++)
  {  
    fv[ia].ndofd = generate_local_dof_ids(grid.ne,grid.nce,grid.nn,fv[ia].ndofn,grid.node,
                                          grid.element,grid.coel,grid.b_elems,mpi_comm,ia);
    com[ia].DomDof[myrank] = generate_global_dof_ids(grid.ne,grid.nce,grid.nn,fv[ia].ndofn,grid.node,
                                                     grid.element,grid.coel,grid.b_elems,mpi_comm,ia);
    // Gather degrees of freedom from all domains
    MPI_Allgather(MPI_IN_PLACE,1,MPI_LONG,com[ia].DomDof,1,MPI_LONG,mpi_comm);

    // Make integer copy of DomDof.  May eventually switch everything to
    // integer since 64 bit
    int *dist = aloc1i(com[ia].nproc+1);
    // build_dist(DomDof,dist,nproc);
    build_distribution(com[ia].DomDof,dist,mpi_comm);    

    for(long ib=0;ib<com[ia].nproc;ib++)
    {
      fv[ia].Gndof += com[ia].DomDof[ib];
      if(ia==0)
        grid.Gne += DomNe[ib];
    }
  
    set_HYPRE_row_col_bounds(sol[ia].PGFEM_hypre, fv[ia].Gndof,com[ia].DomDof,myrank);  
  
    renumber_global_dof_ids(grid.ne,grid.nce,grid.n_be,grid.nn,fv[ia].ndofn,com[ia].DomDof,grid.node,
                            grid.element,grid.coel,grid.b_elems,mpi_comm,ia);
    com[ia].NBN = distribute_global_dof_ids(grid.ne,grid.nce,grid.n_be,grid.nn,fv[ia].ndofn,ndim,grid.node,
                            grid.element,grid.coel,grid.b_elems, com[ia].hints, mpi_comm,ia);
    
    // ALlocate Ap, Ai      
    com[ia].Ap = aloc1i(com[ia].DomDof[myrank]+1);
    com[ia].comm  = (COMMUN) PGFEM_calloc (1,sizeof(COMMUN_1));
    initialize_commun(com[ia].comm );

    com[ia].Ai = Psparse_ApAi(com[ia].nproc,myrank,grid.ne,grid.n_be,grid.nn,fv[ia].ndofn,fv[ia].ndofd,
                              grid.element,grid.b_elems,grid.node,com[ia].Ap,grid.nce,grid.coel,com[ia].DomDof,
                              &(com[ia].GDof),com[ia].comm ,mpi_comm,options.cohesive,ia);
    pgfem_comm_build_fast_maps(com[ia].comm ,fv[ia].ndofd,com[ia].DomDof[myrank],com[ia].GDof);

    // Total number of nonzeros and skyline
    int APP  = 0;
    long sky = 0;
    MPI_Reduce (&(com[ia].Ap[com[ia].DomDof[myrank]]),&APP,1,MPI_INT,MPI_SUM,0,mpi_comm);
    long temp_int = skyline((int) com[ia].DomDof[myrank],com[ia].Ap,com[ia].Ai,dist[myrank]);

    MPI_Reduce (&temp_int,&sky,1,MPI_INT,MPI_SUM,0,mpi_comm);

    //---------------------------------------------------------------------- 
    // print simulation setting info
    //----------------------------------------------------------------------
    //---->    
    if (myrank == 0)
    {
      if(ia==0) // print onece
        err += print_PGFem3D_run_info(argc, argv, &grid, com+ia, &load, gem, &options);

      PGFEM_printf ("---------------------------------------------\n");
      PGFEM_printf ("Physics name: %s\n", mp.physicsname[ia]);      
      PGFEM_printf ("Total number of degrees of freedom       : %ld\n", fv[ia].Gndof);
      PGFEM_printf ("Total number of nonzeros in the matrix   : %d\n",APP);
      PGFEM_printf ("Symmetric skyline (including diagonal)   : %ld\n",sky);
    }
    //<---------------------------------------------------------------------
    free(dist);    
  }

  dealoc1l (DomNe);
  dealoc1l (DomNn);
  
  if (myrank == 0 && options.cohesive == 1)
    PGFEM_printf ("Number of cohesive elements              : %ld\n",grid.Gnce);  
  

  {        
    if((load.sups[mp_id_M])->npd > 0){
      load.sup_defl = aloc1((load.sups[mp_id_M])->npd);
    } else {
      load.sup_defl = NULL;
    }    
    
    /*=== TESTING ===*/
    double *nodal_forces = PGFEM_calloc(fv[mp_id_M].ndofd,sizeof(double));
    int n_feats = 0;
    int n_sur_trac_elem = 0;
    SUR_TRAC_ELEM *ste = NULL;
    {
      int *feat_type = NULL;
      int *feat_id = NULL;
      double *loads = NULL;
      
      char *trac_fname = NULL;
      alloc_sprintf(&trac_fname,"%s/traction.in",options.ipath);
      
      read_applied_surface_tractions_fname(trac_fname,&n_feats,
              &feat_type,&feat_id,&loads);
      
      generate_applied_surface_traction_list(grid.ne,grid.element,
              n_feats,feat_type,
              feat_id,&n_sur_trac_elem,
              &ste);
      
      compute_applied_traction_res(fv[mp_id_M].ndofn,grid.node,grid.element,
              n_sur_trac_elem,ste,
              n_feats,loads,
              nodal_forces, mp_id_M);
      
      double tmp_sum = 0.0;
      for(int i=0; i<fv[mp_id_M].ndofd; i++){
        tmp_sum += nodal_forces[i];
      }
      
      MPI_Allreduce(MPI_IN_PLACE,&tmp_sum,1,MPI_DOUBLE,
              MPI_SUM,mpi_comm);
      
      if(myrank == 0){
        PGFEM_printf("Total load from surface tractions: %.8e\n\n",tmp_sum);
      }
      
      free(feat_type);
      free(feat_id);
      free(loads);
      free(trac_fname);
    }

    //----------------------------------------------------------------------
    // read solver file ( *.in.st)
    // file pointer (solver file) will not be freed and
    // saved in order to read loads increments as time is elapsing.
    //----------------------------------------------------------------------
    //---->
    err += read_solver_file(&time_steps, &mat, fv+mp_id_M, sol+mp_id_M, &load, &arc, crpl, &options, myrank);
    //<---------------------------------------------------------------------    
                         
    
    // Nonlinear solver
    if (myrank == 0) {
      if (sol[0].FNR == 0 || sol[0].FNR == 1) {
        PGFEM_printf ("\nNONLINEAR SOLVER : NEWTON-RAPHSON METHOD\n");
      }
      if ((sol[0].FNR == 2 || sol[0].FNR == 3) && arc.ARC == 0){
        PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Crisfield\n");
      }
      if ((sol[0].FNR == 2 || sol[0].FNR == 3) && arc.ARC == 1) {
        PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Simo\n");
      }
    }
    
    /* HYPRE INITIALIZATION ROUTINES */
    if(options.solverpackage == HYPRE){ /* HYPRE */
      /* Initialize HYPRE */
      for(int ia=0; ia<physicsno; ia++)
      {
        hypre_initialize(com[ia].Ap,
                         com[ia].Ai,
                         com[ia].DomDof[myrank],
                         sol[ia].iter_max_sol,
                         sol[ia].err,
                         sol[ia].PGFEM_hypre,
                         &options,mpi_comm);
      }
    }

    /* alocation of the sigma vector */
    for(int ia=0; ia<mp.physicsno; ia++)
      err += construct_field_varialbe(fv+ia, &grid, com+ia, &options, myrank, mp.physics_ids[ia]);
    
    if (sol[mp_id_M].FNR == 2 || sol[mp_id_M].FNR == 3)
      err += construct_arc_length_variable(&arc, fv+mp_id_M, com+mp_id_M, myrank);

    /* push nodal_forces to s->R */
    vvplus(fv[mp_id_M].R,nodal_forces,fv[mp_id_M].ndofd);
    
    /* alocation of the eps vector */
    initialize_damage(grid.ne,grid.element,mat.hommat,fv[mp_id_M].eps,options.analysis_type);


    
    if (options.analysis_type == CM) {
      /* parameter list and initialize const. model at int points.
       * NOTE: should catch/handle returned error flag...
       */
      char *cm_filename = NULL;
      alloc_sprintf(&cm_filename,"%s/model_params.in",options.ipath);
      FILE *cm_in = PGFEM_fopen(cm_filename, "r");
      read_model_parameters_list(&param_list, mat.nhommat, mat.hommat, cm_in);
      free(cm_filename);
      fclose(cm_in);
      init_all_constitutive_model(fv[mp_id_M].eps,grid.ne,grid.element,mat.nhommat,param_list);
    }
    
    /* alocation of pressure variables */
    int nVol = 1;
    switch(options.analysis_type){
      case TF:
        nVol = 1;
        if(grid.element[0].toe==10 && fv[mp_id_M].ndofn==3)
        {
          fv[mp_id_M].npres = 1;
          nVol = 1;
        }
        break;
      case STABILIZED: case MINI: case MINI_3F:
        if(fv[mp_id_M].npres != 4){
          fv[mp_id_M].npres = 4;
          if(myrank == 0){
            PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 4.\n"
                    "Re-setting to 4 and continuing...\n");
          }
        }
        break;
      case DISP: // intented not to have break
      case CM:
        if(fv[mp_id_M].npres != 0){
          fv[mp_id_M].npres = 0;
          if (myrank == 0) {
            PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 0.\n"
                    "Re-setting to 0 and continuing...\n");
          }
        }
        break;
      default:
        if(fv[mp_id_M].npres != 1){
          fv[mp_id_M].npres = 1;
          if (myrank == 0) {
            PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 1.\n"
                    "Re-setting to 1 and continuing...\n");
          }
        }
        break;
    }/* switch */
    build_pressure_nodes (grid.ne,fv[mp_id_M].npres,grid.element,fv[mp_id_M].sig,fv[mp_id_M].eps,options.analysis_type);
    build_crystal_plast (grid.ne,grid.element,fv[mp_id_M].sig,fv[mp_id_M].eps,crpl,
            options.analysis_type,options.plc);
    
    /* \/ initialized element varialbes */
    if(options.analysis_type==TF)
    {
      for (int e=0;e<grid.ne;e++)
      {
        if(fv[mp_id_M].npres==1)
        {
          fv[mp_id_M].eps[e].d_T   = (double *) PGFEM_calloc(3,sizeof(double));
          for(int a=0; a<3; a++)
            fv[mp_id_M].eps[e].d_T[a] = 0.0;
        }
        
        fv[mp_id_M].eps[e].T   = (double *) PGFEM_calloc(nVol*3,sizeof(double));
        for(int a=0; a<nVol*3; a++)
          fv[mp_id_M].eps[e].T[a] = 1.0;
      }
    }
    /* /\ initialized element varialbes */
        
    /*******************************************************************/
    /* this is for inertia */    
                
    int restart_tim = options.restart;
    
    double tnm1[2] = {-1.0,-1.0};
    err += read_initial_values(&grid, 
                               &mat, 
                               fv+mp_id_M, 
                               sol+mp_id_M, 
                               &load, 
                               &time_steps, 
                               &options, 
                               &restart_tim, 
                               tnm1, 
                               myrank);
    
    for(long idx_a = 0; idx_a<grid.nn; idx_a++)
    {
      for(long idx_b = 0; idx_b<fv[mp_id_M].ndofn; idx_b++)
      {
        long id = grid.node[idx_a].id_map[mp_id_M].id[idx_b];
        if(id>0)
          fv[mp_id_M].u_np1[id-1] = fv[mp_id_M].u_n[idx_a*fv[mp_id_M].ndofn + idx_b];
      }
    }
    
    /*******************************************************************/
    
    /* set finite deformations variables */
    set_fini_def (grid.ne,fv[mp_id_M].npres,grid.element,fv[mp_id_M].eps,fv[mp_id_M].sig,options.analysis_type);
    if (options.analysis_type == FS_CRPL){
      set_fini_def_pl (grid.ne,fv[mp_id_M].npres,grid.element,fv[mp_id_M].eps,fv[mp_id_M].sig,crpl,
              options.analysis_type,options.plc);
    }
    
    /*  NODE (PRESCRIBED DEFLECTION)- SUPPORT COORDINATES generation
     * of the load vector  */
    time_steps.dt_np1 = arc.dt0 = time_steps.times[1] - time_steps.times[0];
    if (time_steps.dt_np1 == 0.0){
      if (myrank == 0){
        PGFEM_printf("Incorrect dt\n");
      }
      PGFEM_Comm_code_abort(mpi_comm,0);
    }
    
    load_vec_node_defl (fv[mp_id_M].f_defl,grid.ne,fv[mp_id_M].ndofn,grid.element,grid.b_elems,grid.node,mat.hommat,
            mat.matgeom,load.sups[mp_id_M],fv[mp_id_M].npres,sol[mp_id_M].nor_min,fv[mp_id_M].sig,fv[mp_id_M].eps,time_steps.dt_np1,
            crpl,options.stab,fv[mp_id_M].u_np1,fv[mp_id_M].u_n,&options,sol[mp_id_M].alpha,mp_id_M);
    
    /*  NODE - generation of the load vector  */
    load_vec_node (fv[mp_id_M].R,load.nln,ndim,load.znod,grid.node,mp_id_M);
    /*  ELEMENT - generation of the load vector  */
    load_vec_elem_sur (fv[mp_id_M].R,load.nle_s,ndim,grid.element,load.zele_s);
    
    /* R   -> Incramental forces
     * RR  -> Total forces for sudivided increment
     * RRn -> Total force after equiblirium */
    
    vvplus  (fv[mp_id_M].f,fv[mp_id_M].R,fv[mp_id_M].ndofd);
    vvplus  (fv[mp_id_M].RR,fv[mp_id_M].f,fv[mp_id_M].ndofd);
    vvminus (fv[mp_id_M].f,fv[mp_id_M].f_defl,fv[mp_id_M].ndofd);
    
    /* Transform LOCAL load vector to GLOBAL */
    if(sol[mp_id_M].FNR == 2 || sol[mp_id_M].FNR == 3)
      LToG (fv[mp_id_M].R,arc.BS_R,myrank,com[mp_id_M].nproc,fv[mp_id_M].ndofd,com[mp_id_M].DomDof,com[mp_id_M].GDof,com[mp_id_M].comm ,mpi_comm);
    
    /* Prescribed deflection */
    for (long i=0;i<(load.sups[mp_id_M])->npd;i++){
      load.sup_defl[i] = (load.sups[mp_id_M])->defl_d[i];
    }
    
    /*=== NO PERIODIC ===*/
    
    arc.lm = arc.DLM = arc.DET0 = 0.0;
    /*PD = 1;*/
    arc.DAL = arc.DLM0 = arc.dAL0;
    long tim;
    tim = arc.AT = arc.ITT = 0;
    
    /* compute un-deformed volume */
    oVolume = 0;
    GVolume = T_VOLUME (grid.ne,ndim,grid.element,grid.node);
    MPI_Allreduce (&GVolume,&oVolume,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    if (myrank == 0){
      PGFEM_printf ("oVolume = %12.12f\n",oVolume);
    }
    VVolume = oVolume;
    
    /*=== BEGIN SOLVE ===*/
    time_steps.dt_np1 = time_steps.times[1] - time_steps.times[0];
    
    ///////////////////////////////////////////////////////////////////
    // start time stepping
    ///////////////////////////////////////////////////////////////////
    while (time_steps.nt > tim)
    {
      time_steps.tim    = tim;
      time_steps.dt_n   = time_steps.dt_np1;
      time_steps.dt_np1 = time_steps.times[tim+1] - time_steps.times[tim];
      if (time_steps.dt_np1 <= 0.0){
        if (myrank == 0) {
          PGFEM_printf("Incorrect dt\n");
        }
        PGFEM_Comm_code_abort(mpi_comm,0);
      }
      
      if (myrank == 0){
        PGFEM_printf("\nFinite deformations time step %ld) "
                " Time %e | dt = %e\n",
                tim,time_steps.times[tim+1],time_steps.dt_np1);
      }
      
      /*=== NEWTON RAPHSON ===*/
      if (sol[mp_id_M].FNR == 0 || sol[mp_id_M].FNR == 1)
      {
        //----------------------------------------------------------------------
        // file pointer (solver file) is active and used to update loads increments        
        // fv[mp_id_M].R   -> Incramental forces
        // fv[mp_id_M].RR  -> Total forces for sudivided increment
        // fv[mp_id_M].RRn -> Total force after equiblirium    
        // push nodal_forces to s->R        
        //----------------------------------------------------------------------
        //---->        
        err += read_and_apply_load_increments(&grid, fv+mp_id_M, &load, &mp, tim, mpi_comm, myrank);
        if (load.tim_load[tim] == 1 && tim != 0)
          vvplus(fv[mp_id_M].R,nodal_forces,fv[mp_id_M].ndofd);
        //<---------------------------------------------------------------------
    
        int n_step = 0;
        sol[mp_id_M].n_step = &n_step;
        sol[mp_id_T].n_step = &n_step;

        //----------------------------------------------------------------------
        // add load increments util time reaches the restart point
        //----------------------------------------------------------------------
        //---->
        if(tim<restart_tim+1)
        {
          for (long i=0;i<(load.sups[mp_id_M])->npd;i++){
            (load.sups[mp_id_M])->defl[i] += (load.sups[mp_id_M])->defl_d[i];
            (load.sups[mp_id_M])->defl_d[i] = 0.0;
          }
          (tim)++;
          continue;
        }
        //<---------------------------------------------------------------------
        
        if(tim==restart_tim+1 && tnm1[1]>0)
        {
          time_steps.times[tim-1] = tnm1[1]; // tnm1[0] = times[tim-2]
                                                        // tnm1[1] = times[tim-1]
                                                        // tnm1[2] = times[tim]
          
          // if restart_tim==0: tim = 1
          if(tim>=2)              // if restart_tim==1: tim = 2
            time_steps.times[tim-2] = tnm1[0];
        }


        //----------------------------------------------------------------------
        // Perform Newton Raphson interation
        //----------------------------------------------------------------------
        //---->        
        fflush(PGFEM_stdout);
        
        hypre_time += Multiphysics_Newton_Raphson(1,&grid,&mat,fv,sol,&load,com,&time_steps,
                                          crpl,mpi_comm,VVolume,&options,&mp);
        
        /* Null global vectors */
        for (long i=0;i<fv[mp_id_M].ndofd;i++){
          fv[mp_id_M].RRn[i] += fv[mp_id_M].R[i];
          fv[mp_id_M].RR[i] = fv[mp_id_M].RRn[i];
          fv[mp_id_M].R[i] = 0.0;
        }
        
        /* null the prescribed displacement increment */
        nulld(load.sup_defl,(load.sups[mp_id_M])->npd);
      }/* end NR */
      
      /*=== ARC LENGTH ===*/
      if(sol[mp_id_M].FNR == 2 || sol[mp_id_M].FNR == 3){
        double dlm = Arc_length_test(&grid,&mat,fv+mp_id_M,sol+mp_id_M,&load,com+mp_id_M,&time_steps,
                                     crpl,&arc,mpi_comm,VVolume,&options, 0);
        
        /* Load multiplier */
        arc.lm += dlm;
        
        /* Total force vector */
        for (long i=0;i<fv[mp_id_M].ndofd;i++){
          fv[mp_id_M].RR[i] = arc.lm*fv[mp_id_M].R[i];
        }
      }/* end AL */
      
      /*=== OUTPUT ===*/
      /* Calculating equvivalent Mises stresses and strains vectors */
      Mises (grid.ne,fv[mp_id_M].sig,fv[mp_id_M].eps,options.analysis_type);
      
      /* update output stuff for CM interface */
      if(options.analysis_type == CM && options.cm!=0){
        constitutive_model_update_output_variables(fv[mp_id_M].sig,fv[mp_id_M].eps,grid.node,grid.element,grid.ne,
                time_steps.dt_np1,&options, sol[mp_id_M].alpha);
      }
      
      /* print tractions on marked features */
      {
        double *sur_forces = NULL;
        if(n_feats > 0){
          sur_forces = PGFEM_calloc(n_feats*ndim,sizeof(double));
          compute_resultant_force(n_feats,n_sur_trac_elem,
                  ste,grid.node,grid.element,
                  fv[mp_id_M].sig,fv[mp_id_M].eps,sur_forces);
          MPI_Allreduce(MPI_IN_PLACE,sur_forces,n_feats*ndim,
                  MPI_DOUBLE,MPI_SUM,mpi_comm);
          if(myrank == 0){
            PGFEM_printf("Forces on marked features:\n");
            print_array_d(PGFEM_stdout,sur_forces,n_feats*ndim,
                    n_feats,ndim);
            fflush(PGFEM_stdout);
          }
        }
        free(sur_forces);
      }
      
      if(options.comp_print_reaction)
      {
        double dts[2];
        if(tim==0)
          dts[DT_N] = time_steps.times[tim+1] - time_steps.times[tim];
        else
          dts[DT_N] = time_steps.times[tim] - time_steps.times[tim-1];
        
        dts[DT_NP1] = time_steps.times[tim+1] - time_steps.times[tim];
        
        fd_res_compute_reactions(fv[mp_id_M].ndofn, fv[mp_id_M].npres, fv[mp_id_M].d_u, fv[mp_id_M].u_np1, grid.element, grid.node,
                mat.matgeom, mat.hommat, load.sups[mp_id_M], fv[mp_id_M].eps, fv[mp_id_M].sig, sol[mp_id_M].nor_min,
                crpl, dts, time_steps.times[tim+1], options.stab, mpi_comm,
                &options, sol[mp_id_M].alpha, fv[mp_id_M].u_n, fv[mp_id_M].u_nm1,mp_id_M);
      }
      
      if(options.comp_print_macro)
      {
        /* Calculate macro deformation gradient */
        double *GF = computeMacroF(grid.element,grid.ne,grid.node,grid.nn,fv[mp_id_M].eps,oVolume,mpi_comm);
        double *GS = computeMacroS(grid.element,grid.ne,grid.node,grid.nn,fv[mp_id_M].sig,oVolume,mpi_comm);
        double *GP = computeMacroP(grid.element,grid.ne,grid.node,grid.nn,fv[mp_id_M].sig,fv[mp_id_M].eps,oVolume,mpi_comm);
        
        /* print GF & GS to file */
        if(myrank==0)
        {
          
          sprintf(filename,"%s_macro.out.%ld",out_dat,tim);
          FILE *out = fopen(filename,"w");
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[0],GF[1],GF[2]);
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[3],GF[4],GF[5]);
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[6],GF[7],GF[8]);
          PGFEM_fprintf(out,"\n");
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\t",GS[0],GS[1],GS[2]);
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GS[3],GS[4],GS[5]);
          PGFEM_fprintf(out,"\n");
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[0],GP[1],GP[2]);
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[3],GP[4],GP[5]);
          PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[6],GP[7],GP[8]);
          fclose(out);
        }
        
        free(GF);
        free(GS);
        free(GP);
      }
      
      if (time_steps.print[tim] == 1 && options.vis_format != VIS_NONE ) {
        if(options.ascii){
          ASCII_output(&options,mpi_comm,tim,time_steps.times,grid.Gnn,grid.nn,grid.ne,grid.nce,fv[mp_id_M].ndofd,
                  com[mp_id_M].DomDof,com[mp_id_M].Ap,sol[mp_id_M].FNR,arc.lm,fv[mp_id_M].pores,VVolume,grid.node,grid.element,load.sups[mp_id_M],
                  fv[mp_id_M].u_np1,fv[mp_id_M].eps,fv[mp_id_M].sig,fv[mp_id_M].sig_n,grid.coel);
        } /* End ASCII output */
        
        switch(options.vis_format){
          case VIS_ELIXIR:/* Print to elix file */
            sprintf (filename,"%s_%d.elx%ld",out_dat,myrank,tim);
            elixir (filename,grid.nn,grid.ne,ndim,grid.node,grid.element,load.sups[mp_id_M],fv[mp_id_M].u_np1,fv[mp_id_M].sig,
                    fv[mp_id_M].sig_n,fv[mp_id_M].eps,options.smoothing,grid.nce,grid.coel,&options);
            break;
          case VIS_ENSIGHT:/* Print to EnSight files */
            sprintf (filename,"%s",out_dat);
            EnSight (filename,tim,time_steps.nt,grid.nn,grid.ne,ndim,grid.node,grid.element,load.sups[mp_id_M],
                    fv[mp_id_M].u_np1,fv[mp_id_M].sig,fv[mp_id_M].sig_n,fv[mp_id_M].eps,options.smoothing,grid.nce,grid.coel,
                    /*nge,geel,ngn,gnod,*/sol[mp_id_M].FNR,arc.lm,ensight,mpi_comm,
                    &options);
            break;
          case VIS_VTK:/* Print to VTK files */
            if(myrank == 0){
              VTK_print_master(options.opath,options.ofname,
                      tim,com[mp_id_M].nproc,&options);
            }
            
            VTK_print_vtu(options.opath,options.ofname,tim,
                    myrank,grid.ne,grid.nn,grid.node,grid.element,load.sups[mp_id_M],fv[mp_id_M].u_np1,fv[mp_id_M].sig,fv[mp_id_M].eps,
                    &options,mp_id_M);
            
            // save restart files
            if(SAVE_RESTART_FILE)
            {
              write_restart(fv[mp_id_M].u_nm1, fv[mp_id_M].u_n, &options,
                      grid.element, grid.node,fv[mp_id_M].sig,fv[mp_id_M].eps,load.sups[mp_id_M],
                      myrank, grid.ne, grid.nn, fv[mp_id_M].ndofn, fv[mp_id_M].ndofd, tim, time_steps.times, fv[mp_id_M].NORM,
                      mp_id_M);
            }
            
            if (options.cohesive == 1){
              if(myrank == 0){
                VTK_print_cohesive_master(options.opath,
                        options.ofname,tim,com[mp_id_M].nproc,
                        &options);
              }
              
              VTK_print_cohesive_vtu(options.opath,options.ofname,
                      tim,myrank,grid.nce,grid.node,grid.coel,load.sups[mp_id_M],fv[mp_id_M].u_np1,ensight,
                      &options,mp_id_M);
            }
            break;
          default: /* no output */ break;
        }/* switch(format) */
        
      }/* end output */
      
      if (myrank == 0){
        PGFEM_printf("\n");
        PGFEM_printf("*********************************************\n");
        PGFEM_printf("*********************************************\n");
      }
      
      tim++;
    }/* end while */
    
    destroy_applied_surface_traction_list(n_sur_trac_elem,ste);
    free(nodal_forces);  
  }

  //----------------------------------------------------------------------
  // deallocate objects
  //----------------------------------------------------------------------
  //---->
  if (options.solverpackage == HYPRE)
  {  
    destroy_PGFEM_HYPRE_solve_info(sol[mp_id_M].PGFEM_hypre);
    destroy_PGFEM_HYPRE_solve_info(sol[mp_id_T].PGFEM_hypre);
  }
    
  err += destruct_time_stepping(&time_steps);  
  
  for(int ia=0; ia<physicsno; ia++)
    err += destruct_field_varialbe(fv+ia, &grid, &options, mp.physics_ids[ia]);
  
  free(fv);  

  err += destruct_loading_steps(&load, &mp);
  err += destroy_model_parameters_list(mat.nhommat,param_list);
  err += destruct_material(&mat, &options);
  err += destruct_grid(&grid, &options, &mp);
  err += destruct_communication_structure(com+mp_id_M);  
  err += destruct_multiphysics(&mp);
  
  com[mp_id_T].bndel = NULL;
  com[mp_id_T].hints = NULL;
  err += destruct_communication_structure(com+mp_id_T); 

  if (sol[mp_id_M].FNR == 2 || sol[mp_id_M].FNR == 3)
    err += destruct_arc_length_variable(&arc);
    
  destroy_ensight(ensight);
  //<---------------------------------------------------------------------  

  total_time += MPI_Wtime(); // measure time spent

  //----------------------------------------------------------------------
  // print time of analysis and finalize
  //----------------------------------------------------------------------
  //---->
  if (myrank == 0)
    err += print_PGFem3D_final(total_time, hypre_time, &usage, myrank);  

  PGFEM_finalize_io();
  
  int flag_MPI_finalized;
  MPI_Finalized(&flag_MPI_finalized);
  if(!flag_MPI_finalized)
  {
    if(myrank==0)
      printf("MPI finalizing\n");
    
    MPI_Finalize();
  }
  //<---------------------------------------------------------------------
  
  return err;
}
