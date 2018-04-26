#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PFEM3d.h"
#include "Arc_length.h"
#include "Newton_Raphson.h"
#include "PGFEM_io.h"
#include "PGFem3D_options.h"
#include "Printing.h"
#include "SetGlobalNodeNumbers.h"
#include "Psparse_ApAi.h"
#include "allocation.h"
#include "applied_traction.h"
#include "bounding_element.h"
#include "bounding_element_utils.h"
#include "build_distribution.h"
#include "comm_hints.h"
#include "computeMacroF.h"
#include "computeMacroS.h"
#include "constitutive_model.h"
#include "dynamics.h"
#include "element.h"
#include "enumerations.h"
#include "fd_residuals.h"
#include "gen_path.h"
#include "generate_dof_ids.h"
#include "homogen.h"
#include "in.h"
#include "incl.h"
#include "initialize_damage.h"
#include "interface_macro.h"
#include "load.h"
#include "matice.h"
#include "matrix_printing.h"
#include "node.h"
#include "out.h"
#include "post_processing.h"
#include "print_dist.h"
#include "profiler.h"
#include "read_cryst_plast.h"
#include "read_input_file.h"
#include "renumber_ID.h"
#include "restart.h"
#include "set_fini_def.h"
#include "skyline.h"
#include "three_field_element.h"
#include "utils.h"
#include "vtk_output.h"
#include "set_initial_plastic_deformation_gradient.h"
#include <cstdlib>
#include <cassert>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "crystal_plasticity_integration.h"

namespace {
using namespace pgfem3d;
const constexpr int periodic = 0;
}

/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/
#define SAVE_RESTART_FILE 1

/// Print PGFem3D running options
/// This will print out mesh info, analysis options
///
/// \param[in] argc number of arguments passed through command line
/// \param[in] argv arguments passed through command line
/// \param[in] grid mesh info (Grid object)
/// \param[in] com commuincation info (CommunicationStructure object)
/// \param[in] load info for loading steps (LoadingSteps object)
/// \param[in] opts structure PGFem3D option
static void print_PGFem3D_run_info(int argc,
                                  char *argv[],
                                  const Grid *grid,
                                  const CommunicationStructure *com,
                                  const LoadingSteps *load,
                                  const PGFem3D_opt *opts)
{
  PrintTitleV1();
  for (int i = 0, e = argc; i < e; ++i) {
    PGFEM_printf("%s ", argv[i]);
  }

  PGFEM_printf("\n\n");

  print_interpreted_options(opts);

  if (opts->multi_scale) {
    if (load->sups[MULTIPHYSICS_MECHANICAL]->npd >= 9) {
      PGFEM_printf("*** BULK Multiscale Modeling ***\n");
    } else {
      PGFEM_printf("*** INTERFACE Multiscale Modeling ***\n");
    }
  }

  PGFEM_printf("\n");
  PGFEM_printf("Number of total nodes                    : %ld\n", grid->Gnn);
  PGFEM_printf("Number of nodes on domain interfaces     : %ld\n", com->NBN);
  PGFEM_printf("Total number of elements                 : %ld\n", grid->Gne);
  PGFEM_printf("Number of elems on the COMM interfaces   : %ld\n", grid->Gnbndel);
  PGFEM_printf("Total number of bounding (surf) elems    : %d\n",  grid->Gn_be);
}

/// Print PGFem3D run time info
/// This will print out how much time takes for the finite element analysis
///
/// \param[in] total_time total simulation time
/// \param[in] startup_time - time untill it starts the time-steps iterations 
/// \param[in] hypre_time time took for linear system solver
/// \param[in] stiffmat_time - A matrix construction time of equ Ax=b 
/// \param[in] residuals_time - b matrix construction time of equ Ax=b
/// \param[in] usage detailed time info
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int print_PGFem3D_final(const Multiphysics& mp,
                        double total_time,
                        double startup_time,
                        std::vector<double> &hypre_time,
                        std::vector<double> &stiffmat_time,
                        std::vector<double> &residuals_time,
                        int myrank)
{
  int err = 0;
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  PGFEM_printf("\n");
  PGFEM_printf("Time of analysis on processor [%d] - "
               " System %ld.%ld, User %ld.%ld.\n\n",
               myrank, usage.ru_stime.tv_sec, usage.ru_stime.tv_usec,
               usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);

  PGFEM_printf("Total time (no MPI_Init()) = %f\n", total_time);
  PGFEM_printf("Startup time               = %f\n\n", startup_time);

  for(int mp_id = 0; mp_id<mp.physicsno; mp_id++)
  {
    PGFEM_printf("%s::Hypre solve time          = %f\n", mp.physicsname[mp_id], hypre_time[mp_id]);
    PGFEM_printf("%s::Stiffmat compute time     = %f\n", mp.physicsname[mp_id], stiffmat_time[mp_id]);
    PGFEM_printf("%s::Residualsmat compute time = %f\n\n", mp.physicsname[mp_id], residuals_time[mp_id]);
  }

  return err;
}

/// Prints exascale metrics
/// Currently, this will output the total number of ODEs operations
///
/// \param[in] local_EXA_metric exascale metric counter for total number of integration iterations
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
void print_EXA_metrics(const MPI_Comm mpi_comm, 
                       const int myrank,
                       std::vector<double> hypre_time,
                       std::vector<double> residuals_time,
                       double total_time,
                       const PGFem3D_opt *opts)
{
  int MPI_process_number;
  MPI_Comm_size(mpi_comm, &MPI_process_number);
  
  if (myrank == 0){
    double total_residual_time = 0.0;   //residual time across all physics
    for (auto& n : residuals_time)
      total_residual_time += n;
      
    double total_hypre_time = 0.0;      //hypre time across all physics
    for (auto& n : hypre_time)
      total_hypre_time += n;
      
    double total_combined_time = total_residual_time + total_hypre_time;
    
    double EXA_Numerator = (total_residual_time / total_combined_time) * total_EXA_metric + 
                              (total_hypre_time / total_combined_time) * dof_EXA_metric;
                              
    if (opts->print_EXA_details){
      PGFEM_printf("Total Residual time: %f\n", total_residual_time);
      PGFEM_printf("Total Hypre time: %f\n", total_hypre_time);
      PGFEM_printf("Total number of DOF computations: %ld\n", dof_EXA_metric);
      PGFEM_printf("Total number of MPI processes: %d\n", MPI_process_number);
      PGFEM_printf("Final EXA metric numerator: %f\n", EXA_Numerator);
    }
    
    double EXA_Denominator = total_time * MPI_process_number;
    PGFEM_printf("\nFinal EXA metric: %f\n\n", EXA_Numerator/EXA_Denominator);
  }
}


/// print simulation results
/// output format is VTK so that this function calls VTK_IO library
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] FV array of field variable object
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
/// \param[in] tim current time step number
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int print_results(Grid *grid,
                  MaterialProperty *mat,
                  FieldVariables *FV,
                  Solver *SOL,
                  LoadingSteps *load,
                  CommunicationStructure *COM,
                  TimeStepping *time_steps,
                  CRPL *crpl,
                  Ensight *ensight,
                  PRINT_MULTIPHYSICS_RESULT *pmr,
                  MPI_Comm mpi_comm,
                  const double oVolume,
                  const double VVolume,
                  const PGFem3D_opt *opts,
                  const Multiphysics& mp,
                  long tim,
                  int myrank)
{
  int err = 0;

  Solver                 *sol = NULL;
  FieldVariables          *fv = NULL;
  CommunicationStructure *com = NULL;
  SUPP sup = NULL;
  int mp_id_M = -1;

  for(int ia = 0; ia<mp.physicsno; ia++)
  {
    if(mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
    {
      mp_id_M = ia;
      sol = SOL + mp_id_M;
      fv  = FV  + mp_id_M;
      com = COM + mp_id_M;
      sup = load->sups[mp_id_M];
    }
  }

  // output file name
  char filename[500],out_dat[500];
  sprintf(out_dat,"%s/%s",opts->opath,opts->ofname);

  if(mp_id_M >= 0)
  {
    if(opts->comp_print_reaction)
    {
      double dts[2];
      if(tim==0)
        dts[DT_N] = time_steps->times[tim+1] - time_steps->times[tim];
      else
        dts[DT_N] = time_steps->times[tim] - time_steps->times[tim-1];

      dts[DT_NP1] = time_steps->times[tim+1] - time_steps->times[tim];

      sol->run_integration_algorithm = 0;
      err += fd_res_compute_reactions_MP(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,
                                         mp_id_M,time_steps->times[tim+1],dts);
      sol->run_integration_algorithm = 1;
    }

    if(opts->comp_print_macro)
    {
      /* Calculate macro deformation gradient */
      double *GF = computeMacroF(grid->element,grid->ne,grid->node,grid->nn,fv->eps,oVolume,mpi_comm);
      double *GS = computeMacroS(grid->element,grid->ne,grid->node,grid->nn,fv->sig,oVolume,mpi_comm);
      double *GP = computeMacroP(grid->element,grid->ne,grid->node,grid->nn,fv->sig,fv->eps,oVolume,mpi_comm);

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
  }

  if (time_steps->print[tim] == 1 && opts->vis_format != VIS_NONE )
  {
    if(opts->ascii && mp_id_M >= 0)
    {
      ASCII_output(opts,mpi_comm,tim,time_steps->times,grid->Gnn,grid->nn,grid->ne,grid->nce,fv->ndofd,
                   com->DomDof,com->Ap,sol->FNR,sol->arc->lm,fv->pores,VVolume,grid->node,grid->element,sup,
                   fv->u_np1,fv->eps,fv->sig,fv->sig_n,grid->coel);
    } /* End ASCII output */

    if(opts->vis_format == VIS_VTK)
    {
      if(myrank == 0)
        err += VTK_write_multiphysics_master(pmr,mp.total_write_no,opts,tim,myrank,COM[0].nproc);

      err += VTK_write_multiphysics_vtu(grid,mat,FV,load,pmr,mp.total_write_no,opts,tim,myrank);

      // print cohesive element results
      if((mp_id_M >= 0) && (opts->cohesive == 1))
      {
        if(myrank == 0)
          VTK_print_cohesive_master(opts->opath,opts->ofname,tim,com->nproc,opts);

        VTK_print_cohesive_vtu(opts->opath,opts->ofname,tim,myrank,
                               grid->nce,grid->node,grid->coel,sup,fv->u_np1,ensight,
                               opts,mp_id_M);
      }
    }
    else
    {
      if(mp_id_M >= 0)
      {
        switch(opts->vis_format)
        {
         case VIS_ELIXIR:/* Print to elix file */
          sprintf (filename,"%s_%d.elx%ld",out_dat,myrank,tim);
          elixir (filename,grid->nn,grid->ne,grid->nsd,grid->node,grid->element,sup,fv->u_np1,fv->sig,
                  fv->sig_n,fv->eps,opts->smoothing,grid->nce,grid->coel,opts);
          break;
         case VIS_ENSIGHT:/* Print to EnSight files */
          sprintf (filename,"%s",out_dat);
          EnSight (filename,tim,time_steps->nt,grid->nn,grid->ne,grid->nsd,grid->node,grid->element,sup,
                   fv->u_np1,fv->sig,fv->sig_n,fv->eps,opts->smoothing,grid->nce,grid->coel,
                   /*nge,geel,ngn,gnod,*/sol->FNR,sol->arc->lm,ensight,mpi_comm,
                   opts);
          break;
         case VIS_VTK:/* Print to VTK files */
         default: /* no output */ break;
        }/* switch(format) */
      }
    }
  }/* end output */

  return err;
}

/// write restart files
///
/// When check point is active, write restart files. Even if check point is not active but
/// close to walltime (from command line), write restart files because the run is about to
/// finish or be killed.
///
/// \param[in] grid a mesh object
/// \param[in] FV array of field variable object
/// \param[in] load object for loading
/// \param[in] time_steps object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] tim current time step number
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] time_step_start time measure when time stepping starts for step tim
/// \param[in] time_0 time measure when simulation starts
/// \return non-zero on internal error
int write_restart_files(Grid *grid,
                        FieldVariables *FV,
                        LoadingSteps *load,
                        TimeStepping *time_steps,
                        PGFem3D_opt *opts,
                        const Multiphysics& mp,
                        long tim,
                        MPI_Comm mpi_comm,
                        int myrank,
                        double time_step_start,
                        double time_0)
{
  int err = 0;

  int write_restart_global = 0;
  if(time_steps->print[tim] == 1)
    write_restart_global = 1;
  else
  {
    if(opts->walltime>0)
    {
      // Make decesion to write restart when time is close to walltime
      // even if tim is not check point
      double time_step_end = MPI_Wtime();
      double time_taken = time_step_end - time_step_start;

      int write_restart_local  = 0;

      if(opts->walltime - time_taken*3.0 < (time_step_end + time_0))
        write_restart_local = 1;

      MPI_Allreduce (&write_restart_local,&write_restart_global,1,MPI_INT,MPI_MAX,mpi_comm);
      if(write_restart_global>0)
      {
        if(myrank==0)
        {
          printf("INFO: write restart file since PGFem3D is about to done.(walltime = %f[s], now = %f[s], time taken = %f[s])\n",
                 opts->walltime, time_step_end + time_0, time_taken);
        }

        opts->walltime = -1.0;
      }
    }
  }

  if(SAVE_RESTART_FILE && write_restart_global>0)
    write_restart(grid,FV,load,time_steps,opts,mp,myrank,tim);

  return err;
}


int single_scale_main(int argc,char *argv[])
{
  int err = 0;

  const constexpr int ndim = 3;

  /* Create MPI communicator. Currently aliased to MPI_COMM_WORLD but
   * may change */
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  /* CRYSTAL PLASTICITY */
  CRPL *crpl = nullptr;

  /* MI */
  double GVolume = 0.0;
  double oVolume = 0.0; // original volume
  double VVolume = 0.0; // deformed volume

  /* ***** Set up debug log ***** */
  // FILE *debug_log = NULL;
  /* debug_log = fopen("debug.log","w"); */
  // debug_log = stdout;
  /* debug_log = stderr; */

  /*=== END INITIALIZATION === */

  int initialized = 1;
  if (MPI_Initialized(&initialized)) {
    PGFEM_Abort();
  }
  if (!initialized and MPI_Init(&argc, &argv)) {
    PGFEM_Abort();
  }

  int myrank = 0;
  if (MPI_Comm_rank(mpi_comm, &myrank)) {
    PGFEM_Abort();
  }

  int nproc = 0;
  if (MPI_Comm_size(mpi_comm, &nproc)) {
    PGFEM_Abort();
  }

  PGFEM_initialize_io(NULL, NULL);

  if (myrank == 0) {
    PGFEM_printf("=== SINGLE SCALE ANALYSIS ===\n\n");
  }

  double total_time = 0.0;
  double startup_time = 0.0;
  /* MPI_Barrier(mpi_comm); */
  total_time -= MPI_Wtime();
  startup_time -= MPI_Wtime();

#if (MPI_VERSION < 2)
# define MPI_Comm_set_errhandler MPI_Errhandler_set
#endif

  MPI_Comm_set_errhandler(mpi_comm, MPI_ERRORS_ARE_FATAL);

  if (myrank == 0) {
    PGFEM_printf("\n\nInitializing PFEM3d\n\n");
  }

  /*=== Parse the command line for options ===*/
  PGFem3D_opt options;
  if (argc <= 2) {
    if (myrank == 0) {
      print_usage(stdout);
    }
    exit(0);
  }
  set_default_options(&options);
  re_parse_command_line(myrank, 2, argc, argv, &options);
  if (myrank == 0) {
    print_options(stdout, &options);
  }

  //----------------------------------------------------------------------
  // create and initialization of PGFem3D objects
  //----------------------------------------------------------------------
  //---->
  // Multiphysics setting
  int mp_id_M = -1;
  Multiphysics mp;
  err += read_multiphysics_settings(mp,&options,myrank);

  std::vector<FieldVariables> fv(mp.physicsno);
  std::vector<CommunicationStructure> com(mp.physicsno);

  for (int ia = 0; ia < mp.physicsno; ++ia) {
    err += field_varialbe_initialization(&fv[ia]);
    fv[ia].ndofn = mp.ndim[ia];

    if (mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL) {
      mp_id_M = ia;
      switch (options.analysis_type) {
       case STABILIZED:                         //intented to over flow.
       case MINI:
       case MINI_3F:
        fv[ia].ndofn = 4;
        break;
       default:
        break; // do nothing
      }
    }

    fv[ia].n_coupled = mp.coupled_ids[ia][0];
    // create memories for saving coupling info
    if (0 < mp.coupled_ids[ia][0]) {
      fv[ia].coupled_physics_ids = PGFEM_malloc<int>(mp.coupled_ids[ia][0]);
      fv[ia].fvs = PGFEM_malloc<FieldVariables*>(mp.coupled_ids[ia][0]);
    }

    // save coupling info
    for (int ib = 0; ib < mp.coupled_ids[ia][0]; ++ib) {
      int mp_cp_id = mp.coupled_ids[ia][ib+1];
      // tells physics e.g.) fv[ia].coupled_physics_ids[ib] == MULTIPHYSICS_MECHANICAL
      // tells physics e.g.) fv[ia].coupled_physics_ids[ib] == MULTIPHYSICS_THERMAL
      //                                      :                    :
      fv[ia].coupled_physics_ids[ib] = mp.physics_ids[mp_cp_id];
      fv[ia].fvs[ib] = &fv[mp_cp_id];
    }

    err += communication_structure_initialization(&com[ia]);

    com[ia].nproc = nproc;
  }

  TimeStepping time_steps;
  err += time_stepping_initialization(&time_steps);

  Grid grid;
  err += grid_initialization(&grid); // grid.nsd = 3 is the default

  MaterialProperty mat;
  err += material_initialization(&mat);

  LoadingSteps load;
  err += loading_steps_initialization(&load);
  err += construct_loading_steps(&load, mp);

  /* for attaching debugger */
  while (options.debug);

  /* visualization */
  /* Ensight */
  Ensight *ensight = nullptr;
  switch (options.vis_format) {
   case VIS_ENSIGHT:
   case VIS_VTK:
    ensight = new Ensight{};
    break;
   case VIS_NONE:
   break;
   default:
    PGFEM_printerr("Unexpected visualization format %d\n", options.vis_format);
    PGFEM_Abort();
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

  if(make_path(options.opath,DIR_MODE) != 0){
    if(myrank == 0){
      PGFEM_printf("Could not create path (%s)!\n"
                   "Please check input and try again.\n\n",options.opath);
      print_usage(stdout);
    }
    PGFEM_Comm_code_abort(mpi_comm,-1);
  }


  //<---------------------------------------------------------------------
  /* set up solver variables */
  std::vector<Solver> sol(mp.physicsno);

  //----------------------------------------------------------------------
  // read main input files ( *.in)
  //----------------------------------------------------------------------
  //---->
  if (read_mesh_file(&grid, &mat, fv.data(), sol.data(), &load, mp, mpi_comm,
                     &options))
  {
    PGFEM_printerr("[%d]ERROR: incorrectly formatted input file!\n", myrank);
    PGFEM_Abort();
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
  for(int ia=1; ia<mp.physicsno; ia++)
    com[ia].hints = com[0].hints;

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    if(mp.physics_ids[ia]!=MULTIPHYSICS_MECHANICAL)
      continue;

    /*=== OVERRIDE PRESCRIBED DISPLACEMENTS ===*/
    if(options.override_pre_disp){
      if(override_prescribed_displacements(load.sups[ia],&options) != 0){
        PGFEM_printerr("[%d]ERROR: an error was encountered when"
                       " reading the displacement override file.\n"
                       "Be sure that there are enough prescribed"
                       " displacements in the file.\n",myrank);
        PGFEM_Abort();
      }
    }

    /*=== MULTISCALE INFORMATION ===*/
    if(options.multi_scale && (mp_id_M >=0))
    {
      (load.sups[mp_id_M])->multi_scale = options.multi_scale;
      int ms_err = read_interface_macro_normal_lc(options.ipath,load.sups[ia]);
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
    snprintf(bnd_file, sizeof(bnd_file), "%s/%s%d.in.bnd", options.ipath, options.ifname, myrank);
    read_bounding_elements_fname(bnd_file, 3, &(grid.n_be), &(grid.b_elems), mpi_comm);
    bounding_element_set_local_ids(grid.n_be, grid.b_elems, grid.element);
    bounding_element_reverse_mapping(grid.n_be, grid.b_elems, grid.element);
  }

  /*==== ADDITIONAL SETUP ===*/

  /* list of elements with prescribed deflection */
  for(int ia=0; ia<mp.physicsno; ia++)
    list_el_prescribed_def(load.sups[ia],grid.node,grid.element,grid.b_elems,grid.ne,grid.n_be,grid.nn);

  /* list of elements on the COMMUNICATION boundary */
  //build for 0
  com[0].nbndel = 0;
  com[0].bndel = list_boundary_el(grid.ne,grid.element,grid.nn,grid.node,myrank,&(com[0].nbndel));

  for(int ia=0; ia<mp.physicsno; ia++)
  {

    if(mp.physics_ids[ia]==MULTIPHYSICS_MECHANICAL)
    {
      //material matrices (Mechanical part) of the phases
      Mat_3D_orthotropic (mat.nmat,mat.mater,options.analysis_type);

      long ***a = NULL;
      a = aloc3l (mat.nmat,mat.nmat,fv[ia].n_concentrations);
      mat.nhommat = list (a,grid.ne,mat.nmat,fv[ia].n_concentrations,grid.element);

      //alocation of the material matrices
      mat.hommat = build_hommat (mat.nhommat);

      //creates material matrices of the homogeneous medium : LOCAL
      // COORDINATE SYSTEM
      hom_matrices (a,grid.ne,mat.nmat,fv[ia].n_concentrations,grid.element,mat.mater,mat.matgeom,
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

  for(int ia=0; ia<mp.physicsno; ia++)
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

    renumber_global_dof_ids(grid.ne,grid.nce,grid.n_be,grid.nn,fv[ia].ndofn,com[ia].DomDof,grid.node,
                            grid.element,grid.coel,grid.b_elems,mpi_comm,ia);
    com[ia].NBN = distribute_global_dof_ids(grid.ne,grid.nce,grid.n_be,grid.nn,fv[ia].ndofn,ndim,grid.node,
                                            grid.element,grid.coel,grid.b_elems, com[ia].hints, mpi_comm,ia);

    // ALlocate Ap, Ai
    com[ia].Ap = aloc1i(com[ia].DomDof[myrank]+1);
    com[ia].comm  = PGFEM_calloc (COMMUN_1, 1);
    initialize_commun(com[ia].comm );

    com[ia].Ai = Psparse_ApAi(com[ia].nproc,myrank,grid.ne,grid.n_be,grid.nn,fv[ia].ndofn,fv[ia].ndofd,
                              grid.element,grid.b_elems,grid.node,com[ia].Ap,grid.nce,grid.coel,com[ia].DomDof,
                              &(com[ia].GDof),com[ia].comm ,mpi_comm,options.cohesive,com[ia].hints,ia);
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
      if (ia==0) // print onece
        print_PGFem3D_run_info(argc, argv, &grid, &com[ia], &load, &options);

      PGFEM_printf ("---------------------------------------------\n");
      PGFEM_printf ("Physics name: %s\n", mp.physicsname[ia]);
      PGFEM_printf ("Total number of degrees of freedom       : %ld\n", fv[ia].Gndof);
      PGFEM_printf ("Total number of nonzeros in the matrix   : %d\n", APP);
      PGFEM_printf ("Symmetric skyline (including diagonal)   : %ld\n",sky);
    }
    //<---------------------------------------------------------------------
    free(dist);
  }

  dealoc1l (DomNe);
  dealoc1l (DomNn);

  if (myrank == 0 && options.cohesive == 1)
    PGFEM_printf ("Number of cohesive elements              : %ld\n",grid.Gnce);

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    if((load.sups[ia])->npd > 0){
      load.sup_defl[ia] = aloc1((load.sups[ia])->npd);
    } else {
      load.sup_defl[ia] = NULL;
    }
  }

  // initialize hypre, stiffmat, residuals time
  std::vector<double> hypre_time(mp.physicsno);
  std::vector<double> stiffmat_time(mp.physicsno);
  std::vector<double> residuals_time(mp.physicsno);


  { 
    // set for surface tractions
    double *nodal_forces = NULL;
    SURFACE_TRACTION_ELEM *ste = NULL;
    int n_feats = 0;
    int n_sur_trac_elem = 0;

    if(mp_id_M >=0)
    {
      nodal_forces = PGFEM_calloc(double, fv[mp_id_M].ndofd);

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
    err += read_solver_file(&time_steps,&mat,fv.data(),sol.data(),&load,crpl,mp,&options,myrank);
    //<---------------------------------------------------------------------

    if(myrank == 0)
    {
      PGFEM_printf ("\n");
      // Nonlinear solver
      for(int ia=0; ia<mp.physicsno; ia++)
      {
        PGFEM_printf ("NONLINEAR SOLVER (%s): ", mp.physicsname[ia]);
        switch(sol[ia].FNR)
        {
         case 0:
         case 1:
          PGFEM_printf ("NEWTON-RAPHSON METHOD");
          if(sol[ia].set_initial_residual)
            PGFEM_printf (" with computing 1st residual by perturbing disp. with %e", sol[ia].du);

          PGFEM_printf ("\n");
          break;
         case 2:
         case 3:
          if(sol[ia].arc->ARC == 0)
            PGFEM_printf ("ARC-LENGTH METHOD - Crisfield\n");

          if(sol[ia].arc->ARC == 1)
            PGFEM_printf ("ARC-LENGTH METHOD - Simo\n");
          break;
        }
      }
      PGFEM_printf ("\n");
    }

    /* Sparse INITIALIZATION ROUTINES */
    for (int ia = 0, e = mp.physicsno; ia < e; ++ia) {
      sol[ia].system = pgfem3d::solvers::SparseSystem::Create(options, mpi_comm,
                                                              com[ia].Ap,
                                                              com[ia].Ai,
                                                              com[ia].DomDof,
                                                              sol[ia].iter_max_sol,
                                                              sol[ia].err);
    }

    /* alocation of the sigma vector */
    for(int ia=0; ia<mp.physicsno; ia++)
    {
      err += construct_field_varialbe(&fv[ia], &grid, &com[ia], &options, mp, myrank, ia);
      if(mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL) // only mechanical part
      {
        /* push nodal_forces to s->R */
        vvplus(fv[ia].R,nodal_forces,fv[ia].ndofd);

        /* alocation of the eps vector */
        initialize_damage(grid.ne,grid.element,mat.hommat,fv[ia].eps,options.analysis_type);

        if (options.analysis_type == CM || options.analysis_type == CM3F) {
          /* parameter list and initialize const. model at int points.
           * NOTE: should catch/handle returned error flag...
           */
          char *cm_filename = NULL;
          alloc_sprintf(&cm_filename,"%s/model_params.in",options.ipath);
          FILE *cm_in = PGFEM_fopen(cm_filename, "r");
          read_model_parameters_list(mat.nhommat, mat.hommat, cm_in);
          free(cm_filename);
          fclose(cm_in);
          init_all_constitutive_model(fv[ia].eps,grid.ne,grid.element,mat.nhommat,mat.hommat);
          err += prepare_temporal_field_varialbes(&fv[ia],&grid,1,&options);
        }
        else
          err += prepare_temporal_field_varialbes(&fv[ia],&grid,0,&options);

        /* alocation of pressure variables */
        switch(options.analysis_type){
         case TF: // intended not to have break
         case CM3F:
          fv[ia].npres = 1;
          fv[ia].nVol = 1;
          /*
            if(fv[ia].ndofn==3) // discontinuous pressure
            {
            switch(grid.element[0].toe)
            {
            case 8: // P2/P0/V0
            fv[ia].npres = 8;
            fv[ia].nVol = 8;
            break;
            case 10: // Q1/P0/V0
            fv[ia].npres = 4;
            fv[ia].nVol = 4;
            break;
            default:
            fv[ia].npres = 1;
            fv[ia].nVol = 1;
            }
            }
            else // continuous pressure
            {
            switch(grid.element[0].toe)
            {
            case 8: // P2/P1/V0
            fv[ia].npres = 0;
            fv[ia].nVol = 8;
            break;
            case 10: // Q1/P0/V0
            fv[ia].npres = 0;
            fv[ia].nVol = 4;
            break;
            default:
            fv[ia].npres = 0;
            fv[ia].nVol = 1;
            }
            }*/
          break;
         case STABILIZED: case MINI: case MINI_3F:
          if(fv[ia].npres != 4){
            fv[ia].npres = 4;
            if(myrank == 0){
              PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 4.\n"
                           "Re-setting to 4 and continuing...\n");
            }
          }
          break;
         case DISP: // intended not to have break
         case CM:
          if(fv[ia].npres != 0){
            fv[ia].npres = 0;
            if (myrank == 0) {
              PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 0.\n"
                           "Re-setting to 0 and continuing...\n");
            }
          }
          break;
         default:
          if(fv[ia].npres != 1){
            fv[ia].npres = 1;
            if (myrank == 0) {
              PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 1.\n"
                           "Re-setting to 1 and continuing...\n");
            }
          }
          break;
        }/* switch */
        build_pressure_nodes (grid.ne,fv[ia].npres,grid.element,fv[ia].sig,fv[ia].eps,options.analysis_type);
        build_crystal_plast (grid.ne,grid.element,fv[ia].sig,fv[ia].eps,crpl,
                             options.analysis_type,options.plc);

        // \/ initialized element varialbes
        if(options.analysis_type==CM3F || options.analysis_type==TF)
          fv[ia].tf.construct(grid.ne,fv[ia].npres,fv[ia].nVol);
                   
        // /\ initialized element varialbes */        
      }
      else
        err += prepare_temporal_field_varialbes(&fv[ia],&grid,0,&options);
    }


    //----------------------------------------------------------------------
    // set writting output options for Multiphysics
    //----------------------------------------------------------------------
    //---->
    PRINT_MULTIPHYSICS_RESULT *pmr = PGFEM_malloc<PRINT_MULTIPHYSICS_RESULT>(mp.total_write_no);
    err += VTK_construct_PMR(&grid, fv.data(), mp, pmr);
    //<---------------------------------------------------------------------


    //----------------------------------------------------------------------
    // read initial conditions
    //----------------------------------------------------------------------
    //---->
    double tnm1[2] = {-1.0,-1.0};
    
    for(int ia=0; ia<mp.physicsno; ia++)
    {
      // set inital plastic deformation
      if(mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL 
        && (options.analysis_type == CM || options.analysis_type == CM3F))
      {
        set_initial_plastic_deformation_gradient(&grid,fv.data()+ia,&mat,sol.data()+ia,&load,com.data()+ia,
                                                 mpi_comm, &options, mp, ia, myrank);
      }
    }    
        
    err += read_initial_values(&grid,&mat,fv.data(),sol.data(),&load,&time_steps,&options,mp,tnm1,myrank);
    for(int ia=0; ia<mp.physicsno; ia++)
    {      
      // set temporal variables
      for(int ib=0; ib<grid.nn*fv[ia].ndofn; ib++)
      {
        fv[ia].temporal->u_n[ib]   = fv[ia].u_n[ib];
        fv[ia].temporal->u_nm1[ib] = fv[ia].u_nm1[ib];
      }
    }
    //<---------------------------------------------------------------------

    // set the first time step size
    time_steps.dt_np1 = time_steps.times[1] - time_steps.times[0];
    if (time_steps.dt_np1 == 0.0)
    {
      if (myrank == 0){
        PGFEM_printf("Incorrect dt\n");
      }
      PGFEM_Comm_code_abort(mpi_comm,0);
    }

    for(int ia=0; ia<mp.physicsno; ia++)
    {
      if(mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      {
        /* set finite deformations variables */
        set_fini_def(grid.ne,fv[ia].npres,grid.element,fv[ia].eps,fv[ia].sig,options.analysis_type);
        if (options.analysis_type == FS_CRPL)
        {
          set_fini_def_pl(grid.ne,fv[ia].npres,grid.element,fv[ia].eps,fv[ia].sig,crpl,
                          options.analysis_type,options.plc);
        }
      }

      //  NODE (PRESCRIBED DEFLECTION)- SUPPORT COORDINATES generation
      // of the load vector
      err += compute_load_vector_for_prescribed_BC(&grid,&mat,&fv[ia],&sol[ia],&load,time_steps.dt_np1,crpl,
                                                   &options,mp,ia,myrank);

      if(mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      {

        /*  NODE - generation of the load vector  */
        load_vec_node(fv[ia].R,load.nln,ndim,load.znod,grid.node,ia);
        /*  ELEMENT - generation of the load vector  */
        load_vec_elem_sur(fv[ia].R,load.nle_s,ndim,grid.element,load.zele_s);
      }
      /* R   -> Incramental forces
       * RR  -> Total forces for sudivided increment
       * RRn -> Total force after equiblirium */

      vvplus  (fv[ia].f, fv[ia].R,     fv[ia].ndofd);
      vvplus  (fv[ia].RR,fv[ia].f,     fv[ia].ndofd);
      vvminus (fv[ia].f, fv[ia].f_defl,fv[ia].ndofd);

      // set extra variables for arc lengh
      if(sol[ia].FNR == 2 || sol[ia].FNR == 3)
      {
        err += construct_arc_length_variable(sol[ia].arc, &fv[ia], &com[ia], myrank);
        // Transform LOCAL load vector to GLOBAL
        LToG (fv[ia].R,sol[ia].arc->BS_R,myrank,com[ia].nproc,fv[ia].ndofd,com[ia].DomDof,com[ia].GDof,com[ia].comm ,mpi_comm);
        sol[ia].arc->dt0 = time_steps.dt_np1;
        sol[ia].arc->DAL = sol[ia].arc->DLM0 = sol[ia].arc->dAL0;
      }

      for (long ib=0;ib<(load.sups[ia])->npd;ib++)
        load.sup_defl[ia][ib] = (load.sups[ia])->defl_d[ib];
    }

    /*=== NO PERIODIC ===*/
    long tim = 0;

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

    // compute startup time
    startup_time += MPI_Wtime();

    ///////////////////////////////////////////////////////////////////
    // start time stepping
    ///////////////////////////////////////////////////////////////////
    while (time_steps.nt > tim)
    {
      double time_step_start = MPI_Wtime();

      if(tim>options.restart)
      {
        time_steps.tim    = tim;
        time_steps.dt_n   = time_steps.dt_np1;
        time_steps.dt_np1 = time_steps.times[tim+1] - time_steps.times[tim];
        if(time_steps.dt_np1 <= 0.0)
        {
          if(myrank == 0)
            PGFEM_printf("Incorrect dt\n");
          PGFEM_Comm_code_abort(mpi_comm,0);
        }

        if(myrank==0)
        {
          PGFEM_printf("\nFinite deformations time step %ld)  Time %e | dt = %e\n",
                       tim,time_steps.times[tim+1],time_steps.dt_np1);
        }
      }

      /*=== NEWTON RAPHSON ===*/
      if(sol[0].FNR == 0 || sol[0].FNR == 1)
      {
        //----------------------------------------------------------------------
        // file pointer (solver file) is active and used to update loads increments
        // fv[mp_id_M].R   -> Incramental forces
        // fv[mp_id_M].RR  -> Total forces for sudivided increment
        // fv[mp_id_M].RRn -> Total force after equiblirium
        // push nodal_forces to s->R
        //----------------------------------------------------------------------
        //---->
        err += read_and_apply_load_increments(&grid, fv.data(), &load, mp, tim, mpi_comm, myrank);

        if(mp_id_M>=0)
        {
          if(load.tim_load[mp_id_M][tim] == 1 && tim != 0)
            vvplus(fv[mp_id_M].R,nodal_forces,fv[mp_id_M].ndofd);
        }
        //<---------------------------------------------------------------------

        //----------------------------------------------------------------------
        // add load increments util time reaches the restart point
        //----------------------------------------------------------------------
        //---->
        if(tim<options.restart+1)
        {
          for(int ia=0; ia<mp.physicsno; ia++)
          {
            for (long i=0;i<(load.sups[ia])->npd;i++){
              (load.sups[ia])->defl[i] += (load.sups[ia])->defl_d[i];
              (load.sups[ia])->defl_d[i] = 0.0;
            }
          }
          (tim)++;
          continue;
        }
        //<---------------------------------------------------------------------

        for(int ia=0; ia<mp.physicsno; ia++)
          sol[ia].n_step = 0;

        if(tim==options.restart+1 && tnm1[1]>0)
        {
          time_steps.times[tim-1] = tnm1[1]; // tnm1[0] = times[tim-2]
                                             // tnm1[1] = times[tim-1]
                                             // tnm1[2] = times[tim]

                                             // if options.restart==0: tim = 1
          if(tim>=2)                         // if options.restart==1: tim = 2
            time_steps.times[tim-2] = tnm1[0];
        }


        //----------------------------------------------------------------------
        // Perform Newton Raphson interation
        //----------------------------------------------------------------------
        //---->
        fflush(PGFEM_stdout);

        Multiphysics_Newton_Raphson(hypre_time, stiffmat_time, residuals_time, &grid, 
                                    &mat, fv.data(), sol.data(), &load, com.data(), 
                                    &time_steps, crpl, mpi_comm, VVolume, &options, mp);


        for(int ia = 0; ia<mp.physicsno; ia++)
        {
          /* null the prescribed BCs increment */
          nulld(load.sup_defl[ia],(load.sups[ia])->npd);

          /* Null global vectors */
          for (long i=0;i<fv[ia].ndofd;i++){
            fv[ia].RRn[i] += fv[ia].R[i];
            fv[ia].RR[i]   = fv[ia].RRn[i];
            fv[ia].R[i]    = 0.0;
          }
        }
      }/* end NR */

      /*=== ARC LENGTH ===*/

      if(mp_id_M >= 0)
      {
        if(sol[mp_id_M].FNR == 2 || sol[mp_id_M].FNR == 3)
        {
          double dlm = Multiphysics_Arc_length(&grid, &mat, &fv[mp_id_M],
                                               &sol[mp_id_M], &load,
                                               &com[mp_id_M], &time_steps,
                                               crpl, mpi_comm, VVolume,
                                               &options, mp, 0);

          /* Load multiplier */
          sol[mp_id_M].arc->lm += dlm;

          /* Total force vector */
          for (long i=0;i<fv[mp_id_M].ndofd;i++){
            fv[mp_id_M].RR[i] = sol[mp_id_M].arc->lm*fv[mp_id_M].R[i];
          }
        }/* end AL */

        /*=== OUTPUT ===*/
        /* update output stuff for CM interface */
        if((options.analysis_type == CM || options.analysis_type == CM3F) && options.cm!=0)
        {
          constitutive_model_update_output_variables(&grid,
                                                     &mat,
                                                     fv.data(),
                                                     &load,
                                                     &options,
                                                     mp,
                                                     mp_id_M,
                                                     time_steps.dt_np1,
                                                     sol[mp_id_M].alpha);
        }

        /* Calculating equvivalent Mises stresses and strains vectors */
        Mises (grid.ne,fv[mp_id_M].sig,fv[mp_id_M].eps,options.analysis_type);

        /* print tractions on marked features */
        {
          double *sur_forces = NULL;
          if(n_feats > 0){
            sur_forces = PGFEM_calloc(double, n_feats*ndim);
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
      }

      // print simulation results
      err += print_results(&grid, &mat, fv.data(), sol.data(), &load,
                           com.data(), &time_steps, crpl, ensight, pmr,
                           mpi_comm, oVolume, VVolume, &options, mp, tim,
                           myrank);

      err += write_restart_files(&grid, fv.data(), &load, &time_steps, &options,
                                 mp, tim, mpi_comm, myrank, time_step_start,
                                 total_time);

      if (myrank == 0){
        total_EXA_metric += perTimestep_EXA_metric;    //accumulate the numerator for the EXA metric equation per timestep
        perTimestep_EXA_metric = 0;                    //reset EXA metric accumulator
        PGFEM_printf("*********************************************\n");
        PGFEM_printf("*********************************************\n");
      }

      tim++;
    }/* end while */

    if(mp_id_M >=0)
    {
      destroy_applied_surface_traction_list(n_sur_trac_elem,ste);
      free(nodal_forces);
    }
    if(pmr!=NULL) free(pmr);
  }

  //----------------------------------------------------------------------
  // deallocate objects
  //----------------------------------------------------------------------
  //---->
  for(int ia=0; ia<mp.physicsno; ia++) {
    delete sol[ia].system;
  }


  err += destruct_time_stepping(&time_steps);

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    if(mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL &&
       (options.analysis_type == CM || options.analysis_type == CM3F))
      err += destory_temporal_field_varialbes(&fv[ia],1,&options);
    else
      err += destory_temporal_field_varialbes(&fv[ia],0,&options);

    err += destruct_field_varialbe(&fv[ia], &grid, &options, mp, ia);
  }

  err += destruct_loading_steps(&load, mp);
  err += destruct_material(&mat, &options);
  err += destruct_grid(&grid, &options, mp);

  for(int ia=0; ia<mp.physicsno; ia++)
  {
    if(ia>0)
    {
      com[ia].bndel = NULL;
      com[ia].hints = NULL;
    }
    err += destruct_communication_structure(&com[ia]);
    if(mp.physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
    {
      if(sol[ia].FNR == 2 || sol[ia].FNR == 3)
      {
        err += destruct_arc_length_variable(sol[ia].arc);
        free(sol[ia].arc);
      }
    }
  }

  delete ensight;
  //<---------------------------------------------------------------------

  total_time += MPI_Wtime(); // measure time spent

  //----------------------------------------------------------------------
  // print time of analysis and finalize
  //----------------------------------------------------------------------
  //---->
  if (myrank == 0)
    err += print_PGFem3D_final(mp, total_time, startup_time, hypre_time, stiffmat_time, 
                               residuals_time, myrank);

  /* print EXA_metrics */
  print_EXA_metrics(mpi_comm, myrank, hypre_time, residuals_time, total_time, &options);

  err += destruct_multiphysics(mp);
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
