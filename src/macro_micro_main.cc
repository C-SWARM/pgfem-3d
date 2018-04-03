/*** This is the main function for the fully-coupled multiscale modeling */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Arc_length.h"
#include "Newton_Raphson.h"
#include "PGFEM_io.h"
#include "PFEM3d.h"
#include "PGFem3D_options.h"
#include "Printing.h"
#include "Psparse_ApAi.h"
#include "SetGlobalNodeNumbers.h"
#include "allocation.h"
#include "applied_traction.h"
#include "build_distribution.h"
#include "computeMacroF.h"
#include "computeMacroS.h"
#include "dynamics.h"
#include "enumerations.h"
#include "fd_residuals.h"
#include "gen_path.h"
#include "homogen.h"
#include "in.h"
#include "incl.h"
#include "interface_macro.h"
#include "load.h"
#include "macro_micro_functions.h"
#include "matice.h"
#include "matrix_printing.h"
#include "microscale_information.h"
#include "ms_cohe_job_list.h"
#include "out.h"
#include "pgf_fe2_macro_client.h"
#include "pgf_fe2_micro_server.h"
#include "pgf_fe2_restart.h"
#include "pgf_fe2_compute_max_n_jobs.h"
#include "print_dist.h"
#include "profiler.h"
#include "read_cryst_plast.h"
#include "renumber_ID.h"
#include "set_fini_def.h"
#include "skyline.h"
#include "solver_file.h"
#include "utils.h"
#include "vtk_output.h"
#include <cassert>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

namespace {
const constexpr int ndim = 3;
const constexpr long ARC = 1;
}

int multi_scale_main(int argc, char* argv[])
{
  int err = 0;
  int mp_id = 0;
  /* intitialize MPI */
  if (MPI_Init(&argc, &argv)) {
    PGFEM_Abort();
  }

  /* initialize PGFEM_io */
  PGFEM_initialize_io(NULL,NULL);

  int nproc_macro = 0;
  int micro_group_size = 0;
  int macro_start = 0;
  int macro_argc = 0;
  char **macro_argv = NULL;
  int micro_start = 0;
  int micro_argc = 0;
  char **micro_argv = NULL;
  int debug = 0;
  int EXA_metric = 0;

  /* get macro and micro parts of the command line */
  int rank_world = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_world);
  get_macro_micro_option_blocks(rank_world, argc, argv,
                                &macro_start,&macro_argc,
                                &micro_start,&micro_argc,
                                &nproc_macro,&micro_group_size,
                                &debug);

  macro_argv = argv + macro_start;

  /* ensure that the microscale gets -ms option. Replace
     '-micro-start' with '-ms' and increment micro_argc */
  micro_argv = argv + micro_start - 1;
  micro_argc++;
  sprintf(micro_argv[1],"-ms");

  while(debug);

  PGFEM_mpi_comm *mpi_comm = PGFEM_calloc(PGFEM_mpi_comm, 1);
  err += initialize_PGFEM_mpi_comm(MPI_COMM_WORLD,mpi_comm);

  if(mpi_comm->rank_world == 0){
    PGFEM_printf("=== COUPLED MULTISCALE ANALYSIS ===\n\n");
  }

  err += PGFEM_mpi_comm_MM_split(nproc_macro,micro_group_size,mpi_comm);


  /*=== READ COMM HINTS ===*/
  //allocate memory
  CommunicationStructure *com = new CommunicationStructure{};
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  PGFem3D_opt options;
  MACROSCALE *macro = NULL;
  MICROSCALE *micro = NULL;
  //load comm hints, macro/micro class, macro/micro filenames

  if (mpi_comm->valid_macro) {
    /*=== MACROSCALE ===*/
    re_parse_command_line(myrank, 1, macro_argc, macro_argv, &options);

    const char* fn = Comm_hints_filename(options.ipath, options.ifname, myrank);
    com[0].hints = Comm_hints_construct();
    int ch_err = Comm_hints_read_filename(com[0].hints, fn);

    initialize_MACROSCALE(&macro);
    build_MACROSCALE(macro,mpi_comm->macro,macro_argc,macro_argv,mp_id,com[mp_id].hints);
    if (ch_err) {
      Comm_hints_destroy(com[0].hints);
      com[0].hints = NULL;
      if (myrank == 0) {
        PGFEM_printerr("WARNING: One or more procs could not load communication hints.\n"
                       "Proceeding using fallback functions.\n");
      }
    }

  }
  else {
    /*====== MICROSCALE =======*/
    re_parse_command_line(myrank, 1, micro_argc, micro_argv, &options);


    int micro_rank;
    MPI_Comm_rank(mpi_comm->micro, &micro_rank);
    const char* fn = Comm_hints_filename(options.ipath, options.ifname, micro_rank);
    com[0].hints = Comm_hints_construct();
    int ch_err = Comm_hints_read_filename(com[0].hints, fn);

    initialize_MICROSCALE(&micro);
    build_MICROSCALE(micro, mpi_comm->micro, micro_argc, micro_argv, mp_id,
                     com[mp_id].hints);

    if (ch_err) {
      Comm_hints_destroy(com[0].hints);
      com[0].hints = NULL;
      if (myrank == 0) {
        PGFEM_printerr("WARNING: One or more procs could not load communication hints.\n"
                       "Proceeding using fallback functions.\n");
      }
    }

  }

  /*=== INITIALIZE SCALES ===*/
  if (mpi_comm->valid_macro) {/*=== MACROSCALE ===*/
    initialize_MACROSCALE(&macro);
    build_MACROSCALE(macro, mpi_comm->macro, macro_argc, macro_argv, mp_id,
                     com[mp_id].hints);
    build_MACROSCALE_solution(macro);
  }
  else if (mpi_comm->valid_micro) {/*=== MICROSCALE ===*/
    PGFEM_redirect_io_micro();
    initialize_MICROSCALE(&micro);

    /*=== REDIRECT MICROSCALE I/O ===*/
    {
      /* no output */
      PGFEM_redirect_io_null();
      parse_command_line(micro_argc, micro_argv, mpi_comm->rank_micro,
                         micro->opts);
      PGFEM_redirect_io_micro();

      /* create the directory for log output and set logging
         filename */
      int nproc_world = 0;
      MPI_Comm_size(mpi_comm->world,&nproc_world);
      int group_id = mpi_comm->rank_micro_all/micro_group_size;
      int dir_len = snprintf(NULL,0,"%s/log",micro->opts->opath)+1;
      char *dir_name = PGFEM_calloc(char, dir_len);
      sprintf(dir_name,"%s/log",micro->opts->opath);
      make_path(dir_name,DIR_MODE);
      dir_len = snprintf(NULL,0,"%s/group_%05d",dir_name,group_id)+1;
      char *fname = PGFEM_calloc(char, dir_len);
      sprintf(fname,"%s/group_%05d",dir_name,group_id);

      /* reinitialize I/O */
      PGFEM_finalize_io();
      PGFEM_initialize_io(NULL,fname);
      PGFEM_redirect_io_micro();
      free(dir_name);
      free(fname);
    }

    /*=== BUILD MICROSCALE ===*/
    build_MICROSCALE(micro, mpi_comm->micro, micro_argc, micro_argv, mp_id,
                     com[mp_id].hints);
  } else {
    PGFEM_printerr("[%d]ERROR: neither macro or microscale!\n%s:%s:%d",
                   mpi_comm->rank_world,__func__,__FILE__,__LINE__);
    PGFEM_Abort();
  }

  int n_jobs_max = -1;
  err += pgf_FE2_compute_max_n_jobs(macro,mpi_comm,&n_jobs_max);

  /*=== build the macroscacle clients ===*/
  pgf_FE2_macro_client *client = NULL;

  /*=== Build MICROSCALE server and solutions ===*/
  if (mpi_comm->valid_micro) {
    /* allocate space for maximum number of jobs to be computed. */
    build_MICROSCALE_solutions(n_jobs_max,micro);

    /* start the microscale servers. This function does not exit until
       a signal is passed from the macroscale via
       pgf_FE2_macro_client_send_exit */
    err += pgf_FE2_micro_server_START(mpi_comm,micro,mp_id,EXA_metric);

    /* destroy the microscale */
    destroy_MICROSCALE(micro);

  } else { /*=== MACROSCALE ===*/
    /* initialize the client */
    pgf_FE2_macro_client_init(&client);

    /* create the list of jobs */
    pgf_FE2_macro_client_create_job_list(client,n_jobs_max,macro,mpi_comm,mp_id);

    /* determine the initial job assignment*/
    pgf_FE2_macro_client_assign_initial_servers(client,mpi_comm);

    COMMON_MACROSCALE *c = macro->common;
    MACROSCALE_SOLUTION *s = macro->sol;
    int nproc_macro = 0;
    char filename[500];
    MPI_Comm_size(mpi_comm->macro,&nproc_macro);

    /* Create a context for passing stuff to Newton Raphson */
    MS_SERVER_CTX *ctx = PGFEM_calloc(MS_SERVER_CTX, 1);
    ctx->client = client;
    ctx->mpi_comm = mpi_comm;
    ctx->macro = macro;

    /*=== COMPUTE APPLIED FORCES ON MARKED SURFACES ===*/
    double *nodal_forces = PGFEM_calloc(double, c->ndofd);
    SURFACE_TRACTION_ELEM *ste = NULL;
    int n_feats = 0;
    int n_sur_trac_elem = 0;
    {
      int *feat_type = NULL;
      int *feat_id = NULL;
      double *loads = NULL;
      char *trac_fname = NULL;
      alloc_sprintf(&trac_fname,"%s/traction.in",macro->opts->ipath);

      read_applied_surface_tractions_fname(trac_fname,&n_feats,
                                           &feat_type,&feat_id,&loads);

      generate_applied_surface_traction_list(c->ne,c->elem,
                                             n_feats,feat_type,
                                             feat_id,&n_sur_trac_elem,
                                             &ste);

      compute_applied_traction_res(c->ndofn,c->node,c->elem,
                                   n_sur_trac_elem,ste,
                                   n_feats,loads,
                                   nodal_forces,mp_id);

      double tmp_sum = 0.0;
      for(int i=0; i<c->ndofd; i++){
        tmp_sum += nodal_forces[i];
      }
      MPI_Allreduce(MPI_IN_PLACE,&tmp_sum,1,MPI_DOUBLE,
                    MPI_SUM,mpi_comm->macro);

      if(mpi_comm->rank_macro == 0){
        PGFEM_printf("Total load from surface tractions: %.8e\n\n",tmp_sum);
      }

      free(feat_type);
      free(feat_id);
      free(loads);
      free(trac_fname);
    }

    /* push nodal_forces to s->R */
    vvplus  (s->R,nodal_forces,c->ndofd);

    /*=== SOLUTION PROCESS ===*/
    /*=== READ SOLVER FILE ===*/
    SOLVER_FILE *solver_file = NULL;
    if (macro->opts->override_solver_file) {
      if (mpi_comm->rank_macro == 0) {
        PGFEM_printf("Overriding the default solver file with:\n%s\n",
                     macro->opts->solver_file);
      }
      solver_file_open(macro->opts->solver_file, &solver_file);
    } else {
      /* use the default file/filename */
      char *filename = NULL;
      alloc_sprintf (&filename,"%s/%s%d.in.st",macro->opts->ipath,
                     macro->opts->ifname,mpi_comm->rank_macro);
      solver_file_open(filename, &solver_file);
      free(filename);
    }
    s->tim = 0;
    solver_file_read_header(solver_file);

    /* allocate macro_solution times and copy from solver_file */
    if(solver_file->n_step > 2){
      free(s->times);
      s->times = PGFEM_malloc<double>(solver_file->n_step + 1);
    }
    memcpy(s->times,solver_file->times, ((solver_file->n_step + 1) *
                                         sizeof(*(s->times))));

    /* Nonlinear solver */
    if (mpi_comm->rank_macro == 0) {
      switch (solver_file->nonlin_method) {
       case NEWTON_METHOD:
        PGFEM_printf ("\nNONLINEAR SOLVER : NEWTON-RAPHSON METHOD\n");
        break;

       case ARC_LENGTH_METHOD:
       case AUX_ARC_LENGTH_METHOD:
        if (ARC == 0)
          PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Crisfield\n");
        else if (ARC == 1)
          PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Simo\n");
        break;

       default:
        PGFEM_printerr("Undefined method! ABORT\n");
        PGFEM_Abort();
        break;
      }
    }

    // initialize hypre time
    std::vector<double> hypre_time(mp_id + 1);

    if (macro->opts->restart >= 0) {
      if (mpi_comm->rank_macro == 0) {
        PGFEM_printf("Restarting from step %d\n\n",macro->opts->restart);
      }

      /* increment load to restart step */
      solver_file_scan_to_step(solver_file,macro->opts->restart,
                               c->supports->npd,c->supports->defl_d);

      /* read restart files and set current equilibrium state */
      pgf_FE2_restart_read_macro(macro,macro->opts->restart,mp_id);
      s->tim = macro->opts->restart;

      /* send a job to compute the first tangent */
      pgf_FE2_macro_client_send_jobs(client,mpi_comm,macro,
                                     JOB_NO_COMPUTE_EQUILIBRIUM);

      /* turn off restart at the macroscale */
      macro->opts->restart = -1;

      /* print tractions on marked features */
      {
        double *sur_forces = NULL;
        if(n_feats > 0){
          sur_forces = PGFEM_calloc(double, n_feats*ndim);
          compute_resultant_force(n_feats,n_sur_trac_elem,
                                  ste,c->node,c->elem,
                                  s->sig_e,s->eps,sur_forces);
          MPI_Allreduce(MPI_IN_PLACE,sur_forces,n_feats*ndim,
                        MPI_DOUBLE,MPI_SUM,mpi_comm->macro);
          if(mpi_comm->rank_macro == 0){
            PGFEM_printf("RESTART: Forces on marked features:\n");
            print_array_d(PGFEM_stdout,sur_forces,n_feats*ndim,
                          n_feats,ndim);
          }
        }
        free(sur_forces);
      }

      /* increment macroscale time */
      s->tim ++;

    } else {
      /* not restarting, need to compute initial load/RHS */

      /* send signal to microscale to compute initial tangent */
      pgf_FE2_macro_client_send_jobs(client,mpi_comm,macro,
                                     JOB_NO_COMPUTE_EQUILIBRIUM);

      /*  NODE (PRESCRIBED DEFLECTION)- SUPPORT COORDINATES generation
          of the load vector  */
      s->dt = s->times[s->tim + 1] - s->times[s->tim];
      if (s->dt == 0.0){
        if (mpi_comm->rank_macro == 0){
          PGFEM_printf("Incorrect dt\n");
        }
        PGFEM_Abort();
      }

      compute_load_vector_for_prescribed_BC_multiscale(c,s,macro->opts,solver_file->nonlin_tol,
                                                       mpi_comm->rank_macro,EXA_metric);

      /*=== do not support node/surf loads ===*/
      /* /\*  NODE - generation of the load vector  *\/ */
      /* load_vec_node (R,nln,ndim,znod,node); */
      /* /\*  ELEMENT - generation of the load vector  *\/ */
      /* load_vec_elem_sur (R,nle_s,ndim,elem,zele_s); */

      /* R   -> Incramental forces
         RR  -> Total forces for sudivided increment
         RRn -> Total force after equiblirium */

      vvplus  (s->f,s->R,c->ndofd);
      vvplus  (s->RR,s->f,c->ndofd);
      vvminus (s->f,s->f_defl,c->ndofd);

      /* Transform LOCAL load vector to GLOBAL */
      LToG (s->R,s->BS_R,mpi_comm->rank_macro,nproc_macro,
            c->ndofd,c->DomDof,c->GDof,c->pgfem_comm,
            c->mpi_comm);
    }

    /* Prescribed deflection */
    double *sup_defl = NULL;
    assert(c->supports->npd > 0);
    sup_defl = PGFEM_calloc(double, c->supports->npd);

    double pores = 0.0;
    double lm = 0.0;
    double dlm = 0.0;
    double DET = 0.0;
    double dt0 = 0.0;
    double dlm0 = 0.0;
    double DLM = 0.0;
    double dAL = 0.0;
    long AT = 0;
    long ITT = 0;
    long init = 0;

    /*=== BEGIN SOLVE ===*/
    assert(s->tim >= 0);
    while (solver_file->n_step > (unsigned)s->tim) {
      s->dt = dt0 = s->times[s->tim+1] - s->times[s->tim];
      if (s->dt <= 0.0){
        if (mpi_comm->rank_macro == 0) {
          PGFEM_printf("Incorrect dt\n");
        }
        PGFEM_Abort();
      }

      if (mpi_comm->rank_macro == 0) {
        PGFEM_printf("\nFinite deformations time step %ld) "
                     " Time %f | dt = %10.10f\n",
                     s->tim,s->times[s->tim+1],s->dt);
      }

      /*=== NEWTON RAPHSON ===*/
      if (solver_file->nonlin_method == NEWTON_METHOD){
        int load_err = 0;

        if(s->tim > 0){ /* get load increment */
          load_err = solver_file_read_load(solver_file,s->tim,
                                           c->supports->npd,
                                           c->supports->defl_d);
        }

        /* copy the load increment */
        memcpy(sup_defl,c->supports->defl_d,c->supports->npd*sizeof(double));

        if (mpi_comm->rank_macro == 0 && load_err){
          PGFEM_printf ("Incorrect load input for Time = 0\n");
          PGFEM_Abort();
        }

        int n_step = 0;

        // perform Newton Raphson iteration
        hypre_time[mp_id] += Newton_Raphson_multiscale(1,c,s,solver_file,ctx,macro->opts,sup_defl,&pores,&n_step,EXA_metric);

        /* Null global vectors */
        for (int i=0;i<c->ndofd;i++){
          s->RRn[i] += s->R[i];
          s->RR[i] = s->RRn[i];
          s->R[i] = 0.0;
        }

        /* null the prescribed displacement increment */
        nulld(sup_defl,c->supports->npd);
      }/* end NR */

      /*=== ARC LENGTH ===*/
      if (solver_file->nonlin_method == ARC_LENGTH_METHOD
          || solver_file->nonlin_method == AUX_ARC_LENGTH_METHOD){

        /* initialize arc-length specific variables */
        if(!init){
          lm = dlm = DLM = DET = 0.0;
          dAL = dlm0 = solver_file->nonlin_method_opts[0];
          AT = ITT = 0;
          init++;
        }
        // char out_dat[500];
        // double tmp_val = ((s->times[s->tim+1]-s->times[s->tim])
        //                   /dt0*solver_file->nonlin_method_opts[1]);

        dlm = Arc_length_multiscale(c,s,solver_file,macro->opts,
                                    &pores,dt0,lm,&DET,&dlm0,&DLM,&AT,
                                    ARC,&ITT,&dAL,sup_defl,EXA_metric);

        /* Load multiplier */
        lm += dlm;

        /* Total force vector */
        for (int i=0;i<c->ndofd;i++){
          s->RR[i] = lm*s->R[i];
        }
      }/* end AL */

      /*=== OUTPUT ===*/
      /* Calculating equvivalent Mises stresses and strains vectors */
      Mises (c->ne,s->sig_e,s->eps,macro->opts->analysis_type);

      /* print tractions on marked features */
      {
        double *sur_forces = NULL;
        if(n_feats > 0){
          sur_forces = PGFEM_calloc(double, n_feats*ndim);
          compute_resultant_force(n_feats,n_sur_trac_elem,
                                  ste,c->node,c->elem,
                                  s->sig_e,s->eps,sur_forces);
          MPI_Allreduce(MPI_IN_PLACE,sur_forces,n_feats*ndim,
                        MPI_DOUBLE,MPI_SUM,mpi_comm->macro);
          if(mpi_comm->rank_macro == 0){
            PGFEM_printf("Forces on marked features:\n");
            print_array_d(PGFEM_stdout,sur_forces,n_feats*ndim,
                          n_feats,ndim);
          }
        }
        free(sur_forces);
      }

      double dts[2];

      if(s->tim==0)
        dts[DT_N] = s->times[s->tim+1] - s->times[s->tim];
      else
        dts[DT_N] = s->times[s->tim] - s->times[s->tim-1];

      dts[DT_NP1] = s->times[s->tim+1] - s->times[s->tim];

      fd_res_compute_reactions_multiscale(c,s,solver_file,macro->opts,dts,EXA_metric);

      if (solver_file->print_steps[s->tim] == 1){

        /* do not transfer data between servers */
        pgf_FE2_macro_client_rebalance_servers(client,mpi_comm,
                                               FE2_REBALANCE_NONE);

        /* Send print jobs */
        pgf_FE2_macro_client_send_jobs(client,mpi_comm,macro,JOB_PRINT);

        if (macro->opts->vis_format != VIS_NONE ) {

          /* *ADDITIONAL* ASCII formatted output.
           * *NOT* VTK ASCII format */
          if(macro->opts->ascii){
            int Gnn = 0;
            ASCII_output(macro->opts,c->mpi_comm,s->tim,s->times,
                         Gnn,c->nn,c->ne,c->nce,c->ndofd,
                         c->DomDof,c->Ap,solver_file->nonlin_method,
                         lm,pores,c->VVolume,
                         c->node,c->elem,c->supports,
                         s->r,s->eps,s->sig_e,s->sig_n,c->coel);
          } /* End ASCII output */

          switch(macro->opts->vis_format){
           case VIS_ELIXIR:/* Print to elix file */
            sprintf (filename,"%s/%s_%d.elx%d",macro->opts->opath,
                     macro->opts->ofname,mpi_comm->rank_macro,s->tim);
            elixir (filename,c->nn,c->ne,ndim,c->node,c->elem,
                    c->supports,s->r,s->sig_e,s->sig_n,s->eps,
                    macro->opts->smoothing,c->nce,c->coel,macro->opts);
            break;
           case VIS_ENSIGHT:/* Print to EnSight files */
            sprintf (filename,"%s/%s",macro->opts->opath,
                     macro->opts->ofname);
            EnSight (filename,s->tim,solver_file->n_step,c->nn,c->ne,ndim,c->node,
                     c->elem,c->supports,s->r,s->sig_e,s->sig_n,s->eps,
                     macro->opts->smoothing,c->nce,c->coel,
                     solver_file->nonlin_method,lm,c->ensight,c->mpi_comm,
                     macro->opts);
            break;
           case VIS_VTK:/* Print to VTK files */
            if(mpi_comm->rank_macro == 0){
              VTK_print_master(macro->opts->opath,macro->opts->ofname,
                               s->tim,nproc_macro,macro->opts);
            }

            VTK_print_vtu(macro->opts->opath,macro->opts->ofname,s->tim,
                          mpi_comm->rank_macro,c->ne,c->nn,c->node,c->elem,
                          c->supports,s->r,NULL,NULL,s->sig_e,s->eps,
                          macro->opts,mp_id);

            if (macro->opts->cohesive == 1){
              if(mpi_comm->rank_macro == 0){
                VTK_print_cohesive_master(macro->opts->opath,
                                          macro->opts->ofname,
                                          s->tim,nproc_macro,macro->opts);
              }

              VTK_print_cohesive_vtu(macro->opts->opath,macro->opts->ofname,
                                     s->tim,mpi_comm->rank_macro,c->nce,c->node,
                                     c->coel,c->supports,s->r,c->ensight,
                                     macro->opts,mp_id);
            }
            break;
           default: /* no output */ break;
          }/* switch(format) */
        } /* if !VIZ_NONE */

        /* dump a restart file */
        int re_err = pgf_FE2_restart_print_macro(macro);
        assert(re_err == 0);

        /* complete communication cycle w/ microscale */
        int junk = 0;
        pgf_FE2_macro_client_recv_jobs(client,macro,&junk);
      }/* end output */

      if (mpi_comm->rank_macro == 0){
        PGFEM_printf("\n");
        PGFEM_printf("*********************************************\n");
        PGFEM_printf("*********************************************\n");
      }

      s->tim++;
    }/* end while */

    /*=== SEND EXIT SIGNAL TO MICROSCALE SERVER ===*/
    pgf_FE2_macro_client_send_exit(client,mpi_comm);

    /* cleanup */
    free(ctx);

    free(sup_defl);
    destroy_applied_surface_traction_list(n_sur_trac_elem,ste);

    /* destroy the macroscale client */
    pgf_FE2_macro_client_destroy(client);

    /* destroy the macroscale */
    destroy_MACROSCALE(macro);
  } /*=== END OF COMPUTATIONS ===*/

  /*=== PRINT TIME OF ANALYSIS ===*/
  if(mpi_comm->rank_world == 0){
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    PGFEM_printf ("\n");
    PGFEM_printf ("Time of analysis on processor [%d] - "
                  " System %ld.%ld, User %ld.%ld.\n",
                  mpi_comm->rank_world,usage.ru_stime.tv_sec,
                  usage.ru_stime.tv_usec,
                  usage.ru_utime.tv_sec,usage.ru_utime.tv_usec);
  }

  /* destroy the PGFEM communicator */
  err += destroy_PGFEM_mpi_comm(mpi_comm);
  free(mpi_comm);
  delete com;
  /* finalize and exit */
  err += PGFEM_finalize_io();
  err += MPI_Finalize();
  return err;
}
