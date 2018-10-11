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
#include "pgfem3d/MultiscaleCommon.hpp"
#include "read_input_file.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

namespace {
const constexpr int ndim = 3;
const constexpr long ARC = 1;
}

int multi_scale_main(int argc, char* argv[])
{
  int err = 0;
  int mp_id = 0;

  // Start the Boot class to get our PMI info
  Boot *boot = new Boot();
  
  /* initialize PGFEM_io */
  PGFEM_initialize_io(NULL,NULL);

  int nproc_macro = 0;
  int full_micro_np = 0;
  int micro_group_size = 0;
  int micro2_group_size = 0;
  int macro_start = 0;
  int macro_argc = 0;
  char **macro_argv = NULL;
  int micro_start = 0;
  int micro_argc = 0;
  char **micro_argv = NULL;
  int micro2_start = 0;
  int micro2_argc = 0;
  char **micro2_argv = NULL;

  int debug = 0;

  PGFem3D_opt options;
  set_default_options(&options);

  /* get macro and micro parts of the command line */
  int myrank = boot->get_rank();
  get_macro_micro_option_blocks(myrank, argc, argv,
                                &macro_start,&macro_argc,
                                &micro_start,&micro2_start,&micro_argc,&micro2_argc,
                                &nproc_macro,&micro_group_size,&micro2_group_size,&full_micro_np,
                                &debug,&options);

  /*=== Parse the command line for global options ===*/
  /*
  printf("argc: %d, macro_argc: %d, micro_argc: %d\n", argc, macro_argc, micro_argc);
  if (macro_argc+micro_argc > argc) {
    if (argc <= 2) {
      if (myrank == 0) {
	print_usage(stdout);
      }
      exit(0);
    }
    re_parse_command_line(myrank, 1, 1, argv, &options);
  }

  if (myrank == 0) {
    print_options(stdout, &options);
  }
  */  

  /* re-adjust for macro/micro options */
  macro_argv = argv + macro_start;
  /* ensure that the microscale gets -ms option. Replace
     '-micro-start' with '-ms' and increment micro_argc */
  micro_argv = argv + micro_start - 1;
  micro_argc++;
  micro2_argv = argv + micro2_start - 1;
  micro2_argc++;
  sprintf(micro_argv[1],"-ms");
  sprintf(micro2_argv[1],"-ms");
  while(debug);
  
  //----------------------------------------------------------------------
  //  // create and initialization of PGFem3D objects
  //    //----------------------------------------------------------------------
  //      //---->
       // Multiphysics setting
  Multiphysics mp;
  err += read_multiphysics_settings(mp,&options,myrank);
              

  // Create the desired network
  Network *net = pgfem3d::net::Network::Create(options);

  // Create the multiscale communicator object
  MultiscaleComm *mscom = new MultiscaleComm(NET_COMM_WORLD, net);

  if (myrank == 0) {
    PGFEM_printf("=== COUPLED MULTISCALE ANALYSIS ===\n\n");
  }

  // split multiscale communicators
  mscom->MM_split(nproc_macro, micro_group_size,micro2_group_size,full_micro_np);
  
  /*=== READ COMM HINTS ===*/
  //allocate memory
  CommunicationStructure *com = new CommunicationStructure{};
  Macroscale *macro = NULL;
  Microscale *micro = NULL;
  Microscale *micro2= NULL;
  
  //load comm hints, macro/micro class, macro/micro filenames
  if (mscom->valid_macro) {
    /*=== MACROSCALE ===*/
    macro = new Macroscale();
    re_parse_command_line(myrank, 1, macro_argc, macro_argv, &options);
    
    int macro_rank;
    net->comm_rank(mscom->macro, &macro_rank);
    com->hints = new CommHints(options.ipath, options.ifname, myrank);
    try {
      (com->hints)->read_filename(NULL);
    } catch(...) {
      if (myrank == 0) {
        PGFEM_printerr("WARNING: One or more procs could not load communication hints.\n"
                       "Proceeding using fallback functions.\n");
      }
    }

    // prepare macro communication
    com->rank = macro_rank;
    com->nproc = nproc_macro;
    com->boot = boot;
    com->net = net;
    com->comm = mscom->macro; // MS communicators in mscom
      PGFEM_redirect_io_null();
    macro->initialize(macro_argc, macro_argv, com, mp_id,mp);
//      PGFEM_redirect_io_macro();

  } else if(mscom->valid_micro_1) {
    /*====== MICROSCALE =======*/
    micro = new Microscale();
    re_parse_command_line(myrank, 1, micro_argc, micro_argv, &options);

    int micro_rank;
    net->comm_rank(mscom->micro, &micro_rank);
    com->hints = new CommHints(options.ipath, options.ifname, micro_rank);
    try {
      (com->hints)->read_filename(NULL);
    } catch(...) {
      if (myrank == 0) {
        PGFEM_printerr("WARNING: One or more procs could not load communication hints.\n"
                       "Proceeding using fallback functions.\n");
      }
    }

    com->rank = micro_rank;
    com->nproc = micro_group_size;
    com->boot = boot;
    com->net = net;
    com->comm = mscom->micro; // MS communicators in mscom
      PGFEM_redirect_io_null();
    micro->initialize(micro_argc, micro_argv, com, mp_id,mp);
//      PGFEM_redirect_io_micro();

  } else if(mscom->valid_micro_2) {
    /*====== MICROSCALE =======*/
    micro2 = new Microscale();
    re_parse_command_line(myrank, 1, micro2_argc, micro2_argv, &options);

    int micro_rank;
    net->comm_rank(mscom->micro2, &micro_rank);
    com->hints = new CommHints(options.ipath, options.ifname, micro_rank);
    try {
      (com->hints)->read_filename(NULL);
    } catch(...) {
      if (myrank == 0) {
        PGFEM_printerr("WARNING: One or more procs could not load communication hints.\n"
                       "Proceeding using fallback functions.\n");
      }
    }

    com->rank = micro_rank;
    com->nproc = micro2_group_size;
    com->boot = boot;
    com->net = net;
    com->comm = mscom->micro; // MS communicators in mscom

    micro2->initialize(micro2_argc, micro2_argv, com, mp_id,mp);
  }
  
  /*=== INITIALIZE SCALES ===*/
  if (mscom->valid_macro) {/*=== MACROSCALE ===*/
    macro->build_solutions(1);  // only 1 for macroscale
    // done twice to undo dependancy loop 
    macro->initialize(macro_argc, macro_argv, com, mp_id,mp); 
  }
  else if (mscom->valid_micro_1) {/*=== MICROSCALE ===*/
    PGFEM_redirect_io_micro();
    /*=== REDIRECT MICROSCALE I/O ===*/
    {
      if (options.custom_micro == 1) {
        char filenameMS[1024];
        char in_dat[1024];
        sprintf(in_dat,"%s/%s",options.ipath,options.ifname);
        //sprintf(filenameMS,"%s.msm",in_dat); //micro simulation method
        read_simulation_methods(filenameMS,micro->opts);
      }
      /* no output */
      PGFEM_redirect_io_null();
      parse_command_line(micro_argc, micro_argv, mscom->rank_micro,
                         micro->opts);
      PGFEM_redirect_io_micro();

      /* create the directory for log output and set logging
         filename */
      int nproc_world;
      net->comm_size(mscom->world, &nproc_world);
      int group_id = mscom->rank_micro_all/micro_group_size;
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
    micro->initialize(micro_argc, micro_argv, com, mp_id,mp);
  } else if (mscom->valid_micro_2) {/*=== MICROSCALE ===*/
    PGFEM_redirect_io_micro();
    /*=== REDIRECT MICROSCALE I/O ===*/
    {
      if (options.custom_micro == 1) {
        char filenameMS[1024];
        char in_dat[1024];
        sprintf(in_dat,"%s/%s",options.ipath,options.ifname);
        sprintf(filenameMS,"%s.msm",in_dat); //micro simulation method
        read_simulation_methods(filenameMS,micro2->opts);
      }
      /* no output */
      PGFEM_redirect_io_null();
      parse_command_line(micro2_argc, micro2_argv, mscom->rank_micro,
                         micro2->opts);
      PGFEM_redirect_io_micro();

      /* create the directory for log output and set logging
 *          filename */
      int nproc_world;
      net->comm_size(mscom->world, &nproc_world);
      int group_id = mscom->rank_micro_all/micro_group_size;
      int dir_len = snprintf(NULL,0,"%s/log",micro2->opts->opath)+1;
      char *dir_name = PGFEM_calloc(char, dir_len);
      sprintf(dir_name,"%s/log",micro2->opts->opath);
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
       micro->initialize(micro2_argc, micro2_argv, com, mp_id,mp);
  } else {
    PGFEM_printerr("[%d]ERROR: neither macro or microscale or microscale 2!\n%s:%s:%d",
                   mscom->rank_world,__func__,__FILE__,__LINE__);
    PGFEM_Abort();
  }

  int n_jobs_max = -1;
  err += pgf_FE2_compute_max_n_jobs(macro, micro, mscom, &n_jobs_max);

  /*=== build the macroscacle clients ===*/
  pgf_FE2_macro_client *client = NULL;

  /*=== Build MICROSCALE server and solutions ===*/
  if (mscom->valid_micro_1) {
    /* allocate space for maximum number of jobs to be computed. */
    micro->build_solutions(n_jobs_max);

    micro->opts->custom_micro = options.custom_micro;
    micro->opts->auto_micro = options.auto_micro;
    /* start the microscale servers. This function does not exit until
       a signal is passed from the macroscale via
       pgf_FE2_macro_client_send_exit */
    err += pgf_FE2_micro_server_START(mscom, micro, mp_id);

    /* destroy the microscale */
    delete micro;

  } else if (mscom->valid_micro_2) {
    /* allocate space for maximum number of jobs to be computed. */
    micro2->build_solutions(n_jobs_max);

    micro2->opts->custom_micro = options.custom_micro;
    micro2->opts->auto_micro = options.auto_micro;
    err += pgf_FE2_micro_server_START(mscom, micro2, mp_id);

    /* destroy the microscale */
    delete micro2;

  } else { /*=== MACROSCALE ===*/
    /* initialize the client */
    pgf_FE2_macro_client_init(&client, net);

    /* create the list of jobs */
    pgf_FE2_macro_client_create_job_list(client, n_jobs_max, macro,
					 mscom, mp_id);
    
    /* determine the initial job assignment*/
    pgf_FE2_macro_client_assign_initial_servers(client, mscom);

    MultiscaleCommon *c = macro;
    MULTISCALE_SOLUTION *s = macro->sol;
    int nproc_macro = 0;
    char filename[500];
    net->comm_size(mscom->macro, &nproc_macro);
    
    /* Create a context for passing stuff to Newton Raphson */
    MS_SERVER_CTX *ctx = PGFEM_calloc(MS_SERVER_CTX, 1);
    ctx->client = client;
    ctx->mscom = mscom;
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
                                           &feat_type,&feat_id,&loads,
					   myrank);

      generate_applied_surface_traction_list(c->ne,c->elem,
                                             n_feats,feat_type,
                                             feat_id,&n_sur_trac_elem,
                                             &ste);

      compute_applied_traction_res(c->ndofn,c->node,c->elem,
                                   n_sur_trac_elem,ste,
                                   n_feats,loads,s->eps,
                                   nodal_forces,mp_id);

      double tmp_sum = 0.0;
      for(int i=0; i<c->ndofd; i++){
        tmp_sum += nodal_forces[i];
      }
      net->allreduce(NET_IN_PLACE, &tmp_sum, 1, NET_DT_DOUBLE,
		     NET_OP_SUM, mscom->macro);

      if(mscom->rank_macro == 0){
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
      if (mscom->rank_macro == 0) {
        PGFEM_printf("Overriding the default solver file with:\n%s\n",
                     macro->opts->solver_file);
      }
      solver_file_open(macro->opts->solver_file, &solver_file);
    } else {
      /* use the default file/filename */
      char *filename = NULL;
      alloc_sprintf (&filename,"%s/%s%d.in.st",macro->opts->ipath,
                     macro->opts->ifname,mscom->rank_macro);
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
    if (mscom->rank_macro == 0) {
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
      if (mscom->rank_macro == 0) {
        PGFEM_printf("Restarting from step %d\n\n",macro->opts->restart);
      }

      /* increment load to restart step */
      solver_file_scan_to_step(solver_file,macro->opts->restart,
                               c->supports->npd,c->supports->defl_d);

      /* read restart files and set current equilibrium state */
      pgf_FE2_restart_read_macro(macro, macro->opts->restart, mp_id);
      s->tim = macro->opts->restart;

      /* send a job to compute the first tangent */
      pgf_FE2_macro_client_send_jobs(client, mscom, macro,
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
          net->allreduce(NET_IN_PLACE, sur_forces, n_feats*ndim,
			 NET_DT_DOUBLE, NET_OP_SUM, mscom->macro);
          if(mscom->rank_macro == 0){
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
      pgf_FE2_macro_client_send_jobs(client, mscom, macro,
                                     JOB_NO_COMPUTE_EQUILIBRIUM);

      /*  NODE (PRESCRIBED DEFLECTION)- SUPPORT COORDINATES generation
          of the load vector  */
      s->dt = s->times[s->tim + 1] - s->times[s->tim];
      if (s->dt == 0.0){
        if (mscom->rank_macro == 0){
          PGFEM_printf("Incorrect dt\n");
        }
        PGFEM_Abort();
      }

      compute_load_vector_for_prescribed_BC_multiscale(c,s,macro->opts,
						       solver_file->nonlin_tol,
                                                       mscom->rank_macro);

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
      LToG (s->R, s->BS_R, c->ndofd, c);
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
        if (mscom->rank_macro == 0) {
          PGFEM_printf("Incorrect dt\n");
        }
        PGFEM_Abort();
      }

      if (mscom->rank_macro == 0) {
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

        if (mscom->rank_macro == 0 && load_err){
          PGFEM_printf ("Incorrect load input for Time = 0\n");
          PGFEM_Abort();
        }

        int n_step = 0;

        // perform Newton Raphson iteration
        hypre_time[mp_id] += Newton_Raphson_multiscale(1,c,s,solver_file,ctx,macro->opts,sup_defl,&pores,&n_step);

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
                                    ARC,&ITT,&dAL,sup_defl);

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
          net->allreduce(NET_IN_PLACE,sur_forces,n_feats*ndim,
			 NET_DT_DOUBLE,NET_OP_SUM,mscom->macro);
          if(mscom->rank_macro == 0){
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

      fd_res_compute_reactions_multiscale(c,s,solver_file,macro->opts,dts);

      if (solver_file->print_steps[s->tim] == 1){

        /* do not transfer data between servers */
        pgf_FE2_macro_client_rebalance_servers(client,mscom,
                                               FE2_REBALANCE_NONE);

        /* Send print jobs */
        pgf_FE2_macro_client_send_jobs(client,mscom,macro,JOB_PRINT);

        if (macro->opts->vis_format != VIS_NONE ) {

          /* *ADDITIONAL* ASCII formatted output.
           * *NOT* VTK ASCII format */
          if(macro->opts->ascii){
            int Gnn = 0;
            ASCII_output(macro->opts,c,s->tim,s->times,
                         Gnn,c->nn,c->ne,c->nce,c->ndofd,
			 solver_file->nonlin_method,
                         lm,pores,c->VVolume,
                         c->node,c->elem,c->supports,
                         s->r,s->eps,s->sig_e,s->sig_n,c->coel);
          } /* End ASCII output */

          switch(macro->opts->vis_format){
           case VIS_ELIXIR:/* Print to elix file */
            sprintf (filename,"%s/%s_%d.elx%d",macro->opts->opath,
                     macro->opts->ofname,mscom->rank_macro,s->tim);
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
                     solver_file->nonlin_method,lm,c->ensight,c,
                     macro->opts);
            break;
           case VIS_VTK:/* Print to VTK files */
            if(mscom->rank_macro == 0){
              VTK_print_master(macro->opts->opath,macro->opts->ofname,
                               s->tim,nproc_macro,macro->opts);
            }

            VTK_print_vtu(macro->opts->opath,macro->opts->ofname,s->tim,
                          mscom->rank_macro,c->ne,c->nn,c->node,c->elem,
                          c->supports,s->r,NULL,NULL,s->sig_e,s->eps,
                          macro->opts,mp_id);

            if (macro->opts->cohesive == 1){
              if(mscom->rank_macro == 0){
                VTK_print_cohesive_master(macro->opts->opath,
                                          macro->opts->ofname,
                                          s->tim,nproc_macro,macro->opts);
              }

              VTK_print_cohesive_vtu(macro->opts->opath,macro->opts->ofname,
                                     s->tim,mscom->rank_macro,c->nce,c->node,
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

      if (mscom->rank_macro == 0){
        PGFEM_printf("\n");
        PGFEM_printf("*********************************************\n");
        PGFEM_printf("*********************************************\n");
      }

      s->tim++;
    }/* end while */

    /*=== SEND EXIT SIGNAL TO MICROSCALE SERVER ===*/
    pgf_FE2_macro_client_send_exit(client,mscom);

    /* cleanup */
    free(ctx);

    free(sup_defl);
    destroy_applied_surface_traction_list(n_sur_trac_elem,ste);

    /* destroy the macroscale client */
    pgf_FE2_macro_client_destroy(client);

    /* destroy the macroscale */
    delete macro;
  } /*=== END OF COMPUTATIONS ===*/

  /*=== PRINT TIME OF ANALYSIS ===*/
  if (mscom->rank_world == 0){
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    PGFEM_printf ("\n");
    PGFEM_printf ("Time of analysis on processor [%d] - "
                  " System %ld.%ld, User %ld.%ld.\n",
                  mscom->rank_world,usage.ru_stime.tv_sec,
                  usage.ru_stime.tv_usec,
                  usage.ru_utime.tv_sec,usage.ru_utime.tv_usec);
  }

  err += PGFEM_finalize_io();
  
  delete mscom;
  delete com;

  net->finalize();
  
  delete net;
  delete boot;
  /* finalize and exit */
  return err;
}
