/*** This is the main function for the fully-coupled multiscale modeling */


#include "PFEM3d.h"
#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

/* Standard headers/libs */
#include <time.h>
#include <stdlib.h>
#include <sys/time.h> 
#include <sys/resource.h>


/* Extra libs */
#include "renumbering.h"

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
#include "RNPsparse_ApAi.h"
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
#include "microscale_information.h"
#include "ms_cohe_job_list.h"
#include "applied_traction.h"
#include "macro_micro_functions.h"

#include "pgf_fe2_macro_client.h"
#include "pgf_fe2_micro_server.h"

#include "solver_file.h"

static const int ndim = 3;
static const long ARC = 1;

int multi_scale_main(int argc, char **argv)
{
  int err = 0;
  /* intitialize MPI */
  err += MPI_Init(&argc,&argv);
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

  /* get macro and micro parts of the command line */
  {
    int rank_world = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_world);
    get_macro_micro_option_blocks(rank_world,argc,argv,
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
  }

  while(debug);

  PGFEM_mpi_comm *mpi_comm = PGFEM_calloc(1,sizeof(PGFEM_mpi_comm));
  err += initialize_PGFEM_mpi_comm(MPI_COMM_WORLD,mpi_comm);

  if(mpi_comm->rank_world == 0){
    PGFEM_printf("=== COUPLED MULTISCALE ANALYSIS ===\n\n");
  }

  err += PGFEM_mpi_comm_MM_split(nproc_macro,micro_group_size,mpi_comm);

  /*=== INITIALIZE SCALES ===*/
  MACROSCALE *macro = NULL;
  MICROSCALE *micro = NULL;
  if(mpi_comm->valid_macro){/*=== MACROSCALE ===*/
    initialize_MACROSCALE(&macro);
    build_MACROSCALE(macro,mpi_comm->macro,macro_argc,macro_argv);
    build_MACROSCALE_solution(macro);
  } else if(mpi_comm->valid_micro){/*=== MICROSCALE ===*/
    PGFEM_redirect_io_micro();
    initialize_MICROSCALE(&micro);

    /*=== REDIRECT MICROSCALE I/O ===*/
    {
      /* no output */
      PGFEM_redirect_io_null();
      parse_command_line(micro_argc,micro_argv,
			 mpi_comm->rank_micro,micro->opts);
      PGFEM_redirect_io_micro();

      /* create the directory for log output and set logging
	 filename */
      int nproc_world = 0;
      MPI_Comm_size(mpi_comm->world,&nproc_world);
      int group_id = mpi_comm->rank_micro_all/micro_group_size;
      int dir_len = snprintf(NULL,0,"%s/log",micro->opts->opath)+1;
      char *dir_name = PGFEM_calloc(dir_len,sizeof(char));
      sprintf(dir_name,"%s/log",micro->opts->opath);
      make_path(dir_name,DIR_MODE);
      dir_len = snprintf(NULL,0,"%s/group_%05d",dir_name,group_id)+1;
      char *fname = PGFEM_calloc(dir_len,sizeof(char));
      sprintf(fname,"%s/group_%05d",dir_name,group_id);

      /* reinitialize I/O */
      PGFEM_finalize_io();
      PGFEM_initialize_io(NULL,fname);
      PGFEM_redirect_io_micro();
      free(dir_name);
      free(fname);
    }

    /*=== BUILD MICROSCALE ===*/
    build_MICROSCALE(micro,mpi_comm->micro,micro_argc,micro_argv);
  } else {
    PGFEM_printerr("[%d]ERROR: neither macro or microscale!\n%s:%s:%d",
		   mpi_comm->rank_world,__func__,__FILE__,__LINE__);
    PGFEM_Abort();
  }

  /*=== build the macroscacle clients ===*/
  pgf_FE2_macro_client *client = NULL;
  /* hard-code n_jobs_max, needs to be command line opt and/or scaled
     variable */
  const int n_jobs_max = 20;

  /*=== Build MICROSCALE server and solutions ===*/
  if(mpi_comm->valid_micro){
    /* allocate space for maximum number of jobs to be computed. */
    build_MICROSCALE_solutions(n_jobs_max,micro);

    /* start the microscale servers. This function does not exit until
       a signal is passed from the macroscale via
       pgf_FE2_macro_client_send_exit */
    err += pgf_FE2_micro_server_START(mpi_comm,micro);

    /* destroy the microscale */
    destroy_MICROSCALE(micro);

  } else { /*=== MACROSCALE ===*/
    /* initialize the client */
    pgf_FE2_macro_client_init(&client);

    /* create the list of jobs */
    pgf_FE2_macro_client_create_job_list(client,n_jobs_max,macro,mpi_comm);

    /* determine the initial job assignment*/
    pgf_FE2_macro_client_assign_initial_servers(client,mpi_comm);

    /* send the first set of jobs */
    pgf_FE2_macro_client_send_jobs(client,mpi_comm,macro,
				   JOB_NO_COMPUTE_EQUILIBRIUM);

    COMMON_MACROSCALE *c = macro->common;
    MACROSCALE_SOLUTION *s = macro->sol;
    int nproc_macro = 0;
    char filename[500];
    MPI_Comm_size(mpi_comm->macro,&nproc_macro);

    /* Create a context for passing stuff to Newton Raphson */
    MS_SERVER_CTX *ctx = PGFEM_calloc(1,sizeof(MS_SERVER_CTX));
    ctx->client = client;
    ctx->mpi_comm = mpi_comm;
    ctx->macro = macro;

    /*=== COMPUTE APPLIED FORCES ON MARKED SURFACES ===*/
    double *nodal_forces = PGFEM_calloc(c->ndofd,sizeof(double));
    SUR_TRAC_ELEM *ste = NULL;
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
				   nodal_forces);

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
    if(macro->opts->override_solver_file){
      if(mpi_comm->rank_macro == 0){
	PGFEM_printf("Overriding the default solver file with:\n%s\n",
		     macro->opts->solver_file);
      }
      solver_file_open(macro->opts->solver_file,&solver_file);
    } else {
      /* use the default file/filename */
      char *filename = NULL;
      alloc_sprintf (&filename,"%s/%s%d.in.st",macro->opts->ipath,
		     macro->opts->ifname,mpi_comm->rank_macro);
      solver_file_open(filename,&solver_file);
      free(filename);
    }
    s->tim = 0;
    solver_file_read_header(solver_file);

    /* Nonlinear solver */
    if (mpi_comm->rank_macro == 0) {
      switch(solver_file->nonlin_method){
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

    double hypre_time = 0.0;


    /*  NODE (PRESCRIBED DEFLECTION)- SUPPORT COORDINATES generation
	of the load vector  */
    s->dt = solver_file->times[s->tim + 1] - solver_file->times[s->tim];
    if (s->dt == 0.0){
      if (mpi_comm->rank_macro == 0){
	PGFEM_printf("Incorrect dt\n");
      }
      PGFEM_Abort();
    }

    load_vec_node_defl (s->f_defl,c->ne,c->ndofn,c->elem,
			NULL,c->node,c->hommat,
			c->matgeom,c->supports,c->npres,
			solver_file->nonlin_tol,s->sig_e,s->eps,s->dt,
			s->crpl,macro->opts->stab,
			s->r,macro->opts);
    
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

    /* Prescribed deflection */
    double *sup_defl = NULL;
    if(c->supports->npd > 0){
      sup_defl = PGFEM_calloc(c->supports->npd,sizeof(double));
    }

    double pores = 0.0;
    double gama = 0.0;
    double GNOR = 0.0;
    double nor1 = 0.0;
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
    while (solver_file->n_step > s->tim){
      s->dt = dt0 = solver_file->times[s->tim+1] - solver_file->times[s->tim];
      if (s->dt <= 0.0){
	if (mpi_comm->rank_macro == 0) {
	  PGFEM_printf("Incorrect dt\n");
	}
	PGFEM_Abort();
      }

      if (mpi_comm->rank_macro == 0){
	PGFEM_printf("\nFinite deformations time step %ld) "
		     " Time %f | dt = %10.10f\n",
		     s->tim,solver_file->times[s->tim+1],s->dt);
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
	hypre_time += Newton_Raphson ( 1,&n_step,c->ne,0,c->nn,
				       c->ndofn,c->ndofd,c->npres,s->tim,
				       solver_file->times,
				       solver_file->nonlin_tol,s->dt,c->elem,
				       NULL,c->node,c->supports,sup_defl,
				       c->hommat,c->matgeom,s->sig_e,s->eps,
				       c->Ap,c->Ai,s->r,s->f,
				       s->d_r,s->rr,s->R,s->f_defl,
				       s->RR,s->f_u,s->RRn,s->crpl,
				       macro->opts->stab,c->nce,c->coel,
				       solver_file->nonlin_method,
				       &pores,c->SOLVER,s->BS_x,s->BS_f,
				       s->BS_RR,gama,GNOR,nor1,
				       c->lin_err,s->BS_f_u,c->DomDof,
				       c->pgfem_comm,c->GDof,
				       solver_file->n_step,
				       solver_file->max_nonlin_iter,
				       &(s->NORM),c->nbndel,c->bndel,
				       c->mpi_comm,c->VVolume,macro->opts,ctx,0,
				       NULL,NULL);

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
	char out_dat[500];
	double tmp_val = ((s->times[s->tim+1]-s->times[s->tim])
			  /dt0*solver_file->nonlin_method_opts[1]);
	dlm = Arc_length ( c->ne,0,c->nn,c->ndofn,
			   c->ndofd,c->npres,
			   solver_file->n_step,s->tim,
			   solver_file->times,solver_file->nonlin_tol,
			   solver_file->max_nonlin_iter,s->dt,
			   dt0,c->elem,NULL,c->nbndel,
			   c->bndel,c->node,c->supports,sup_defl,
			   c->hommat,c->matgeom,s->sig_e,s->eps,
			   c->Ap,c->Ai,c->SOLVER,
			   s->RRn,s->f_defl,s->crpl,macro->opts->stab,
			   c->nce,c->coel,s->r,s->f,
			   s->d_r,s->D_R,s->rr,s->R,
			   s->RR,s->f_u,s->U,s->DK,
			   s->dR,s->BS_f,s->BS_d_r,s->BS_D_R,
			   s->BS_rr,s->BS_R,s->BS_RR,s->BS_f_u,
			   s->BS_U,s->BS_DK,s->BS_dR,solver_file->nonlin_method,
			   lm,solver_file->nonlin_method_opts[0],&DET,&dlm0,
			   &DLM,macro->opts->vis_format,
			   macro->opts->smoothing,
			   s->sig_n,out_dat,
			   (long*) (solver_file->print_steps),&AT,
			   ARC,tmp_val,&ITT,&dAL,
			   &pores,c->DomDof,c->GDof,c->pgfem_comm,
			   c->lim_zero,&s->NORM,c->mpi_comm,
			   c->VVolume,macro->opts);

	/* Load multiplier */
	lm += dlm;
	dlm = 0.0;
	
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
	  sur_forces = PGFEM_calloc(n_feats*ndim,sizeof(double));
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

      if (solver_file->print_steps[s->tim] == 1
	  && macro->opts->vis_format != VIS_NONE ) {
	/* NOTE: null d_r is sent for print job because jump is
	   computed from updated coordinates! */

	/* do not transfer data between servers */
	pgf_FE2_macro_client_rebalance_servers(client,mpi_comm,
					       FE2_REBALANCE_NONE);

	/* Send print jobs */
	pgf_FE2_macro_client_send_jobs(client,mpi_comm,macro,JOB_PRINT);

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
			c->supports,s->r,s->sig_e,s->eps,
			macro->opts);

	  if (macro->opts->cohesive == 1){
	    if(mpi_comm->rank_macro == 0){
	      VTK_print_cohesive_master(macro->opts->opath,
					macro->opts->ofname,
					s->tim,nproc_macro,macro->opts);
	    }

	    VTK_print_cohesive_vtu(macro->opts->opath,macro->opts->ofname,
				   s->tim,mpi_comm->rank_macro,c->nce,c->node,
				   c->coel,c->supports,s->r,c->ensight,
				   macro->opts);
	  }
	  break;
	default: /* no output */ break;
	}/* switch(format) */

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

  /* finalize and exit */
  err += PGFEM_finalize_io();
  err += MPI_Finalize();
  return err;
}
