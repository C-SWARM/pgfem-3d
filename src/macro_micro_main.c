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
#include "BSprivate.h"


/*=== PFEM3d headers ===*/
#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef ARC_LENGTH_H
#include "Arc_length.h"
#endif

#ifndef BUILD_DISTRIBUTION_H
#include "build_distribution.h"
#endif

#ifndef HOMOGEN_H
#include "homogen.h"
#endif

#ifndef HYPRE_GLOBAL_H
#include "hypre_global.h"
#endif

#ifndef IN_H
#include "in.h"
#endif

#ifndef LOAD_H
#include "load.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

#ifndef MATRIX_PRINTING_H
#include "matrix_printing.h"
#endif

#ifndef NEWTON_RAPHSON_H
#include "Newton_Raphson.h"
#endif

#ifndef OUT_H
#include "out.h"
#endif

#ifndef PRINTING_H
#include "Printing.h"
#endif

#ifndef PRINT_DIST_H
#include "print_dist.h"
#endif

#ifndef PROFILER_H
#include "profiler.h"
#endif

#ifndef PSPARSE_APAI_H
#include "Psparse_ApAi.h"
#endif

#ifndef READ_CRYST_PLAST_H
#include "read_cryst_plast.h"
#endif

#ifndef RENUMBER_ID_H
#include "renumber_ID.h"
#endif

#ifndef RNPSPARSE_APAI_H
#include "RNPsparse_ApAi.h"
#endif

#ifndef SET_FINI_DEF_H
#include "set_fini_def.h"
#endif

#ifndef SKYLINE_H
#include "skyline.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef SETGLOBALNODENUMBERS_H
#include "SetGlobalNodeNumbers.h"
#endif

#ifndef INTERFACE_MACRO_H
#include "interface_macro.h"
#endif

#ifndef COMPUTE_MACRO_F_H
#include "computeMacroF.h"
#endif

#ifndef COMPUTE_MACRO_S_H
#include "computeMacroS.h"
#endif

#ifndef VTK_OUTPUT_H
#include "vtk_output.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifndef GEN_PATH_H
#include "gen_path.h"
#endif

#ifndef MICROSCALE_INFORMATION_H
#include "microscale_information.h"
#endif

#ifndef MS_COHE_JOB_LIST_H
#include "ms_cohe_job_list.h"
#endif

#ifndef APPLIED_TRACTION_H
#include "applied_traction.h"
#endif

#ifndef MACRO_MICRO_FUNCTIONS_H
#include "macro_micro_functions.h"
#endif

/* #define str(a) xstr(a) */
/* #define xstr(a) #a */
/* #define mgs 4 */

static const int ndim = 0;

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
      dir_len = snprintf(NULL,0,"%s/group_%0.5d",dir_name,group_id)+1;
      char *fname = PGFEM_calloc(dir_len,sizeof(char));
      sprintf(fname,"%s/group_%0.5d",dir_name,group_id);

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

  /*=== BUILD INTERCOMMUNICATOR NETWORK ===*/
  PGFEM_ms_job_intercomm *intercomm = NULL;
  int n_jobs = 0;
  if(mpi_comm->valid_mm_inter){ /* MACRO + MICRO-MASTERS */
    int *job_buff_sizes = NULL;
    if(mpi_comm->valid_macro){ /*=== MACROSCALE ===*/
      /* compute number of jobs and job_buff_sizes */
      COMMON_MACROSCALE *c = macro->common;
      err += compute_n_job_and_job_sizes(c,&n_jobs,&job_buff_sizes);
    }
    /* create the intercommunicator on all processes in the
       communicator. Microscale only talks to macroscale and
       vice-versa */
    err += create_PGFEM_ms_job_intercomm(nproc_macro,
					 mpi_comm,n_jobs,
					 job_buff_sizes,
					 &intercomm);
    free(job_buff_sizes);
  }


  if(mpi_comm->valid_micro){/*=== MICROSCALE ===*/
    /* allocate the microscale solutions based on the number of jobs I
       will get from macroscale */
    if(mpi_comm->valid_mm_inter){
      err += PGFEM_comm_info_get_n_comms(intercomm->recv_info,&n_jobs);
    }
    err += MPI_Bcast(&n_jobs,1,MPI_INT,0,mpi_comm->micro);
    build_MICROSCALE_solutions(n_jobs,micro);

    /* The microscale communicator starts a server and computes jobs
       as they are received until it recieves a job tagged with
       JOB_EXIT. */
    err += start_microscale_server(mpi_comm,intercomm,micro);

    /* proceed to cleanup and exit */
  } else { /*=== MACROSCALE ===*/
    /*=== BUILD LIST OF JOBS ===*/
    /* NOTE: only the MACROSCALE stores the list of jobs */
    MS_COHE_JOB_INFO *job_list = NULL;
    COMMON_MACROSCALE *c = macro->common;
    MACROSCALE_SOLUTION *s = macro->sol;
    int nproc_macro = 0;
    char filename[500];
    MPI_Comm_size(mpi_comm->macro,&nproc_macro);
    {
      long Gn_jobs = 0;
      long *n_job_dom = NULL;
    /* use create group but pass MPI_COMM_SELF for ms_comm and
       mpi_comm->macro for macro_mpi_comm */
      err += create_group_ms_cohe_job_list(c->nce,c->coel,c->node,
					   mpi_comm->macro,MPI_COMM_SELF,
					   0,&Gn_jobs,&n_job_dom,&job_list);
      /* error check */
      if(Gn_jobs != n_jobs){
	PGFEM_printerr("ERROR: got different number "
		       "of jobs on macroscale!\n");
	PGFEM_Abort();
      }
      free(n_job_dom);
    }

    PGFEM_server_ctx *send = PGFEM_calloc(1,sizeof(PGFEM_server_ctx));
    PGFEM_server_ctx *recv = PGFEM_calloc(1,sizeof(PGFEM_server_ctx));
    err += initialize_PGFEM_server_ctx(send);
    err += initialize_PGFEM_server_ctx(recv);
    err += build_PGFEM_server_ctx_from_PGFEM_comm_info(intercomm->send_info,send);
    err += build_PGFEM_server_ctx_from_PGFEM_comm_info(intercomm->recv_info,recv);


    /* compute the inital tangent from the microscale */
    err += start_macroscale_compute_jobs(intercomm,macro,
					 JOB_NO_COMPUTE_EQUILIBRIUM,
					 s->r,job_list,send,recv);

    /* err += finish_macroscale_compute_jobs(job_list,macro,send,recv); */

    /* Create a context for passing stuff to Newton Raphson */
    MS_SERVER_CTX *ctx = PGFEM_calloc(1,sizeof(MS_SERVER_CTX));
    ctx->n_jobs = n_jobs;
    ctx->job_list = job_list;
    ctx->send = send;
    ctx->recv = recv;
    ctx->intercomm = intercomm;
    ctx->macro = macro;


    /*=== COMPUTE APPLIED FORCES ON MARKED SURFACES ===*/
    double *nodal_forces = PGFEM_calloc(c->ndofd,sizeof(double));
    { 
      int n_feats = 0;
      int n_sur_trac_elem = 0;
      int *feat_type = NULL;
      int *feat_id = NULL;
      double *loads = NULL;
      SUR_TRAC_ELEM *ste = NULL;

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
      destroy_applied_surface_traction_list(n_sur_trac_elem,ste);
    }

    /* push nodal_forces to s->R */
    vvplus  (s->R,nodal_forces,c->ndofd);

    /*=== SOLUTION PROCESS ===*/
    /*=== READ SOLVER FILE ===*/
    /* override the default solver file with one specified
       at commandline */
    FILE *in1 = NULL;
    if(macro->opts->override_solver_file){
      if(mpi_comm->rank_macro == 0){
	PGFEM_printf("Overriding the default solver file with:\n%s\n",
		     macro->opts->solver_file);
      }
      in1 = PGFEM_fopen(macro->opts->solver_file,"r");
    } else {
      /* use the default file/filename */
      char filename[500];
      sprintf (filename,"%s/%s%d.in.st",macro->opts->ipath,
	       macro->opts->ifname,mpi_comm->rank_macro);
      in1 = PGFEM_fopen(filename,"r");
    }

    double nor_min = 0.0;
    long iter_max = 0;
    long FNR = 0;
    long ARC = 1;
    double dAL0 = 0.0;
    double dALMAX = 0.0;
    double hypre_time = 0.0;
    {
      long npres = 0;
      fscanf (in1,"%lf %ld %ld %ld",&nor_min,&iter_max,&npres,&FNR);
      if (FNR == 2 || FNR == 3){
	fscanf (in1,"%lf %lf",&dAL0,&dALMAX);
      }
    }
    
    /* Nonlinear solver */
    if (mpi_comm->rank_macro == 0) {
      if (FNR == 0 || FNR == 1) {
	PGFEM_printf ("\nNONLINEAR SOLVER : NEWTON-RAPHSON METHOD\n");
      }

      if ((FNR == 2 || FNR == 3) && ARC == 0){
	PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Crisfield\n");
      }
      if ((FNR == 2 || FNR == 3) && ARC == 1) {
	PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Simo\n");
      }
    }

    /* read number of computational times */
    long nt = 0;
    fscanf (in1,"%ld",&nt);
    
    /* Compute times */
    free(s->times);
    s->times = aloc1 (nt+1);
    for (int i=0;i<nt+1;i++){
      fscanf (in1,"%lf",&(s->times[i]));
    }
    
    /* read times for output */
    long n_p = 0;
    fscanf (in1,"%ld",&n_p);
    
    /* Times for printing */
    long *print = times_print (in1,nt,n_p);
    
    long nlod_tim = 0;
    fscanf (in1,"%ld",&nlod_tim);
    
    /* read times dependent load */
    long *tim_load = compute_times_load (in1,nt,nlod_tim);
    

    /*  NODE (PRESCRIBED DEFLECTION)- SUPPORT COORDINATES generation
	of the load vector  */
    s->dt = s->times[1] - s->times[0];
    if (s->dt == 0.0){
      if (mpi_comm->rank_macro == 0){
	PGFEM_printf("Incorrect dt\n");
      }
      PGFEM_Abort();
    }

    load_vec_node_defl (s->f_defl,c->ne,c->ndofn,c->elem,
			NULL,c->node,c->hommat,
			c->matgeom,c->supports,c->npres,
			nor_min,s->sig_e,s->eps,s->dt,
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
      for (int i=0;i<c->supports->npd;i++){
	sup_defl[i] = c->supports->defl_d[i];
      }
    }

    /*=== MODEL ENTITIES ===*/
    MODEL_ENTITY *entities = NULL;
    double *forces = NULL;
    if(macro->opts->me){
      char *me_fname = NULL;
      alloc_sprintf(&me_fname,"%s/entities.in",macro->opts->ipath);
      read_model_entity_list(me_fname,&entities,c->nn,c->ne);
      model_entity_mark_nodes_elems(c->nn,c->node,c->ne,c->elem,entities);
      if(entities->n_entities > 0){
	forces = PGFEM_calloc(entities->n_entities
			      *entities->n_dim,sizeof(double));
      }
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

    /*=== BEGIN SOLVE ===*/
    lm = dlm = DLM = DET = 0.0;
    dAL = dlm0 = dAL0;
    s->tim = AT = ITT = 0;

    while (nt > s->tim){
      s->dt = dt0 = s->times[s->tim+1] - s->times[s->tim];
      if (s->dt <= 0.0){
	if (mpi_comm->rank_macro == 0) {
	  PGFEM_printf("Incorrect dt\n");
	}
	PGFEM_Abort();
      }

      if (mpi_comm->rank_macro == 0){
	PGFEM_printf("\nFinite deformations time step %ld) "
		     " Time %f | dt = %10.10f\n",
		     s->tim,s->times[s->tim+1],s->dt);
      }

      /*=== NEWTON RAPHSON ===*/
      if (FNR == 0 || FNR == 1){
	if (tim_load[s->tim] == 1 && s->tim == 0) {
	  if (mpi_comm->rank_macro == 0){
	    PGFEM_printf ("Incorrect load input for Time = 0\n");
	  }
	  PGFEM_Abort();
	}
	if (tim_load[s->tim] == 1 && s->tim != 0) {
	  /*  read nodal prescribed deflection */
	  for (int i=0;i<c->supports->npd;i++){
	    fscanf (in1,"%lf",&c->supports->defl_d[i]);
	    sup_defl[i] = c->supports->defl_d[i];
	  }

	  /*=== do not support node/surf loads ===*/
	  /* /\* read nodal load in the subdomain *\/ */
	  /* read_nodal_load (in1,nln,ndim,znod); */
	  /* /\* read elem surface load *\/ */
	  /* read_elem_surface_load (in1,nle_s,ndim,elem,zele_s); */
	  /* /\*  NODE - generation of the load vector  *\/ */
	  /* load_vec_node (R,nln,ndim,znod,node); */
	  /* /\*  ELEMENT - generation of the load vector  *\/ */
	  /* load_vec_elem_sur (R,nle_s,ndim,elem,zele_s); */

	  /*
	    R   -> Incramental forces 
	    RR  -> Total forces for sudivided increment 
	    RRn -> Total force after equiblirium
	  */

	  /* push nodal_forces to s->R */
	  vvplus  (s->R,nodal_forces,c->ndofd);

	} /* end load increment */



	hypre_time += Newton_Raphson ( 1,c->ne,0,c->nn,
				       c->ndofn,c->ndofd,c->npres,s->tim,
				       s->times,nor_min,s->dt,c->elem,
				       NULL,c->node,c->supports,sup_defl,
				       c->hommat,c->matgeom,s->sig_e,s->eps,
				       c->Ap,c->Ai,s->r,s->f,
				       s->d_r,s->rr,s->R,s->f_defl,
				       s->RR,s->f_u,s->RRn,s->crpl,
				       macro->opts->stab,c->nce,c->coel,FNR,
				       &pores,c->SOLVER,NULL,NULL,
				       NULL,NULL,s->BS_x,s->BS_f,
				       s->BS_RR,gama,GNOR,nor1,
				       c->lin_err,s->BS_f_u,c->DomDof,
				       c->pgfem_comm,c->GDof,nt,iter_max,
				       &(s->NORM),c->nbndel,c->bndel,
				       c->mpi_comm,c->VVolume,macro->opts,
				       entities,forces,NULL,ctx);

	/* Null global vectors */
	for (int i=0;i<c->ndofd;i++){
	  s->RRn[i] += s->R[i];
	  s->RR[i] = s->RRn[i];
	  s->R[i] = 0.0;
	}

      }/* end NR */

      /*=== ARC LENGTH ===*/
      if (FNR == 2 || FNR == 3){
	char out_dat[500];
	double tmp_val = (s->times[s->tim+1]-s->times[s->tim])/dt0*dALMAX;
	dlm = Arc_length ( c->ne,0,c->nn,c->ndofn,
			   c->ndofd,c->npres,nt,s->tim,
			   s->times,nor_min,iter_max,s->dt,
			   dt0,c->elem,NULL,c->nbndel,
			   c->bndel,c->node,c->supports,sup_defl,
			   c->hommat,c->matgeom,s->sig_e,s->eps,
			   c->Ap,c->Ai,NULL,c->SOLVER,
			   s->RRn,s->f_defl,s->crpl,macro->opts->stab,
			   c->nce,c->coel,s->r,s->f,
			   s->d_r,s->D_R,s->rr,s->R,
			   s->RR,s->f_u,s->U,s->DK,
			   s->dR,s->BS_f,s->BS_d_r,s->BS_D_R,
			   s->BS_rr,s->BS_R,s->BS_RR,s->BS_f_u,
			   s->BS_U,s->BS_DK,s->BS_dR,FNR,
			   lm,dAL0,&DET,&dlm0,
			   &DLM,macro->opts->vis_format,
			   macro->opts->smoothing,
			   s->sig_n,out_dat,print,&AT,
			   ARC,tmp_val,&ITT,&dAL,
			   &pores,c->DomDof,c->GDof,c->pgfem_comm,
			   NULL,NULL,NULL,c->lim_zero,
			   &s->NORM,c->mpi_comm,c->VVolume,macro->opts);

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

      /*=== fully coupled, so why do this? ===*/
      /* /\* Calculate macro deformation gradient *\/ */
      /* double *GF = computeMacroF(c->elem,c->ne,c->node,c->nn, */
      /* 				 s->eps,c->VVolume,c->mpi_comm); */
      /* double *GS = computeMacroS(c->elem,c->ne,c->node,c->nn, */
      /* 				 s->sig_e,c->VVolume,c->mpi_comm); */
      /* double *GP = computeMacroP(c->elem,c->ne,c->node,c->nn, */
      /* 				 s->sig_e,s->eps,c->VVolume,c->mpi_comm); */

      /* /\* print GF & GS to file *\/ */
      /* if(mpi_comm->rank_macro == 0){ */
      /* 	FILE *out = NULL; */
      /* 	sprintf(filename,"%s/%s_macro.out.%d",macro->opts->opath, */
      /* 		macro->opts->ofname,s->tim); */
      /* 	out = PGFEM_fopen(filename,"w"); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[0],GF[1],GF[2]); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[3],GF[4],GF[5]); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[6],GF[7],GF[8]); */
      /* 	PGFEM_fprintf(out,"\n"); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\t",GS[0],GS[1],GS[2]); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GS[3],GS[4],GS[5]); */
      /* 	PGFEM_fprintf(out,"\n"); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[0],GP[1],GP[2]); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[3],GP[4],GP[5]); */
      /* 	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[6],GP[7],GP[8]); */
      /* 	fclose(out); */
      /* } */

      /* free(GF); */
      /* free(GS); */
      /* free(GP); */

      if (print[s->tim] == 1 && macro->opts->vis_format != VIS_NONE ) {
	/* NOTE: null d_r is sent for print job because jump is
	   computed from updated coordinates! */
	err += start_macroscale_compute_jobs(intercomm,macro,
					     JOB_PRINT,
					     s->d_r,
					     job_list,send,recv);
	if(macro->opts->ascii){
	  int Gnn = 0;
	  ASCII_output(macro->opts,c->mpi_comm,s->tim,s->times,
		       Gnn,c->nn,c->ne,c->nce,c->ndofd,
		       c->DomDof,c->Ap,FNR,lm,pores,c->VVolume,
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
	  EnSight (filename,s->tim,nt,c->nn,c->ne,ndim,c->node,
		   c->elem,c->supports,s->r,s->sig_e,s->sig_n,s->eps,
		   macro->opts->smoothing,c->nce,c->coel,
		   /*nge,geel,ngn,gnod,*/FNR,lm,c->ensight,c->mpi_comm,
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
					s->tim,nproc_macro);
	    }

	    VTK_print_cohesive_vtu(macro->opts->opath,macro->opts->ofname,
				   s->tim,mpi_comm->rank_macro,c->nce,c->node,
				   c->coel,c->supports,s->r,c->ensight,
				   macro->opts);
	  }
	  break;
	default: /* no output */ break;
	}/* switch(format) */

	err += finish_macroscale_compute_jobs(job_list,macro,send,recv);
      }/* end output */
 
      if (mpi_comm->rank_macro == 0){
	PGFEM_printf("\n");
	PGFEM_printf("*********************************************\n");
	PGFEM_printf("*********************************************\n");
      }
      
      s->tim++;
    }/* end while */

    /*=== SEND EXIT SIGNAL TO MICROSCALE SERVER ===*/
    err += start_macroscale_compute_jobs(intercomm,macro,
					 JOB_EXIT,s->r,
					 job_list,send,recv);

    err += finish_macroscale_compute_jobs(job_list,macro,send,recv);

    /* cleanup */
    free(ctx);
    err += destroy_PGFEM_server_ctx(send);
    err += destroy_PGFEM_server_ctx(recv);
    free(send);
    free(recv);
    for(int i=0; i<n_jobs; i++){
      destroy_MS_COHE_JOB_INFO(&job_list[i]);
    }
    free(job_list);
    free(sup_defl);
    free(tim_load);
    free(print);
    free(forces);
    destroy_model_entity(entities);
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

  /* destroy the scale information */
  destroy_MACROSCALE(macro);
  destroy_MICROSCALE(micro);

  /* destroy the intercommunicator */
  destroy_PGFEM_ms_job_intercomm(intercomm);
  free(intercomm);

  /* destroy the communicator */
  err += destroy_PGFEM_mpi_comm(mpi_comm);
  free(mpi_comm);

  /* finalize and exit */
  err += PGFEM_finalize_io();
  err += MPI_Finalize();
  return err;
}
