/* HEADER */

/** This file contains functions to parse the command line and print
    helpful information */

#include "PGFem3D_options.h"
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include "enumerations.h"
#include "allocation.h"

/* Generated at http://patorjk.com/software/taag/ */
static const char *prog_name=
  " _______    ______   ________                        ______   _______  \n"
  "/       \\  /      \\ /        |                      /      \\ /       \\ \n"
  "$$$$$$$  |/$$$$$$  |$$$$$$$$/______   _____  ____  /$$$$$$  |$$$$$$$  |\n"
  "$$ |__$$ |$$ | _$$/ $$ |__  /      \\ /     \\/    \\ $$ ___$$ |$$ |  $$ |\n"
  "$$    $$/ $$ |/    |$$    |/$$$$$$  |$$$$$$ $$$$  |  /   $$< $$ |  $$ |\n"
  "$$$$$$$/  $$ |$$$$ |$$$$$/ $$    $$ |$$ | $$ | $$ | _$$$$$  |$$ |  $$ |\n"
  "$$ |      $$ \\__$$ |$$ |   $$$$$$$$/ $$ | $$ | $$ |/  \\__$$ |$$ |__$$ |\n"
  "$$ |      $$    $$/ $$ |   $$       |$$ | $$ | $$ |$$    $$/ $$    $$/ \n"
  "$$/        $$$$$$/  $$/     $$$$$$$/ $$/  $$/  $$/  $$$$$$/  $$$$$$$/  \n";

static char cur_dir = '.';

typedef struct long_opt_descr{
  struct option opt;
  char *descr;
  int sc;
} long_opt_descr;

static const char *opt_list = "AVlh";

static const long_opt_descr analysis_opts[] = {
  /* Analysis options */
  {{"cp",no_argument,NULL,'c'},"\tFinite strain crystal plasticity",0},
  {{"fd",no_argument,NULL,'f'},"\tFinite strain elasticity (Quadradic tetras or hexas)",0},
  {{"cpp",no_argument,NULL,0},"\tFinite strain Portevin-Le Chatelier effect",0},
  {{"st",required_argument,NULL,1},"Stabilized finite strains analysis",0},
  {{"he",no_argument,NULL,2},"\tFinite strains analysis using MINI element",0},
  {{"he3",no_argument,NULL,2},"\tFinite strains analysis using bubble-enhanced 3-Field element",0},
  {{"disp",no_argument,NULL,2},"\tTOTAL Lagrangian displacement-based finite strains analysis",0},
  {{"tf",no_argument,NULL,2},"\tTOTAL Lagrangian displacement-based 3 field finite strains analysis",0},
  {{"cm",required_argument,NULL,2},"Use of constitutive model interface\n"
                                   "\t\t arg = 0: Hyperelasticity\n"
                                   "\t\t arg = 1: Crystal plasticity\n"
                                   "\t\t arg = 2: BPA_plasticity",0},  
  {{"coh",no_argument,NULL,1},"\tCohesive elements",0},
  {{"ms",no_argument,NULL,'m'},("\tINTERFACE or BULK multiscale modeling.\n"
				"\t\tRequires six (6) or nine (9), respectively, prescribed displacements\n"
				"\t\tand a file named \"normal.in\" in the specified\n"
				"\t\tinput directory containing the macroscopic normal,\n"
				"\t\tbounding volume, and cell thickness. Format: [V lc Nx Ny Nz]."),0}
};

static const long_opt_descr solver_opts[] = {
  /* Solver options */
  {{"gmres",no_argument,NULL,3},"\tUse the HYPRE GMRES solver (default)",0},
  {{"boomer",no_argument,NULL,3},"\tUse the HYPRE BoomerAMG solver (no precond)",0},
  {{"bcgstab",no_argument,NULL,3},"Use the HYPRE BiCGSTAB solver",0},
  {{"flex-gmres",no_argument,NULL,3},"Use the HYPRE FlexGMRES solver",0},
  {{"hybrid",no_argument,NULL,3},"\tUse the HYPRE Hybrid (GMRES) solver\n"
   "\t\t(adaptively switches preconditioner)",0},
  {{"kdim",required_argument,NULL,3},"Set the Krylov dimension",0},
  {{"maxit",required_argument,NULL,3},"Set the maximum number of iterations",0}
};

static const long_opt_descr precond_opts[] = {
  /* Preconditioner options */
  {{"pre-euclid",no_argument,NULL,4},"Use the HYPRE Euclid preconditioner (default)",0},
  {{"pre-boomer",no_argument,NULL,4},"Use the HYPRE BoomerAMG preconditioner",0},
  {{"pre-pilut",no_argument,NULL,4},"Use the HYPRE PILUT preconditioner",0},
  {{"pre-sails",no_argument,NULL,4},"Use the HYPRE ParaSAILS preconditioner",0},
  {{"pre-diag",no_argument,NULL,4},"Use the custom diagonal scaling preconditioner",0},
  {{"pre-jacobi",no_argument,NULL,4},"Use the Jacobi scaling preconditioner",0},
  {{"pre-none",no_argument,NULL,4},"Do not use a preconditioner",0}
};

static const long_opt_descr vis_opts[] = {
  /* Visualization options */
  {{"vtk",no_argument,NULL,'V'},"Output in VTK format",1},
  {{"ascii",no_argument,NULL,'A'},"Additional output in ASCII format",1},
};

static const long_opt_descr other_opts[] = {
  /* Other options */
  {{"ipath",required_argument,NULL,'i'},"Path to input files parent directory",0},
  {{"opath",required_argument,NULL,'o'},"Path to output files parent directory",0},
  {{"override-pre-disp",required_argument,NULL,'O'},("\n\t\tOverride the prescribed displacements in *.in\n"
						   "\t\twith those provided in the given file."),0},
  {{"override-solver-file",required_argument,NULL,'O'},("\n\t\tOverride the default solver filename with\n"
						      "\t\tthe provided filename."),0},
  {{"override-material-props",required_argument,NULL,'O'},("\n\t\tOverride the material properties in *.in\n"
                                                           "\t\twith those provided in the given file."),0},
  {{"restart",required_argument,NULL,'r'},("Restart from specified step (FE2 only). Requires original\n"
					   "\t\tinput files and dumped restart files for specified step."),0},
  {{"max-server-jobs",required_argument,NULL,'S'},("\n\t\tSet the maximum number of jobs allowed on a server (FE2)."),0},
  {{"no-migrate",no_argument,NULL,'N'},("Do not migrate cells between servers (FE2)."),0},
  {{"legacy",no_argument,NULL,'l'},"Read files from legacy format",1},
  {{"debug",no_argument,NULL,9999},"\tSend into infinite loop to attach debugger",0},
  {{"help",no_argument,NULL,'h'},"Print this message and exit",1}

};

/* these options may no longer be supported/functional. They are kept
   for documentation purposes only and are ignored if used */
static const long_opt_descr depricated_opts[] = {
  {{"elixir",no_argument,NULL,'X'},"Output in Elixir format [unsupported, use -V]",1},
  {{"ensight",no_argument,NULL,'E'},"Output in EnSight format [outdated, use -V]",1},
  {{"sm",no_argument,NULL,'s'},"Smooth stress field [unsupported]",1},
  {{"me",no_argument,NULL,'m'},"\tCompute nodal forces on model entities"
   " (requires entities.in file) [unsupported]",0},
  {{"pr",no_argument,NULL,'p'},"Periodic domain [outdated, use -ms]",1},
  {{"rn",no_argument,NULL,'r'},"\tRenumber degrees of freedom [unsupported]",0},
};

static const struct option last_opt = {NULL,0,NULL,0};

static const int n_analysis = sizeof(analysis_opts)/sizeof(long_opt_descr);
static const int n_solver = sizeof(solver_opts)/sizeof(long_opt_descr);
static const int n_precond = sizeof(precond_opts)/sizeof(long_opt_descr);
static const int n_vis = sizeof(vis_opts)/sizeof(long_opt_descr);
static const int n_other = sizeof(other_opts)/sizeof(long_opt_descr);
static const int n_depricated = sizeof(depricated_opts)/sizeof(long_opt_descr);

static void print_long_options(FILE* out,
			       const long_opt_descr *opts,
			       int n_opt)
{
  for(int i=0; i<n_opt; i++){
    if(opts[i].sc == 1 && opts[i].opt.has_arg == no_argument){
      /* single character alias & no arg */
      PGFEM_fprintf(out,"-%s [-%c]\t%s\n",opts[i].opt.name,opts[i].opt.val,
	      opts[i].descr);
    } else if(opts[i].sc == 1 && opts[i].opt.has_arg == required_argument){
      /* single character alias & requires arg */
      PGFEM_fprintf(out,"-%s [-%c] (arg)\t%s\n",opts[i].opt.name,opts[i].opt.val,
	      opts[i].descr);
    } else if(opts[i].sc == 0 && opts[i].opt.has_arg == no_argument){
      /* no alias no arg */
      PGFEM_fprintf(out,"-%s\t%s\n",opts[i].opt.name,opts[i].descr);
    } else if(opts[i].sc == 0 && opts[i].opt.has_arg == required_argument){
      PGFEM_fprintf(out,"-%s (arg)\t%s\n",opts[i].opt.name,opts[i].descr);
    }
  }
}

static void copy_options(struct option *A)
{
  const int n_options = n_analysis + n_solver + n_precond + n_vis + n_other + n_depricated;
  int idx = 0;
  for(int i=0; i<n_analysis; i++){
    A[i+idx] = analysis_opts[i].opt;
  }
  idx += n_analysis;
  for(int i=0; i<n_solver; i++){
    A[i+idx] = solver_opts[i].opt;
  }
  idx += n_solver;
  for(int i=0; i<n_precond; i++){
    A[i+idx] = precond_opts[i].opt;
  }
  idx += n_precond;
  for(int i=0; i<n_vis; i++){
    A[i+idx] = vis_opts[i].opt;
  }
  idx += n_vis;
  for(int i=0; i<n_other; i++){
    A[i+idx] = other_opts[i].opt;
  }
  A[n_options] = last_opt;
}

/*====================================================*/
/*                   END LOCALS                       */
/*====================================================*/

void set_default_options(PGFem3D_opt *options)
{
  /* solver options */
  options->solverpackage = 1; /* HYPRE */
  options->solver = HYPRE_GMRES;
  options->precond = EUCLID;
  options->kdim = 500;
  options->maxit = 1000;

  /* analysis options */
  options->analysis_type = -1;
  options->stab = 0;
  options->cohesive = 0;
  options->plc = 0;
  options->multi_scale = 0;


  /* visualization options */
  options->vis_format = VIS_NONE; /* no output */
  options->ascii = 0; /* no ASCII output */
  options->smoothing = -1; /* no smoothing (undefined value) */

  /* other options */
  options->periodic = -1; /* non-periodic (undefined value)*/
  options->renumber = 0;
  options->legacy = 0;
  options->debug = 0;
  options->me = 0;
  options->restart = -1; /* flag >= 0 used to specify both restart and
			    step to start from */
  options->max_n_jobs = 0;
  options->no_migrate = 0;

  /* input overrides */
  options->override_pre_disp = 0;
  options->pre_disp_file = NULL;
  options->override_solver_file = 0;
  options->solver_file = NULL;
  options->override_material_props = NULL;

  /* I/O file names */
  options->ipath = NULL;
  options->opath = NULL;
  options->ifname = NULL;
  options->ofname = NULL;
}

void print_options(FILE *out,
		   const PGFem3D_opt *options)
{
  PGFEM_fprintf(out,"OPTION VALUES:\n");
  PGFEM_fprintf(out,"=== SOLVER OPTIONS ===\n");
  PGFEM_fprintf(out,"Solver package: %d\n",options->solverpackage);
  PGFEM_fprintf(out,"Solver:         %d\n",options->solver);
  PGFEM_fprintf(out,"Preconditioner: %d\n",options->precond);
  PGFEM_fprintf(out,"Kdim:           %d\n",options->kdim);
  PGFEM_fprintf(out,"Max It:         %d\n",options->maxit);

  PGFEM_fprintf(out,"\n=== ANALYSIS OPTIONS ===\n");
  PGFEM_fprintf(out,"Analysis type:      %d\n",options->analysis_type);
  PGFEM_fprintf(out,"Stab parameter:     %.12e\n",options->stab);
  PGFEM_fprintf(out,"Cohesive elems:     %d\n",options->cohesive);
  PGFEM_fprintf(out,"Multi-scale:        %d\n",options->multi_scale);
  PGFEM_fprintf(out,"Constitutive Model: %d\n",options->cm);

  PGFEM_fprintf(out,"\n=== VISUALIZATION OPTIONS ===\n");
  PGFEM_fprintf(out,"Visualization:  %d\n",options->vis_format);
  PGFEM_fprintf(out,"Smoothing:      %d\n",options->smoothing);

  PGFEM_fprintf(out,"\n=== OTHER OPTIONS ===\n");
  PGFEM_fprintf(out,"Periodic:       %d\n",options->periodic);
  PGFEM_fprintf(out,"Renumber:       %d\n",options->renumber);
  PGFEM_fprintf(out,"Legacy format:  %d\n",options->legacy);
  PGFEM_fprintf(out,"Debug:          %d\n",options->debug);
  PGFEM_fprintf(out,"Restart:        %d\n",options->restart);

  PGFEM_fprintf(out,"\n=== FILE OPTIONS ===\n");
  PGFEM_fprintf(out,"IPath:          %s\n",options->ipath);
  PGFEM_fprintf(out,"OPath:          %s\n",options->opath);
  PGFEM_fprintf(out,"Input: %s\n",options->ifname);
  PGFEM_fprintf(out,"Output: %s\n",options->ofname);

  PGFEM_fprintf(out,"\n");
}

void print_usage(FILE* out)
{
  PGFEM_fprintf(out,"%s\n",prog_name);
  PGFEM_fprintf(out,"SS_USAGE: mpirun -np [NP] PGFem3D -SS [options] input output\n");
  PGFEM_fprintf(out,"MS_USAGE: mpirun -np [NP] PGFem3D -MS "
		"-macro-np [P] -micro-group-size [S] "
		"[macro OPTION_BLK] [micro OPTION_BLK]\n"
		"OPTION_BLK: -[scale]-start [options] "
		"input output -[scale]-end\n");
  PGFEM_fprintf(out,"\nAnalysis Options:\n");
  print_long_options(out,analysis_opts,n_analysis);
  PGFEM_fprintf(out,"\nSolver Options:\n");
  print_long_options(out,solver_opts,n_solver);
  PGFEM_fprintf(out,"\nPreconditioner Options:\n");
  print_long_options(out,precond_opts,n_precond);
  PGFEM_fprintf(out,"\nVisuzlization Options:\n");
  print_long_options(out,vis_opts,n_vis);
  PGFEM_fprintf(out,"\nOther Options:\n");
  print_long_options(out,other_opts,n_other);
  PGFEM_fprintf(out,"\nDepricated (ignored) Options:\n");
  print_long_options(out,depricated_opts,n_depricated);
} /* print_usage() */

void parse_command_line(const int argc,
			char **argv,
			const int myrank,
			PGFem3D_opt *options)
{
  re_parse_command_line(myrank,1,argc,argv,options);
}/* parse_command_line() */

void print_interpreted_options(const PGFem3D_opt *opts)
{
  switch(opts->analysis_type){
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
    if (opts->cohesive == 0) {
      PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION : stb = %12.5e\n",
		    opts->stab);
    } else if( opts->cohesive == 1) {
      PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION"
		    " WITH COHESIVE FRACTURE : stb = %12.5e\n",
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
    PGFEM_printf("THREE FIELD MIXED METHOD:\n"
		 "TOTAL LAGRANGIAN DISPLACEMENT, PRESSURE, AND VOLUME BASED ELEMENT\n");
    break; 
  case CM:
    PGFEM_printf("USE CONSTITUTIVE MODEL INTERFACE:\n"
		 "HYPERELASTICITY, CRYSTAL PLASTICITY, AND BPA_PLASTICITY\n");
    break;    
       

  default:
    PGFEM_printerr("ERROR: unrecognized analysis type!\n");
    abort();
    break;
  }
    PGFEM_printf ("\n");
    PGFEM_printf ("SolverPackage: ");
    if (opts->solverpackage == BLOCKSOLVE){
      PGFEM_printf ("BlockSolve95 no longer supported !!!\n");
      abort();
    } else if(opts->solverpackage == HYPRE){
      switch(opts->solver){
      case HYPRE_GMRES: PGFEM_printf ("HYPRE - GMRES\n"); break;
      case HYPRE_BCG_STAB: PGFEM_printf ("HYPRE - BiCGSTAB\n"); break;
      case HYPRE_AMG: PGFEM_printf ("HYPRE - BoomerAMG\n"); break;
      case HYPRE_FLEX: PGFEM_printf ("HYPRE - FlexGMRES\n"); break;
      case HYPRE_HYBRID: PGFEM_printf ("HYPRE - Hybrid (GMRES)\n"); break;
      default:
	PGFEM_printerr("Unrecognized solver package!\n");
	abort();
	break;
      }
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
}

void re_parse_command_line(const int myrank,
			   const int start_idx,
			   const int argc,
			   char **argv,
			   PGFem3D_opt *options)
{
  const int n_options = n_analysis + n_solver + n_precond + n_vis + n_other + n_depricated;
  struct option *opts;/* [n_options+1]; */
  opts = (struct option*) PGFEM_calloc(n_options+1,sizeof(struct option));
  copy_options(opts);
  opterr = 0;

  /* print command line to parse */
  if(myrank == 0){
    PGFEM_printf("*** Parsing options from: ");
    for(int i=start_idx; i<argc; i++){
      PGFEM_printf("%s ",argv[i]);
    }
    PGFEM_printf("***\n");
  }

  int ipath, opath;
  ipath = opath = 0;

  /* reset the external variable optind to 1 (do not process the
     initial command) */
  optind = start_idx;

  /* parse command line */
  int opts_idx = 0;
  int opt = getopt_long_only(argc,argv,opt_list,opts,&opts_idx);
  while(opt != -1){
    switch(opt){

      /* HELP OPTIONS */
    case '?':
      if(myrank == 0){
	PGFEM_printf("Skipping unrecognized option :%s\n",argv[optind-1]);
      }
      break;
    case 'h':
      if(myrank == 0){
  	print_usage(stdout);
      }
      exit(0);

      /* ANALYSIS OPTIONS */
      /* Crystal plasticity */
    case 'c': options->analysis_type = FS_CRPL; break;
      /* finite strain elasticity */
    case 'f': options->analysis_type = FINITE_STRAIN; break;
      /* finite strain Portevin-Le Chatelier */
    case 0: options->analysis_type = FS_CRPL; options->plc = 1; break;

    case 1: /* Stabilized */
      if(strcmp("st",opts[opts_idx].name) == 0){
	options->analysis_type = STABILIZED;
	options->stab = atof(optarg);
      } 
      /* Cohesive */
      else if(strcmp("coh",opts[opts_idx].name) == 0){
	options->cohesive = 1;
      }
      break;

    case 2: /* New finite strain formulations */
      if(strcmp("he",opts[opts_idx].name) == 0){
	options->analysis_type = MINI;
      } else if(strcmp("he3",opts[opts_idx].name) == 0){
	options->analysis_type = MINI_3F;
      } else if(strcmp("disp",opts[opts_idx].name) == 0){
	options->analysis_type = DISP;
      } else if(strcmp("tf",opts[opts_idx].name) == 0){
	options->analysis_type = TF;	
      } else if(strcmp("cm",opts[opts_idx].name) == 0){
	options->analysis_type = CM;
	options->cm = atof(optarg);	
      }
      
      break;

      /* SOLVER OPTIONS */
    case 3:
      if(strcmp("gmres",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->solver = HYPRE_GMRES;
      } else if(strcmp("bcgstab",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->solver = HYPRE_BCG_STAB;
      } else if(strcmp("boomer",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->solver = HYPRE_AMG;
      } else if(strcmp("flex-gmres",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->solver = HYPRE_FLEX;
      } else if(strcmp("hybrid",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->solver = HYPRE_HYBRID;
      } else if(strcmp("kdim",opts[opts_idx].name) == 0){
	options->kdim = (int) atof(optarg);
      } else if(strcmp("maxit",opts[opts_idx].name) == 0){
	options->maxit = (int) atof(optarg);
      }
      break;

      /* PRECOND OPTIONS */
    case 4:
      if(strcmp("pre-pilut",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->precond = PILUT;
      } else if(strcmp("pre-euclid",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->precond = EUCLID;
      } else if(strcmp("pre-boomer",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->precond = BOOMER;
      } else if(strcmp("pre-sails",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->precond = PARA_SAILS;
      } else if(strcmp("pre-diag",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->precond = DIAG_SCALE;
      } else if(strcmp("pre-jacobi",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->precond = JACOBI;
      }else if(strcmp("pre-none",opts[opts_idx].name) == 0){
	options->solverpackage = HYPRE;
	options->precond = NONE;
      }

      break;

      /* VISUALIZATION OPTIONS */
    case 'V': options->vis_format = VIS_VTK; break;
    case 'A': options->ascii = 1; break;

      /* OTHER OPTIONS */
    case 'i':
      options->ipath = optarg;
      ipath = 1;
      break;

    case 'o':
      options->opath = optarg;
      opath = 1;
      break;

    case 'O':
      if(strcmp("override-pre-disp",opts[opts_idx].name) == 0){
	options->override_pre_disp = 1;
	options->pre_disp_file = optarg;
      } else  if(strcmp("override-solver-file",opts[opts_idx].name) == 0){
	options->override_solver_file = 1;
	options->solver_file = optarg;
      } else  if(strcmp("override-material-props",opts[opts_idx].name) == 0){
        options->override_material_props = optarg;
      }
      break;

    case 'm':
      if(strcmp("ms",opts[opts_idx].name) == 0){
	options->multi_scale = 1;
      }
      break;

    case 'S':
      options->max_n_jobs = atoi(optarg);
      break;

    case 'N':
      options->no_migrate = 1;
      break;

    case 'l':
      options->legacy = 1;
      break;

    case 'r':
      options->restart = atoi(optarg);
      break;

    case 9999: /* debug mode */
      options->debug = 1;
      break;

    default:
      PGFEM_printf("How did I get here???\n");
      break;
    }

    opt = getopt_long_only(argc,argv,opt_list,opts,&opts_idx);
  }

  if (optind <= argc - 2){
    options->ifname = argv[argc-2];
    options->ofname = argv[argc-1];

    if(!ipath){/* ipath not given, so get from filename */
      char *sp;
      sp = strrchr(argv[argc-2],'/');
      if(sp != NULL){ /* input files not in current directory */
   	*sp = '\0';
   	options->ipath = argv[argc-2];
   	options->ifname = sp + 1;
      } else {
   	options->ifname = argv[argc-2];
   	options->ipath = &cur_dir;
      }
    }

   if(!opath){/* ipath not given, so get from filename */
      char *sp;
      sp = strrchr(argv[argc-1],'/');
      if(sp != NULL){ /* input files not in current directory */
   	*sp = '\0';
   	options->opath = argv[argc-1];
   	options->ofname = sp + 1;
      } else {
   	options->ofname = argv[argc-1];
   	options->opath = &cur_dir;
      }
    }
   
  } else {
    if(myrank == 0){
      PGFEM_printf("ERROR: must provide input and output file names!\n"
		   "Please see usage information below:\n\n");
      print_usage(PGFEM_stdout);
    }
    exit(1);
  }

  free(opts);
}

void get_macro_micro_option_blocks(const int myrank,
				   const int argc,
				   /* const */ char **argv,
				   int *macro_start,
				   int *macro_argc,
				   int *micro_start,
				   int *micro_argc,
				   int *macro_nproc,
				   int *micro_group_size,
				   int *debug)
{
  /* loop through command line and look for markers */
  enum{MACRO_START,
       MACRO_ARGC,
       MICRO_START,
       MICRO_ARGC,
       MACRO_NPROC,
       MICRO_SIZE,
       N_OPT};
  int *got_opt = calloc(N_OPT,sizeof(int));
  for(int i=0; i<argc; i++){
    const char *arg_str = argv[i];

    /* -debug */
    if(strcmp(arg_str,"-debug") == 0){ *debug = 1; continue;}

    /* -help */
    if(strcmp(arg_str,"-help") == 0){
      if(!myrank) print_usage(PGFEM_stdout);
      exit(0);
    }

    /* -macro-np */
    if(!strcmp(arg_str,"-macro-np") && !got_opt[MACRO_NPROC]){
      sscanf(argv[++i],"%d",macro_nproc);
      got_opt[MACRO_NPROC] = 1;
      continue;
    }

    /* -micro-group-size */
    if(!strcmp(arg_str,"-micro-group-size") && !got_opt[MICRO_SIZE]){
      sscanf(argv[++i],"%d",micro_group_size);
      got_opt[MICRO_SIZE] = 1;
      continue;
    }

    /* -macro-start */
    if(!strcmp(arg_str,"-macro-start")
       && !got_opt[MACRO_START] /* have not set option */
       && (got_opt[MICRO_START] == got_opt[MICRO_ARGC]) /* not in micro */
       ){
      got_opt[MACRO_START] = 1;
      *macro_start = i;
      continue;
    }

    /* -macro-end */
    if(!strcmp(arg_str,"-macro-end")
       && (!got_opt[MACRO_ARGC] && got_opt[MACRO_START]) /* in macro */
       && (got_opt[MICRO_START] == got_opt[MICRO_ARGC]) /* not in micro */
       ){
      got_opt[MACRO_ARGC] = 1;
      *macro_argc = i - *macro_start;
      continue;
    }

    /* -micro-start */
    if(!strcmp(arg_str,"-micro-start")
       && !got_opt[MICRO_START] /* have not set option */
       && (got_opt[MACRO_START] == got_opt[MACRO_ARGC]) /* not in macro */
       ){
      got_opt[MICRO_START] = 1;
      *micro_start = i;
      continue;
    }

    /*  -micro-end */
    if(!strcmp(arg_str,"-micro-end")
       && (!got_opt[MICRO_ARGC] && got_opt[MICRO_START]) /* in micro */
       &&(got_opt[MACRO_START] == got_opt[MACRO_ARGC]) /* not in macro */
       ){
      got_opt[MICRO_ARGC] = 1;
      *micro_argc = i - *micro_start;
      continue;
    }

  } /* end parse command line */

  for(int i=0; i<N_OPT; i++){
    if(!got_opt[i] && !myrank){
      PGFEM_printf("ERROR parsing macro/micro option blocks!\n");
      print_usage(PGFEM_stdout);
      abort();
    }
  }

  /* exit function */
  free(got_opt);
}
