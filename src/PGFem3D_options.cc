#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

/** This file contains functions to parse the command line and print
    helpful information */
#include "pgfem3d/Communication.hpp"
#include "PGFem3D_options.h"
#include "PGFEM_io.h"
#include "allocation.h"
#include "constitutive_model.h"
#include "enumerations.h"
#include "utils.h"
#include <getopt.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>

using namespace pgfem3d;

namespace {
/* Generated at http://patorjk.com/software/taag/ */
const char *prog_name =
 " _______    ______   ________                        ______   _______  \n"
 "/       \\  /      \\ /        |                      /      \\ /       \\ \n"
 "$$$$$$$  |/$$$$$$  |$$$$$$$$/______   _____  ____  /$$$$$$  |$$$$$$$  |\n"
 "$$ |__$$ |$$ | _$$/ $$ |__  /      \\ /     \\/    \\ $$ ___$$ |$$ |  $$ |\n"
 "$$    $$/ $$ |/    |$$    |/$$$$$$  |$$$$$$ $$$$  |  /   $$< $$ |  $$ |\n"
 "$$$$$$$/  $$ |$$$$ |$$$$$/ $$    $$ |$$ | $$ | $$ | _$$$$$  |$$ |  $$ |\n"
 "$$ |      $$ \\__$$ |$$ |   $$$$$$$$/ $$ | $$ | $$ |/  \\__$$ |$$ |__$$ |\n"
 "$$ |      $$    $$/ $$ |   $$       |$$ | $$ | $$ |$$    $$/ $$    $$/ \n"
 "$$/        $$$$$$/  $$/     $$$$$$$/ $$/  $$/  $$/  $$$$$$/  $$$$$$$/  \n";

const char cur_dir = '.';

using LongOption = option;

/// An option structure that allows us to associate descriptions with getopt
/// options.
struct Option {
  LongOption opt;
  const char *descr;
  int sc;

  void print(FILE* out) const {
    if (sc == 1 && opt.has_arg == no_argument) {
      /* single character alias & no arg */
      PGFEM_fprintf(out, "-%s [-%c]\t%s\n", opt.name, opt.val, descr);
    }
    else if (sc == 1 && opt.has_arg == required_argument) {
      /* single character alias & requires arg */
      PGFEM_fprintf(out, "-%s [-%c] (arg)\t%s\n", opt.name, opt.val, descr);
    }
    else if (sc == 0 && opt.has_arg == no_argument) {
      /* no alias no arg */
      PGFEM_fprintf(out,"-%s\t%s\n", opt.name, descr);
    }
    else if (sc == 0 && opt.has_arg == required_argument) {
      PGFEM_fprintf(out,"-%s (arg)\t%s\n", opt.name, descr);
    }
  }
};

const Option analysis_opts[] = {
  /* Analysis options */
  {{"cp",no_argument,NULL,'c'},"\tFinite strain crystal plasticity",0},
  {{"fd",no_argument,NULL,'f'},"\tFinite strain elasticity (Quadradic tetras or hexas)",0},
  {{"cpp",no_argument,NULL,0},"\tFinite strain Portevin-Le Chatelier effect",0},
  {{"st",required_argument,NULL,1},"Stabilized finite strains analysis",0},
  {{"he",no_argument,NULL,2},"\tFinite strains analysis using MINI element",0},
  {{"he3",no_argument,NULL,2},"\tFinite strains analysis using bubble-enhanced 3-Field element",0},
  {{"disp",no_argument,NULL,2},"\tTOTAL Lagrangian displacement-based finite strains analysis",0},
  {{"disp-cm",no_argument,NULL,2},"(BETA) TOTAL Lagrangian displacement-based finite strains analysis\n"
   "\t\t Provides access to CM models",0},
  {{"tf",no_argument,NULL,2},"\tTOTAL Lagrangian displacement-based 3 field finite strains analysis",0},
  {{"cm",required_argument,NULL,2},"Use of constitutive model interface\n"
   "\t\t arg = 0: Updated Lagrangian\n"
   "\t\t arg = 1: Total Lagrangian\n"
   "\t\t arg = 2: Mixed analysis mode",0},
  {{"cm3f",required_argument,NULL,2},"Three-field mixed method with constitutive model interface\n"
   "\t\t arg = 0: Updated Lagrangian\n"
   "\t\t arg = 1: Total Lagrangian\n"
   "\t\t arg = 2: Mixed analysis mode",0},
  {{"coh",no_argument,NULL,1},"\tCohesive elements",0},
  {{"ms",no_argument,NULL,'m'},("\tINTERFACE or BULK multiscale modeling.\n"
                                "\t\tRequires six (6) or nine (9), respectively, prescribed displacements\n"
                                "\t\tand a file named \"normal.in\" in the specified\n"
                                "\t\tinput directory containing the macroscopic normal,\n"
                                "\t\tbounding volume, and cell thickness. Format: [V lc Nx Ny Nz]."),0}
};

const Option solver_opts[] = {
  /* Solver options */
  {{"gmres",no_argument,NULL,3},"\tUse the HYPRE GMRES solver (default)",0},
  {{"boomer",no_argument,NULL,3},"\tUse the HYPRE BoomerAMG solver (no precond)",0},
  {{"bcgstab",no_argument,NULL,3},"Use the HYPRE BiCGSTAB solver",0},
  {{"flex-gmres",no_argument,NULL,3},"Use the HYPRE FlexGMRES solver",0},
  {{"hybrid",no_argument,NULL,3},"\tUse the HYPRE Hybrid (GMRES) solver\n"
   "\t\t(adaptively switches preconditioner)",0},
  {{"kdim",required_argument,NULL,3},"Set the Krylov dimension",0},
  {{"maxit",required_argument,NULL,3},"Set the maximum number of iterations",0},
  {{"noLS",no_argument,NULL,3},"\tNo use line search, default = Yes",0},
  {{"at",no_argument,NULL,3},"\tUse adaptive time stepping, default = No",0},
  {{"noCCE",no_argument,NULL,3},"\tNo converge check on energy norm, default = Yes",0}
};

const Option precond_opts[] = {
  /* Preconditioner options */
  {{"pre-euclid",no_argument,NULL,4},"Use the HYPRE Euclid preconditioner (default)",0},
  {{"pre-boomer",no_argument,NULL,4},"Use the HYPRE BoomerAMG preconditioner",0},
  {{"pre-pilut",no_argument,NULL,4},"Use the HYPRE PILUT preconditioner",0},
  {{"pre-sails",no_argument,NULL,4},"Use the HYPRE ParaSAILS preconditioner",0},
  {{"pre-diag",no_argument,NULL,4},"Use the custom diagonal scaling preconditioner",0},
  {{"pre-jacobi",no_argument,NULL,4},"Use the Jacobi scaling preconditioner",0},
  {{"pre-none",no_argument,NULL,4},"Do not use a preconditioner",0}
};

const Option network_opts[] = {
  /* Network options */
  {{"isir",no_argument,NULL,5},"Use Isend/Irecv network (MPI) (default)",0},
  {{"pwc",no_argument,NULL,5},"Use put-with-completion (PWC) network (Photon)",0},
};
  
const Option vis_opts[] = {
  /* Visualization options */
  {{"vtk",no_argument,NULL,'V'},"Output in VTK format",1},
  {{"ascii",no_argument,NULL,'A'},"Additional output in ASCII format",1},
};

const Option other_opts[] = {
  /* Other options */
  {{"ipath",required_argument,NULL,'i'},"Path to input files parent directory",0},
  {{"opath",required_argument,NULL,'o'},"Path to output files parent directory",0},
  {{"override-pre-disp",required_argument,NULL,'O'},("\n\t\tOverride the prescribed displacements in *.in\n"
                                                     "\t\twith those provided in the given file."),0},
  {{"override-solver-file",required_argument,NULL,'O'},("\n\t\tOverride the default solver filename with\n"
                                                        "\t\tthe provided filename."),0},
  {{"override-material-props",required_argument,NULL,'O'},("\n\t\tOverride the material properties in *.in\n"
                                                           "\t\twith those provided in the given file."),0},
  {{"restart",required_argument,NULL,'r'},("Restart from specified step. Requires original\n"
                                           "\t\tinput files and dumped restart files for specified step."),0},
  {{"max-server-jobs",required_argument,NULL,'S'},("\n\t\tSet the maximum number of jobs allowed on a server (FE2)."),0},
  {{"no-migrate",no_argument,NULL,'N'},("Do not migrate cells between servers (FE2)."),0},
  {{"legacy",no_argument,NULL,'l'},"Read files from legacy format",1},
  {{"debug",no_argument,NULL,9999},"\tSend into infinite loop to attach debugger",0},
  {{"help",no_argument,NULL,'h'},"Print this message and exit",1},
  {{"no-compute-reactions",no_argument,NULL,'R'},"\n\t\tNo compute and print reaction forces",0},
  {{"no-compute-macro",no_argument,NULL,'M'},"\n\t\tNo compute and print macro values (GF,GS,GP)",0},
  {{"walltime",required_argument,NULL,'w'},("\n\t\tSet Walltime[s] and write restart files nearby this walltime.\n"
                                            "\t\tDefault is -1.0 (no actions)"),0},
};

/* these options may no longer be supported/functional. They are kept
   for documentation purposes only and are ignored if used */
const Option depricated_opts[] = {
  {{"elixir",no_argument,NULL,'X'},"Output in Elixir format [unsupported, use -V]",1},
  {{"ensight",no_argument,NULL,'E'},"Output in EnSight format [outdated, use -V]",1},
  {{"sm",no_argument,NULL,'s'},"Smooth stress field [unsupported]",1},
  {{"me",no_argument,NULL,'m'},"\tCompute nodal forces on model entities"
   " (requires entities.in file) [unsupported]",0},
  {{"pr",no_argument,NULL,'p'},"Periodic domain [outdated, use -ms]",1},
  {{"rn",no_argument,NULL,'r'},"\tRenumber degrees of freedom [unsupported]",0},
};

const Option null_opts[] = {
  {{NULL,0,NULL,0},nullptr,0}
};

/// The total number of options that we have defined.
constexpr const int N_OPTS = {
  size(analysis_opts) + size(network_opts) + size(solver_opts) + size(precond_opts) +
  size(vis_opts) + size(other_opts) + size(depricated_opts) + size(null_opts)
};

void initialize(LongOption *A) {
  for (const auto& o : analysis_opts)   { *A++ = o.opt; }
  for (const auto& o : network_opts)    { *A++ = o.opt; }
  for (const auto& o : solver_opts)     { *A++ = o.opt; }
  for (const auto& o : precond_opts)    { *A++ = o.opt; }
  for (const auto& o : vis_opts)        { *A++ = o.opt; }
  for (const auto& o : other_opts)      { *A++ = o.opt; }
  for (const auto& o : depricated_opts) { *A++ = o.opt; }
  for (const auto& o : null_opts)       { *A++ = o.opt; }
}
} // namespace

void set_default_options(PGFem3D_opt *options)
{
  /* network option */
  options->network = NETWORK_ISIR;
  
  /* solver options */
  options->solverpackage = HYPRE;
  options->solver = SOLVER_GMRES;
  options->precond = PRECOND_EUCLID;
  options->kdim = 500;
  options->maxit = 1000;
  options->solution_scheme_opt[LINE_SEARCH]              = 1;
  options->solution_scheme_opt[ADAPTIVE_TIME_STEPPING]   = 0;
  options->solution_scheme_opt[CVG_CHECK_ON_ENERGY_NORM] = 1;

  /* analysis options */
  options->analysis_type = -1;
  options->stab = 0;
  options->cohesive = 0;
  options->plc = 0;
  options->multi_scale = 0;
  options->cm = -1;

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
  options->comp_print_reaction = 1;
  options->comp_print_macro = 1;
  options->print_EXA_details = false;

  /* I/O file names */
  options->ipath = NULL;
  options->opath = NULL;
  options->ifname = NULL;
  options->ofname = NULL;
  options->walltime = -1.0;
}

void print_options(FILE *out, const PGFem3D_opt *options)
{
  PGFEM_fprintf(out,"OPTION VALUES:\n");
  PGFEM_fprintf(out,"=== SOLVER OPTIONS ===\n");
  PGFEM_fprintf(out,"Solver package: %d\n", options->solverpackage);
  PGFEM_fprintf(out,"Solver:         %d\n", options->solver);
  PGFEM_fprintf(out,"Preconditioner: %d\n", options->precond);
  PGFEM_fprintf(out,"Kdim:           %d\n", options->kdim);
  PGFEM_fprintf(out,"Max It:         %d\n", options->maxit);
  if(options->solution_scheme_opt[LINE_SEARCH])
    PGFEM_fprintf(out,"LINE SEARCH is enabled\n");
  else
    PGFEM_fprintf(out,"LINE SEARCH is disabled\n");

  if(options->solution_scheme_opt[ADAPTIVE_TIME_STEPPING])
    PGFEM_fprintf(out,"ADAPTIVE TIME STEPPING is enabled\n");
  else
    PGFEM_fprintf(out,"ADAPTIVE TIME STEPPING is disabled\n");

  if(options->solution_scheme_opt[CVG_CHECK_ON_ENERGY_NORM])
    PGFEM_fprintf(out,"EXIT WITH ENERGY NORM IS CONVERGED is enabled\n");
  else
    PGFEM_fprintf(out,"EXIT WITH ENERGY NORM IS CONVERGED is disabled\n");

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
  PGFEM_fprintf(out,"Walltime:       %f[s]\n", options->walltime);
  PGFEM_fprintf(out,"Network:        %d\n",options->network);
  
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
  PGFEM_fprintf(out,"MS_USAGE: mpirun -np [NP] PGFem3D -MS [network] "
                "-macro-np [P] -micro-group-size [S] "
                "[macro OPTION_BLK] [micro OPTION_BLK]\n"
                "OPTION_BLK: -[scale]-start [options] "
                "input output -[scale]-end\n");
  PGFEM_fprintf(out,"\nAnalysis Options:\n");
  for (const auto& opt : analysis_opts) { opt.print(out); }
  PGFEM_fprintf(out,"\nNetwork Options:\n");
  for (const auto& opt : network_opts) { opt.print(out); }
  PGFEM_fprintf(out,"\nSolver Options:\n");
  for (const auto& opt : solver_opts) { opt.print(out); }
  PGFEM_fprintf(out,"\nPreconditioner Options:\n");
  for (const auto& opt : precond_opts) { opt.print(out); }
  PGFEM_fprintf(out,"\nVisualization Options:\n");
  for (const auto& opt : vis_opts) { opt.print(out); }
  PGFEM_fprintf(out,"\nOther Options:\n");
  for (const auto& opt : other_opts) { opt.print(out); }
  PGFEM_fprintf(out,"\nDepricated (ignored) Options:\n");
  for (const auto& opt : depricated_opts) { opt.print(out); }
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
  switch (opts->analysis_type) {
   case ELASTIC:
    PGFEM_printf("ELASTIC ANALYSIS\n");
    break;
   case TP_ELASTO_PLASTIC:
    PGFEM_printf("TWO PHASE COMPOSITE SYSTEM : ELASTO-PLASTIC ANALYSIS\n");
    break;
   case FS_CRPL:
    PGFEM_printf("FINITE STRAIN CRYSTAL ELASTO-PLASTICITY\n");
    break;
   case FINITE_STRAIN:
    if (opts->cohesive == 0) {
      PGFEM_printf("FINITE STRAIN ELASTICITY\n");
    } else {
      PGFEM_printf("FINITE STRAIN ELASTICITY WITH COHESIVE FRACTURE\n");
    }
    break;
   case STABILIZED:
    if (opts->cohesive == 0 && opts->gem == 0) {
      PGFEM_printf("FINITE STRAIN STABILIZED FORMULATION : stb = %12.5e\n",
                   opts->stab);
    }
    else if (opts->cohesive == 1) {
      PGFEM_printf("FINITE STRAIN STABILIZED FORMULATION"
                   " WITH COHESIVE FRACTURE : stb = %12.5e\n",
                   opts->stab);
    }
    else if (opts->gem == 1) {
      PGFEM_printf("GENERALIZED FINITE ELEMENT METHOD\n");
      PGFEM_printf("FINITE STRAIN STABILIZED FORMULATION : stb = %12.5e\n",
                    opts->stab);
    }
    else {
      PGFEM_Abort();
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
    // @todo this text was in main.cc... don't know which one is correct
    // PGFEM_printf("FINITE STRAIN TREE FIELDS HYPERELASTICITY:\n"
    //              "TOTAL LAGRANGIAN TREE FIELDS-BASED ELEMENT\n");

    PGFEM_printf("THREE FIELD MIXED METHOD:\n"
                 "TOTAL LAGRANGIAN DISPLACEMENT, PRESSURE, AND VOLUME BASED ELEMENT\n");
    break;
   case CM:
     PGFEM_printf("USE CONSTITUTIVE MODEL INTERFACE:\n"
                  "UPDATED LAGRANGIAN, TOTAL LAGRANGIAN, AND MIXED ANALYSIS MODE\n");
     break;
   case CM3F:
     PGFEM_printf("USE CONSTITUTIVE MODEL INTERFACE: ");
     switch (opts->cm) {
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
   default:
    PGFEM_printerr("ERROR: unrecognized analysis type!\n");
    PGFEM_Abort();
    break;
  }
  PGFEM_printf("\n");

  if (opts->solverpackage == BLOCKSOLVE) {
    PGFEM_printerr("BlockSolve95 no longer supported !!!\n");
    PGFEM_Abort();
  }

  if (opts->solver >= SOLVER_MAX) {
    PGFEM_printerr("Unrecognized solver package!\n");
    PGFEM_Abort();
  }

  PGFEM_printf ("Network: %s\n",
                NETWORK_OPTS[opts->network]);
  
  PGFEM_printf ("SolverPackage: %s - %s\n",
                SOLVER_PACKAGE_OPTS[opts->solverpackage],
                SOLVER_OPTS[opts->solver]);

  
  PGFEM_printf("Preconditioner: ");
  switch (opts->precond) {
   case PRECOND_PARA_SAILS: PGFEM_printf ("HYPRE - PARASAILS\n"); break;
   case PRECOND_PILUT: PGFEM_printf ("HYPRE - PILUT\n"); break;
   case PRECOND_EUCLID: PGFEM_printf ("HYPRE - EUCLID\n"); break;
   case PRECOND_BOOMER: PGFEM_printf ("HYPRE - BoomerAMG\n"); break;
   case PRECOND_NONE: PGFEM_printf ("PGFEM3D - NONE\n"); break;
   case PRECOND_DIAG_SCALE: PGFEM_printf ("PGFEM3D - DIAGONAL SCALE\n"); break;
   case PRECOND_JACOBI: PGFEM_printf ("PGFEM3D - JACOBI\n"); break;
  }
}

void re_parse_command_line(const int myrank,
                           const int start_idx,
                           const int argc,
                           char **argv,
                           PGFem3D_opt *options)
{
  // Build the array of options that we need for the getopt API by extracting
  // them from our Option arrays.
  LongOption opts[N_OPTS];
  initialize(opts);
  opterr = 0;

  options->cm = -1; //default: no use of constitutive model

  /* print command line to parse */
  if (myrank == 0) {
    PGFEM_printf("*** Parsing options from: ");
    for(int i = start_idx, e = argc; i < e; ++i) {
      PGFEM_printf("%s ", argv[i]);
    }
    PGFEM_printf("***\n");
  }

  int ipath = 0, opath = 0;

  /* reset the external variable optind to 1 (do not process the
     initial command) */
  optind = start_idx;

  /* parse command line */
  int opt, opts_idx;
  while ((opt = getopt_long_only(argc, argv, "AVlh", opts, &opts_idx)) != -1)
  {
    switch (opt) {
      /* HELP OPTIONS */
     case '?':
      if (myrank == 0) {
        PGFEM_printf("Skipping unrecognized option :%s\n", argv[optind-1]);
      }
      break;
     case 'h':
      if (myrank == 0) {
        print_usage(stdout);
      }
      exit(0);

      /* ANALYSIS OPTIONS */
     case 'c':
      /* Crystal plasticity */
      options->analysis_type = FS_CRPL;
      break;
     case 'f':
      /* finite strain elasticity */
      options->analysis_type = FINITE_STRAIN;
      break;
     case 0:
      /* finite strain Portevin-Le Chatelier */
      options->analysis_type = FS_CRPL;
      options->plc = 1;
      break;
     case 1:
      /* Stabilized */
      if (strcmp("st", opts[opts_idx].name) == 0) {
        options->analysis_type = STABILIZED;
        options->stab = strtod(optarg, nullptr);
      }
      /* Cohesive */
      else if (strcmp("coh", opts[opts_idx].name) == 0) {
        options->cohesive = 1;
      }
      break;

     case 2: /* New finite strain formulations */
      if (strcmp("he", opts[opts_idx].name) == 0) {
        options->analysis_type = MINI;
      }
      else if (strcmp("he3", opts[opts_idx].name) == 0) {
        options->analysis_type = MINI_3F;
      }
      else if(strcmp("disp", opts[opts_idx].name) == 0) {
        options->analysis_type = DISP;
      }
      else if(strcmp("tf", opts[opts_idx].name) == 0) {
        options->analysis_type = TF;
      }
      else if(strcmp("cm", opts[opts_idx].name) == 0) {
        options->analysis_type = CM;
        options->cm = strtof(optarg, nullptr);
      }
      else if(strcmp("cm3f", opts[opts_idx].name) == 0) {
        options->analysis_type = CM3F;
        options->cm = strtof(optarg, nullptr);
      }
      else if(strcmp("disp-cm", opts[opts_idx].name) == 0) {
        options->analysis_type = CM;
        options->cm = DISP;
      }

      break;

      /* SOLVER OPTIONS */
     case 3:
      if (strcmp("gmres", opts[opts_idx].name) == 0) {
        options->solver = SOLVER_GMRES;
      }
      else if (strcmp("bcgstab", opts[opts_idx].name) == 0) {
        options->solver = SOLVER_BCG_STAB;
      }
      else if (strcmp("boomer", opts[opts_idx].name) == 0) {
        options->solver = SOLVER_AMG;
      }
      else if (strcmp("flex-gmres", opts[opts_idx].name) == 0) {
        options->solver = SOLVER_FLEX;
      }
      else if (strcmp("hybrid", opts[opts_idx].name) == 0) {
        options->solver = SOLVER_HYBRID;
      }
      else if (strcmp("kdim", opts[opts_idx].name) == 0) {
        options->kdim = strtol(optarg, nullptr, 10);
      }
      else if (strcmp("maxit", opts[opts_idx].name) == 0) {
        options->maxit = strtol(optarg, nullptr, 10);
      }
      else if (strcmp("noLS", opts[opts_idx].name) == 0) {
        options->solution_scheme_opt[LINE_SEARCH] = 0;
      }
      else if (strcmp("at", opts[opts_idx].name) == 0) {
        options->solution_scheme_opt[ADAPTIVE_TIME_STEPPING] = 1;
      }
      else if (strcmp("noCCE", opts[opts_idx].name) == 0) {
        options->solution_scheme_opt[CVG_CHECK_ON_ENERGY_NORM] = 0;
      }
      break;

      /* PRECOND OPTIONS */
     case 4:
      if (strcmp("pre-pilut", opts[opts_idx].name) == 0) {
        options->precond = PRECOND_PILUT;
      }
      else if(strcmp("pre-euclid", opts[opts_idx].name) == 0) {
        options->precond = PRECOND_EUCLID;
      }
      else if (strcmp("pre-boomer", opts[opts_idx].name) == 0) {
        options->precond = PRECOND_BOOMER;
      }
      else if (strcmp("pre-sails", opts[opts_idx].name) == 0) {
        options->precond = PRECOND_PARA_SAILS;
      }
      else if (strcmp("pre-diag", opts[opts_idx].name) == 0) {
        options->precond = PRECOND_DIAG_SCALE;
      }
      else if (strcmp("pre-jacobi", opts[opts_idx].name) == 0) {
        options->precond = PRECOND_JACOBI;
      }
      else if (strcmp("pre-none", opts[opts_idx].name) == 0) {
        options->precond = PRECOND_NONE;
      }
      break;

      /* NETWORK OPTIONS */
    case 5:
      if (strcmp("isir", opts[opts_idx].name) == 0) {
	options->network = NETWORK_ISIR;
      }
      else if(strcmp("pwc", opts[opts_idx].name) == 0) {
        options->network = NETWORK_PWC;
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
      if (strcmp("override-pre-disp", opts[opts_idx].name) == 0) {
        options->override_pre_disp = 1;
        options->pre_disp_file = optarg;
      }
      else if (strcmp("override-solver-file", opts[opts_idx].name) == 0) {
        options->override_solver_file = 1;
        options->solver_file = optarg;
      }
      else if (strcmp("override-material-props", opts[opts_idx].name) == 0) {
        options->override_material_props = optarg;
      }
      break;

     case 'm':
      if (strcmp("ms", opts[opts_idx].name) == 0) {
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

     case 'R':
      options->comp_print_reaction = 0;
      break;

     case 'M':
      options->comp_print_macro = 0;
      break;

     case 'w':
      options->walltime = atoi(optarg);
      break;

     default:
      PGFEM_printf("How did I get here???\n");
      PGFEM_Abort();
      break;
    }
  }

  if (optind <= argc - 2) {
    options->ifname = argv[argc-2];
    options->ofname = argv[argc-1];

    if (!ipath) { /* ipath not given, so get from filename */
      if (char* sp = strrchr(argv[argc-2],'/')) { /* input files not in current directory */
        *sp = '\0';
        options->ipath = argv[argc-2];
        options->ifname = sp + 1;
      } else {
        options->ifname = argv[argc-2];
        options->ipath = &cur_dir;
      }
    }

    if (!opath) { /* ipath not given, so get from filename */
      if (char* sp = strrchr(argv[argc-1],'/')) { /* input files not in current directory */
        *sp = '\0';
        options->opath = argv[argc-1];
        options->ofname = sp + 1;
      }
      else {
        options->ofname = argv[argc-1];
        options->opath = &cur_dir;
      }
    }
  }
  else {
    if (myrank == 0) {
      PGFEM_printf("ERROR: must provide input and output file names!\n"
                   "Please see usage information below:\n\n");
      print_usage(PGFEM_stdout);
    }
    exit(1);
  }
}

/* loop through command line and look for markers */
void get_macro_micro_option_blocks(int myrank,
                                   int argc,
                                   char *argv[],
                                   int *macro_start,
                                   int *macro_argc,
                                   int *micro_start,
                                   int *micro_argc,
                                   int *macro_nproc,
                                   int *micro_group_size,
                                   int *debug)
{
  enum {
    MACRO_START,
    MACRO_ARGC,
    MICRO_START,
    MICRO_ARGC,
    MACRO_NPROC,
    MICRO_SIZE,
    N_OPT
  };

  int got_opt[N_OPT] = {};

  for (int i = 0, e = argc; i < e; ++i) {
    const char *arg = argv[i];

    /* -debug */
    if (strcmp(arg,"-debug") == 0) {
      *debug = 1;
      continue;
    }

    /* -help */
    if (strcmp(arg,"-h") == 0 and myrank == 0) {
      print_usage(PGFEM_stdout);
      exit(0);
    }

    if (strcmp(arg,"-help") == 0 and myrank == 0) {
      print_usage(PGFEM_stdout);
      exit(0);
    }

    if (strcmp(arg,"--help") == 0 and myrank == 0) {
      print_usage(PGFEM_stdout);
      exit(0);
    }

    /* -macro-np */
    if (!strcmp(arg,"-macro-np") and !got_opt[MACRO_NPROC]) {
      sscanf(argv[++i],"%d",macro_nproc);
      got_opt[MACRO_NPROC] = 1;
      continue;
    }

    /* -micro-group-size */
    if (!strcmp(arg,"-micro-group-size") and !got_opt[MICRO_SIZE]) {
      sscanf(argv[++i],"%d",micro_group_size);
      got_opt[MICRO_SIZE] = 1;
      continue;
    }

    /* -macro-start */
    if (!strcmp(arg,"-macro-start")
        and !got_opt[MACRO_START] /* have not set option */
        and (got_opt[MICRO_START] == got_opt[MICRO_ARGC])) /* not in micro */
    {
      got_opt[MACRO_START] = 1;
      *macro_start = i;
      continue;
    }

    /* -macro-end */
    if (!strcmp(arg,"-macro-end")
        and (!got_opt[MACRO_ARGC] && got_opt[MACRO_START]) /* in macro */
        and (got_opt[MICRO_START] == got_opt[MICRO_ARGC])) /* not in micro */
    {
      got_opt[MACRO_ARGC] = 1;
      *macro_argc = i - *macro_start;
      continue;
    }

    /* -micro-start */
    if (!strcmp(arg,"-micro-start")
        and !got_opt[MICRO_START] /* have not set option */
        and (got_opt[MACRO_START] == got_opt[MACRO_ARGC])) /* not in macro */
    {
      got_opt[MICRO_START] = 1;
      *micro_start = i;
      continue;
    }

    /*  -micro-end */
    if (!strcmp(arg,"-micro-end")
        and (!got_opt[MICRO_ARGC] && got_opt[MICRO_START]) /* in micro */
        and (got_opt[MACRO_START] == got_opt[MACRO_ARGC])) /* not in macro */
    {
      got_opt[MICRO_ARGC] = 1;
      *micro_argc = i - *micro_start;
      continue;
    }
  } /* end parse command line */

  for (int i = 0, e = N_OPT; i < e; ++i) {
    if (!got_opt[i] && !myrank) {
      PGFEM_printf("ERROR parsing macro/micro option blocks!\n");
      print_usage(PGFEM_stdout);
      PGFEM_Abort();
    }
  }
}
