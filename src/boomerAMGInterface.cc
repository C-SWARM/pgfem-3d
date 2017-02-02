#include "boomerAMGInterface.h"
#include "PGFEM_mpi.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

typedef enum coarsen_type{CLJP=0,
			  RS=1,/* not recommended */
			  RS3=3,
			  FALGOUT=6, /* default */
			  CLJP_DEBUG=7,
			  PMIS=8, /* lower complexities */
			  PMIS_DEBUG=9,
			  HMIS=10, /* lower complexities */
			  RS1=11, /* not recommended */
			  CGC=21,
			  CGC_E=22} coarsen_type;

typedef enum cycle_type{V_CYCLE=1, W_CYCLE=2} cycle_type;
typedef enum cycle_stage{FINEST, DOWN, UP, COARSEST} cycle_stage;

typedef enum relax_type{RELAX_JACOBI,
			GAUSS_SEIDEL_S, /* very slow */
			GAUSS_SEIDEL_P, /*slow */
			F_SOR, /* default */
			B_SOR,
			S_SOR,
			GAUSS_ELIM /* coarsest level only */} relax_type;

typedef enum relax_order{NATURAL, COARSE_FINE /* default */} relax_order;

typedef enum interpolation_type{CLASS_MOD, /*default */
				LS_INTERP, /* for Geo. Smooth MG */
				DIRECT_SEP_WT,
				MULTI_PASS,
				MULTI_PASS_SEP_WT,
				E_CLASS_MOD, /* recommend use w/ PMIS/HMIS */
				E_NO_C_CLASS_MOD, /* recommend use w/ PMIS/HMIS */
				STANDARD,
				STANDARD_S_WT,
				CLASS_BLOC, /* nodal systems only */
				CLASS_BLOCK_DIAG, /* nodal systems only */
				FF_INTERP,
				FF1_INTERP} interpolation_type;

typedef enum smooth_type{SM_SCHWARZ=6, /* (default) Extra option routines:
					  HYPRE_BoomerAMG...
					  ...SetDomainType,
					  ...SetOverlap,
					  ...SetVariant,
					  ...SetSchwarzRlxWeight,
					  ... SetSchwarzUseNonSymm */
			 SM_PILUT=7, /* ...SetDropTol, ...SetMaxNzPerRow */
			 SM_SAILS=8, /* ...SetSym, ...SetLevel, ...SetFilter,
					...SetThreshold */
			 SM_EUCLID=9 /* ...SetEuclidFile */} smooth_type;

typedef enum print_level{PRINT_NONE,
			 PRINT_SETUP,
			 PRINT_SOLVE,
			 PRINT_ALL} print_level;

typedef enum schwarz_variant{SV_H_MULT, /* (default) no proc bnd overlap */
			     SV_H_ADD, /* no proc bnd overlap */
			     SV_ADDITIVE,
			     SV_H_MULT_BND} schwarz_variant;

typedef enum schwarz_overlap{SOL_NONE,
			     SOL_MIN, /* (default) */
			     SOL_NEIGHBORS} schwarz_overlap;

typedef enum schwarz_domain{SD_POINT,
			    SD_NODE,
			    SD_GEN /* (default) */} schwarz_domain;

typedef enum schwarz_sym{S_SPD /* (default) */,S_NS} schwarz_sym;



int initializeBoomerAMG(HYPRE_Solver *hypre_pc,
			 boomerAMGOptions *options,
			 const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);
  err += HYPRE_BoomerAMGCreate(hypre_pc);

  if (options->useDefault){

    err += HYPRE_BoomerAMGSetCoarsenType(*hypre_pc, options->coarsenType);
    err += HYPRE_BoomerAMGSetStrongThreshold(*hypre_pc, options->strongThreshold);
    err += HYPRE_BoomerAMGSetMaxIter(*hypre_pc, options->maxIterations);
    err += HYPRE_BoomerAMGSetPrintLevel(*hypre_pc, options->printLevel);
    err += HYPRE_BoomerAMGSetMaxLevels(*hypre_pc,options->maxMGLevels);
    err += HYPRE_BoomerAMGSetNumSweeps(*hypre_pc,options->numSweeps);
    err += HYPRE_BoomerAMGSetSmoothNumLevels(*hypre_pc,options->smoothNumLevels);

    err += HYPRE_BoomerAMGSetSchwarzUseNonSymm(*hypre_pc,S_NS);

    err += HYPRE_BoomerAMGSetCycleType(*hypre_pc,options->cycleType);
    err += HYPRE_BoomerAMGSetRelaxType(*hypre_pc,options->relaxType);
    err += HYPRE_BoomerAMGSetAggNumLevels(*hypre_pc,options->coarsenNumAggLevels);

    if(myrank == 0){
      PGFEM_printf("***BoomerAMG Options:***\n");
      PGFEM_printf("\tCoarsen type:     %d\n",options->coarsenType);
      PGFEM_printf("\tMax Levels:       %d\n",options->maxMGLevels);
      PGFEM_printf("\tNum Sweeps:       %d\n",options->numSweeps);
      PGFEM_printf("\tCycle Type:       %d\n",options->cycleType);
      PGFEM_printf("\tStrong Threshold: %f\n",options->strongThreshold);
    }

  }

  return err;
}

void setBoomerAMGOptions(boomerAMGOptions *options)
{
  /* Currently must change option values here.  In future, will allow
     reading of a file for option parameters */

  options->useDefault = 1;

  /* These are the default options->  If useDefault != 0, they may be
     overwritten. */
  options->maxIterations = 1;
  options->strongThreshold = 0.9;
  options->printLevel = PRINT_NONE;
  options->coarsenType = FALGOUT;
  options->maxMGLevels = 25;
  options->cycleType = V_CYCLE;
  options->numSweeps = 3;
  options->relaxType = F_SOR;
  options->coarsenNumAggLevels = 0; 
  options->smoothNumLevels = 10;

/*   if(!options->useDefault){ */
/*   /\* Options for tuning solver/preconditioner to specific problem *\/ */
/*   /\* General AMG method options *\/ */
/*   options->maxMGLevels = ; */
/*   options->maxRowSum = ; */
/*   options->measureType = ; */
/*   options->cycleType = ; */

/*   /\* sweep options *\/ */
/*   options->numSweeps = ; */
/*   options->numSweepsAtCycle = ; */
/*   options->sweepCycleLevel = ; */

/*   /\* relaxation options *\/ */
/*   options->relaxType = ; */
/*   options->cycleRelaxType = ; */
/*   options->relaxCycleLevel = ; */
/*   options->relaxOrder = ; */
/*   options->relaxWeight = ; */
/*   options->levelRelaxWeight = ; */
/*   options->relaxLevel = ; */
/*   options->omega = ; */
/*   options->levelOmega = ; */
/*   options->omegaLevel = ; */

/*   /\* interpolation options *\/ */
/*   options->truncationFactor = ; */
/*   options->interpMaxElem = ; */
/*   options->interpType = ; */

/*   /\* smoothing options *\/ */
/*   options->smoothType = ; */
/*   options->smoothNumLevels = ; */
/*   options->smoothNumSweeps = ; */

/*   /\* Schwarz smoothing options *\/ */
/*   options->schwarzVariant = ; */
/*   options->schwarzOverlap = ; */
/*   options->schwarzDomainType = ; */
/*   options->schwarzRelaxWeight = ; */

/*   /\* ParaSAILS smoothing options *\/ */
/*   options->psSymmetry = ; */
/*   options->psNumLevels = ; */
/*   options->psThreshold = ; */
/*   options->psFilter = ; */

/*   /\* PILUT smoothing options *\/ */
/*   options->pilutDropTol = ; */
/*   options->pilutMaxNZRow = ; */

/*   /\* coarsening options *\/ */
/*   options->coarsenNumAggLevels = ; */
/*   options->coarsenAggNumPaths = ; */

/*   /\* GSMG and LS smoothing options *\/ */
/*   options->useGSMG = ; */
/*   options->numSamples = ; */

/*   /\* Output options *\/ */
/*   options->boomerDebug = ; */
/*   options->logging = ; */
/*   } */
}
