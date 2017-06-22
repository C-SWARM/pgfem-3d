#ifndef BOOMER_AMG_INTERFACE_H
#define BOOMER_AMG_INTERFACE_H

#include "data_structure.h"
#include "PGFEM_mpi.h"

#ifndef HYPRE_GLOBAL_H
#include "hypre_global.h"
#endif


#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */


typedef struct boomerAMGOptions{

  /* Required options to set */
  int useDefault;

  /* Options for tuning solver/precondiditioner to specific problem */
  /* General AMG method options */
  /* double tolerance; meaningless when used as precon */
  int maxIterations;
  int maxMGLevels;
  double strongThreshold;
  double maxRowSum;
  int coarsenType;
  int measureType;
  int cycleType;

  /* sweep options */
  int numSweeps;
  int numSweepsAtCycle;
  int sweepCycleLevel;

  /* relaxation options */
  int relaxType;
  int cycleRelaxType;
  int relaxCycleLevel;
  int relaxOrder;
  double relaxWeight;
  double levelRelaxWeight;
  int relaxLevel;
  double omega;
  double levelOmega;
  int omegaLevel;

  /* interpolation options */
  double truncationFactor;
  int interpMaxElem;
  int interpType;

  /* smoothing options */
  int smoothType;
  int smoothNumLevels;
  int smoothNumSweeps;

  /* Schwarz smoothing options */
  int schwarzVariant;
  int schwarzOverlap;
  double schwarzDomainType;
  double schwarzRelaxWeight;

  /* ParaSAILS smoothing options */
  int psSymmetry;
  int psNumLevels;
  int psThreshold;
  double psFilter;

  /* PILUT smoothing options */
  double pilutDropTol;
  int pilutMaxNZRow;

  /* coarsening options */
  int coarsenNumAggLevels;
  int coarsenAggNumPaths;

  /* GSMG and LS smoothing options */
  int useGSMG;
  int numSamples;

  /* Output options */
  int boomerDebug;
  int printLevel;
  int logging;

} boomerAMGOptions;

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */


int initializeBoomerAMG(HYPRE_Solver *hypre_pc,
			boomerAMGOptions *options,
			const MPI_Comm mpi_comm);

void setBoomerAMGOptions(boomerAMGOptions *options);



#endif
