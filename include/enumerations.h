/* HEADER */

#pragma once
#ifndef ENUMERATIONS_H
#define ENUMERATIONS_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

typedef  enum {ELASTIC=0,
	TP_ELASTO_PLASTIC=1,
	FS_CRPL=2,
	FINITE_STRAIN=3,
	STABILIZED=4,
	MINI=5,
	MINI_3F=6,
	DISP=7,
  TF=8,
  CM=9,
	ANALYSIS_MAX} ANALYSIS_TYPE;

typedef  enum {PARA_SAILS=0,
	PILUT=1,
	EUCLID=2,
	BOOMER=3,
	NONE=4,
	DIAG_SCALE,
	JACOBI} PRECOND_TYPE;

typedef  enum {BLOCKSOLVE=0,
	HYPRE=1} SOLVER_PACKAGE;

typedef  enum {HYPRE_GMRES=0,
	HYPRE_BCG_STAB=1,
	HYPRE_AMG=2,
	HYPRE_FLEX,
	HYPRE_HYBRID} SOLVER_TYPE;

typedef  enum {VIS_ELIXIR,VIS_ENSIGHT,VIS_VTK,VIS_NONE,VIS_OTHER} VIS_TYPE;

typedef  enum {LINE_SEARCH, 
        ADAPTIVE_TIME_STEPPING,
        CVG_CHECK_ON_ENERGY_NORM, 
        SOLUTION_SCHEME_OPT_NO} SOLUTION_SCHEME_OPT;

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
