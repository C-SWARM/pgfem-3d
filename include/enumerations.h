/* HEADER */

#pragma once
#ifndef ENUMERATIONS_H
#define ENUMERATIONS_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  enum {ELASTIC=0,
	TP_ELASTO_PLASTIC=1,
	FS_CRPL=2,
	FINITE_STRAIN=3,
	STABILIZED=4,
	MINI=5,
	MINI_3F=6,
	DISP=7,
	ANALYSIS_MAX} ANALYSIS_TYPE;

  enum {PARA_SAILS=0,
	PILUT=1,
	EUCLID=2,
	BOOMER=3,
	NONE=4,
	DIAG_SCALE,
	JACOBI} PRECOND_TYPE;

  enum {BLOCKSOLVE=0,
	HYPRE=1} SOLVER_PACKAGE;

  enum {HYPRE_GMRES=0,
	HYPRE_BCG_STAB=1,
	HYPRE_AMG=2,
	HYPRE_FLEX,
	HYPRE_HYBRID} SOLVER_TYPE;

  enum {VIS_ELIXIR,VIS_ENSIGHT,VIS_VTK,VIS_NONE,VIS_OTHER} VIS_TYPE;

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
