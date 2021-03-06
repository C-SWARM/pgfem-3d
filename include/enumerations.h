/* HEADER */
#ifndef PGFEM3D_ENUMERATIONS_H
#define PGFEM3D_ENUMERATIONS_H

enum {
  ELASTIC,
  TP_ELASTO_PLASTIC,
  FS_CRPL,
  FINITE_STRAIN,
  STABILIZED,
  MINI,
  MINI_3F,
  DISP,
  TF,
  CM,
  CM3F,
  ANALYSIS_MAX
};

constexpr const char* ANALYSIS_OPTS[ANALYSIS_MAX] = {
  /*[ELASTIC]           =*/ "ELASTIC",
  /*[TP_ELASTO_PLASTIC] =*/ "TP_ELASTO_PLASTIC",
  /*[FS_CRPL]           =*/ "FS_CRPL",
  /*[FINITE_STRAIN]     =*/ "FINITE_STRAIN",
  /*[STABILIZED]        =*/ "STABILIZED",
  /*[MINI]              =*/ "MINI",
  /*[MINI_3F]           =*/ "MINI_3F",
  /*[DISP]              =*/ "DISP",
  /*[TF]                =*/ "TF",
  /*[CM]                =*/ "CM",
  /*[CM3F]              =*/ "CM3F"
};

enum {
  PRECOND_PARA_SAILS,
  PRECOND_PILUT,
  PRECOND_EUCLID,
  PRECOND_BOOMER,
  PRECOND_NONE,
  PRECOND_DIAG_SCALE,
  PRECOND_JACOBI,
  PRECOND_IFPACK2,
  PRECOND_MUELU,
  PRECOND_MAX
};

constexpr const char* PRECOND_OPTS[PRECOND_MAX] = {
  /*[PRECOND_PARA_SAILS] =*/ "PARA_SAILS",
  /*[PRECOND_PILUT]      =*/ "PILUT",
  /*[PRECOND_EUCLID]     =*/ "EUCLID",
  /*[PRECOND_BOOMER]     =*/ "BOOMER",
  /*[PRECOND_NONE]       =*/ "NONE",
  /*[PRECOND_DIAG_SCALE] =*/ "DIAG_SCALE",
  /*[PRECOND_JACOBI]     =*/ "JACOBI",
  /*[PRECOND_IFPACK2]    =*/ "IFPACK2",
  /*[PRECOND_MUELU]      =*/ "MUELU"
};

enum {
  NETWORK_ISIR,
  NETWORK_PWC,
  NETWORK_ENV,
  NETWORK_MAX
};

constexpr const char* NETWORK_OPTS[NETWORK_MAX] = {
  /*[NETWORK_ISIR] =*/ "ISIR",
  /*[NETWORK_PWC] =*/ "PWC",
  /*[NETWORK_ENV] =*/ "ENV"
};

enum {
  BLOCKSOLVE,
  HYPRE,
  MTL,
  TRILINOS,
  SOLVER_PACKAGE_MAX
};

constexpr const char* SOLVER_PACKAGE_OPTS[SOLVER_PACKAGE_MAX] = {
  /*[BLOCKSOLVE] =*/ "BLOCKSOLVE",
  /*[HYPRE] =*/ "HYPRE",
  /*[MTL] =*/ "MTL",
  /*[TRILINOS]=*/ "TRILINOS"
};

enum {
  SOLVER_GMRES,
  SOLVER_BCG_STAB,
  SOLVER_AMG,
  SOLVER_FLEX,
  SOLVER_HYBRID,
  SOLVER_MAX
};

constexpr const char* SOLVER_OPTS[SOLVER_MAX] = {
  /*[SOLVER_GMRES]    =*/ "GMRES",
  /*[SOLVER_BCG_STAB] =*/ "BiCGSTAB",
  /*[SOLVER_AMG]      =*/ "BoomerAMG",
  /*[SOLVER_FLEX]     =*/ "FlexGMRES",
  /*[SOLVER_HYBRID]   =*/ "Hybrid (GMRES)"
};

enum {
  VIS_ELIXIR,
  VIS_ENSIGHT,
  VIS_VTK,
  VIS_NONE,
  VIS_OTHER,
  VIS_MAX
};

constexpr const char* VIS_OPTS[VIS_MAX] = {
  /*[VIS_ELIXIR]  =*/ "ELIXIR",
  /*[VIS_ENSIGHT] =*/ "ENSIGHT",
  /*[VIS_VTK]     =*/ "VTK",
  /*[VIS_NONE]    =*/ "NONE",
  /*[VIS_OTHER]   =*/ "OTHER"
};

enum {
  LINE_SEARCH,
  ADAPTIVE_TIME_STEPPING,
  CVG_CHECK_ON_ENERGY_NORM,
  NO_SUBDIVISION_LIMITS,
  SOLUTION_RECOVERY,
  SOLUTION_SCHEME_OPT_NO
};

constexpr const char* SOLUTION_SCHEME_OPTS[SOLUTION_SCHEME_OPT_NO] = {
  /*[LINE_SEARCH]              =*/ "LINE_SEARCH",
  /*[ADAPTIVE_TIME_STEPPING]   =*/ "ADAPTIVE_TIME_STEPPING",
  /*[CVG_CHECK_ON_ENERGY_NORM] =*/ "CVG_CHECK_ON_ENERGY_NORM"
};

#endif // #define PGFEM3D_ENUMERATIONS_H
