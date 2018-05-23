/* HEADER */
#pragma once
#ifndef PGFEM_OPTIONS_H
#define PGFEM_OPTIONS_H

#include "enumerations.h"
#include <cstdio>

/** Structure containing option parameters to set from the command
    line. */
struct PGFem3D_opt {
  /* network option */
  int network;
  
  /* solver options */
  int solverpackage; /* HYPRE or BS95 (BS95 is obsolete) */
  int solver;
  int precond;
  int kdim;
  int maxit;
  int solution_scheme_opt[SOLUTION_SCHEME_OPT_NO]; // set numerical solutions scheme options
  // solution_scheme_opt[0]: for line search
  // solution_scheme_opt[1]: for adaptive time stepping
  // solution_scheme_opt[2]: check convergence on energy norm

  /* analysis options */
  int analysis_type;
  double stab;
  int cohesive;
  static constexpr const int gem = 0;           // generalized fem
  int plc;
  int multi_scale;
  int cm;

  /* visualization options */
  int vis_format; /* gr2 */
  int ascii;
  int smoothing; /* gr4 */

  /* other options */
  int periodic; /* gr5, OBSOLETE */
  int renumber; /* unused */
  int legacy;
  int debug;
  int me; /* model entity */
  int restart;
  int max_n_jobs;
  int no_migrate;
  int comp_print_reaction; // compute and print reaction forces
  int comp_print_macro;    // compute and print macro values (GF, GP, GS)
  int print_EXA;           // output information about the EXA metric
  bool comp_print_max_pressure; // compute and print maximum element pressure together with coordinate

  /* input overrides */
  int override_pre_disp;
  char *pre_disp_file;
  int override_solver_file;
  char *solver_file;
  char *override_material_props;

  /* I/O file names */
  const char *ipath;
  const char *opath;
  const char *ifname;
  const char *ofname;
  double walltime;
};

/** Set the default option flags. */
void set_default_options(PGFem3D_opt *options);

/** Print the PGFem3D_opt structure. */
void print_options(FILE *out, const PGFem3D_opt *options);

/** Print usage message outlining the available options. */
void print_usage(FILE *out);

/** Parse the command line and set option values */
void parse_command_line(const int argc,
                        char **argv,
                        const int myrank,
                        PGFem3D_opt *options);

/** parse the command line starting from a argv[start_idx] */
void re_parse_command_line(const int myrank,
                           const int start_idx,
                           const int argc,
                           char **argv,
                           PGFem3D_opt *options);

/** get start_idx and argc to pass to re_parse_command_line for
    micro and macro option blocks. Set the macro/micro sizes and
    return debug flag if found in any option block */
void get_macro_micro_option_blocks(int myrank,
                                   int argc,
                                   char *argv[],
                                   int *macro_start,
                                   int *macro_argc,
                                   int *micro_start,
                                   int *micro_argc,
                                   int *macro_nproc,
                                   int *micro_group_size,
                                   int *debug);

/** Print a summary of the solution method */
void print_interpreted_options(const PGFem3D_opt *opts);

#endif /* #ifndef PFEM_OPTIONS_H */
