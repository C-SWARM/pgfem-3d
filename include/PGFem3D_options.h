/* HEADER */
#pragma once
#ifndef PGFEM_OPTIONS_H
#define PGFEM_OPTIONS_H

#include "PGFEM_io.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure containing option parameters to set from the command
      line. */
  typedef struct PGFEM3D_OPTIONS{
    /* solver options */
    int solverpackage; /* HYPRE or BS95 (BS95 is obsolete) */
    int solver;
    int precond;
    int kdim;
    int maxit;

    /* analysis options */
    int analysis_type;
    double stab;
    int cohesive;
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

    /* input overrides */
    int override_pre_disp;
    char *pre_disp_file;
    int override_solver_file;
    char *solver_file;
    char *override_material_props;

    /* I/O file names */
    char *ipath;
    char *opath;
    char *ifname;
    char *ofname;
  } PGFem3D_opt;

  /** Set the default option flags. */
  void set_default_options(PGFem3D_opt *options);

  /** Print the PGFem3D_opt structure. */
  void print_options(FILE *out,
		     const PGFem3D_opt *options);

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
  void get_macro_micro_option_blocks(const int myrank,
				     const int argc,
				     /* const */ char **argv,
				     int *macro_start,
				     int *macro_argc,
				     int *micro_start,
				     int *micro_argc,
				     int *macro_nproc,
				     int *micro_group_size,
				     int *debug);

  /** Print a summary of the solution method */
  void print_interpreted_options(const PGFem3D_opt *opts);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PFEM_OPTIONS_H */
