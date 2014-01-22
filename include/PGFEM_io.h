/* HEADER */
/* This is a short wrapper to handle Micro-Macro output to stdout and
   stderr. Other GENERAL I/O functions may be put here as well. */
#ifndef PGFEM_IO_H
#define PGFEM_IO_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef INCLUDED_STDIO_H
#include <stdio.h>
#define INCLUDED_STDIO_H
#endif

  /* FILE pointers that get switched about */
  extern FILE *PGFEM_stdout;
  extern FILE *PGFEM_stderr;
  extern FILE *PGFEM_null;

  /** Initializes ability to redirect stdout/stderr streams */
  int PGFEM_initialize_io(const char *macro_filename,
			  const char *micro_filename);

  /** Finalizes/removes ability to redirect stdout/stderr streams */
  int PGFEM_finalize_io();

  /** Redirect output for microscale */
  int PGFEM_redirect_io_micro();

  /** Redirect output for macroscale */
  int PGFEM_redirect_io_macro();

  /** Redirect to /dev/null */
  int PGFEM_redirect_io_null();

  /** Wrapper for fprintf */
  int PGFEM_fprintf(FILE *stream, const char *format,...);

  /** Wrapper for printf. Prints to PGFEM_stdout */
  int PGFEM_printf(const char *format,...);

  /** Print to PGFEM_stderr */
  int PGFEM_printerr(const char *format,...);

  /** This function is typically called from a macro to pass function,
      file and line information from invocation point. */
  FILE* PGFEM_FOPEN(const char *filename,
		    const char *mode,
		    const char *function,
		    const char *file,
		    const long line);


#define PGFEM_fopen(filename,mode)			\
  PGFEM_FOPEN(filename,mode,__func__,__FILE__,__LINE__)
#define PGFEM_fclose(a) fclose(a)

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */
#endif /* #ifndef  */

/* include block

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

 */
