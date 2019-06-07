/* HEADER */

/**
 * \file This is a short wrapper to handle Micro-Macro output to stdout and
 *  stderr. Other GENERAL I/O functions may be put here as well.
 *
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef PGFEM_IO_H
#define PGFEM_IO_H

#include <stdio.h>
#include <string>

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

#define PGFEM_fopen(filename,mode)                      \
  PGFEM_FOPEN(filename,mode,__func__,__FILE__,__LINE__)

#define PGFEM_fclose(a) fclose(a)

namespace pgfem3d {
class ScopedFile {
 public:
  ScopedFile(const char* fn, const char* mode) : fp_(PGFEM_fopen(fn, mode)) {
  }

  ~ScopedFile() {
    if (fp_) {
      PGFEM_fclose(fp_);
    }
  }

  operator bool() const {
    return fp_ != nullptr;
  }

  operator FILE*() {
    return fp_;
  }

 private:
  FILE *fp_;
};

static inline ScopedFile scoped_fopen(const char* fn, const char* mode) {
  return { fn, mode };
}

static inline ScopedFile scoped_fopen(std::string fn, const char* mode) {
  return { fn.c_str(), mode };
}
}
#endif /* #ifndef  */
