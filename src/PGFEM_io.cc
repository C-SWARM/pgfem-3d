#include "PGFEM_io.h"
#include <string.h>
#include <stdarg.h>
#include "PGFEM_mpi.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

/* FILE pointers for macro and micro output */
FILE *PGFEM_macro_stdout = NULL;
FILE *PGFEM_macro_stderr = NULL;
FILE *PGFEM_micro_stdout = NULL;
FILE *PGFEM_micro_stderr = NULL;

/* main stdout/err file pointers */
FILE *PGFEM_stdout = NULL;
FILE *PGFEM_stderr = NULL;
FILE *PGFEM_null = NULL;

static int initialized = 0;

/* initialize PGFEM_stdout/err if PGFEM_initialize_io has not been
   called */
static inline int init_null()
{
  if(!initialized) return PGFEM_initialize_io(NULL,NULL);
  else return 0;
}

/** allocate enough space for the  filename and populate */
static void create_filename(const char *base,
                const char *ext,
                char **dest)
{
  const int len = strlen(base) + strlen(ext) + 1 /* for terminating char */;
  *dest = PGFEM_calloc(char, len);
  sprintf(*dest,"%s%s",base,ext);
}

/* points the PGFEM_* streams to out and err respectively */
static inline int redirect_output(FILE *out,
                  FILE *err)
{

  if(out == NULL) PGFEM_stdout = stdout; else PGFEM_stdout = out;
  if(err == NULL) PGFEM_stderr = stderr; else PGFEM_stderr = err;
  return 0;
}

int PGFEM_initialize_io(const char *macro_filename,
            const char *micro_filename)
{
  int err = 0;
  initialized = 1;
  const char mode[] = "a";

  if(PGFEM_null == NULL){
    PGFEM_null = PGFEM_fopen("/dev/null",mode);
  }

  /* set the macroscale output filenames */
  if(macro_filename != NULL){
    char *macro_stdout = NULL;
    char *macro_stderr = NULL;
    create_filename(macro_filename,".stdout.log",&macro_stdout);
    create_filename(macro_filename,".stderr.log",&macro_stderr);
    PGFEM_macro_stdout = PGFEM_fopen(macro_stdout,mode);
    PGFEM_macro_stderr = PGFEM_fopen(macro_stderr,mode);
    setvbuf(PGFEM_macro_stdout,NULL,_IOLBF,0);
    setvbuf(PGFEM_macro_stderr,NULL,_IONBF,0);
    free(macro_stdout);
    free(macro_stderr);
  }

  /* set the microscale output filenames */
  if(micro_filename != NULL){
    char *micro_stdout = NULL;
    char *micro_stderr = NULL;
    create_filename(micro_filename,".stdout.log",&micro_stdout);
    create_filename(micro_filename,".stderr.log",&micro_stderr);
    PGFEM_micro_stdout = PGFEM_fopen(micro_stdout,mode);
    PGFEM_micro_stderr = PGFEM_fopen(micro_stderr,mode);
    setvbuf(PGFEM_micro_stdout,NULL,_IOLBF,0);
    setvbuf(PGFEM_micro_stderr,NULL,_IONBF,0);
    free(micro_stdout);
    free(micro_stderr);
  }

  /* redirect output for macroscale */
  err += PGFEM_redirect_io_macro();

  return err;
}

int PGFEM_finalize_io()
{
  /* if the files do not point to stdout/stderr, free them */
  if(PGFEM_macro_stdout != stdout
     && PGFEM_macro_stdout != NULL){
    fclose(PGFEM_macro_stdout);
    PGFEM_macro_stdout = NULL;
  }

  if(PGFEM_macro_stderr != stderr
     && PGFEM_macro_stderr != NULL){
    fclose(PGFEM_macro_stderr);
    PGFEM_macro_stderr = NULL;
  }

  if(PGFEM_micro_stdout != stdout
     && PGFEM_micro_stdout != NULL){
    fclose(PGFEM_micro_stdout);
    PGFEM_micro_stdout = NULL;
  }

  if(PGFEM_micro_stderr != stderr
     && PGFEM_micro_stderr != NULL){
    fclose(PGFEM_micro_stderr);
    PGFEM_micro_stderr = NULL;
  }

  if(PGFEM_null != NULL){
    fclose(PGFEM_null);
    PGFEM_null = NULL;
  }

  /* reset the output to macro (console) */
  return PGFEM_redirect_io_macro();
}

int PGFEM_redirect_io_micro()
{
  return redirect_output(PGFEM_micro_stdout,PGFEM_micro_stderr);
}

int PGFEM_redirect_io_macro()
{
  return redirect_output(PGFEM_macro_stdout,PGFEM_macro_stderr);
}

int PGFEM_redirect_io_null()
{
  return redirect_output(PGFEM_null,PGFEM_null);
}
int PGFEM_fprintf(FILE * stream, const char * format,...)
{
  init_null();
  int result = 0;
  va_list args;
  va_start(args,format);
  result += vfprintf(stream,format,args);
  va_end(args);
  return result;
}

int PGFEM_printf(const char *format,...)
{
  init_null();
  int result = 0;
  va_list args;
  va_start(args,format);
  result += vfprintf(PGFEM_stdout,format,args);
  va_end(args);
  return result;
}

int PGFEM_printerr(const char *format,...)
{
  init_null();
  int result = 0;
  va_list args;
  va_start(args,format);
  result += vfprintf(PGFEM_stderr,format,args);
  va_end(args);
  return result;
}

/** This function is typically called from a macro to pass function,
    file and line information from invocation point. */
FILE* PGFEM_FOPEN(const char *filename,
          const char *mode,
          const char *function,
          const char *file,
          const long line)
{
  FILE *stream = fopen(filename,mode);
  if(stream == NULL){
    int err_rank = 0;
    PGFEM_Error_rank(&err_rank);
    PGFEM_printerr("[%d]ERROR: cannot open file (%s)! %s:%s:%ld\n",
           err_rank,filename,function,file,line);
    PGFEM_Abort();
  }
  return stream;
}
