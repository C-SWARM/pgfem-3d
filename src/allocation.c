/* HEADER */

/* This file defines the allocation functions which include more
   descriptive debugging messages. Note that long is used in lew of
   size_t so that signs may be checked for error reporting */
#include "allocation.h"
#include "PGFEM_mpi.h"
#include <string.h>
#include <unistd.h>

#include "PGFEM_io.h"

#ifdef NDEBUG
#define DEBUG_PGFEM_ALLOC 0
#define DEBUG_PGFEM_ALLOC_LOG 0
#else
#define DEBUG_PGFEM_ALLOC 1
#define DEBUG_PGFEM_ALLOC_LOG 0
#endif

static const size_t ptr_size = sizeof(void*);

static inline int get_err_rank(int *rank) { return PGFEM_Error_rank(rank); }

/** Write memory state (kilobytes) to debug file. THIS DRAMATICALLY
    REDUCES PERFORMACE and seems to be a poor measure of memory
    usage... */
static void log_mem_state_kb()
{
  static int called = 0;
  static char log_name[100];
  static char mem_name[100];
  if(!called){
    int myrank = 0; 
    get_err_rank(&myrank);
    const size_t pid = getpid();
    sprintf(log_name,"mem_%.5d.log",myrank);
    sprintf(mem_name,"/proc/%ld/statm",pid);
    called = 1;
  }
  FILE *log = fopen(log_name,"a"); /* append to log */
  if(log == NULL){
    PGFEM_printerr("WARNING: could not open (%s). No output written.\n",
		   log_name);
    return;
  } else {
    FILE *mem_state = PGFEM_fopen(mem_name,"r");
    if(mem_state == NULL){
      PGFEM_printerr("WARNING: could not open (%s)."
		     " No output written.\n",mem_name);
      PGFEM_fclose(log);
      return;
    }
    /* set the mem_state file to NO BUFFERING */
    setvbuf(mem_state,NULL,_IONBF,0);
    fseek(mem_state,0,SEEK_SET);
    int size = 0;
    int resident = 0;
    int shared = 0;
    fscanf(mem_state,"%d %d %d",&size,&resident,&shared);
    fclose(mem_state);
    PGFEM_fprintf(log,"%d %d\n",4*resident,4*size);
    PGFEM_fclose(log);
    return;
  }
}

/** Note this function is called typically called through the macro
    PGFEM_calloc so that __func__ __FILE__ and __LINE__ can be passed
    from the point of invocation without as much code */
void* PGFEM_CALLOC(const long nelem,
		   const long size,
		   const char *function,
		   const char *file,
		   const long line)
{
  int err = 0;
  int myrank = 0;
  void *val = NULL;
  if(DEBUG_PGFEM_ALLOC){
    if(nelem < 0){
      err = 1;
    } else if(size < 0){
      err = 2;
    } else if(nelem == 0){
      err = 3;
    }

    switch(err){
    case 0: default: val = calloc(nelem,size); break;
    case 3:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]Warning: attempted to allocate 0 elements. %s:%s:%ld\n",
		     myrank,function,file,line);
      val = calloc(1,size); break;
    case 1:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    case 2:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative sized elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    }
  } /* No debugging messages */
  else if (nelem == 0){
    val = calloc(1,size);
  } else {
    val = calloc(nelem,size);
  }

  /* Check that valid pointer was returned */
  if(val == NULL){
    PGFEM_printerr("[%d]ERROR: calloc returned NULL! %s:%s:%ld\n",
		   myrank,function,file,line);
    PGFEM_Abort();
  }

  if(DEBUG_PGFEM_ALLOC_LOG){
    log_mem_state_kb();
  }
  return val;
}

void** PGFEM_CALLOC_2(const long nelem1,
		      const long nelem2,
		      const long size,
		      const char *function,
		      const char *file,
		      const long line)
{
  int err = 0;
  int myrank = 0;
  void **val = NULL;
  long n = nelem1;
  long m = nelem2;
  if(nelem1 == 0) n = 1;
  if(nelem2 == 0) m = 1;

  if(DEBUG_PGFEM_ALLOC){
    if(nelem1 < 0 || nelem2 < 0){
      err = 1;
    } else if(size < 0){
      err = 2;
    } else if (nelem1 == 0 || nelem2 == 0){
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]Warning: attempted to allocate 0 elements. %s:%s:%ld\n",
		     myrank,function,file,line);
    }

    switch(err){
    case 0: default:
      val = PGFEM_CALLOC(n,ptr_size,function,file,line);
      for(long i=0; i<n; i++){
	val[i] = PGFEM_CALLOC(m,size,function,__FILE__,__LINE__);
      }
      break;
    case 1:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    case 2:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative sized elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    }
  } /* No debugging messages */
  else {
    val = PGFEM_CALLOC(n,ptr_size,function,file,line);
    for(long i=0; i<n; i++){
      val[i] = PGFEM_CALLOC(m,size,function,__FILE__,__LINE__);
    }
  }
  return val;
}

void*** PGFEM_CALLOC_3(const long nelem1,
		       const long nelem2,
		       const long nelem3,
		       const long size,
		       const char *function,
		       const char *file,
		       const long line)
{
  int err = 0;
  int myrank = 0;
  void ***val = NULL;
  long n = nelem1;
  long m = nelem2;
  long p = nelem3;
  if(nelem1 == 0) n = 1;
  if(nelem2 == 0) m = 1;
  if(nelem3 == 0) p = 1;

  if(DEBUG_PGFEM_ALLOC){
    if(nelem1 < 0 || nelem2 < 0 || nelem3 < 0){
      err = 1;
    } else if(size < 0){
      err = 2;
    } else if (nelem1 == 0 || nelem2 == 0 || nelem3 == 0){
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]Warning: attempted to allocate 0 elements. %s:%s:%ld\n",
		     myrank,function,file,line);
    }

    switch(err){
    case 0: default:
      val = PGFEM_CALLOC(n,ptr_size,function,file,line);
      for(long i=0; i<n; i++){
	val[i] = PGFEM_CALLOC(m,ptr_size,function,__FILE__,__LINE__);
	for(long j=0; j<m; j++){
	  val[i][j] = PGFEM_CALLOC(p,size,function,__FILE__,__LINE__);
	}
      }
      break;
    case 1:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    case 2:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative sized elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    }
  } /* No debugging messages */
  else {
    val = PGFEM_CALLOC(n,ptr_size,function,file,line);
    for(long i=0; i<n; i++){
      val[i] = PGFEM_CALLOC(m,ptr_size,function,__FILE__,__LINE__);
      for(long j=0; j<m; j++){
	val[i][j] = PGFEM_CALLOC(p,size,function,__FILE__,__LINE__);
      }
    }
  }

  return val;
}

void**** PGFEM_CALLOC_4(const long nelem1,
			const long nelem2,
			const long nelem3,
			const long nelem4,
			const long size,
			const char *function,
			const char *file,
			const long line)
{
  int err = 0;
  int myrank = 0;
  void ****val = NULL;
  long n = nelem1;
  long m = nelem2;
  long p = nelem3;
  long r = nelem4;
  if(nelem1 == 0) n = 1;
  if(nelem2 == 0) m = 1;
  if(nelem3 == 0) p = 1;
  if(nelem4 == 0) r = 1;

  if(DEBUG_PGFEM_ALLOC){
    if(nelem1 < 0 || nelem2 < 0 || nelem3 < 0 || nelem4 < 0){
      err = 1;
    } else if(size < 0){
      err = 2;
    } else if (nelem1 == 0 || nelem2 == 0
	       || nelem3 == 0 || nelem4 == 0){
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]Warning: attempted to allocate 0 elements. %s:%s:%ld\n",
		     myrank,function,file,line);
    }

    switch(err){
    case 0: default:
      val = PGFEM_CALLOC(n,ptr_size,function,file,line);
      for(long i=0; i<n; i++){
	val[i] = PGFEM_CALLOC(m,ptr_size,function,__FILE__,__LINE__);
	for(long j=0; j<m; j++){
	  val[i][j] = PGFEM_CALLOC(p,ptr_size,function,__FILE__,__LINE__);
	  for(long k=0; k<p; k++){
	    val[i][j][k] = PGFEM_CALLOC(r,size,function,__FILE__,__LINE__);
	  }
	}
      }
      break;
    case 1:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    case 2:
      get_err_rank(&myrank);
      PGFEM_printerr("[%d]ERROR: cannot allocate negative sized elements! %s:%s:%ld\n",
		     myrank,function,file,line);
      PGFEM_Abort();
      break;
    }
  } /* No debugging messages */
  else {
    val = PGFEM_CALLOC(n,ptr_size,function,file,line);
    for(long i=0; i<n; i++){
      val[i] = PGFEM_CALLOC(m,ptr_size,function,__FILE__,__LINE__);
      for(long j=0; j<m; j++){
	val[i][j] = PGFEM_CALLOC(p,ptr_size,function,__FILE__,__LINE__);
	for(long k=0; k<p; k++){
	  val[i][j][k] = PGFEM_CALLOC(r,size,function,__FILE__,__LINE__);
	}
      }
    }
  }

  return val;
}

void PGFEM_free2(void **a,
		 const long nelem)
{
  for(long i=0; i<nelem; i++){
    PGFEM_free(a[i]);
  }
  PGFEM_free(a);
}

void PGFEM_free3(void ***a,
		 const long nelem1,
		 const long nelem2)
{
  for(long i=0; i<nelem1; i++){
    for(long j=0; j<nelem2; j++){
      PGFEM_free(a[i][j]);
    }
    PGFEM_free(a[i]);
  }
  PGFEM_free(a);
}

void PGFEM_free4(void ****a,
		 const long nelem1,
		 const long nelem2,
		 const long nelem3)
{
  for(long i=0; i<nelem1; i++){
    for(long j=0; j<nelem2; j++){
      for(long k=0; k<nelem3; k++){
	PGFEM_free(a[i][j][k]);
      }
      PGFEM_free(a[i][j]);
    }
    PGFEM_free(a[i]);
  }
  PGFEM_free(a);
}
