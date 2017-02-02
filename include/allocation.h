/* HEADER */
#pragma once
#ifndef _PGFEM_ALLOCATION_H_
#define _PGFEM_ALLOCATION_H_
#define ALLOCATION_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Note this function is called typically called through the macro
      PGFEM_calloc so that __func__ __FILE__ and __LINE__ can be passed
      from the point of invocation without as much code */
  void* PGFEM_CALLOC(const long nelem,
		     const long size,
		     const char *function,
		     const char *file,
		     const long line);

  void** PGFEM_CALLOC_2(const long nelem1,
			const long nelem2,
			const long size,
			const char *function,
			const char *file,
			const long line);

  void*** PGFEM_CALLOC_3(const long nelem1,
			 const long nelem2,
			 const long nelem3,
			 const long size,
			 const char *function,
			 const char *file,
			 const long line);

  void**** PGFEM_CALLOC_4(const long nelem1,
			  const long nelem2,
			  const long nelem3,
			  const long nelem4,
			  const long size,
			  const char *function,
			  const char *file,
			  const long line);

  void PGFEM_free2(void **a,
		   const long nelem);

  void PGFEM_free3(void ***a,
		   const long nelem1,
		   const long nelem2);

  void PGFEM_free4(void ****a,
		   const long nelem1,
		   const long nelem2,
		   const long nelem3);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

/* define macro wrappers to pass in function/file/line info */
#define PGFEM_calloc(nelem,size) \
  PGFEM_CALLOC(nelem,size,__func__,__FILE__,__LINE__)
#define PGFEM_calloc2(nelem1,nelem2,size) \
  PGFEM_CALLOC_2(nelem1,nelem2,size,__func__,__FILE__,__LINE__)
#define PGFEM_calloc3(nelem1,nelem2,nelem3,size) \
  PGFEM_CALLOC_3(nelem1,nelem2,nelem3,size,__func__,__FILE__,__LINE__)
#define PGFEM_calloc4(nelem1,nelem2,nelem3,nelem4,size) \
  PGFEM_CALLOC_4(nelem1,nelem2,nelem3,nelem4,size,__func__,__FILE__,__LINE__)
#define PGFEM_free(a) free(a)

/* define macro wrappers to replace incl.c allocation/deallocation
   functions */
#ifndef PGFEM_MACRO_ALLOCATION
#define PGFEM_MACRO_ALLOCATION
#define aloc1i(m) PGFEM_calloc(m,sizeof(int))
#define aloc1l(m) PGFEM_calloc(m,sizeof(long))
#define aloc2l(m,n) (long**) PGFEM_calloc2(m,n,sizeof(long))
#define aloc3l(m,n,p) (long***) PGFEM_calloc3(m,n,p,sizeof(long))
#define aloc1(m) PGFEM_calloc(m,sizeof(double))
#define aloc2(m,n) (double**) PGFEM_calloc2(m,n,sizeof(double))
#define aloc3(m,n,p) (double***) PGFEM_calloc3(m,n,p,sizeof(double))
#define aloc4(m,n,p,r) (double****) PGFEM_calloc4(m,n,p,r,sizeof(double))
#define dealoc1i(a) PGFEM_free(a)
#define dealoc1l(a) PGFEM_free(a)
#define dealoc2l(a,m) PGFEM_free2((void**) a,m)
#define dealoc3l(a,m,n) PGFEM_free3((void***) a,m,n)
#define dealoc1(a) PGFEM_free(a)
#define dealoc2(a,m) PGFEM_free2((void**) a,m)
#define dealoc3(a,m,n) PGFEM_free3((void***) a,m,n)
#define dealoc4(a,m,n,p) PGFEM_free4((void****) a,m,n,p)
#endif /* #ifndef PGFEM_MACRO_ALLOCATION */

#endif /* #ifndef _PGFEM_ALLOCATION_H_ */
