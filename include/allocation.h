/* HEADER */
#pragma once
#ifndef _PGFEM_ALLOCATION_H_
#define _PGFEM_ALLOCATION_H_
#define ALLOCATION_H

#include "pgfem3d/Communication.hpp"
#include <stdlib.h>

// Help with legacy malloc.
template <class T>
static inline T* PGFEM_malloc(size_t n = 1) {
  return static_cast<T*>(malloc(n * sizeof(T)));
}

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

void* PGFEM_CALLOC_PIN(const long nelem,
    const long size,
    multiscale::net::Network *net,
    multiscale::net::Key *key,
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

void PGFEM_free_unpin(void *a, long bytes, multiscale::net::Network *net);

/* define macro wrappers to pass in function/file/line info */
/* define macro wrappers to pass in function/file/line info */
#define PGFEM_calloc(T, nelem)                      \
  static_cast<T*>(PGFEM_CALLOC((nelem), sizeof(T),  \
          __func__, __FILE__, __LINE__))

#define PGFEM_calloc2(T, nelem1, nelem2)                                \
  reinterpret_cast<T**>(PGFEM_CALLOC_2((nelem1), (nelem2), sizeof(T),   \
          __func__, __FILE__, __LINE__))

#define PGFEM_calloc3(T, nelem1, nelem2, nelem3)                        \
  reinterpret_cast<T***>(PGFEM_CALLOC_3((nelem1), (nelem2), (nelem3),   \
          sizeof(T), __func__, __FILE__, __LINE__))

#define PGFEM_calloc4(T, nelem1, nelem2, nelem3, nelem4)                \
  reinterpret_cast<T****>(PGFEM_CALLOC_4((nelem1), (nelem2), (nelem3),  \
          (nelem4), sizeof(T), __func__, __FILE__, __LINE__))

#define PGFEM_calloc_pin(T, nelem, net, key)                            \
  static_cast<T*>(PGFEM_CALLOC_PIN((nelem), sizeof(T), net, key,	\
          __func__, __FILE__, __LINE__))

#define PGFEM_free(a) free(a)

#define aloc1i(m)      PGFEM_calloc(int, (m))
#define aloc1l(m)      PGFEM_calloc(long, (m))
#define aloc2l(m,n)    PGFEM_calloc2(long, (m), (n))
#define aloc3l(m,n,p)  PGFEM_calloc3(long, (m), (n), (p))
#define aloc1(m)       PGFEM_calloc(double, (m))
#define aloc2(m,n)     PGFEM_calloc2(double, (m), (n))
#define aloc3(m,n,p)   PGFEM_calloc3(double, (m), (n), (p))
#define aloc4(m,n,p,r) PGFEM_calloc4(double, (m), (n), (p), r)

#define dealoc1i(a)      PGFEM_free(a)
#define dealoc1l(a)      PGFEM_free(a)
#define dealoc2l(a,m)    PGFEM_free2(reinterpret_cast<void**>(a), (m))
#define dealoc3l(a,m,n)  PGFEM_free3(reinterpret_cast<void***>(a), (m), (n))
#define dealoc1(a)       PGFEM_free(a)
#define dealoc2(a,m)     PGFEM_free2(reinterpret_cast<void**>(a), (m))
#define dealoc3(a,m,n)   PGFEM_free3(reinterpret_cast<void***>(a), (m), (n))
#define dealoc4(a,m,n,p) PGFEM_free4(reinterpret_cast<void****>(a), (m), (n), (p))

#endif /* #ifndef _PGFEM_ALLOCATION_H_ */
