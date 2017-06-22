/* HEADER */
#pragma once
#ifndef INCL_H
#define INCL_H

#include "data_structure.h"
#include "element.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "allocation.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

void build_elem_inelas (long ne,
			ELEMENT *elem);

void build_pressure_nodes (long ne,
			   long npres,
			   ELEMENT *elem,
			   SIG *sig,
			   EPS *eps,
			   const int analysis);

void build_crystal_plast (long ne,
			  ELEMENT *elem,
			  SIG *sig,
			  EPS *eps,
			  CRPL *crpl,
			  const int analysis,
			  const int plc);

void nulld (double *a,
	    long n);

void nulld2 (double **a,
	     long m,
	     long n);

/*=== Replaced in allocation.h ===*/
#ifndef PGFEM_MACRO_ALLOCATION
int* aloc1i (int m);

void dealoc1i (int *a);

long* aloc1l (long m);

void dealoc1l (long *a);

long** aloc2l (long m,long n);

void dealoc2l (long **a,long m);

long*** aloc3l (long m,long n,long p);

void dealoc3l (long ***a,long m,long n);

double* aloc1 (long m);

void dealoc1 (double *a);

double** aloc2 (long m,long n);

void dealoc2 (double **a,long m);

double*** aloc3 (long m,long n,long p);

void dealoc3 (double ***a,long m,long n);

double**** aloc4 (long m,long n,long p,long q);

void dealoc4 (double ****a,long m,long n,long p);
#endif /* #ifndef PGFEM_MACRO_ALLOCATION */

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */


#endif
