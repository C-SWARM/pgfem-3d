/* HEADER */
/**
 * AUTHORS:
 * Karel Matous, UNiversity of Notre Dame, kmatous [at] nd.edu
 */

#pragma once
#ifndef MATICE_H
#define MATICE_H

#include "PGFEM_mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

struct VAL_W_KEY {
  int val;
  int key;
};

typedef struct VAL_W_KEY val_key;

/** sort structures by value keeping key in tact */
int compare_val_w_key(const void *a,
		      const void *b);

int compare_long (const void *a,
		  const void *b);

int compare_int (const void *a,
		 const void *b);

long round1 (double a);

void nas_AB (double **A,
	     double **B,
	     double **C,
	     long m,
	     long n,
	     long p);

void nas_ATB (double **AT,
	      double **B,
	      double **C,
	      long m,
	      long n,
	      long p);

void nas_ABT (double **AT,
	      double **B,
	      double **C,
	      long m,
	      long n,
	      long p);

void inv_I(double **A,
	   double **I,
	   long m);

void vvplus (double *a,
	     double *b,
	     long n);

void vvminus (double *a,
	      double *b,
	      long n);

void mtv (double *a,
	  double *b,
	  double *c,
	  long m,
	  long n);

void mtvc (double *a,
	   double *b,
	   double *c,
	   long m,
	   long n);

void mv (double *a,
	 double *b,
	 double *c,
	 long m,
	 long n);

void mvc (double *a,
	  double *b,
	  double *c,
	  long m,
	  long n);

void mtmccr (double *a,
	     double *b,
	     double *c,
	     long m,
	     long n,
	     long p);

void copyi (long *a,
	    long *b,
	    long n);

void copyd (double *a,
	    double *b,
	    long n);

void normalization (double *a,
		    double nor,
		    long n);

void mv_sky (double *a,
	     double *b,
	     double *c,
	     long *adr,
	     long n);

void utv (double *a,
	  double *b,
	  double *c,
	  long *adr,
	  long n);

void ltv (double *a,
	  double *b,
	  double *c,
	  long *adr,
	  long n);

void mmt (double *a,
	  double *b,
	  double *c,
	  long ra,
	  long ca,
	  long rb);

void mm (double *a,
	 double *b,
	 double *c,
	 long l,
	 long m,
	 long n);

void mtm (double *a,
	  double *b,
	  double *c,
	  long ra,
	  long ca,
	  long cb);

double ss (const double *a,
	   const double *b,
	   const long n);

void nor_vec (double *BS_a,
	      long n,
	      MPI_Comm mpi_comm);

void nor_vec_serial (double *a,
		     long n,
		     int myrank);

void per_sym (double ***e);

void cross_product (double *a,
		    double *b,
		    double *c);

void tred2 (double **a,
	    int n,
	    double *d,
	    double *e);

void tqli (double *d,
	   double *e,
	   int n,
	   double **z);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
