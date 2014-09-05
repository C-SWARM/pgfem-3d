/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#pragma once
#ifndef LOAD_H
#define LOAD_H

#include <stdlib.h>
#include <stdio.h>

/**
 * Internal ideintifying structure for a LOAD object.
 */
struct LOAD_ID{
  size_t proc;
  size_t elem;
};
#ifndef TYPEDEF_LOAD_ID
#define TYPEDEF_LOAD_ID
typedef struct LOAD_ID LOAD_ID;
#endif

void LOAD_ID_copy(LOAD_ID *dest,
		  const LOAD_ID *src);

/**
 * Compare two LOAD_ID objects.
 * 
 * Return -1 if a < b, 0 if a == b, 1 if a > b. Suitable for use in
 * sorting/searching functions.
 */
int LOAD_ID_compare(const void *a,
		    const void *b);

/**
 * Print a _single_ LOAD_ID object to a file.
 */
void LOAD_ID_print(FILE *out,
		   const LOAD_ID *id);

struct LOAD{
  double load;
  size_t part_id;
  LOAD_ID id;
};
#ifndef TYPEDEF_LOAD
#define TYPEDEF_LOAD
typedef struct LOAD LOAD;
#endif

/**
 * Copy _single_ LOAD object from src to dest. Does nothing if src ==
 * dest;
 */ 
void LOAD_copy(LOAD *dest,
	       const LOAD *src);


/**
 * Swap _single_ LOAD object A with _single_ LOAD object B. Does nothing if A == B.
 */
void LOAD_swap(LOAD *A,
	       LOAD *B);

/**
 * See compare_LOAD_ID, but called on LOAD object.
 */
int LOAD_compare_id(const void *a,
		    const void *b);

/**
 * Compare the internal load value of two LOAD objects.
 *
 * Return -1 if a < b, 0 if a == b, 1 if a > b. Suitable for use in
 * sorting/searching functions.
 */
int LOAD_compare_load(const void *a,
		      const void *b);

/**
 * See LOAD_compare_load, but gives reverse ordering.
 */
int LOAD_compare_r_load(const void *a,
			const void *b);

/**
 * Compare the part_id of two LOAD objects. See also size_t_comp.
 */
int LOAD_compare_part_id(const void *a,
			 const void *b);

/**
 * Compare by part_id, then by load.
 */
int LOAD_compare_part_id_load(const void *a,
			      const void *b);
/**
 * Determine whether two LOAD objects are equal within some absolute
 * tolerance (tol >= 0).
 *
 * Return |a - b| >= tol.
 */
int LOAD_approx_equal(const LOAD *a,
		      const LOAD *b,
		      const double tol);

/**
 * Print a _single_ LOAD object to a file.
 */
void LOAD_print(FILE *out,
		const LOAD *load);
#endif
