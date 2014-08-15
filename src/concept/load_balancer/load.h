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
 * Print a _single_ LOAD object to a file.
 */
void LOAD_print(FILE *out,
		const LOAD *load);
#endif
