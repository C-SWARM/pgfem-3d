/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#pragma once
#ifndef LB_TRANSFER_H
#define LB_TRANSFER_H

#include "load.h"
#include <stdlib.h>
#include <stdio.h>

struct PARTITION_LIST;
#ifndef TYPEDEF_PARTITION_LIST
#define TYPEDEF_PARTITION_LIST
typedef struct PARTITION_LIST PARTITION_LIST;
#endif

struct TRANSFER{
  size_t to;
  size_t from;
  LOAD_ID id;
};
#ifndef TYPEDEF_TRANSFER
#define TYPEDEF_TRANSFER
typedef struct TRANSFER TRANSFER;
#endif

/**
 * Copy _single_ TRANSFER object from src to dest.
 *
 * Does not perform copy if src == dest.
 */
void TRANSFER_copy(TRANSFER *dest,
		   const TRANSFER *src);

void TRANSFER_print(FILE *out,
		    const TRANSFER *T);

struct TRANSFER_LIST{
  size_t size;
  size_t max_size;
  TRANSFER *transfers;
};
#ifndef TYPEDEF_TRANSFER_LIST
#define TYPEDEF_TRANSFER_LIST
typedef struct TRANSFER_LIST TRANSFER_LIST;
#endif

void TRANSFER_LIST_build(TRANSFER_LIST *TL,
			 const size_t max_size);

void TRANSFER_LIST_destroy(TRANSFER_LIST *TL);

/**
 * Compute the list of required loads to transfer.
 *
 * Calls TRANSFER_LIST_destroy/build internally to create fresh list.
 */
void TRANSFER_LIST_compute(const PARTITION_LIST *PL,
			   TRANSFER_LIST *TL);
/**
 * Append a TRANSFER object to the list.
 *
 * Returns non-zero if TRANSFER_LIST is full and does not append
 * TRANSFER.
 */
int TRANSFER_LIST_push(TRANSFER_LIST *TL,
		       const TRANSFER *T);

/**
 * T/F (1/0) if the list is full.
 */
int TRANSFER_LIST_is_full(const TRANSFER_LIST *TL);

/**
 * T/F (1/0) if the list is empty.
 */
int TRANSFER_LIST_is_empty(const TRANSFER_LIST *TL);

/**
 * Print the TRANSFER_LIST object.
 */
void TRANSFER_LIST_print(FILE *out,
			 const TRANSFER_LIST *TL);

#endif
