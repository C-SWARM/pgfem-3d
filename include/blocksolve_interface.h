/* HEADER */

/**
 * @file This is a temporary file to remove some dependency on
 * BlockSolve95.
 */

#pragma once
#ifndef BLOCKSOLVE_INTERFACE_H
#define BLOCKSOLVE_INTERFACE_H

#include "BSprivate.h"

/** Write a BlockSolve95 matrix to a file in MATLAB sparse format */
void write_mat_matlab(char *str,
		      BSspmat *A,
		      BSprocinfo  *procinfo);

/** Set all entries in a BlockSolve95 matrix to 0.0 */
void null_BSspmat (BSspmat *K);

/**
 * Allocate a Blocksolve Matrix.
 */
BSspmat *BSalloc_A (int start_num,
		    int n,
		    int *rp,
		    int *cval,
		    BSprocinfo *procinfo);
#endif
