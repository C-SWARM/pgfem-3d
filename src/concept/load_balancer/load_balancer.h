/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef LOAD_BALANCER_H
#define LOAD_BALANCER_H

/* pre-declare PARTITION_LIST and typedef */
struct PARTITION_LIST;
#ifndef TYPEDEF_PARTITION_LIST
#define TYPEDEF_PARTITION_LIST
typedef struct PARTITION_LIST PARTITION_LIST;
#endif

/**
 * Re-balance the PARTITION_LIST using the greedy heuristic.
 *
 * Given an initial partition list PL, generate the S = union of P_i
 * in PL. Re-assign the elements of S to P_i in PL according to the
 * greedy heuristic. Returns 0 on success. May call abort() on
 * internal error.
 */
int load_balancer_greedy(PARTITION_LIST *PL);

#endif
