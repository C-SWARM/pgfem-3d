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

/**
 * Swap jobs that require communication for those that don't if they
 * take approximately the same time (within some _absolute_ tolerance,
 * tol > 0).
 *
 * Modifies the PARTITION_LIST, indices to internal data may become
 * invalid. Returns 0 on success. May call abort() on internal error.
 */
int load_balancer_reduce_comm(PARTITION_LIST *PL,
			      const double tol);
#endif
