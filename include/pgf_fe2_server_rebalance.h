/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#pragma once
#ifndef PGF_FE2_SERVER_REBALANCE_H
#define PGF_FE2_SERVER_REBALANCE_H

#include <stdlib.h>
#include "PGFEM_mpi.h"
#include "microscale_information.h"

/**
 * Structure for handling information about how a server is to be
 * rebalanced. Describes what jobs are to be kept, what jobs should be
 * sent (and to whom), and what jobs are to be received (and from
 * whom). The total number of jobs on the server after the rebalancing
 * is n_keep + n_recv.
 */
/* typedef int* pgf_FE2_server_rebalance; */
struct pgf_FE2_server_rebalance;
typedef struct pgf_FE2_server_rebalance pgf_FE2_server_rebalance;


/**
 * Return the size fo the rebalance object in bytes.
 */
size_t pgf_FE2_server_rebalance_n_bytes(const pgf_FE2_server_rebalance *t);

/**
 * Return a pointer to the internal bufeer for communication purposes.
 */
void* pgf_FE2_server_rebalance_buff(const pgf_FE2_server_rebalance *t);

/**
 * Build the data structure and set metadata.
 */
void pgf_FE2_server_rebalance_build(pgf_FE2_server_rebalance **t,
				    const size_t n_keep,
				    const size_t n_send,
				    const size_t n_recv);

/**
 * Initialize a pgf_FE2_server_rebalance from a buffer. Invalidates
 * the handle 'buffer' on exit. Any memory pointed to by buffer is
 * destroyed on a call to ..._destroy(t);
 */
void pgf_FE2_server_rebalance_build_from_buffer(pgf_FE2_server_rebalance **t,
						void **buffer);
/**
 * Destroy the data strucuture.
 */
void pgf_FE2_server_rebalance_destroy(pgf_FE2_server_rebalance *t);

/**
 * Get the number of jobs to keep.
 */
int pgf_FE2_server_rebalance_n_keep(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer containing job ids to keep.
 */
int* pgf_FE2_server_rebalance_keep_buf(const pgf_FE2_server_rebalance *t);

/**
 * Get the number of jobs to send.
 */
int pgf_FE2_server_rebalance_n_send(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer containing job ids to send.
 */
int* pgf_FE2_server_rebalance_send_buf(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer contining destinations of
 * communication.
 */
int* pgf_FE2_server_rebalance_send_dest(const pgf_FE2_server_rebalance *t);

/**
 * Get the number of jobs to receive.
 */
int pgf_FE2_server_rebalance_n_recv(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer containing job ids to recv.
 */
int* pgf_FE2_server_rebalance_recv_buf(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer contining sources of
 * communication.
 */
int* pgf_FE2_server_rebalance_recv_src(const pgf_FE2_server_rebalance *t);

/**
 * Posts non-blocking communication for exchanging all of the
 * jobs. Allocates buffers as necisarry for maintaining consistant
 * send information and assigns locations for microstructures to
 * receive.
 */
int pgf_FE2_server_rebalance_post_exchange(const pgf_FE2_server_rebalance *t,
					   const PGFEM_mpi_comm *mpi_comm,
					   MICROSCALE *micro);
/**
 * Finalizes exchange communications and releases internal
 * buffers. Need to extend this functionality to overlab with further
 * computation.
 */
int pgf_FE2_server_rebalance_finalize_exchange(const PGFEM_mpi_comm *mpi_comm);

#endif
