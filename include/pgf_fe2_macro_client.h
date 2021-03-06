/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */
#pragma once
#ifndef PGF_FE2_MACRO_CLIENT_H
#define PGF_FE2_MACRO_CLIENT_H

#include "pgfem3d/Communication.hpp"
#include "pgf_fe2_rebalancer.h"
#include <stdlib.h>

/* fully encapsulate the client */
struct pgf_FE2_macro_client;

/**
 * Initialize a handle to a client.
 */
void pgf_FE2_macro_client_init(pgf_FE2_macro_client **client,
			       multiscale::net::Network *n);

/**
 * Destroy a client object. Returns NULL for the client handle.
 */
void pgf_FE2_macro_client_destroy(pgf_FE2_macro_client *client);

/**
 * Get/create list of jobs owned by the client and to be computed by
 * one of the servers.
 *
 * Does NOT communicate w/ microscale servers
 */
void pgf_FE2_macro_client_create_job_list(pgf_FE2_macro_client *client,
					  const int n_jobs_max,
					  const pgfem3d::Macroscale *macro,
					  const multiscale::MultiscaleCommunicator *mscom,
					  const int mp_id);

/**
 * Generate initial partitioning of jobs to compute on servers.
 *
 * Communicates w/ microscale servers
 */
void pgf_FE2_macro_client_assign_initial_servers(pgf_FE2_macro_client *client,
					const multiscale::MultiscaleCommunicator *mscom);

/**
 * Reassign jobs to balance load on servers. Send new assignment
 * information to servers to allow data migration while macroscale
 * does other stuff. Either this function or
 * pgf_FE2_macro_client_send_ exit MUST be called after each call to
 * pgf_FE2_macro_client_recv_jobs before
 * pgf_FE2_macro_client_send_jobs can be executed.
 */
void pgf_FE2_macro_client_rebalance_servers(pgf_FE2_macro_client *client,
			const multiscale::MultiscaleCommunicator *mscom,
                        const int heuristic);

/**
 * Send jobs to servers to be computed.
 */
void pgf_FE2_macro_client_send_jobs(pgf_FE2_macro_client *client,
		    const multiscale::MultiscaleCommunicator *mscom,
		    const pgfem3d::Macroscale *macro,
                    const int job_type);

/**
 * Receive finished jobs from the servers. Returns the maximum number
 * of sub-steps taken at the microscale.
 */
void pgf_FE2_macro_client_recv_jobs(pgf_FE2_macro_client *client,
		    pgfem3d::Macroscale *macro,
                    int *max_micro_sub_step);

/**
 * Send signal to servers to exit.
 */
void pgf_FE2_macro_client_send_exit(pgf_FE2_macro_client *client,
			   const multiscale::MultiscaleCommunicator *mscom);

#endif
