/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#pragma once
#ifndef MICRO_SERVER_H
#define MICRO_SERVER_H

#include <stdlib.h>
#include "PGFEM_mpi.h"
#include "pgf_fe2_job.h"

struct pgf_FE2_server_stats{
  double avg;
  double std;
  double min;
  double max;
};
typedef struct pgf_FE2_server_stats pgf_FE2_server_stats;

struct pgf_FE2_server_job_list{
  size_t n_jobs;
  pgf_FE2_job *jobs;
  pgf_FE2_server_stats stats;
};
typedef struct pgf_FE2_server_job_list pgf_FE2_server_job_list;

typedef int* pgf_FE2_server_rebalance;

void pgf_FE2_server_rebalance_build(pgf_FE2_server_rebalance *t,
				    const size_t n_keep,
				    const size_t n_send,
				    const size_t n_recv);

void pgf_FE2_server_rebalance_destroy(pgf_FE2_server_rebalance *t);

/**
 * Get the number of jobs to keep.
 */
int pgf_FE2_server_rebalance_n_keep(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer containing job ids to keep.
 */
int* pgf_FE2_server_rebalance_keep_buff(const pgf_FE2_server_rebalance *t);

/**
 * Get the number of jobs to send.
 */
int pgf_FE2_server_rebalance_n_send(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer containing job ids to send.
 */
int* pgf_FE2_server_rebalance_send_buff(const pgf_FE2_server_rebalance *t);

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
int* pgf_FE2_server_rebalance_recv_buff(const pgf_FE2_server_rebalance *t);

/**
 * Get pointer to beginning of buffer contining sources of
 * communication.
 */
int* pgf_FE2_server_rebalance_recv_src(const pgf_FE2_server_rebalance *t);


#endif
