/* HEADER */
/**
 * Define Comm_hints interface.
 *
 * AUTHORS:
 *  Matthew Mosby, University of Notre Dame <mmosby1@nd.edu>
 */
#pragma once
#ifndef COMM_HINTS_H
#define COMM_HINTS_H

#include <stdio.h>

/**
 * Typedef handle to opaque Comm_hints object.
 */
typedef struct COMM_HINTS Comm_hints;

/**
 * Construct a Comm_hints object.
 */
Comm_hints* Comm_hints_construct();

/**
 * Destroy a Comm_hints object.
 *
 * \return non-zero on internal error.
 */
int Comm_hints_destroy(Comm_hints *hints);

/**
 * Read Comm_hints from a file.
 *
 * \return non-zero on internal error.
 */
int Comm_hints_read(Comm_hints *hints,
                    FILE *in);

/**
 * Read Comm_hints from a file given the filename.
 *
 * \return non-zero on internal error.
 */
int Comm_hints_read_filename(Comm_hints *hints,
                             const char *fn);

/**
 * Write the Comm_hints to a file.
 *
 * \return non-zero on internal error.
 */
int Comm_hints_write(const Comm_hints *hints,
                     FILE *out);

/**
 * Write Comm_hints to a file given the filename.
 *
 * \return non-zero on internal error.
 */
int Comm_hints_write_filename(const Comm_hints *hints,
                              const char *fn);


/**
 * \return the number of processes to send to during the assembly
 * process.
 */
int Comm_hints_nsend(const Comm_hints *hints);

/**
 * \return the number of processes to receive from during the assembly
 * process.
 */
int Comm_hints_nrecv(const Comm_hints *hints);

/**
 * Provide const access to the list of processors to send to during
 * the assembly process.
 *
 * \return const pointer to list of processors **DO NOT FREE THIS
 * POINTER**
 */
const int* Comm_hints_send_list(const Comm_hints *hints);

/**
 * Provide const access to the list of processors to receive from
 * during the assembly process.
 *
 * \return const pointer to list of processors **DO NOT FREE THIS
 * POINTER**
 */
const int* Comm_hints_recv_list(const Comm_hints *hints);

/**
 * Utility function to construct the Comm_hints filename.
 *
 * \return string containing filename. User responsible for
 * deallocating memory.
 */
char* Comm_hints_filename(const char *ipath,
                          const char *base_filename,
                          const int rank);

#endif /* #ifndef  */
