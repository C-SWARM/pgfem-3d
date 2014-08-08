/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef LOAD_BALANCER_H
#define LOAD_BALANCER_H

#include "load_list.h"

int load_balancer(LOAD_LIST *servers,
		  const size_t n_servers);

#endif
