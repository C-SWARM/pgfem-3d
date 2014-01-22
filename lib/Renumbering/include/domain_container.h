#ifndef DOMAIN_CONTAINER_H
#define DOMAIN_CONTAINER_H

#include "sort_container.h"

typedef struct Domain_container{
  int owner,    /**< The rank of the processor which owns the domain */
    preference, /**< The preference for this domain on the owner */
    value,      /**< The number of entries for this domain on the owner */
    covered;    /**< Boolean for whether the domain is covered */
} Domain_container;

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Initializes Domain container (resolves the domain conflicts as well). */
  void initialize_domain_container(int n, sort_container *sort, Domain_container *Domain);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef DOMAIN_CONTAINER_H */
