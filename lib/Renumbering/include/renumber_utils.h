#ifndef RENUMBER_UTILS_H
#define RENUMBER_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Comparator for integers */
  int compare_integer (const void * a, const void * b);

  /** Comparator for sort container which retains original index. */
  int compare_sort_container (const void * a, const void * b);

  /** Same as compare_sort_container, but sorts from high to low. */
  int compare_sort_container_reverse (const void * a, const void * b);

  /** Comparator to revert sort container to original configuration. */
  int compare_sort_container_revert (const void * a, const void * b);

#ifdef __cplusplus
{
#endif /* #ifdef __cplusplus */

#endif /* #ifndef RENUMBER_UTILS_H */
