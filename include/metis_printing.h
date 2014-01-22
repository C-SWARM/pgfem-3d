#ifndef METIS_PRINTING_H
#define METIS_PRINTING_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Prints local reorder array to file local::global::new global */
  void print_loc_reorder(char *name, int *order, int lndof, int offset);

  /** Prints global reorder array to file global::new global */
  void print_reorder(char *name, int *order, int gndof);

  /** Prints sizes vector */
  void print_sizes(char *name, int *sizes, int nproc);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef METIS_PRINTING_H */
