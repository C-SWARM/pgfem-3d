#ifndef ADDRESSES_H
#define ADDRESSES_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Computes addresses of the diagonal matrix elements. */
  void addresses (long *adr,long ndof);

  /** Computes addresses of the diagonal matrix elements. */
  void addresses_nonsym (long *adr,long ndof);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ADDRESSES_H */
