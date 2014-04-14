/* HEADER */
#pragma once
#ifndef MESH_LOAD_H
#define MESH_LOAD_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of load on nodes */
  struct ZATNODE{
    long nod;
    double *load;
  };
  typedef struct ZATNODE ZATNODE;

  /** Structure of load on elements */
  struct ZATELEM{
    /** Number of loaded element */
    long elem;
    /** Number of nodes on th loaded surface */
    long *sur;
    /** Vector of load */
    double *load;
  };
  typedef struct ZATELEM ZATELEM;

  ZATNODE* build_zatnode(long ndofn,long nln);

  void destroy_zatnode(ZATNODE* zn, long nln);

  ZATELEM* build_zatelem(long ndofn,long nle);

  void destroy_zatelem(ZATELEM* ze, long nle);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
