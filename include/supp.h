/* HEADER */
#ifndef SUPP_H
#define SUPP_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of support nodes */
  struct SUPP_1{
    /** 
	nfn - number of free nodes
	nsn - number of supported nodes
	ndn - number of nodes with prescribed deflection
	npd - number of magnidudes of prescribed deflection
	nde - number of elements with prescribed deflection
	nd_be - number of bounding elements w/ pre. def.
    */
    long nfn,nsn,ndn,npd,nde,nd_be;
    /** 
	free - list of free nodes
	supp - list of supported nodes
	lnpd - list of nodes with prescribed deflection
	lepd - list of elements w/ pre. def
	lbepd - -- " -- bounding elements ---"---
    */
    long *free,*supp,*lnpd,*lepd,*lbepd;
    /** defl - value of deflection at the node. */
    double *defl;
    double *defl_d;

    /** multiscale information. For interfaces, the displacement jump
	vector is computed from the first six (6) prescribed
	displacements */
    double *F0; /**< macroscale deformation gradient contribution */
    double *N0; /**< macroscale unit normal */
    double lc; /**< multiscale scaling parameter. For interfaces, is the
		 thickness of the unit cell/interface layer. */
    double v0; /**< Volume of microstructure V = B U H [Miehe 2002]*/
    int multi_scale; /*< flag for whether or not multiscale routines
		       should be executed. */
  };
  typedef struct SUPP_1 SUPP_1;
  typedef struct SUPP_1 *SUPP;

  void destroy_supp(SUPP sup);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
