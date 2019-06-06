/* HEADER */
#pragma once
#ifndef SUPP_H
#define SUPP_H

/** Structure of support nodes */
struct SUPP_1 {
  long       nfn;             //!< number of free nodes
  long       nsn;             //!< number of supported nodes
  long       ndn;             //!< number of nodes with prescribed deflection
  long       npd;             //!< number of magnidudes of prescribed deflection
  long       nde;             //!< number of elements with prescribed deflection
  long     nd_be;             //!< number of bounding elements w/ pre. def.
  long     *free;             //!< list of free nodes
  long     *supp;             //!< list of supported nodes
  long     *lnpd;             //!< list of nodes with prescribed deflection
  long     *lepd;             //!< list of elements w/ pre. def
  long    *lbepd;             //!< list of bounding elements
  double   *defl;             //!< value of deflection at the node.
  double *defl_d;

  /// Multiscale information. For interfaces, the displacement jump vector is
  /// computed from the first six (6) prescribed displacements
  /// @{
  double *F0;          //!< macroscale deformation gradient contribution
  double *N0;          //!< macroscale unit normal
  double  lc;          //!< multiscale scaling parameter. For interfaces, is the
                       //   thickness of the unit cell/interface layer.
  double  v0;          //!< Volume of microstructure V = B U H [Miehe 2002]
  int multi_scale;     //!< flag for whether or not multiscale routines should
                       //   be executed.
  /// @}
};

using SUPP = SUPP_1*;

void destroy_supp(SUPP sup);

#endif /* #ifndef  */
