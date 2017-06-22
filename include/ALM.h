#ifndef ALM_H
#define ALM_H

#include "data_structure.h"

#include "PGFEM_mpi.h"

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef PGFEM_COMM_H
#include "pgfem_comm.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#include "PGFem3D_data_structure.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  double D_lam_ALM (long ndofd,
		    double *BS_rr,
		    double *BS_d_r,
		    double *BS_D_R,
		    double *BS_R,
		    double *BS_DK,
		    double dlm,
		    double dAL,
		    long *DomDof,
		    MPI_Comm mpi_comm);

  /** Returns 0. */
  double d_ALM2 (long ndofd,
		 double *rr,
		 double *R,
		 double *DK,
		 double d_lm);

  /** Returns 0. */
  double d_lam_ALM2 (long ndofd,
		     double *rr,
		     double *R,
		     double *DK,
		     double dAL,
		     double DET,
		     double DET0,
		     double dlm0,
		     double nor_min,
		     double *dR);

/// D_lam_ALM2_MP
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv array of field variable object
/// \param[in] sol object array for solution scheme
/// \param[in] load object for loading
/// \param[in] com object array for communications
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dlm Arc_length parameter
/// \param[in] dAL Arc_length parameter
/// \param[in] dt times step size
/// \return computed DLM
double D_lam_ALM2_MP(GRID *grid,
                     MATERIAL_PROPERTY *mat,
                     FIELD_VARIABLES *fv,
                     SOLVER_OPTIONS *sol,
                     LOADING_STEPS *load,
                     COMMUNICATION_STRUCTURE *com,
                     CRPL *crpl,
                     MPI_Comm mpi_comm,
                     const PGFem3D_opt *opts,
                     MULTIPHYSICS *mp,
                     int mp_id,
                     double dlm,
                     double dAL,
                     double dt);		     

  double d_ALM4 (long ndofd,
		 double *BS_rr,
		 double *BS_DK,
		 double dlm,
		 long *DomDof,
		 MPI_Comm mpi_comm);

  double d_lam_ALM4 (long ndofd,
		     double *BS_rr,
		     double *BS_DK,
		     double *BS_dR,
		     double dAL,
		     long *DomDof,
		     MPI_Comm mpi_comm);

  double D_lam_ALM4 (long ndofd,
		     double *BS_rr,
		     double *BS_d_r,
		     double *BS_D_R,
		     double *BS_DK,
		     double dlm,
		     double dAL,
		     long *DomDof,
		     MPI_Comm mpi_comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ALM_H */
