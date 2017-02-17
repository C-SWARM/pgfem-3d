/* HEADER */

#pragma once
#ifndef ARC_LENGTH_H
#define ARC_LENGTH_H

#include "PGFem3D_data_structure.h"
#include "PGFEM_mpi.h"
#include "hypre_global.h"
#include "element.h"
#include "supp.h"
#include "hommat.h"
#include "matgeom.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "pgfem_comm.h"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Arc length procedure.
   */


#ifndef TYPE_ARC_LENGTH_VARIABLES
#define TYPE_ARC_LENGTH_VARIABLES
typedef struct ARC_LENGTH_VARIABLES ARC_LENGTH_VARIABLES;
#endif

/// Arc length struct. It has additional variables for the Arc length scheme.
/// Member variables are taken from lagcy code and not yet been fully
/// identified how the member variables are used in arc length analysis.
struct ARC_LENGTH_VARIABLES
{
  double dt0;
  double *D_R;
  double *U;
  double *DK;
  double *dR;
  double *BS_d_r;
  double *BS_D_R;
  double *BS_rr;
  double *BS_R;
  double *BS_U;
  double *BS_DK;
  double *BS_dR;
  double lm;
  double dAL0;
  double DET0;
  double DLM0;
  double DLM;
  long AT;
  long ARC;
  double dALMAX;
  long ITT;
  double DAL;  
};

/// initialize arc length variable object
/// 
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \return non-zero on internal error
int arc_length_variable_initialization(ARC_LENGTH_VARIABLES *arc);


/// construct arc length variable object
/// 
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \param[in] fv an object containing all field variables
/// \param[in] com an object for communication
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int construct_arc_length_variable(ARC_LENGTH_VARIABLES *arc,
                                  FIELD_VARIABLES *fv,
                                  COMMUNICATION_STRUCTURE *com,
                                  int myrank);

int destruct_arc_length_variable(ARC_LENGTH_VARIABLES *arc);

/// perform arc length analysis
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] time_steps object for time stepping
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \return load multiplier
double Arc_length_test(GRID *grid,
                       MATERIAL_PROPERTY *mat,
                       FIELD_VARIABLES *variables,
                       SOLVER_OPTIONS *sol,
                       LOADING_STEPS *load,
                       COMMUNICATION_STRUCTURE *com,
                       PGFem3D_TIME_STEPPING *time_steps, 
                       CRPL *crpl,
                       MPI_Comm mpi_comm,
                       const double VVolume,
                       const PGFem3D_opt *opts,
                       const int mp_id);

  double Arc_length (long ne,
		     int n_be,
		     long nn,
		     long ndofn,
		     long ndofd,
		     long npres,
		     long nt,
		     long tim,
		     double *times,
		     double nor_min,
		     long iter_max,
		     double dt,
		     double dt0,
		     ELEMENT *elem,
		     BOUNDING_ELEMENT *b_elems,
		     long nbndel,
		     long *bndel,
		     NODE *node,
		     SUPP sup,
		     double *sup_defl,
		     HOMMAT *hommat,
		     MATGEOM matgeom,
		     SIG *sig_e,
		     EPS *eps,
		     int *Ap,
		     int *Ai,
		     PGFEM_HYPRE_solve_info *PGFEM_hypre,
		     double *RRn,
		     double *f_defl,
		     CRPL *crpl,
		     double stab,
		     long nce,
		     COEL *coel,
		     double *r,
		     double *f,
		     double *d_r,
		     double *D_R,
		     double *rr,
		     double *R,
		     double *RR,
		     double *f_u,
		     double *U,
		     double *DK,
		     double *dR,
		     double *BS_f,
		     double *BS_d_r,
		     double *BS_D_R,
		     double *BS_rr,
		     double *BS_R,
		     double *BS_RR,
		     double *BS_f_u,
		     double *BS_U,
		     double *BS_DK,
		     double *BS_dR,
		     long FNR,
		     double lm,
		     double dAL0,
		     double *DET0,
		     double *DLM0,
		     double *dlmdlm,
		     long gr2,
		     long gr4,
		     SIG *sig_n,
		     char *out_dat,
		     long *print,
		     long *AT,
		     long ARC,
		     double dALMAX,
		     long *ITT,
		     double *DAL,
		     double *pores,
		     /*long nge,
		       GEEL *geel,
		       long ngn,
		       GNOD *gnod*/
		     long *DomDof,
		     int GDof,
		     COMMUN comm,
		     double err,
		     double *NORM,
		     MPI_Comm mpi_comm,
		     const double VVolume,
		     const PGFem3D_opt *opts,
		     const int mp_id);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ARC_LENGTH_H */
