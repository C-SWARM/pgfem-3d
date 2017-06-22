/* HEADER */

#pragma once
#ifndef ARC_LENGTH_H
#define ARC_LENGTH_H
#include "data_structure.h"
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
#include "macro_micro_functions.h"
#include "solver_file.h"

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
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return load multiplier
double Multiphysics_Arc_length(GRID *grid,
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
                               MULTIPHYSICS *mp,
                               const int mp_id);

/// Multiscale simulation interface to perform Newton Raphson iteration
///
/// data structures have defined for multiscale simulation differently. 
/// In order to use multiphysics data sturcture for multiscale simulation,
/// data should be reformed from data structure from multiscale to data structure for multiphysics
///
/// \param[in] c structure of macroscale information
/// \param[in,out] s contains the information for the history-dependent solution
/// \param[in] solver_file structure for storing/updating the data
/// \param[in] opts structure PGFem3D option
/// \param[in,out] pores opening volume of failed cohesive interfaces
/// \param[in] dt0 time step size before subdivision, t(n+1) - t(n)
/// \param[in] lm Load multiplier level
/// \param[in,out] DET0 Arc Lengh parameter
/// \param[in,out] DLM0 Arc Lengh parameter  
/// \param[in,out] DLM  Arc Lengh parameter
/// \param[in,out] AT Arc Lengh parameter      
/// \param[in] ARC if ARC=0 : Arc Length method - Crisfield
///                   ARC=1 : Arc Length method - Simo
/// \param[in,out] ITT Arc Lengh parameter  
/// \param[in,out] DAL Arc Lengh parameter
/// \param[in] sup_defl Prescribed deflection
/// \return load multiplier
double Arc_length_multiscale(COMMON_MACROSCALE *c,
                             MACROSCALE_SOLUTION *s,
                             SOLVER_FILE *solver_file,
                             const PGFem3D_opt *opts,
                             double *pores,
                             double dt0,
                             double lm,
                             double *DET0,  
                             double *DLM0,  
                             double *DLM,
                             long *AT,      
                             long ARC,
                             long *ITT,     
                             double *DAL,
                             double *sup_defl);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ARC_LENGTH_H */
