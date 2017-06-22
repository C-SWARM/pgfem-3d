/* HEADER */
#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include "data_structure.h"
#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
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

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef ENSIGHT_H
#include "ensight.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#include "PGFem3D_data_structure.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Print the master VTK file (call on only 1 CPU)*/
  void VTK_print_master(char *path,
			char *base_name,
			int time,
			int nproc,
			const PGFem3D_opt *opts);

  /** Print master VTK file for cohesive elements (call on only 1 CPU)*/
  void VTK_print_cohesive_master(char *path,
				 char *base_name,
				 int time,
				 int nproc,
				 const PGFem3D_opt *opts);

  /** Print the individual vtu files */
  void VTK_print_vtu(char *path,
		     char *base_name,
		     int time,
		     int myrank,
		     long ne,
		     long nn,
		     NODE *node,
		     ELEMENT *elem,
		     SUPP sup,
		     double *r,
		     SIG *sig,
		     EPS *eps,
		     const PGFem3D_opt *opts,
		     const int mp_id);

  /** Print the individual vtu files for the cohesive elements */
  void VTK_print_cohesive_vtu(char *path,
			      char *base_name,
			      int time,
			      int myrank,
			      long nce,
			      NODE *node,
			      COEL *coel,
			      SUPP sup,
			      double *r,
			      ENSIGHT ensight,
			      const PGFem3D_opt *opts,
			      const int mp_id);

struct PRINT_MULTIPHYSICS_RESULT;
#ifndef TYPE_PRINT_MULTIPHYSICS_RESULT
#define TYPE_PRINT_MULTIPHYSICS_RESULT
typedef struct PRINT_MULTIPHYSICS_RESULT PRINT_MULTIPHYSICS_RESULT;
#endif

/// Mechanical part of index of output vailabes
/// for writing simulation resuls 
typedef enum{MECHANICAL_Var_Displacement,
             MECHANICAL_Var_MacroDisplacement,
             MECHANICAL_Var_NodalPressure,
             MECHANICAL_Var_CauchyStress,
             MECHANICAL_Var_EulerStrain,
             MECHANICAL_Var_EffectiveStrain,
             MECHANICAL_Var_EffectiveStress,
             MECHANICAL_Var_CellProperty,
             MECHANICAL_Var_Damage,
             MECHANICAL_Var_Chi,
             MECHANICAL_Var_F,
             MECHANICAL_Var_P,
             MECHANICAL_Var_W,
             MECHANICAL_Var_ElementPressure,
             MECHANICAL_Var_ElementVolume,
             MECHANICAL_Var_Density,
             MECHANICAL_Var_HydrostaticStress,
             MECHANICAL_Var_NO} MECHANICAL_Var;

/// Thermal part of index of output vailabes
/// for writing simulation resuls       
typedef enum{THERMAL_Var_Temperature,
             THERMAL_Var_HeatFlux,
             THERMAL_Var_HeatGenerations,     
             Thermal_Var_NO} THERMAL_Var;

/// Chemical part of index of output vailabes
/// for writing simulation resuls       
typedef enum{CHEMICAL_VAR_SPECIES,
             CHEMICAL_Var_NO} CHEMICAL_Var;     

/// function pointer for generalizing writing simulation results
/// Using this function pointer, different variables can be written
/// in a one function call using PRINT_MULTIPHYSICS_RESULT.
/// e.g. PRINT_MULTIPHYSICS_RESULT *pmr;
///      pmr->write_vtk = [function name to write any variables]     
typedef int (*write_vtk) (FILE *out,
                          GRID *grid,
                          const MATERIAL_PROPERTY *mat,
                          FIELD_VARIABLES *fv,                          
                          LOADING_STEPS *load,
                          PRINT_MULTIPHYSICS_RESULT *pmr,
                          const PGFem3D_opt *opts); 

/// structure for writing multiphysics simulation results
/// In writing vtk outputs, this object is constructed and used
/// to pass vaiable types, number of components, and a manner (function)
/// to write variables based on physics
struct PRINT_MULTIPHYSICS_RESULT
{
  int mp_id;                /// multiphyiscs id, start from 0 and + integer
  int physics_id;           /// physics id, Mechanical=0, Thermal = 1, Chemical = 2, ...
  char variable_name[1024]; /// variable name to be witten
  int is_point_data;        /// identifer of data (1=point, others cell data)
  void *p_data;             /// data pointer containing variblas to be written
                            /// When VTK_write_data_double function is used, 
                            /// p_data should be allocated in advance.
  int m_row;                /// number of rows
  int m_col;                /// number of columns
  char data_type[1024];     /// data type (Int64, Float64, Float32, ...)
  write_vtk write_vtk;      /// function pointer for generalizing writing simulation results
};


/// write vtk master file based on physics
///
/// \param[in] pD a PRINT_MULTIPHYSICS_RESULT struct for writing results based on physics
/// \param[in] datano number data(vaialbes) to be written
/// \param[in] opts structure PGFem3D option
/// \param[in] time time step number
/// \param[in] myrank current process rank
/// \param[in] nproc number of MPI processes
/// \return non-zero on internal error
int VTK_write_multiphysics_master(PRINT_MULTIPHYSICS_RESULT *pD,
                                  int datano,
                                  const PGFem3D_opt *opts,
                                  int time,
                                  int myrank,
                                  int nproc);

/// write simulation results in vtk format based on physics
///
/// \param[in] grid an object containing all mesh data
/// \param[in] mat a material object
/// \param[in] FV array of field variables
/// \param[in] load object for loading
/// \param[in] pD a PRINT_MULTIPHYSICS_RESULT struct for writing results based on physics
/// \param[in] datano number data(vaialbes) to be written
/// \param[in] opts structure PGFem3D option
/// \param[in] time time step number
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int VTK_write_multiphysics_vtu(GRID *grid,
                               const MATERIAL_PROPERTY *mat,
                               FIELD_VARIABLES *FV,
                               LOADING_STEPS *load,
                               PRINT_MULTIPHYSICS_RESULT *pD,
                               int datano,
                               const PGFem3D_opt *opts,
                               int time,
                               int myrank);

/// construct PRINT_MULTIPHYSICS_RESULT array based on physics
///
/// \param[in] grid an object containing all mesh data
/// \param[in] FV array of field variables
/// \param[in] mp an object for multiphysics stepping
/// \param[in] pmr a PRINT_MULTIPHYSICS_RESULT struct for writing results based on physics
/// \return non-zero on internal error                               
int VTK_construct_PMR(GRID *grid,
                      FIELD_VARIABLES *FV,
                      MULTIPHYSICS *mp,
                      PRINT_MULTIPHYSICS_RESULT *pmr);
                  
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  VTK_OUTPUT_H */
