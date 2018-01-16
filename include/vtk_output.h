#ifndef PGFEM3D_VTK_OUTPUT_H
#define PGFEM3D_VTK_OUTPUT_H

#include "PGFem3D_data_structure.h"
#include "PGFem3D_options.h"
#include "cohesive_element.h"
#include "data_structure.h"
#include "element.h"
#include "ensight.h"
#include "eps.h"
#include "node.h"
#include "sig.h"
#include "supp.h"

/** Print the master VTK file (call on only 1 CPU)*/
void VTK_print_master(const char *path,
                      const char *base_name,
                      int time,
                      int nproc,
                      const PGFem3D_opt *opts);

/** Print master VTK file for cohesive elements (call on only 1 CPU)*/
void VTK_print_cohesive_master(const char *path,
                               const char *base_name,
                               int time,
                               int nproc,
                               const PGFem3D_opt *opts);

/** Print the individual vtu files */
void VTK_print_vtu(const char *path,
                   const char *base_name,
                   int time,
                   int myrank,
                   long ne,
                   long nn,
                   Node *node,
                   Element *elem,
                   SUPP sup,
                   double *r,
                   double *P,
                   double *V,
                   SIG *sig,
                   EPS *eps,
                   const PGFem3D_opt *opts,
                   const int mp_id);

/** Print the individual vtu files for the cohesive elements */
void VTK_print_cohesive_vtu(const char *path,
                            const char *base_name,
                            int time,
                            int myrank,
                            long nce,
                            Node *node,
                            COEL *coel,
                            SUPP sup,
                            double *r,
                            Ensight *ensight,
                            const PGFem3D_opt *opts,
                            const int mp_id);

struct PRINT_MULTIPHYSICS_RESULT;

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
             MECHANICAL_Var_PrincipalStress,
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
typedef int (*write_vtk_t) (FILE *out,
                            Grid *grid,
                            const MaterialProperty *mat,
                            FieldVariables *fv,
                            LoadingSteps *load,
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
  write_vtk_t write_vtk;    /// function pointer for generalizing writing simulation results
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
int VTK_write_multiphysics_vtu(Grid *grid,
                               const MaterialProperty *mat,
                               FieldVariables *FV,
                               LoadingSteps *load,
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
int VTK_construct_PMR(Grid *grid,
                      FieldVariables *FV,
                      Multiphysics *mp,
                      PRINT_MULTIPHYSICS_RESULT *pmr);

#endif /* #define PGFEM3D_VTK_OUTPUT_H */
