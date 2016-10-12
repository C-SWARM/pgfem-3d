/* HEADER */
#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

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

enum{MECHANICAL_Var_Displacement,
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
     MECHANICAL_Var_NO} MECHANICAL_Var;
      
enum{THERMAL_Var_Temperature,
     THERMAL_Var_HeatFlux,
     Thermal_Var_NO} THERMAL_Var;
     
typedef int (*write_vtk) (FILE *out,
                          GRID *grid,
                          FIELD_VARIABLES *fv,                          
                          LOADING_STEPS *load,
                          PRINT_MULTIPHYSICS_RESULT *pmr,
                          const PGFem3D_opt *opts); 

struct PRINT_MULTIPHYSICS_RESULT
{
  int mp_id;
  int physics_id;
  char variable_name[1024];
  int is_point_data;  
  void *p_data;
  int m_row;
  int m_col;
  char data_type[1024];
  write_vtk write_vtk;  
};

int VTK_write_multiphysics_master(PRINT_MULTIPHYSICS_RESULT *pD,
                                  int datano,
                                  const PGFem3D_opt *opts,
                                  int time,
                                  int myrank,
                                  int nproc);

int VTK_write_multiphysics_vtu(GRID *grid,
                               FIELD_VARIABLES *FV,
                               LOADING_STEPS *load,
                               PRINT_MULTIPHYSICS_RESULT *pD,
                               int datano,
                               const PGFem3D_opt *opts,
                               int time,
                               int myrank);
                               
int VTK_construct_PMR(GRID *grid,
                      FIELD_VARIABLES *FV,
                      MULTIPHYSICS *mp,
                      PRINT_MULTIPHYSICS_RESULT *pmr);
                  
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  VTK_OUTPUT_H */
