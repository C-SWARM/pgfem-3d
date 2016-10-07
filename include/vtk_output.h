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

typedef int (*write_vtk) (PRINT_MULTIPHYSICS_RESULT *pmr,
                          GRID *grid,
                          LOADING_STEPS *load,
                          FILE *out); 

struct PRINT_MULTIPHYSICS_RESULT
{
  int variable_id;
  char *variable_name;
  char *data_format;
  int data_type;
  void *p_data;
  int m_row;
  int m_col;
  write_vtk write_vtk;  
};


int write_VTK_out(GRID *grid,
                  PRINT_MULTIPHYSICS_RESULT *pD);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  VTK_OUTPUT_H */
