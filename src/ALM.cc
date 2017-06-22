#include "ALM.h"
#include <math.h>

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef FD_RESIDUALS_H
#include "fd_residuals.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#include "Arc_length.h"

#ifndef ALM_DEBUG
#define ALM_DEBUG 0
#endif

static const double PII = 3.141592653589793238462643;

double D_lam_ALM (long ndofd,
		  double *BS_rr,
		  double *BS_d_r,
		  double *BS_D_R,
		  double *BS_R,
		  double *BS_DK,
		  double dlm,
		  double dAL,
		  long *DomDof,
		  MPI_Comm mpi_comm)
{
  return 0.0;
}

/***********************************************************************/
/***********************************************************************/

double d_ALM2 (long ndofd,
	       double *rr,
	       double *R,
	       double *DK,
	       double d_lm)
{
  return 0.0;
}

double d_lam_ALM2 (long ndofd,
		   double *rr,
		   double *R,
		   double *DK,
		   double dAL,
		   double DET,
		   double DET0,
		   double dlm0,
		   double nor_min,
		   double *dR)
{
  return 0.0;
}

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
                     double dt)                                        
{
  return 0.0;
}


/*********************************************************************/
/*********************************************************************/

double d_ALM4 (long ndofd,
	       double *BS_rr,
	       double *BS_DK,
	       double dlm,
	       long *DomDof,
	       MPI_Comm mpi_comm)
/*
  SIMO
*/
{
  return 0.0;
}

double d_lam_ALM4 (long ndofd,
		   double *BS_rr,
		   double *BS_DK,
		   double *BS_dR,
		   double dAL,
		   long *DomDof,
		   MPI_Comm mpi_comm)
/*
  SIMO
*/
{
  return 0.0;
}

double D_lam_ALM4 (long ndofd,
		   double *BS_rr,
		   double *BS_d_r,
		   double *BS_D_R,
		   double *BS_DK,
		   double dlm,
		   double dAL,
		   long *DomDof,
		   MPI_Comm mpi_comm)
/*
  SIMO  
*/
{
  return 0.0;
}
