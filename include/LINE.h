#ifndef PGFEM3D_LINE_H
#define PGFEM3D_LINE_H

#include "PGFEM_mpi.h"
#include "PGFem3D_data_structure.h"
#include "PGFem3D_options.h"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "pgfem_comm.h"
#include "sig.h"
#include "supp.h"

/// Line search algorith for multiphysics mode
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv array of field variable object
/// \param[in] sol object array for solution scheme
/// \param[in] load object for loading
/// \param[in] com object array for communications
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] dts time step sizes from t(n-1) to t(n), and t(n) to t(n+1)
/// \param[in] t t(n+1)
/// \param[in] mp_id mutiphysics id
/// \param[in] nor Normalize norm
/// \param[in] nor2 Normalize norm
/// \param[in] nor1 norm
/// \param[in] LS1 1/2*f'*f (f=residual)
/// \param[in] iter Newton Raphson iteration number
/// \param[in] max_damage physics based evolution
/// \param[in] dissipation volume weighted dissipation
/// \param[in] tim current time step number
/// \param[in] STEP subdivision number
/// \return info id about convergence
long LINE_S3_MP(Grid *grid,
                MaterialProperty *mat,
                FieldVariables *fv,
                pgfem3d::Solver *sol,
                LoadingSteps *load,
                CommunicationStructure *com,
                CRPL *crpl,
                MPI_Comm mpi_comm,
                const PGFem3D_opt *opts,
                Multiphysics *mp,
                double *dts,
                double t,
                int mp_id,
                double *nor,
                double *nor2,
                double nor1,
                double LS1,
                long iter,
                double *max_damage,
                double *dissipation,
                long tim,
                long STEP);

/// Line search algorithm for Arc Length
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv array of field variable object
/// \param[in] sol object array for solution scheme
/// \param[in] load object for loading
/// \param[in] com object array for communications
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] dts time step sizes from t(n-1) to t(n), and t(n) to t(n+1)
/// \param[in] mp_id mutiphysics id
/// \param[in,out] nor Normalize norm
/// \param[in,out] nor2 Normalize norm
/// \param[in] nor1 norm
/// \param[in] LS1 1/2*f'*f (f=residual)
/// \param[in] iter Newton Raphson iteration number
/// \param[in] max_damage physics based evolution
/// \param[in] dissipation volume weighted dissipation
/// \param[in] STEP subdivision number
/// \param[in] tim current time step number
/// \param[in,out] DLM Arc Length parameter
/// \param[out] gama line search parameter
/// \param[in] dlm Arc Length parameter
/// \param[in] dAL Arc Length parameter
/// \return info id about convergence
long ALINE_S3_MP(Grid *grid,
                 MaterialProperty *mat,
                 FieldVariables *fv,
                 pgfem3d::Solver *sol,
                 LoadingSteps *load,
                 CommunicationStructure *com,
                 CRPL *crpl,
                 MPI_Comm mpi_comm,
                 const PGFem3D_opt *opts,
                 Multiphysics *mp,
                 double *dts,
                 int mp_id,
                 double *nor,
                 double *nor2,
                 double nor1,
                 double LS1,
                 long iter,
                 double *max_damage,
                 double *dissipation,
                 long tim,
                 long STEP,
                 double *DLM,
                 double *gama,
                 double dlm,
                 double dAL);

#endif /* #define PGFEM3D_LINE_H */
