#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "LINE.h"
#include "ALM.h"
#include "Arc_length.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "Newton_Raphson.h"
#include "PGFEM_io.h"
#include "bounding_element_utils.h"
#include "dynamics.h"
#include "enumerations.h"
#include "fd_residuals.h"
#include "integration.h"
#include "macroscopic_load_AL.h"
#include "matice.h"
#include "utils.h"
#include "vol_damage_int_alg.h"
#include <cmath>
#include <cstring>

namespace {
const constexpr int periodic = 0;
using pgfem3d::Solver;
}

/// Line search algorithm for multiphysics mode
///
/// \param[in,out] residuals_loc_time - residuals matrix compute time
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
long LINE_S3_MP(double *residuals_loc_time, 
                Grid *grid,
                MaterialProperty *mat,
                FieldVariables *fv,
                Solver *sol,
                LoadingSteps *load,
                CommunicationStructure *com,
                CRPL *crpl,
                MPI_Comm mpi_comm,
                const PGFem3D_opt *opts,
                const Multiphysics& mp,
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
                long STEP)
{
  long i,N,M,INFO,GInfo;
  double LS2,slope,tmplam,rhs1,rhs2,AL{},a,b,f2{},disc,scale,nor3;
  const char *error[]={"inf","-inf","nan"};
  char str1[500];
  double dt = dts[DT_NP1];

  int myrank,nproc;
  MPI_Comm_rank(mpi_comm,&myrank);
  MPI_Comm_size(mpi_comm,&nproc);

  scale = LS1; LS1 /= scale;

  slope = -2.*LS1;

  /* Compute norm LS2 = 1./2. * ss (f,f,ndofd)/scale; */
  nor3 = ss (fv->BS_f,fv->BS_f,com->DomDof[myrank]);
  /* Gather residual nor from Domains */
  MPI_Allreduce (&nor3,&LS2,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  LS2 *= 1./(2.*scale);

  sol->gama = 1.0;
  INFO = 0;
  if (iter == 0)
    return (0);

  if (myrank == 0)
    PGFEM_printf("nor = %12.8e :: slope = %12.8e || LS1 = %12.8e"
                 "|| LS2 = %12.8e || HU = %12.8e\n",
                 *nor,slope,LS1,LS2,LS1 + 1.e-4*(sol->gama)*slope);

  while (LS2 > LS1 + 1.e-4*(sol->gama)*slope){

    if (sol->gama == 1.0){
      tmplam = -slope/(2.*(LS2-LS1-slope));
    }
    else{
      rhs1 = LS2 - LS1 - (sol->gama)*slope;
      rhs2 = f2 - LS1 - AL*slope;
      a = (rhs1/((sol->gama)*(sol->gama))-rhs2/(AL*AL))/((sol->gama)-AL);
      b = (-AL*rhs1/((sol->gama)*(sol->gama))+(sol->gama)*rhs2/(AL*AL))/((sol->gama)-AL);
      if (a == 0.0) tmplam = -slope/(2.*b);
      else{
        disc = b*b - 3.*a*slope;
        if (disc < 0.0) tmplam = 0.5*(sol->gama);
        else
          if (b <= 0.0) tmplam = (-b+sqrt(disc))/(3.0*a);
          else tmplam = -slope/(b+sqrt(disc));
      }/* else a == 0 */
      if (tmplam > 0.5*(sol->gama)) tmplam = 0.5*(sol->gama);
    }/* end (sol->gama) == 1*/

    /* Gama */
    AL = (sol->gama);
    if (tmplam >= 0.1*(sol->gama))
      sol->gama = tmplam;
    else
      sol->gama *= 0.1;

    if (sol->gama < 1.e-1*(sol->nor_min))
      return (1);

    /* Decrease increment of displacements */
    for (i=0;i<fv->ndofd;i++)
      fv->f[i] = (sol->gama)*fv->dd_u[i];

    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL) {
      INFO = integration_alg (grid->ne,fv->ndofn,fv->ndofd,fv->npres,crpl,grid->element,
                              grid->node,fv->d_u,fv->f,load->sups[mp_id],mat->matgeom,mat->hommat,fv->eps,
                              fv->sig,tim,iter,dt,sol->nor_min,STEP,1,opts,mp_id);

      /* Gather INFO from all domains */
      MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
      if (GInfo == 1) {
        INFO = 1;
        return (1);
      }
    }

    /* Update deformations */
    for (i=0;i<fv->ndofd;i++) {
      fv->f[i] += fv->d_u[i];
      fv->f_u[i] = 0.0;
    }

    INFO = vol_damage_int_alg(grid->ne,fv->ndofn,fv->f,fv->u_np1,grid->element,grid->node,
                              mat->hommat,load->sups[mp_id],dt,iter,mpi_comm,
                              fv->eps,fv->sig,max_damage,dissipation,
                              opts->analysis_type,mp_id);

    bounding_element_communicate_damage(grid->n_be,grid->b_elems,grid->ne,fv->eps,mpi_comm);

    MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
    if (GInfo == 1) {
      if(myrank == 0){
        PGFEM_printf("Inverted element detected (vol_damage_int_alg).\n"
                     "Subdividing load.\n");
      }
      INFO = 1;
      return 1;
    }

    /* Residuals */
    INFO = 0;
    *residuals_loc_time += compute_residuals_for_NR(&INFO,grid,mat,fv,sol,load,crpl,mpi_comm,
                                                    opts,mp,mp_id,t,dts, 1);

    /* Compute Euclidian norm */
    for (i=0;i<fv->ndofd;i++)
      fv->f[i] = fv->RR[i] - fv->f_u[i];

    /* Local to global transformation */
    LToG(fv->f,fv->BS_f,myrank,nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);

    nor3 = ss (fv->BS_f,fv->BS_f,com->DomDof[myrank]);

    /* Gather residual nor from Domains */
    MPI_Allreduce(&nor3,nor,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    nor3 = *nor;
    *nor = *nor2 = sqrt(nor3);
    *nor /= nor1;

    sprintf (str1,"%f",*nor);
    for (N=0;N<3;N++){
      M = strcmp(error[N],str1);

      if (M == 0){
        if (myrank == 0)
          PGFEM_printf("ERROR in the algorithm LS2 : nor = %s\n",error[N]);
        return (1);
      }
    }

    /* New LS and gama */
    f2 = LS2;
    LS2 = 1./2. * nor3/scale;

    if (myrank ==0)
      PGFEM_printf ("gama = %12.12f || tmplam = %12.8f ||"
                    "NOR = %12.8e  || LS2 = %12.8e\n",
                    (sol->gama),tmplam,*nor,LS2);

  }/* end LINE SEARCH */

  return (INFO);
}

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
                 Solver *sol,
                 LoadingSteps *load,
                 CommunicationStructure *com,
                 CRPL *crpl,
                 MPI_Comm mpi_comm,
                 const PGFem3D_opt *opts,
                 const Multiphysics& mp,
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
                 double dAL)
{
  ARC_LENGTH_VARIABLES *arc = sol->arc;
  double dt = dts[DT_NP1];
  double t = 0.0;

  long i,N,M,INFO,GInfo;
  double LS2,slope,tmplam,rhs1,rhs2,AL{},a,b,f2{},disc,scale,tmp;
  const char *error[]={"inf","-inf","nan"};
  char str1[500];

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  scale = LS1;
  LS1 /= scale;
  slope = -2.*LS1;
  tmp = ss(fv->BS_f,fv->BS_f,com->DomDof[myrank]);
  MPI_Allreduce(&tmp,&LS2,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  LS2 *= 1./(2.*scale);

  *gama = 1.0;
  INFO = 0;
  if (iter == 0) return (0);

  if (myrank == 0)
    PGFEM_printf("nor = %12.8e :: slope = %12.8e || LS1 = %12.8e "
                 "|| LS2 = %12.8e || HU = %12.8e\n",
                 *nor,slope,LS1,LS2,LS1 + 1.e-4**gama*slope);

  while (LS2 > LS1 + 1.e-4**gama*slope){

    if (*gama == 1.0){
      tmplam = -slope/(2.*(LS2-LS1-slope));
    }
    else{
      rhs1 = LS2 - LS1 - *gama*slope;
      rhs2 = f2 - LS1 - AL*slope;
      a = (rhs1/(*gama**gama)-rhs2/(AL*AL))/(*gama-AL);
      b = (-AL*rhs1/(*gama**gama)+*gama*rhs2/(AL*AL))/(*gama-AL);
      if (a == 0.0) tmplam = -slope/(2.*b);
      else{
        disc = b*b - 3.*a*slope;
        if (disc < 0.0) tmplam = 0.5**gama;
        else
          if (b <= 0.0) tmplam = (-b+sqrt(disc))/(3.0*a);
          else tmplam = -slope/(b+sqrt(disc));
      }/* else a == 0 */
      if (tmplam > 0.5**gama) tmplam = 0.5**gama;
    }/* end *gama == 1*/

    /* Gama */
    AL = *gama;
    if (tmplam >= 0.1**gama)
      *gama = tmplam;
    else
      *gama *= 0.1;
    if (*gama < 1.e-1*(sol->nor_min))
      return (1);

    /* Decrease increment of displacements */
    for (i=0;i<com->DomDof[myrank];i++) {
      arc->BS_D_R[i] = *gama*(arc->BS_U[i]);
      fv->BS_f[i] = *gama*(arc->BS_rr[i]);
    }

    /* dlam */
    if (arc->ARC == 0){
      *DLM = D_lam_ALM (fv->ndofd,fv->BS_f,arc->BS_d_r,arc->BS_D_R,arc->BS_R,arc->BS_DK,
                        dlm,dAL,com->DomDof,mpi_comm);
    }
    if (arc->ARC == 1)
      *DLM = D_lam_ALM4 (fv->ndofd,fv->BS_f,arc->BS_d_r,arc->BS_D_R,arc->BS_DK,dlm,
                         dAL,com->DomDof,mpi_comm);

    sprintf (str1,"%f",*DLM);
    for (N=0;N<3;N++){
      M = strcmp(error[N],str1);

      if (M == 0) {
        if (myrank == 0) PGFEM_printf("Complex root in ARC-LENGTH method\n");
        return (1);
      }
    }

    /* Displacement update */
    for (i=0;i<com->DomDof[myrank];i++)
      arc->BS_D_R[i] = *gama * (arc->BS_U[i] + *DLM*(arc->BS_rr[i]));

    /* macroscopic load */
    if (periodic == 1) macroscopic_load_AL (fv->eps[0].type,arc->lm  +dlm+*DLM,fv->eps);

    /* Global to local transformation */
    GToL(arc->BS_D_R,arc->D_R,myrank,nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);

    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL) {
      INFO=integration_alg (grid->ne,fv->ndofn,fv->ndofd,fv->npres,crpl,grid->element,grid->node,
                            fv->d_u,arc->D_R,load->sups[mp_id] ,mat->matgeom ,mat->hommat,fv->eps,fv->sig ,
                            tim,iter,dt,sol->nor_min,STEP,0,opts,mp_id);

      /* Gather INFO from all domains */
      MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
      if (GInfo == 1) {
        INFO = 1;
        return (1);
      }
    }

    /* Update deformations */
    for (i=0;i<fv->ndofd;i++) {
      fv->f[i] = fv->d_u[i] + arc->D_R[i];
      fv->f_u[i] = 0.0;
    }

    INFO = vol_damage_int_alg(grid->ne,fv->ndofn,fv->f,fv->u_np1,grid->element,grid->node,
                              mat->hommat,load->sups[mp_id] ,dt,iter,mpi_comm,
                              fv->eps,fv->sig ,max_damage,dissipation,
                              opts->analysis_type,mp_id);

    MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
    if (GInfo == 1) {
      if(myrank == 0){
        PGFEM_printf("Inverted element detected (vol_damage_int_alg).\n"
                     "Subdividing load.\n");
      }
      INFO = 1;
      return 1;
    }

    /* Residuals */
    fd_residuals_MP(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,mp_id,t,dts,1);


    /* Compute Euclidean norm */
    for (i=0;i<fv->ndofd;i++)
      fv->f[i] = (arc->lm   + dlm + *DLM)*(fv->R[i]) - fv->f_u[i];

    /* L -> G : f */
    LToG(fv->f,fv->BS_f,myrank,nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);

    *nor = ss(fv->BS_f,fv->BS_f,com->DomDof[myrank]);
    MPI_Allreduce(nor,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

    *nor = *nor2 = sqrt (tmp);
    *nor /= nor1;
    sprintf (str1,"%f",*nor);
    for (N=0;N<3;N++){
      M = strcmp(error[N],str1);
      if (M == 0){
        if (myrank == 0)
          PGFEM_printf("ERROR in the algorithm LS2 : nor = %s\n",error[N]);
        return (1);
      }
    }

    /* New LS and gama */
    f2 = LS2;
    LS2 = 1./2. * tmp/scale;

    if (myrank == 0)
      PGFEM_printf ("gama = %12.12f || tmplam = %12.8f || NOR = %12.8e "
                    " || LS2 = %12.8e\n",
                    *gama,tmplam,*nor,LS2);

  }/* end LINE SEARCH */

  return (INFO);
}
