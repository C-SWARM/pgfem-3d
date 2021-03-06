#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Arc_length.h"
#include "ALM.h"
#include "LINE.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "PGFEM_io.h"
#include "PGFem3D_data_structure.h"
#include "displacement_based_element.h"
#include "dynamics.h"
#include "enumerations.h"
#include "fd_increment.h"
#include "fd_residuals.h"
#include "incl.h"
#include "integration.h"
#include "macro_micro_functions.h"
#include "macroscopic_load_AL.h"
#include "matice.h"
#include "out.h"
#include "press_theta.h"
#include "res_fini_def.h"
#include "solver_file.h"
#include "stabilized.h"
#include "stiffmat_fd.h"
#include "subdivision.h"
#include "utils.h"
#include "vol_damage_int_alg.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>

#ifndef ARC_DEBUG
#define ARC_DEBUG 0
#endif

#ifndef ARC_UPDATE
#define ARC_UPDATE 0
#endif

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#ifndef PFEM_DEBUG_ALL
#define PFEM_DEBUG_ALL 0
#endif

static constexpr int periodic = 0;

using namespace multiscale::net;

namespace {
  using pgfem3d::CommunicationStructure;
  using pgfem3d::MultiscaleCommon;
  using pgfem3d::MULTISCALE_SOLUTION;
}

/// initialize arc length variable object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs
///                  except ARC=1)
///
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \return non-zero on internal error
int arc_length_variable_initialization(ARC_LENGTH_VARIABLES *arc)
{
  int err = 0;
  arc->dt0    = 0.0;
  arc->D_R    = NULL;
  arc->U      = NULL;
  arc->DK     = NULL;
  arc->dR     = NULL;
  arc->BS_d_r = NULL;
  arc->BS_D_R = NULL;
  arc->BS_rr  = NULL;
  arc->BS_R   = NULL;
  arc->BS_U   = NULL;
  arc->BS_DK  = NULL;
  arc->BS_dR  = NULL;
  arc->lm     = 0.0;
  arc->dAL0   = 0.0;
  arc->DET0   = 0.0;
  arc->DLM0   = 0.0;
  arc->DLM    = 0.0;
  arc->AT     = 0;
  arc->ARC    = 1; // ARC LENGTH METHOD
                   // ARC == 0 :: Crisfield
                   // ARC == 1 :: Simo ONLY ARC LENG METHOD FOR PARALELL
  arc->dALMAX = 0.0;
  arc->ITT    = 0;
  arc->DAL    = 0.0;
  return err;
}

/// construct arc length variable object
/// create memory space for member arrays and structs
///
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \param[in] fv an object containing all field variables
/// \param[in] com an object for communication
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int construct_arc_length_variable(ARC_LENGTH_VARIABLES *arc,
                                  FieldVariables *fv,
                                  const CommunicationStructure *com) {
  int err = 0;
  int myrank  = com->rank;
  arc->D_R    = aloc1(fv->ndofd);
  arc->U      = aloc1(fv->ndofd);
  arc->DK     = aloc1(fv->ndofd);
  arc->dR     = aloc1(fv->ndofd);
  arc->BS_d_r = aloc1(com->DomDof[myrank]);
  arc->BS_D_R = aloc1(com->DomDof[myrank]);
  arc->BS_rr  = aloc1(com->DomDof[myrank]);
  arc->BS_R   = aloc1(com->DomDof[myrank]);
  arc->BS_U   = aloc1(com->DomDof[myrank]);
  arc->BS_DK  = aloc1(com->DomDof[myrank]);
  arc->BS_dR  = aloc1(com->DomDof[myrank]);
  return err;
}

/// destruct arc length variable object
/// free memory spaces for member arrays and structs
///
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \return non-zero on internal error
int destruct_arc_length_variable(ARC_LENGTH_VARIABLES *arc)
{
  int err = 0;
  if(NULL != arc->D_R)    free(arc->D_R);
  if(NULL != arc->U)      free(arc->U);
  if(NULL != arc->DK)     free(arc->DK);
  if(NULL != arc->dR)     free(arc->dR);
  if(NULL != arc->BS_d_r) free(arc->BS_d_r);
  if(NULL != arc->BS_D_R) free(arc->BS_D_R);
  if(NULL != arc->BS_rr)  free(arc->BS_rr);
  if(NULL != arc->BS_R)   free(arc->BS_R);
  if(NULL != arc->BS_U)   free(arc->BS_U);
  if(NULL != arc->BS_DK)  free(arc->BS_DK);
  if(NULL != arc->BS_dR)  free(arc->BS_dR);
  err += arc_length_variable_initialization(arc);
  return err;
}

/// Perform arc length analysis.
/// It calls the lagacy Arc_length function which requires many function arguments.
/// Later, this function needs to be combined with Arc_length function.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] time_steps object for time stepping
/// \param[in] com object for communication
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in, out] arc an object for Arc length scheme, cotains variables related to Arc length
/// \param[in] com object for communication
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return load multiplier
double Multiphysics_Arc_length(Grid *grid,
                               MaterialProperty *mat,
                               FieldVariables *fv,
                               pgfem3d::Solver *sol,
                               LoadingSteps *load,
                               const CommunicationStructure *com,
                               TimeStepping *time_steps,
                               CRPL *crpl,
                               const double VVolume,
                               const PGFem3D_opt *opts,
                               const Multiphysics& mp,
                               const int mp_id)
{
  ARC_LENGTH_VARIABLES *arc = sol->arc;
  double dALMAX = (time_steps->dt_np1) / (arc->dt0) * (arc->dALMAX);
  char out_dat[500];
  sprintf(out_dat, "%s/%s", opts->opath, opts->ofname);

  double dts[2];
  dts[0] = dts[1] = time_steps->dt_np1;
  double t = 0.0;
  double nor, nor1, nor2, dlm, dlm0{}, DLM, DET=0.0, dAL;
  double DT, DDLM, ddlm, ERROR, LS1, gama, pdt, tmp, nor3;
  long iter, INFO, STEP, DIV, ST, GAMA, OME, FI, gam, TYPE, GInfo;
  int ART;
  int myrank = com->rank;
  const char *error[] = {
    "inf",
    "-inf",
    "nan"
  };

  char str1[500], jmeno[600];
  FILE *out;
  struct rusage usage;

  /* damage substep criteria */
  const double max_damage_per_step = 0.05;
  double max_damage = 0.0;
  double alpha = 0.0;

  /* dissipation */
  double dissipation = 0.0;

  /* BlockSolve 95 */
  double BS_nor = 0.0;
  int BS_iter;

  /* TEMPORARY bounding element compile testing */
  /* int grid->n_be = 0; */
  /* BOUNDING_ELEMENT *b_elems = NULL; */

  TYPE = 1;
  pdt = time_steps->dt_np1;

  switch (opts->analysis_type) {
   case STABILIZED:
   case MINI:
   case MINI_3F:
    fv->ndofn = 4;
    break;
   default:
    break;
  }

  /* SUBDIVISION */
  DT = DDLM = ddlm = dlm = 0.0;
  DIV = ST = GAMA = OME = INFO = ART = gam = FI = iter = 0;
  STEP = 1;

 rest:
  dlm0 = subdiv_arc(INFO, &(time_steps->dt_np1), arc->dt0, &STEP, &DIV,
                    time_steps->tim, time_steps->times, &ST, grid->ne,
                    fv->ndofd, fv->npres, grid->element, crpl, fv->eps, fv->sig,
                    load->sups[mp_id], load->sup_defl[mp_id], fv->dd_u, fv->d_u,
                    arc->D_R, fv->f_defl, fv->f, &GAMA, &DT, &OME, opts->stab,
                    arc->dAL0, arc->DAL, dALMAX, sol->nor_min, dlm0,
                    &(arc->ITT), iter, sol->iter_max, TYPE, com,
                    opts->analysis_type);

  dts[DT_NP1] = time_steps->dt_np1;

  ART = (periodic == 1) ? 1 : 0;
  if (TYPE == 1) dAL = dlm0;
  /* } */

  dlm0 *= (arc->DLM0) / fabs(arc->DLM0);

  INFO = 0;
  while (STEP > DIV) {
    if ((STEP > 1 || ST == 1) && com->rank == 0) {
      PGFEM_printf("\nSTEP = %ld :: NS =  %ld || Time %f | dt = %10.10f\n",
                   DIV, STEP, time_steps->times[time_steps->tim+1],
                   time_steps->dt_np1);
    }

    if (myrank == 0) {
      PGFEM_printf("dlm0 = %12.12e\n",dlm0);
    }

    iter = 0;
    if (periodic == 1) {
      if (TYPE == 1) {
        nor = dlm0;
        dlm0 = 0.0;
      }
      /* macroscopic load */
      macroscopic_load_AL(fv->eps[0].type, arc->lm + dlm0, fv->eps);
    }/* end periodic */

    if (periodic == 1) {
      nulld(fv->R, fv->ndofd);
    }

    assert(opts->solverpackage == HYPRE || opts->solverpackage == TRILINOS);

    /* Null the matrix */
    sol->system->zero();

    stiffmat_fd_MP(grid, mat, fv, sol, load, com, crpl, opts, mp,
                   mp_id, time_steps->dt_np1, iter);

    /* Assemble the matrix */
    sol->system->assemble();

    if (periodic == 1 && TYPE == 1) {
      dlm0 = nor;
    }
    /* Transform LOCAL load vector to GLOBAL */
    if (periodic == 1) {
      LToG(fv->R, arc->BS_R, fv->ndofd, com);
    }

    /*** SOLVE THE SYSTEM FOR DIRECTION ***/
    {
      SOLVER_INFO s_info;
      sol->system->solveSystem(opts, arc->BS_R, arc->BS_rr, time_steps->tim,
                               iter, com->DomDof, &s_info);
      if(myrank == 0) {
        solve_system_check_error(stdout,s_info);
      }
      BS_nor = s_info.res_norm;
      BS_iter = s_info.n_iter;
    }

    if (ARC_DEBUG && myrank == 0) {
      PGFEM_printf("Completed solve I!\n");
    }

    if (BS_nor > 100.*(sol->err) || BS_iter < 0) {
      INFO = 1;
      if (myrank == 0) {
        PGFEM_printf("ERROR in the BSpar_solve : nor = %8.8e || iter = %d\n",
                     BS_nor, BS_iter);
      }
      goto rest;
    }

    sprintf (str1, "%f", BS_nor);
    for (int N = 0; N < 3; N++) {
      int M = 10;
      M = strcmp(error[N], str1);
      if (M == 0) {
        if (myrank == 0) {
          PGFEM_printf("ERROR in the BSpar_solve : nor = %s\n", error[N]);
        }
        INFO = 1;
        goto rest;
      }
    }

    if (ARC_DEBUG && myrank == 0) {
      PGFEM_printf("Starting GToL... ");
    }

    /* Transform GLOBAL displacement vector to LOCAL */
    GToL(arc->BS_rr, fv->dd_u, fv->ndofd, com);

    /* Null R for periodic */
    if (periodic == 1) {
      nulld(fv->R, fv->ndofd);
      nulld(arc->BS_R, com->DomDof[myrank]);
    }

    /* Scaling parameter DK */
    if (time_steps->tim == 0 && iter == 0) {
      if (arc->ARC == 0) { /* diag_KK (fv->ndofd,k,com->Ap,com->Ai,arc->DK); */
        for (int i = 0, e = com->DomDof[myrank]; i < e; ++i) {
          arc->BS_DK[i] = 1.0;
        }
      }
      if (arc->ARC == 1){ /* diag_KK (fv->ndofd,k,com->Ap,com->Ai,arc->DK); */
        for (int i = 0, e = com->DomDof[myrank]; i < e; ++i) {
          arc->BS_DK[i] = 1.0;
        }
      }
    }

    /* First load multiplier */
    if (arc->ARC == 0) {
      if (TYPE == 0) {
        dAL = d_ALM2(fv->ndofd, arc->BS_rr, arc->BS_R, arc->BS_DK, dlm0);
      }
      dlm = d_lam_ALM2(fv->ndofd, arc->BS_rr, arc->BS_R, arc->BS_DK, dAL, DET,
                       arc->DET0, arc->DLM0, sol->nor_min, arc->dR);
    }
    if (arc->ARC == 1) {
      if (TYPE == 0) {
        dAL = d_ALM4(fv->ndofd, arc->BS_rr, arc->BS_DK, dlm0, com->DomDof, com);
      }
      dlm = d_lam_ALM4(fv->ndofd, arc->BS_rr, arc->BS_DK, arc->BS_dR, dAL,
                       com->DomDof, com);
    }

    /* First load multiplier */
    dlm0 = dlm;
    if (myrank == 0) {
      PGFEM_printf("dAL = %12.15e :: dALmax = %12.15e || dlm0  = %12.15e\n\n",
                   dAL, dALMAX, dlm0);
    }

    /* Limiting the arc length */
    if (dAL > dALMAX) {
      if (myrank == 0) {
        PGFEM_printf("*** Arc length too large: Restart with smaller arc ***\n");
      }
      INFO = 1;
      goto rest;
    }

    /* macroscopic load */
    if (periodic == 1) {
      macroscopic_load_AL(fv->eps[0].type, arc->lm+dlm, fv->eps);
    }

    /* First displacement increment dr_0 */
    for (int i = 0, e = fv->ndofd; i < e; ++i) {
      fv->d_u[i] = dlm * fv->dd_u[i];
    }

    /* INCREMENT: L -> G */
    LToG(fv->d_u, arc->BS_d_r, fv->ndofd, com);

    /* Pressure and volume change THETA */
    nulld(fv->f_u, fv->ndofd);
    switch (opts->analysis_type) {
     case FS_CRPL:
     case FINITE_STRAIN:
      press_theta(grid->ne, fv->ndofn, fv->npres, grid->element, grid->node,
                  fv->f_u, fv->d_u, load->sups[mp_id], mat->matgeom,
                  mat->hommat, fv->eps, fv->sig, iter, sol->nor_min,
                  time_steps->dt_np1, crpl, opts, mp_id);
      break;
     case MINI:
      MINI_update_bubble(grid->element, grid->ne, grid->node, fv->ndofn,
                         load->sups[mp_id], fv->eps, fv->sig, mat->hommat,
                         fv->f_u, fv->d_u, iter, mp_id);
      break;
     case MINI_3F:
      MINI_3f_update_bubble(grid->element, grid->ne, grid->node, fv->ndofn,
                            load->sups[mp_id], fv->eps, fv->sig, mat->hommat,
                            fv->f_u, fv->d_u, iter, mp_id);
      break;
     default:
      break;
    }

    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL) {
      INFO = integration_alg(grid->ne, fv->ndofn, fv->ndofd, fv->npres, crpl,
                             grid->element, grid->node, fv->f_u, fv->d_u,
                             load->sups[mp_id], mat->matgeom, mat->hommat,
                             fv->eps, fv->sig, time_steps->tim, iter,
                             time_steps->dt_np1, sol->nor_min, STEP, 0, opts,
                             mp_id);

      /* Gather INFO from all domains */
      com->net->allreduce(&INFO, &GInfo, 1, NET_DT_LONG, NET_OP_BOR, com->comm);

      if (GInfo == 1) {
        INFO = 1;
        goto rest;
      }
    }

    vol_damage_int_alg(grid->ne, fv->ndofn, fv->d_u, fv->u_np1, grid->element,
                       grid->node, mat->hommat, load->sups[mp_id],
                       time_steps->dt_np1, iter, com, fv->eps, fv->sig,
                       &max_damage, &dissipation, opts->analysis_type, mp_id);

    /* Residuals */
    fd_residuals_MP(grid, mat, fv, sol, load, crpl, com, opts, mp, mp_id,
                    t, dts, 0);

    /* Compute Euclidian norm */
    for (int i = 0, e = fv->ndofd; i < e; ++i) {
      fv->f[i] = (arc->lm + dlm) * fv->R[i] - fv->f_u[i];
    }

    /* fv->BS_f : L->G */
    LToG(fv->f, fv->BS_f, fv->ndofd, com);

    nor = ss(fv->BS_f, fv->BS_f, com->DomDof[myrank]);

    /* Gather nor from each domain */
    com->net->allreduce(&nor, &nor3, 1, NET_DT_DOUBLE, NET_OP_SUM, com->comm);
    nor = nor2 = sqrt(nor3);
    nor1 = nor;

    sprintf (str1,"%f",nor);
    for (int N = 0; N < 3; N++) {
      int M = 10;
      M = strcmp(error[N], str1);
      if (M == 0) {
        if (myrank == 0) {
          PGFEM_printf("ERROR in the algorithm : nor = %s\n", error[N]);
        }
        INFO = 1;
        goto rest;
      }
    }

    if (time_steps->tim == 0 && iter == 0) {
      fv->NORM = nor1;
    }

    /* THIS IS THE GLOBAL-LOCAL TOLERANCE */
    nor1 = fv->NORM;

    /* Normalize norm */
    if (nor < 1.e-10) {
      if (myrank == 0) {
        PGFEM_printf ("Scale your units : Close to the computer precision\n");
      }
      nor = 0.0;
    }
    else {
      nor /= nor1;
    }

    /* For very bad convergence problems */
    if (fabs(dlm) < sol->nor_min) {
      nor = nor2;
      ERROR = 10.0 * (sol->nor_min);
      if (myrank == 0) {
        PGFEM_printf ("Reducing and changing error || ERROR = %12.12e\n",ERROR);
      }
    }
    else {
      ERROR = sol->nor_min;
    }

    if (myrank == 0) {
      getrusage(RUSAGE_SELF, &usage);
      PGFEM_printf("(%ld) IT = %d : R = %8.8e :: ||f||/||f0|| ="
                   " [%8.8e] || [%8.8e] :: S %ld.%ld, U %ld.%ld\n",
                   iter, BS_iter, BS_nor, nor, nor2, usage.ru_stime.tv_sec,
                   usage.ru_stime.tv_usec, usage.ru_utime.tv_sec,
                   usage.ru_utime.tv_usec);
    }

    /*************/
    /* iter =  1 */
    /*************/

    iter = 1;
    DLM = 0.0;
    while (nor > ERROR) {

      /* Null the residual vector */
      nulld(fv->f_u, fv->ndofd);

      assert(opts->solverpackage == HYPRE || opts->solverpackage == TRILINOS);
      /* Null the matrix */
      sol->system->zero();

      stiffmat_fd_MP(grid, mat, fv, sol, load, com, crpl, opts, mp,
                     mp_id, time_steps->dt_np1, iter);

      /* Assemble the matrix */
      sol->system->assemble();

      /* Transform unequibriated force: L -> G */
      if (periodic != 1){
        for (int i = 0, e = com->DomDof[myrank]; i < e; ++i) {
          fv->BS_f_u[i] = arc->BS_R[i];
        }
      }
      else {
        LToG(fv->f_u, fv->BS_f_u, fv->ndofd, com);
      }

      /*** SOLVE THE SYSTEM EQUATIONS ***/
      {
        SOLVER_INFO s_info;
        sol->system->solveSystem(opts, fv->BS_f_u, arc->BS_rr, time_steps->tim,
                                 iter, com->DomDof, &s_info);
        if (myrank == 0) {
          solve_system_check_error(stdout,s_info);
        }
        BS_nor = s_info.res_norm;
        BS_iter = s_info.n_iter;
      }

      /* Check for correct solution */
      if (BS_nor > 100.*(sol->err) || BS_iter < 0) {
        INFO = 1;
        if (myrank == 0) {
          PGFEM_printf("ERROR in the BSpar_solve : nor = %8.8e || iter = %d\n",
                       BS_nor, BS_iter);
        }
        goto rest;
      }

      sprintf (str1, "%f", BS_nor);
      for (int N = 0; N < 3; N++) {
        int M = 10;
        M = strcmp(error[N], str1);
        if (M == 0) {
          if (myrank == 0) {
            PGFEM_printf("ERROR in the BSpar_solve : nor = %s\n", error[N]);
          }
          INFO = 1;
          goto rest;
        }
      }


      if (myrank == 0) {
        getrusage(RUSAGE_SELF, &usage);
        PGFEM_printf("(%ld.5) Intermediate Solve Info: IT = %d : R = %8.8e ::"
                     " S %ld.%ld, U %ld.%ld\n", (iter-1), BS_iter, BS_nor,
                     usage.ru_stime.tv_sec, usage.ru_stime.tv_usec,
                     usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);
      }

      /* initial disp. G-> L */
      GToL(arc->BS_rr, fv->dd_u, fv->ndofd, com);

      /* Solve i-th incremnt */
      {
        SOLVER_INFO s_info;
        sol->system->solveSystem(opts, fv->BS_f, arc->BS_U, time_steps->tim,
                                 iter, com->DomDof, &s_info);
        if(myrank == 0) {
          solve_system_check_error(stdout, s_info);
        }
        BS_nor = s_info.res_norm;
        BS_iter = s_info.n_iter;
      }

      /* Check for correct solution */
      if (BS_nor > 100.*(sol->err) || BS_iter < 0) {
        INFO = 1;
        if (myrank == 0) {
          PGFEM_printf("ERROR in the BSpar_solve : nor = %8.8e || iter = %d\n",
                       BS_nor, BS_iter);
        }
        goto rest;
      }

      sprintf(str1, "%f", BS_nor);
      for (int N = 0; N < 3; N++) {
        int M = 10;
        M = strcmp(error[N], str1);
        if (M == 0) {
          if (myrank == 0) {
            PGFEM_printf("ERROR in the BSpar_solve : nor = %s\n", error[N]);
          }
          INFO = 1;
          goto rest;
        }
      }

      /* increment disp. G-> L */
      GToL(arc->BS_U, arc->U, fv->ndofd, com);

      /* dlam */
      if (arc->ARC == 0) {
        if (arc->AT == 1 || ART == 0) {
          DLM = D_lam_ALM(fv->ndofd, arc->BS_rr, arc->BS_d_r, arc->BS_U,
                          arc->BS_R, arc->BS_DK, dlm, dAL, com->DomDof,
                          com);
        }
        else {
          DLM = D_lam_ALM2_MP(grid, mat, fv, sol, load, com, crpl,
                              opts, mp, mp_id, dlm, dAL, time_steps->dt_np1);
        }
      }
      if (arc->ARC == 1) {
        DLM = D_lam_ALM4(fv->ndofd, arc->BS_rr, arc->BS_d_r, arc->BS_U,
                         arc->BS_DK, dlm, dAL, com->DomDof, com);
      }
      sprintf(str1, "%f", DLM);
      for (int N = 0; N < 3; N++) {
        int M = 10;
        M = strcmp(error[N], str1);
        if (M == 0) {
          if (myrank == 0) {
            PGFEM_printf("Complex root in ARC-LENGHT method\n");
          }
          INFO = 1;
          goto rest;
        }
      }

      /* Periodic loading */
      if (periodic == 1) {
        macroscopic_load_AL(fv->eps[0].type, arc->lm + dlm + DLM, fv->eps);
      }

      /* Update deformation */
      for (int i = 0, e = fv->ndofd; i < e; ++i) {
        arc->D_R[i] = arc->U[i] + DLM * fv->dd_u[i];
      }

      /* LINE SEARCH */
      tmp = ss(fv->BS_f, fv->BS_f, com->DomDof[myrank]);
      com->net->allreduce(&tmp, &LS1, 1, NET_DT_DOUBLE, NET_OP_SUM, com->comm);
      LS1 *= 1./2;

      /* Pressure and volume change THETA */
      switch (opts->analysis_type) {
       case FS_CRPL:
       case FINITE_STRAIN:
        press_theta(grid->ne, fv->ndofn, fv->npres, grid->element, grid->node,
                    fv->d_u, arc->D_R, load->sups[mp_id], mat->matgeom,
                    mat->hommat, fv->eps, fv->sig, iter, sol->nor_min,
                    time_steps->dt_np1, crpl, opts, mp_id);
        break;
       case MINI:
        MINI_update_bubble(grid->element, grid->ne, grid->node, fv->ndofn,
                           load->sups[mp_id], fv->eps, fv->sig, mat->hommat,
                           fv->d_u, arc->D_R, iter, mp_id);
        break;
       case MINI_3F:
        MINI_3f_update_bubble(grid->element, grid->ne, grid->node, fv->ndofn,
                              load->sups[mp_id], fv->eps, fv->sig, mat->hommat,
                              fv->d_u, arc->D_R, iter, mp_id);
        break;
       default:
        break;
      }

      /*************************/
      /* INTEGRATION ALGORITHM */
      /*************************/
      if (opts->analysis_type == FS_CRPL) {
        INFO = integration_alg(grid->ne, fv->ndofn, fv->ndofd, fv->npres, crpl,
                               grid->element, grid->node, fv->d_u, arc->D_R,
                               load->sups[mp_id], mat->matgeom, mat->hommat,
                               fv->eps, fv->sig, time_steps->tim, iter,
                               time_steps->dt_np1, sol->nor_min, STEP, 0, opts, mp_id);

        /* Gather INFO from all domains */
        com->net->allreduce(&INFO, &GInfo, 1, NET_DT_LONG, NET_OP_BOR, com->comm);
        if (GInfo == 1) {
          INFO = 1;
          goto rest;
        }
      }

      /* Update deformations */
      for (int i = 0, e = fv->ndofd; i < e; ++i) {
        fv->f[i] = fv->d_u[i] + arc->D_R[i];
        fv->f_u[i] = 0.0;
      }

      vol_damage_int_alg(grid->ne, fv->ndofn, fv->f, fv->u_np1, grid->element,
                         grid->node, mat->hommat, load->sups[mp_id],
                         time_steps->dt_np1, iter, com, fv->eps, fv->sig,
                         &max_damage, &dissipation, opts->analysis_type, mp_id);

      /* Residuals */
      fd_residuals_MP(grid, mat, fv, sol, load, crpl, com, opts, mp, mp_id,
                      t, dts, 1);

      /* Compute Euclidean norm */
      for (int i = 0, e = fv->ndofd; i < e; ++i) {
        fv->f[i] = (arc->lm + dlm + DLM) * fv->R[i] - fv->f_u[i];
      }

      /* Residuals L -> G */
      LToG(fv->f, fv->BS_f, fv->ndofd, com);

      nor = ss(fv->BS_f, fv->BS_f, com->DomDof[myrank]);
      com->net->allreduce(&nor, &tmp, 1, NET_DT_DOUBLE, NET_OP_SUM, com->comm);
      nor = nor2 = sqrt (tmp);
      nor /= nor1;

      sprintf (str1,"%f",nor);
      for (int N = 0; N < 3; N++) {
        int M = 10;
        M = strcmp(error[N], str1);
        if (M == 0) {
          if (myrank == 0) {
            PGFEM_printf("ERROR in the algorithm : nor = %s\n", error[N]);
          }
          INFO = 1;
          goto rest;
        }
      }

      /* LINE SEARCH */
      if (ART == 0) {
        INFO = ALINE_S3_MP(grid, mat, fv, sol, load, com, crpl, opts,
                           mp, dts, mp_id, &nor, &nor2, nor1, LS1, iter,
                           &max_damage, &dissipation, time_steps->tim, STEP,
                           &DLM, &gama, dlm, dAL);

        /* Gather INFO from all domains */
        com->net->allreduce(&INFO, &GInfo, 1, NET_DT_LONG, NET_OP_BOR, com->comm);
        if (GInfo == 1) {
          if (myrank == 0) {
            PGFEM_printf ("Error in the Line Search 3 algorithm\n");
          }
          INFO = 1;
          goto rest;
        }
      }

      /* For very bad convergence problems */
      if (fabs(dlm0) < sol->nor_min) {
        nor = nor2;
      }

      /* Load multiplied update */
      dlm += DLM;

      /* Deformation update */
      for (int i = 0, e = fv->ndofd; i < e; ++i) {
        fv->d_u[i] += arc->D_R[i];
      }

      /* Total increment: L -> G */
      LToG(fv->d_u, arc->BS_d_r, fv->ndofd, com);

      if (myrank == 0) {
        getrusage(RUSAGE_SELF,&usage);
        PGFEM_printf("(%ld) IT = %d : R = %8.8e :: "
                     "||f||/||f0|| = [%8.8e] || [%8.8e] ::"
                     " S %ld.%ld, U %ld.%ld\n",
                     iter, BS_iter, BS_nor, nor, nor2,
                     usage.ru_stime.tv_sec, usage.ru_stime.tv_usec,
                     usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);
      }

      /* Max number of iterations restart */
      if (iter > (sol->iter_max-1) && nor > ERROR) {
        if (nor > 20.*ERROR || iter > (sol->iter_max + 2)) {
          if (nor < 5.*ERROR) {
            if (myrank == 0) {
              PGFEM_printf("I will take it\n");
            }
            nor = 0.0;
          }
          else {
            if (myrank == 0) {
              PGFEM_printf("Error in the iteration : iter > iter_max\n");
            }
            INFO = 1;
            goto rest;
          }
        }
      } /* end iter > iter_max */
      ++iter;
    }/* while nor > sol->nor_min */

    /* before increment after convergence, check max damage */
    alpha = max_damage / max_damage_per_step;
    com->net->allreduce(NET_IN_PLACE, &alpha, 1, NET_DT_DOUBLE, NET_OP_MAX, com->comm);
    if (alpha > 1.0) {
      if (myrank == 0) {
        PGFEM_printf("Damage value (%f) is greater than max. damage/step (%f).\n"
                     "Subdividing to maintain accuracy of damage law.\n",
                     alpha * max_damage_per_step, max_damage_per_step);
      }
      INFO = 1;
      goto rest;
    }

    ST = GAMA = gam = 0;
    ART = (periodic == 1) ? 1 : 0;

    /* macroscopic load */
    if (periodic == 1) {
      macroscopic_load_AL(fv->eps[0].type, arc->lm + dlm, fv->eps);
    }

    /* increment coheisve elements */
    if (opts->cohesive) {
      increment_cohesive_elements(grid->nce, grid->coel, &(fv->pores),
                                  grid->node, load->sups[mp_id], fv->d_u,
                                  mp_id);
    }

    /* Finite deformations increment */
    switch (opts->analysis_type) {
     case FS_CRPL:
     case FINITE_STRAIN:
      fd_increment(grid->ne, grid->nn, fv->ndofn, fv->npres, mat->matgeom,
                   mat->hommat, grid->element, grid->node, load->sups[mp_id],
                   fv->eps, fv->sig, fv->d_u, fv->u_np1, sol->nor_min, crpl,
                   time_steps->dt_np1, grid->nce, grid->coel, &(fv->pores),
                   com, VVolume, opts, mp_id);
      break;
     case STABILIZED:
      st_increment(grid->ne, grid->nn, fv->ndofn, fv->ndofd, mat->matgeom,
                   mat->hommat, grid->element, grid->node, load->sups[mp_id],
                   fv->eps, fv->sig, fv->d_u, fv->u_np1, sol->nor_min,
                   opts->stab, time_steps->dt_np1, grid->nce, grid->coel,
                   &(fv->pores), com, opts->cohesive, mp_id);
      break;
     case MINI:
      MINI_increment(grid->element, grid->ne, grid->node, grid->nn, fv->ndofn,
                     load->sups[mp_id], fv->eps, fv->sig, mat->hommat, fv->d_u,
                     com, mp_id);
      break;
     case MINI_3F:
      MINI_3f_increment(grid->element, grid->ne, grid->node, grid->nn,
                        fv->ndofn, load->sups[mp_id], fv->eps, fv->sig,
                        mat->hommat, fv->d_u, com, mp_id);
      break;
     case DISP:
      DISP_increment(grid->element, grid->ne, grid->node, grid->nn, fv->ndofn,
                     load->sups[mp_id], fv->eps, fv->sig, mat->hommat, fv->d_u,
                     fv->u_np1, com, mp_id);
      break;
     default:
      break;
    }

    /* Add deformation increment into displacement vector */
    vvplus(fv->u_np1, fv->d_u, fv->ndofd);

    /* Null prescribed increment deformation */
    for (int i = 0, e = load->sups[mp_id]->npd; i < e; ++i) {
      load->sups[mp_id]->defl_d[i] = 0.0;
    }

    for (int i = 0, e = fv->ndofd; i < e; ++i) {
      arc->dR[i] = fv->d_u[i];
      fv->d_u[i] = arc->D_R[i] = fv->dd_u[i] = fv->f_defl[i] = fv->f[i] = arc->U[i] = 0.0;
    }

    /* Tranform arc->dR : L -> G */
    LToG(arc->dR, arc->BS_dR, fv->ndofd, com);

    /* Prevent circling */
    nor = (time_steps->tim == 0) ? 1.0 : (sqrt((dlm - fabs(arc->DLM)) *
                                               (dlm - fabs(arc->DLM))) /
                                          sqrt((arc->DLM) * (arc->DLM)));
    arc->AT = (nor < 1.e-2 && (arc->DLM) < 0.0 && dlm > 0.0) ? 1 : 0;

    /* Load multiplicator and Determinant */
    DDLM += dlm;
    arc->lm += dlm;
    (arc->DLM) = dlm;
    ddlm += fabs (dlm);
    arc->DLM0 = dlm0;
    arc->DET0 = DET;
    arc->ITT = iter;
    arc->DAL = dAL;

    if (myrank == 0) {
      PGFEM_printf ("Inelastic step [%ld] :: STEP = %ld :: NS =  %ld"
                    " || lm = %12.12f : DDLM = %12.12f : dlm = %12.12f\n\n",
                    time_steps->tim, DIV, STEP, arc->lm, DDLM, dlm);
    }

    dlm = 0.0;
    if (periodic == 1 && myrank == 0) {
      PGFEM_printf("\nThe total F\n");
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          PGFEM_printf("%12.12f  ", fv->eps[0].F[i][j]);
        }
        PGFEM_printf("\n");
      }
      PGFEM_printf("\n");
    }
    if (arc->AT == 1 && myrank == 0) {
      PGFEM_printf ("*** Trying to prevent cycling ***\n\n");
    }

    /************* TEST THE UPDATE FROM N TO N+1  *************/
    if(PFEM_DEBUG || PFEM_DEBUG_ALL || ARC_UPDATE) {
      for (int i = 0, e = fv->ndofd; i < e; ++i) {
        fv->f_u[i] = 0.0;
        fv->d_u[i] = 0.0;
      }
      fd_residuals_MP(grid, mat, fv, sol, load, crpl, com, opts, mp, mp_id,
                      t, dts, 0);
      for (int i = 0, e = fv->ndofd; i < e; ++i) {
        fv->f[i] = (arc->lm)*(fv->R[i]) - fv->f_u[i];
      }

      LToG(fv->f, fv->BS_f, fv->ndofd, com);
      nor = ss(fv->BS_f, fv->BS_f, com->DomDof[myrank]);
      com->net->allreduce(&nor, &tmp, 1, NET_DT_DOUBLE, NET_OP_SUM, com->comm);
      nor = sqrt(tmp);

      if (myrank == 0) {
        PGFEM_printf("NORM NORM = %12.12e || %12.12f\n", nor, arc->lm);
      }
    } /* End debug */
    /************* TEST THE UPDATE FROM N TO N+1  *************/

    /* Calculating equivalent Mises stresses and strains vectors */
    Mises(grid->ne, fv->sig, fv->eps, opts->analysis_type);

    /********************************************************************/
    /***  PRINT INTERMIDIATE OUTPUT FILES FOR CONSTANT LOAD LEVEL     ***/
    /********************************************************************/

    if (ddlm > pdt * (arc->dAL0) / (arc->dt0) && STEP > (DIV + 1)) {

      if (time_steps->print[time_steps->tim] == 1 &&
          opts->vis_format == VIS_ELIXIR) {
        sprintf(jmeno, "%s_%d.out%ld_%ld", out_dat, myrank, time_steps->tim, FI);

        if ((out = fopen(jmeno,"w")) == nullptr ) {
          PGFEM_printf("Output file is not possible to open\n");
          PGFEM_printf("Check the output file and run program again\n");
          return (0);
        }

        /*logo (out);*/

        PGFEM_fprintf(out, "\n");
        PGFEM_fprintf(out, "FINITE DEFORMA. + CRYSTAL PLASTICITY ANALYSIS  : Step - %ld, Time - %f\n", time_steps->tim, time_steps->times[time_steps->tim+1]);
        PGFEM_fprintf(out, "Number of nodes                                : %ld\n", grid->nn);
        PGFEM_fprintf(out, "Number of elements                             : %ld\n", grid->ne);
        PGFEM_fprintf(out, "Number of equations                            : %ld\n", fv->ndofd);
        PGFEM_fprintf(out, "Number of elements in the matrix - SPARSE      : %d\n", com->Ap[fv->ndofd]);
        PGFEM_fprintf(out, "Load multiplier level                          : %f\n", arc->lm);
        PGFEM_fprintf(out, "\n");
        getrusage(RUSAGE_SELF, &usage);
        PGFEM_fprintf(out, "Time of solution of the system of equations  -  System %ld.%ld, User %ld.%ld\n",
                      usage.ru_stime.tv_sec, usage.ru_stime.tv_usec,
                      usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);

        /* Print MACRO fileds */
        if (periodic == 1) {
          macro_fields_out(out, fv->eps, opts);
        }

        /* Print deformations to output file */
        deform(out, grid->node, grid->element, grid->nn, grid->ne, fv->ndofd,
               load->sups[mp_id], fv->u_np1);
        /* Print stress to output file */
        stress_out(out, grid->ne, grid->nn, grid->element, fv->sig, fv->sig_n,
                   opts->smoothing);
        /* Print strain to output file */
        strain_out(out, grid->ne, grid->element, fv->eps, opts);
        /* Print Fe */
        deform_grad_out(out, grid->ne, grid->element, fv->eps);
        /* Print cohesive elements */
        if (opts->cohesive == 1) {
          cohesive_out(out, grid->nce, grid->coel);
        }

        fclose(out);

        sprintf(jmeno, "%s_.elx%ld_%ld", out_dat, time_steps->tim, FI);
        FI++;

        /* Print to elix file */
        elixir(jmeno, grid->nn, grid->ne, fv->ndofn, grid->node, grid->element,
               load->sups[mp_id], fv->u_np1, fv->sig, fv->sig_n, fv->eps,
               opts->smoothing, grid->nce, grid->coel,
               /*,nge,geel,ngn,gnod*/
               opts);
      }
      ddlm = 0.0;
    }/* end ddlm > dAL0 */

    DIV++;
    if (STEP > 2 && DIV == 2)
      goto rest;
  }/*end SUBDIVISION */

  return (DDLM);
}


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
double Arc_length_multiscale(MultiscaleCommon *c,
                             MULTISCALE_SOLUTION *s,
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
                             double *sup_defl)
{
  // initialize and define multiphysics
  Multiphysics mp;
  int id = MULTIPHYSICS_MECHANICAL;
  int ndim = c->ndofn;
  int write_no = 0;

  vector<int> coupled_ids;
  char *physicsname = (char *) malloc(sizeof(char)*1024);
  {
    coupled_ids.push_back(0);
    sprintf(physicsname, "Mechanical");

    mp.physicsno      = 1;
    mp.physicsname    = &physicsname;
    mp.physics_ids    = &id;
    mp.ndim           = &ndim;
    mp.write_no       = &write_no;
    mp.coupled_ids.push_back(coupled_ids);
    mp.total_write_no = 0;
  }

  // initialize and define mesh object
  Grid grid;
  grid_initialization(&grid);
  {
    grid.ne          = c->ne;
    grid.nn          = c->nn;
    grid.element     = c->elem;
    grid.b_elems     = NULL;
    grid.node        = c->node;
    grid.nce         = c->nce;
    grid.coel        = c->coel;
    grid.n_be        = 0;
  }

  // initialize and define field variables
  FieldVariables fv;
  {
    field_varialbe_initialization(&fv);
    fv.ndofn  = c->ndofn;
    fv.ndofd  = c->ndofd;
    fv.npres  = c->npres;
    fv.sig    = s->sig_e;
    fv.eps    = s->eps;
    fv.u_np1  = s->r;
    fv.f      = s->f;
    fv.d_u    = s->d_r;
    fv.dd_u   = s->rr;
    fv.R      = s->R;
    fv.f_defl = s->f_defl;
    fv.RR     = s->RR;
    fv.f_u    = s->f_u;
    fv.RRn    = s->RRn;
    fv.pores  = *pores;
    //fv.BS_x   = s->BS_x;
    fv.BS_f   = s->BS_f;
    fv.BS_RR  = s->BS_RR;
    fv.BS_f_u = s->BS_f_u;
    fv.sig_n  = s->sig_n;
    fv.NORM   = s->NORM;
  }

  ARC_LENGTH_VARIABLES arc;
  {
    arc_length_variable_initialization(&arc);
    arc.dt0     = dt0;
    arc.D_R     = s->D_R;
    arc.U       = s->U;
    arc.DK      = s->DK;
    arc.dR      = s->dR;
    arc.BS_d_r  = s->BS_d_r;
    arc.BS_D_R  = s->BS_D_R;
    arc.BS_rr   = s->BS_rr;
    arc.BS_R    = s->BS_R;
    arc.BS_U    = s->BS_U;
    arc.BS_DK   = s->BS_DK;
    arc.BS_dR   = s->BS_dR;
    arc.lm      = lm;
    arc.dAL0    = solver_file->nonlin_method_opts[0];
    arc.DET0    = *DET0;
    arc.DLM0    = *DLM0;
    arc.DLM     = *DLM;
    arc.AT      = *AT;
    arc.ARC     = ARC;
    arc.ITT     = *ITT;
    arc.DAL     = *DAL;
  }

  /// initialize and define iterative solver object
  pgfem3d::Solver sol{};
  {
    sol.nor_min    = solver_file->nonlin_tol;
    sol.FNR        = solver_file->nonlin_method;
    sol.iter_max   = solver_file->max_nonlin_iter;
    sol.system     = c->SOLVER;
    sol.err        = c->lim_zero;
    sol.microscale = nullptr;
  }

  sol.arc = &arc;

  // initialize and define loading steps object
  LoadingSteps load;
  {
    loading_steps_initialization(&load);
    load.sups     = &(c->supports);
    load.sup_defl = &sup_defl;
  }

  // initialize and define material properties
  MaterialProperty mat;
  {
    material_initialization(&mat);
    mat.hommat  = c->hommat;
    mat.matgeom = c->matgeom;
  }

  /// initialize and define communication structures
  CommunicationStructure com;
  {
    communication_structure_initialization(&com);
    com.Ap     = c->Ap;
    com.Ai     = c->Ai;
    com.DomDof = c->DomDof;
    com.GDof   = c->GDof;
    com.nbndel = c->nbndel;
    com.bndel  = c->bndel;
    com.boot   = c->boot;
    com.net    = c->net;
    com.comm   = c->comm;
    com.rank   = c->rank;
    com.nproc  = c->nproc;
  }

  /// initialize and define time stepping variable
  TimeStepping ts;
  {
    time_stepping_initialization(&ts);
    ts.nt  = solver_file->n_step;
    ts.tim = s->tim;
    ts.times = s->times;
    ts.dt_n   = 0.0;
    ts.dt_np1 = s->dt;
    ts.print  = (long*) (solver_file->print_steps);
  }

  double dlm = Multiphysics_Arc_length(&grid, &mat, &fv, &sol, &load, &com, &ts,
                                       s->crpl, c->VVolume, opts, mp, 0);
  
  s->NORM = fv.NORM;
  *pores  = fv.pores;

  *DET0 = arc.DET0;
  *DLM0 = arc.DLM0;
  *DLM  = arc.DLM;
  *AT   = arc.AT;
  *ITT  = arc.ITT;
  *DAL  = arc.DAL;

  free(physicsname);

  return dlm;
}
