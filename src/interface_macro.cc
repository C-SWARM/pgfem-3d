#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "interface_macro.h"
#include "allocation.h"
#include "enumerations.h"
#include "utils.h"
#include <mkl_cblas.h>
#include <cmath>
#include <cstring>

int read_interface_macro_normal_lc(const char *in_dir, SUPP sup) {
  std::string fn = in_dir;
  fn += "/normal.in";

  /* open file and read */
  if (FILE *in = PGFEM_fopen(fn.c_str(), "r")) {
    CHECK_SCANF(in, "%lf %lf %lf %lf %lf",
                &sup->v0, &sup->lc,
                &sup->N0[0], &sup->N0[1], &sup->N0[2]);
    fclose(in);
  }
  else {
    return 1;
  }

  if (sup->v0 == 0.0) {
    PGFEM_printerr("ERROR: specified 0.0 volume! (%s)\n", fn.c_str());
    return 1;
  }

  /* no problems reading the normal, normalize it */
  double mag = sqrt(sup->N0[0] * sup->N0[0] +
                    sup->N0[1] * sup->N0[1] +
                    sup->N0[2] * sup->N0[2]);
  sup->N0[0] /= mag;
  sup->N0[1] /= mag;
  sup->N0[2] /= mag;
  return 0;
}

int compute_interface_macro_jump_u(double *jump_u,
                                   const SUPP sup,
                                   const int analysis)
{
  int err = 0;
  if(sup->npd < 6){ /* incorrect number of prescribed displacements */
    memset(jump_u,0,3*sizeof(double));
    err = 1;
  } else {
    /* if(analysis == DISP){ /\* get total jump *\/ */
    jump_u[0] = ((sup->defl[0] + sup->defl_d[0])
                 -(sup->defl[3] + sup->defl_d[3]));
    jump_u[1] = ((sup->defl[1] + sup->defl_d[1])
                 -(sup->defl[4] + sup->defl_d[4]));
    jump_u[2] = ((sup->defl[2] + sup->defl_d[2])
                 -(sup->defl[5] + sup->defl_d[5]));
    /* } else { /\* get jump increment *\/ */
    /*   jump_u[0] = (sup->defl_d[0] - sup->defl_d[3]); */
    /*   jump_u[1] = (sup->defl_d[1] - sup->defl_d[4]); */
    /*   jump_u[2] = (sup->defl_d[2] - sup->defl_d[5]); */
    /* } */
  }

  return err;
} /* compute_interface_macro_jump_u() */

int compute_interface_macro_grad_u(double *F_0,
                                   const double lc,
                                   const double *jump_u,
                                   const double *normal)
{
  int err = 0;
  /* F_0 = 1/lc [[u]] ox N */
  memset(F_0,0,9*sizeof(double));
  cblas_dger(CblasRowMajor,3,3,1./lc,jump_u,1,normal,1,F_0,3);
  return err;
} /* compute_interface_macro_grad_u() */

int compute_interface_macro_disp_at_node(double *u_0,
                                         const Node *ptrNode,
                                         const double *F_0,
                                         const int analysis)
{
  int err = 0;
  double coord[3];

  /* if(analysis == DISP){ /\* get reference coords *\/ */
  coord[0] = ptrNode->x1_fd;
  coord[1] = ptrNode->x2_fd;
  coord[2] = ptrNode->x3_fd;
  /* } else { /\* get current coords *\/ */
  /*   coord[0] = ptrNode->x1; */
  /*   coord[1] = ptrNode->x2; */
  /*   coord[2] = ptrNode->x3; */
  /* } */

  /* u_0 = F_0{=grad(u_0)} * X */
  cblas_dgemv(CblasRowMajor,CblasNoTrans,
              3,3,1.0,F_0,3,coord,1,0.0,u_0,1);

  return err;
} /* compute_interface_macro_disp_at_node() */

int compute_macro_grad_u(double *F0,
                         const SUPP sup,
                         const int analysis)
{
  int err = 0;
  if(sup->npd >= 9){ /* Bulk grad(u0) */
    switch(analysis){
     default:
      /* case DISP: */
      F0[0] = sup->defl[0] + sup->defl_d[0];
      F0[1] = sup->defl[1] + sup->defl_d[1];
      F0[2] = sup->defl[2] + sup->defl_d[2];

      F0[3] = sup->defl[3] + sup->defl_d[3];
      F0[4] = sup->defl[4] + sup->defl_d[4];
      F0[5] = sup->defl[5] + sup->defl_d[5];

      F0[6] = sup->defl[6] + sup->defl_d[6];
      F0[7] = sup->defl[7] + sup->defl_d[7];
      F0[8] = sup->defl[8] + sup->defl_d[8];
      break;
      /* default: */
      /*   F0[0] = sup->defl_d[0]; */
      /*   F0[1] = sup->defl_d[1]; */
      /*   F0[2] = sup->defl_d[2]; */

      /*   F0[3] = sup->defl_d[3]; */
      /*   F0[4] = sup->defl_d[4]; */
      /*   F0[5] = sup->defl_d[5]; */

      /*   F0[6] = sup->defl_d[6]; */
      /*   F0[7] = sup->defl_d[7]; */
      /*   F0[8] = sup->defl_d[8]; */
      /*   break; */
    }
  } else if(sup->npd >= 6){ /* interface grad(u0) */
    double *ju = PGFEM_calloc(double, 3);
    err += compute_interface_macro_jump_u(ju,sup,analysis);
    err += compute_interface_macro_grad_u(F0,sup->lc,ju,sup->N0);
    free(ju);
  } else { /* not enough pre-disp */
    memset(F0,0,9*sizeof(double));
    err = 1;
  }
  return err;
}
