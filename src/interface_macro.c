/* HEADER */
#include "interface_macro.h"
#include "mkl_cblas.h"
#include <string.h>
#include <string.h>
#include <math.h>

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

int read_interface_macro_normal_lc(char *in_dir,
				   SUPP sup)
{
  int err = 0;
  const char *filename = "normal.in";
  int s_len = strlen(in_dir) + strlen(filename) + 3 /* 2 terminantion chars + / */;
  char *in_name = PGFEM_calloc (s_len,sizeof(char));
  if(sprintf(in_name,"%s/%s",in_dir,filename) >= s_len){
    PGFEM_printerr("Need to allocate more space for filename in %s\n",__func__);
    err++;
  }

  /* open file and read */
  FILE *in = fopen(in_name,"r");
  if(in == NULL){
    err++;
  } else {
    int n_matched = fscanf(in,"%lf %lf %lf %lf",&sup->lc,
			   &sup->N0[0],&sup->N0[1],&sup->N0[2]);
    if(n_matched != 4){
      PGFEM_printerr("Error reading file! (%s)\n",in_name);
      err++;
    }
  }

  if(err == 0){/* no problems reading the normal, normalize it */
    double mag = sqrt(sup->N0[0]*sup->N0[0]
		      + sup->N0[1]*sup->N0[1]
		      + sup->N0[2]*sup->N0[2]);
    sup->N0[0] /= mag;
    sup->N0[1] /= mag;
    sup->N0[2] /= mag;
  }

  free(in_name);
  fclose(in);
  return err;
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
    if(analysis == DISP){ /* get total jump */
      jump_u[0] = ((sup->defl[0] + sup->defl_d[0])
		   -(sup->defl[3] + sup->defl_d[3]));
      jump_u[1] = ((sup->defl[1] + sup->defl_d[1])
		   -(sup->defl[4] + sup->defl_d[4]));
      jump_u[2] = ((sup->defl[2] + sup->defl_d[2])
		   -(sup->defl[5] + sup->defl_d[5]));
    } else { /* get jump increment */
      jump_u[0] = (sup->defl_d[0] - sup->defl_d[3]);
      jump_u[1] = (sup->defl_d[1] - sup->defl_d[4]);
      jump_u[2] = (sup->defl_d[2] - sup->defl_d[5]);
    }
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
					 const NODE *ptrNode,
					 const double *F_0,
					 const int analysis)
{
  int err = 0;
  double coord[3];

  if(analysis == DISP){ /* get reference coords */
    coord[0] = ptrNode->x1_fd;
    coord[1] = ptrNode->x2_fd;
    coord[2] = ptrNode->x3_fd;
  } else { /* get current coords */
    coord[0] = ptrNode->x1;
    coord[1] = ptrNode->x2;
    coord[2] = ptrNode->x3;
  }

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
    case DISP:
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
    default:
      F0[0] = sup->defl_d[0];
      F0[1] = sup->defl_d[1];
      F0[2] = sup->defl_d[2];
      
      F0[3] = sup->defl_d[3];
      F0[4] = sup->defl_d[4];
      F0[5] = sup->defl_d[5];

      F0[6] = sup->defl_d[6];
      F0[7] = sup->defl_d[7];
      F0[8] = sup->defl_d[8];
      break;
    }
  } else if(sup->npd >= 6){ /* interface grad(u0) */
    double *ju = PGFEM_calloc(3,sizeof(double));
    err += compute_interface_macro_jump_u(ju,sup,analysis);
    err += compute_interface_macro_grad_u(F0,sup->lc,ju,sup->N0);
    free(ju);
  } else { /* not enough pre-disp */
    memset(F0,0,9*sizeof(double));
    err = 1;
  }
  return err;
 }
