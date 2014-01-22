/* HEADER */
/*** This is the PGFem3D main function. It has a branch between the
     macro-micro or single scale execution models. */

/** main function prototype for single scale */
int single_scale_main(int argc, char **argv);

/** main function prototype for multi scale */
int multi_scale_main(int argc, char **argv);


/** define switching variables */
enum{SINGLE_SCALE,MULTI_SCALE};

#ifndef PGFEM_FRAMEWORK
#define PGFEM_FRAMEWORK MULTI_SCALE
#endif

#include <stdio.h>

int main(int argc, char **argv)
{
  int err = 0;

  switch(PGFEM_FRAMEWORK){
  case SINGLE_SCALE: err = single_scale_main(argc,argv); break;
  case MULTI_SCALE: err = multi_scale_main(argc,argv); break;
  default:
    fprintf(stderr,"ERROR: Undefined PGFEM_FRAMEORK (%d)!\n\n",
	    PGFEM_FRAMEWORK);
    err = 1;
    break;
  }

  return err;
}

