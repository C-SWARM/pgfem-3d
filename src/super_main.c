/* HEADER */
/*** This is the PGFem3D main function. It has a branch between the
     macro-micro or single scale execution models. */

/** main function prototype for single scale */
int single_scale_main(int argc, char **argv);

/** main function prototype for multi scale */
int multi_scale_main(int argc, char **argv);


/** define switching variables */
enum{SINGLE_SCALE,MULTI_SCALE,N_FRAMEWORK};
static const char *framework_opts[] = {"-SS","-MS"};

#include <stdio.h>
#include <string.h>

/**
 * get the execution framework from the very first passed option.
 *
 * It is an error to not have any options or to pass an undefined
 * option. \return non-zero on error.
 */
static void get_framework(int argc,
			  char **argv,
			  int *PGFEM_FRAMEWORK)
{
  *PGFEM_FRAMEWORK = N_FRAMEWORK; /* poison value */
  if(argc < 2){
    fprintf(stderr,"ERROR: must provide at least one option!\n");
    return;
  } else {
    /* check for match agains the framework option strings */
    for(int i=0; i<N_FRAMEWORK; i++){
      if(!strcmp(argv[1],framework_opts[i])){
	*PGFEM_FRAMEWORK = i;
	break;
      }
    }
  }
}

static void print_valid_usage(FILE *out)
{
  fprintf(out,"Required usage: PGFem3D [FRAMEWORK] {options}\n");
  fprintf(out,"\tValid FRAMEWORK options: -SS (single scale) or -MS (multi-scale)\n\n");
  fprintf(out,"Run \'PGFem3D [FRAMEWORK] -h\' for other option information.\n");
}

int main(int argc, char **argv)
{
  int err = 0;
  int PGFEM_FRAMEWORK = 0;
  get_framework(argc, argv,&PGFEM_FRAMEWORK);
  switch(PGFEM_FRAMEWORK){
  case SINGLE_SCALE: err = single_scale_main(argc,argv); break;
  case MULTI_SCALE: err = multi_scale_main(argc,argv); break;
  default:
    print_valid_usage(stderr);
    err = 1;
    break;
  }

  return err;
}

