/**
 * This driver program is a general material point simulator that
 * calls the available models in the Constitutive Modeling interface
 * (-cm for PGFem3D). The material point is loaded by a prescribed
 * macroscopic deformation gradient, in this case the average
 * deformation gradient defined for layers (F = 1 + [[u]] / lc ox N).
 *
 * AUTHORS:
 *   Matthew Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */

#include "constitutive_model.h"
#include "hommat.h"
#include "data_structure_c.h"
#include "index_macros.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>

Define_Matrix(double);

#define dim 3
#define tensor 9
//static const double eye[tensor] = {[0] = 1.0, [4] = 1.0, [8] = 1.0};
static const double eye[tensor] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; 
static const double LC = 0.20;
static const double JUMPU = 0.001;
//static const double DIR[dim] = {[2] = 1.0};
static const double DIR[dim] = {0.0, 0.0, 1.0};
//static const double NORMAL[dim] = {[2] = 1.0};
static const double NORMAL[dim] = {0.0, 0.0, 1.0};

/**
 * Compute the deformation gradient.
 */
void get_F0(const double lc,
            const double * restrict jumpu,
            const double * restrict N,
            double *F0)
{
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      F0[idx_2(i,j)] = eye[idx_2(i,j)] + jumpu[i] * N[j] / lc;
    }
  }
}

/**
 * Compute the total jump for the given step [0, inf).
 */
void compute_jumpu(const int step,
                   const double du,
                   const double *restrict dir,
                   double *restrict jumpu)
{
  for (int i =0; i < dim; i++) {
    jumpu[i] = (step + 1) * du * dir[i];
  }
}

/**
 * Read the material properties
 */
int ch_read_props(Constitutive_model *m,
                  Model_parameters *p,
                  HOMMAT *hm,
                  const char *filename)
{
  /* FILE FORMAT:
     # <- comment
     {
     # HOMMAT Properties
     E mu10 mu01 G nu dev vol
     }

     CM_TYPE
     {
     CM MODEL Properties
     }
  */
  int err = 0; 

  FILE *in = fopen(filename, "r");
  if (!in) {
    printf("ERROR: could not open \"%s\"\n", filename);
    return 1;
  }

  err += scan_for_valid_line(in);
  assert( !feof(in) && "Premature EOF");
  int brace = fgetc(in);
  assert(brace == '{' && "Expect opening brace before HOMMAT properties");

  err += scan_for_valid_line(in);
  assert( !feof(in) && "Premature EOF");

  /* read in HOMMAT (i.e., hyper-elastic) properties */
  int nread = fscanf(in, "%lf %lf %lf %lf %lf %d %d",
                     &(hm->E), &(hm->m10), &(hm->m01),&(hm->G), &(hm->nu),
                     &(hm->devPotFlag), &(hm->volPotFlag));
  assert(nread == 7 && "Expected 7 values for HOMMAT properties");

  err += scan_for_valid_line(in);
  assert( !feof(in) && "Premature EOF");
  brace = fgetc(in);
  assert(brace == '}' && "Expect closing brace after HOMMAT properties");

  err += scan_for_valid_line(in);
  assert( !feof(in) && "Premature EOF");
  int model_type = -1;
  nread = fscanf(in, "%d", &model_type);
  assert(nread == 1 && "Expected model type identifier");
  
  /* initialize the Model_parameters object */
  err += model_parameters_initialize(p, hm, model_type);

  /* read in the model parameters */
  err += scan_for_valid_line(in);
  assert( !feof(in) && "Premature EOF");
  brace = fgetc(in);
  assert(brace == '{' && "Expect opening brace before CM properties");
  err += p->read_param(p, in);

  /* initialize cm object */
  err += constitutive_model_initialize(m, p);

  fclose(in);
  return err;
}

static const char *optstring = "";

typedef struct {
  struct option opt;
  char *descr;
} opts_help;

static const opts_help opt_help[] = {
  {{"help", no_argument, NULL, 'h'},"Print this help message"},
  {{"lc", required_argument, NULL, 'l'},"Thickness of layer (default = 0.2)"},
  {{"du", required_argument, NULL, 'u'},"Load increment (default = 0.001)"},
  {{"dt", required_argument, NULL, 't'},"Time increment (default = 0.05)"},
  {{"dir", required_argument, NULL, 'd'},"Load loading direction (default = \"0,0,1\")"},
  {{"normal", required_argument, NULL, 'n'},"Normal to the interface (default = \"0,0,1\")"},
  {{"prop-file", required_argument, NULL, 'i'},"File to read properties from (default = props.in)"},
  {{"output", required_argument, NULL, 'o'},"File to write data to (default = out.dat)"}
}; 

typedef struct {
  double lc;
  double du;
  double dt;
  double dir[dim];
  double normal[dim];
  char *ifname;
  char *ofname;
} opts;

int chd_print_help()
{
  printf("CH_layer_driver: General MPS for a multi-scale cohesive layer.\n");
  printf("USAGE: CH_layer_driver [options]\n");
  printf("OPTIONS:\n");
  const int n_opt = sizeof(opt_help) / sizeof(opts_help);
  for (int i = 0; i < n_opt; i++) {
    if (opt_help[i].opt.has_arg != 0)
      printf("\t--%-9s (arg)  %s\n", opt_help[i].opt.name, opt_help[i].descr);
    else
      printf("\t--%-9s        %s\n", opt_help[i].opt.name, opt_help[i].descr);
  }
  return 1;
}

void nrm_vec(double *a)
{
  double nrm = 0.0;
  for(int i = 0; i < dim; i++) nrm += a[i] * a[i];
  nrm = sqrt(nrm);
  for (int i = 0; i < dim; i++) a[i] /= nrm;
}

int chd_parse_command_line(int argc,
                           char **argv,
                           opts *opt)
{
  int err = 0;

  /* copy option structure from help descriptions */
  const int n_opt = sizeof(opt_help) / sizeof(opts_help);
  struct option *OPTS = calloc(n_opt + 1, sizeof(*OPTS));
  for (int i = 0; i < n_opt; i++){
    OPTS[i] = opt_help[i].opt;
  }

  /* reset external getopt variables */
  opterr = 0;
  optind = 1;

  /* parse the command line */
  int print_help = 0;
  int opts_idx = 0;
  int opt_val = getopt_long_only(argc, argv, optstring, OPTS, &opts_idx);
  char *tmp_arg = NULL;
  while (opt_val != -1) {
    switch (opt_val) {
    case '?':
      printf("Unrecognized option: %s\n", argv[optind]);
      break;
    case 'h':
      print_help = 1;
      break;
    case 'l':
      opt->lc = atof(optarg);
      break;
    case 'u':
      opt->du = atof(optarg);
      break;
    case 't':
      opt->dt = atof(optarg);
      break;
    case 'd':
      /* parse the triple arg formatted as arg,arg,arg */
      tmp_arg = strdup(optarg);
      opt->dir[0] = atof(strtok(tmp_arg, ", "));
      opt->dir[1] = atof(strtok(NULL, ", "));
      opt->dir[2] = atof(strtok(NULL, ", "));
      free(tmp_arg);
      tmp_arg = NULL;
      break;
    case 'n':
      /* parse the triple arg formatted as arg,arg,arg */
      tmp_arg = strdup(optarg);
      opt->normal[0] = atof(strtok(tmp_arg, ", "));
      opt->normal[1] = atof(strtok(NULL, ", "));
      opt->normal[2] = atof(strtok(NULL, ", "));
      free(tmp_arg);
      tmp_arg = NULL;
      break;
    case 'i':
      opt->ifname = strdup(optarg);
      break;
    case 'o':
      opt->ofname = strdup(optarg);
      break;
    default:
      err++;
      break;
    }

    opt_val = getopt_long_only(argc, argv, optstring, OPTS, &opts_idx);
  }

  if (!opt->ifname) opt->ifname = strdup("props.in");
  if (!opt->ofname) opt->ofname = strdup("out.dat");
  if (print_help) err += chd_print_help();

  nrm_vec(opt->dir);
  nrm_vec(opt->normal);

  free(OPTS);
  return err;
}

opts* opts_construct()
{
  opts *opt = malloc(sizeof(*opt));
  /* set default options */
  opt->lc = 0.200;
  opt->du = 0.001;
  opt->dt = 0.05;
  opt->dir[0] = 0.0;
  opt->dir[1] = 0.0;
  opt->dir[2] = 1.0;
  opt->normal[0] = 0.0;
  opt->normal[1] = 0.0;
  opt->normal[2] = 1.0;
  opt->ifname = NULL;
  opt->ofname = NULL;
  return opt;
}

void opts_destroy(opts *opt)
{
  free(opt->ifname);
  free(opt->ofname);
  free(opt);
}

void opts_print(opts *opt)
{
  printf("opt->lc = %g\nopt->du = %g\nopt->dt = %g\n"
         "opt->dir = {%g, %g, %g}\n"
         "opt->normal = {%g, %g, %g}\n"
         "opt->ifname = %s\n"
         "opt->ofname = %s\n",
         opt->lc, opt->du, opt->dt,
         opt->dir[0], opt->dir[1], opt->dir[2],
         opt->normal[0], opt->normal[1], opt->normal[2],
         opt->ifname, opt->ofname);
}

int ch_output_ts(const Constitutive_model *m,
                 const double *jumpu,
                 const double *normal,
                 const double *F0,
                 const double dt,
                 FILE *out)
{
  int err = 0;

  Matrix_double F, S, P;
  Matrix_construct_init(double, S, dim, dim, 0.0);
  Matrix_construct(double, P);
  Matrix_construct(double, F);
  Matrix_init_w_array(F, dim, dim, F0);

  /* compute the total SPK stress tensor */
  err += constitutive_model_defaut_update_elasticity(m, &F,  NULL, &S, 0);

  /* compute the FPK stress tensor */
  Matrix_AxB(P, 1.0, 0.0, F, 0, S, 0);

  /* compute the traction vector */
  double trac[dim] = {0};
  cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim,
              1.0, P.m_pdata, dim, normal, 1, 0.0, trac, 1);

  fprintf(out,"%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",
          jumpu[0], jumpu[1], jumpu[2],
          trac[0],   trac[1],  trac[2]);

  Matrix_cleanup(F);
  Matrix_cleanup(S);
  Matrix_cleanup(P);
  return err;
}

int main(int argc, char **argv)
{
  int err = 0;
  const double alpha = 0.5;
  const int nstep = 20;

  /* parse command line */
  opts *opt = opts_construct();
  err += chd_parse_command_line(argc, argv, opt);
  opts_print(opt);

  FILE *out = fopen(opt->ofname, "w");
  if (!out) {
    err ++;
    printf("ERROR: coupld not open \"%s\"\n", opt->ofname);
  }

  /* construct/initialize objects */
  Constitutive_model *m = calloc(1, sizeof(*m));
  Model_parameters *p = calloc(1, sizeof(*p));
  HOMMAT *hm = calloc(1, sizeof(*hm));

  err += constitutive_model_construct(m);
  err += model_parameters_construct(p);
  err += ch_read_props(m, p, hm, opt->ifname);

  double F0[tensor] = {0};
  double jumpu[tensor] = {0};
  void *ctx = NULL;
  if (err) goto exit_main;

  for (int i = 0; i < nstep; i++) {
    /* compute the current deformation */
    compute_jumpu(i, opt->du, opt->dir, jumpu);
    get_F0(opt->lc, jumpu, opt->normal, F0);

    /* run the integration algorithm */
    err += construct_model_context(&ctx, p->type, F0, opt->dt, alpha,NULL);
    err += p->integration_algorithm(m, ctx);
    err += p->update_state_vars(m);
    err += ch_output_ts(m, jumpu, opt->normal, F0, opt->dt, out);
    err += p->destroy_ctx(&ctx);
    if (err) break;
  }

 exit_main:
  if (err) printf("Caught error, exiting early!\n");
  fclose(out);
  opts_destroy(opt);
  err += constitutive_model_destroy(m);
  err += model_parameters_destroy(p);
  free(m);
  free(p);
  free(hm);
  return err;
}
