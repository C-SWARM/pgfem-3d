#ifndef COHESIVE_POTENTIALS_H
#define COHESIVE_POTENTIALS_H

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

typedef int (*cohesive_fptr)(double *val,
			       const double *jump,
			       const double *normal,
			       const double *props,
			       const double *vars);

typedef int (*cohesive_fptr2)(double *val,
				double *val2,
				const double *jump,
				const double *normal,
				const double *props,
				const double *vars);


typedef struct COHESIVE_PROPS {
int type;
long nprops;
double *props;
cohesive_fptr get_potential;
cohesive_fptr get_traction;
cohesive_fptr2 get_tangents;
} cohesive_props;

struct PPR_PROP_KEY{
unsigned int phi_n,phi_t,sig_max,tau_max,
  alpha,beta,lam_n,lam_t;
};
static const struct PPR_PROP_KEY PPR_pkey = {0,1,2,3,4,5,6,7};

struct NEEDLEMAN_PROP_KEY{
unsigned int sig_c,chi_c,beta,alpha;
};
static const struct NEEDLEMAN_PROP_KEY Needleman_pkey = {0,1,2,3};

struct MA_PROP_KEY{
unsigned int sig_c,chi_c,beta,kappa,alpha;
};
static const struct MA_PROP_KEY MA_pkey = {0,1,2,3,4};

/* van der Waals potential */
struct VDW_PROP_KEY{
unsigned int A,beta;
};
static const struct VDW_PROP_KEY VDW_pkey = {0,1};

/* Leonard-Jones potential */
struct LJ_PROP_KEY{
unsigned int A,X0,beta;
};
static const struct LJ_PROP_KEY LJ_pkey = {0,1,2};

/* enumeration for cohesive model types */
enum {CO_MOD_MS=-1,
	CO_MOD_NEEDLEMAN,
	CO_MOD_MATOUS_ARCINIEGA,
	CO_MOD_PPR,
	CO_MOD_VDW,
	CO_MOD_LJ};

/** Reads/allocates the list of coheisive properties */
int read_cohesive_properties(FILE *in,
			       int *n_mat,
			       cohesive_props **props,
			       MPI_Comm mpi_comm);

/** Properly frees up memory from a list of coheive properties */ 
void destroy_cohesive_props(const int n_mat,
			      cohesive_props *props);

int Needleman_potential(double *p_pot,
			  const double *jump,
			  const double *normal,
			  const double *props,
			  const double *vars);

int Needleman_traction(double *traction,
			 const double *jump,
			 const double *normal,
			 const double *props,
			 const double *vars);

int Needleman_tangents(double *mat_tan,
			 double *geo_tan,
			 const double *jump,
			 const double *normal,
			 const double *props,
			 const double *vars);

/*=== NOTE: The atomistic potentials are one dimensional and do not
  depend on the normal to the surface, but rather the direction of
  the flat plate. The *smooth & rigid* flat plate is assumed to be
  orientd with the normal to its surface in the {0,0,1}'
  direction. The normal passed to the function (normal to the mean
  surface) is hard coded. Note also that the normal is not dependednt
  on the deformation and therefore ther is no geometric tangent. */

int vdWaals_potential(double *p_pot,
			const double *jump,
			const double *normal,
			const double *props,
			const double *vars);

int vdWaals_traction(double *traction,
		       const double *jump,
		       const double *normal,
		       const double *props,
		       const double *vars);

int vdWaals_tangents(double *mat_tan,
		       double *geo_tan,
		       const double *jump,
		       const double *normal,
		       const double *props,
		       const double *vars);

int LJ_potential(double *p_pot,
		   const double *jump,
		   const double *normal,
		   const double *props,
		   const double *vars);

int LJ_traction(double *traction,
		  const double *jump,
		  const double *normal,
		  const double *props,
		  const double *vars);

int LJ_tangents(double *mat_tan,
		  double *geo_tan,
		  const double *jump,
		  const double *normal,
		  const double *props,
		  const double *vars);

/** The multiscale functions return 0 for everything. The tangents
    and traction must be computed/assembled by other routines */
int Multiscale_potential(double *p_pot,
			   const double *jump,
			   const double *normal,
			   const double *props,
			   const double *vars);

int Multiscale_traction(double *traction,
			  const double *jump,
			  const double *normal,
			  const double *props,
			  const double *vars);

int Multiscale_tangents(double *mat_tan,
			  double *geo_tan,
			  const double *jump,
			  const double *normal,
			  const double *props,
			  const double *vars);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
