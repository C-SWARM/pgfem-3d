/* HEADER */
/* This file contains the cohesive potentials */

#include "cohesive_potentials.h"
#include <math.h>
#include <string.h>

#include "allocation.h"
#include "index_macros.h"
#include "cohesive_element_utils.h"

static const double max_hy_term = 25;
static const int ndim = 3;
static const double pi = 3.14159265358979323846;

int read_cohesive_properties(FILE *in,
                 int *n_mat,
                 cohesive_props **props,
                 MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  MPI_Comm_rank(mpi_comm,&myrank);

  int len = 500;
  char *line = PGFEM_calloc(char, len);
  cohesive_props *p = NULL;

  /* get number of materials allowing for blank line(s) */
  int match = 0;
  do{
    if (fgets(line,len,in) == nullptr) {
      abort();
    }
    match = sscanf(line,"%d",n_mat);
  }while(match != 1 && !feof(in));
  if(feof(in)){ err ++; goto exit_function;}

  /* allocate properties */
  (*props) = PGFEM_calloc(cohesive_props, *n_mat);
  /* alias for simpler access/readability */
  p = (*props);

  for(int i=0; i<*n_mat; i++){
    do{
      if (fgets(line,len,in) == nullptr) {
        abort();
      }
      match = sscanf(line,"%d",&(p[i].type));
    }while(match != 1 && !ferror(in));
    if(ferror(in)){ err ++; goto exit_function;}
    switch(p[i].type){
    case CO_MOD_MS:
      {
    if(myrank == 0)
      PGFEM_printerr("Cohesive material %d follows MULTISCALE potential.\n",i);
    p[i].nprops = 0;
    p[i].props = NULL;
    p[i].get_potential = Multiscale_potential;
    p[i].get_traction = Multiscale_traction;
    p[i].get_tangents = Multiscale_tangents;
    break;
      }
    case CO_MOD_NEEDLEMAN:
      {
    if(myrank == 0)
      PGFEM_printerr("Cohesive material %d follows NEEDLEMAN potential.\n",i);

    const struct NEEDLEMAN_PROP_KEY *pk = &Needleman_pkey;
    p[i].nprops = 4;
    p[i].props = PGFEM_calloc(double, p[i].nprops);
    sscanf(line,"%*d %lf %lf %lf %lf",
           &p[i].props[pk->sig_c],
           &p[i].props[pk->chi_c],
           &p[i].props[pk->beta],
           &p[i].props[pk->alpha]);
    p[i].get_potential = Needleman_potential;
    p[i].get_traction = Needleman_traction;
    p[i].get_tangents = Needleman_tangents;
    break;
      }
    case CO_MOD_VDW: /* van der Waals */
      {
    if(myrank == 0)
      PGFEM_printerr("Cohesive material %d follows VDW (van der Waals) potential.\n",i);

    const struct VDW_PROP_KEY *pk = &VDW_pkey;
    p[i].nprops = 2;
    p[i].props = PGFEM_calloc(double, p[i].nprops);
    sscanf(line,"%*d %lf",&p[i].props[pk->A]);
    /* beta is set to 0.0 -> normal component only */
    p[i].props[pk->beta] = 0.0;
    p[i].get_potential = vdWaals_potential;
    p[i].get_traction = vdWaals_traction;
    p[i].get_tangents = vdWaals_tangents;
    break;
      }
    case CO_MOD_LJ: /* Leonard-Jones */
      {
    if(myrank == 0)
      PGFEM_printerr("Cohesive material %d follows LJ (Leonard-Jones) potential.\n",i);

    const struct LJ_PROP_KEY *pk = &LJ_pkey;
    p[i].nprops = 2;
    p[i].props = PGFEM_calloc(double, p[i].nprops);
    sscanf(line,"%*d %lf %lf",&p[i].props[pk->A],
           &p[i].props[pk->X0]);
    /* beta is set to 0.0 -> normal component only */
    p[i].props[pk->beta] = 0.0;
    p[i].get_potential = LJ_potential;
    p[i].get_traction = LJ_traction;
    p[i].get_tangents = LJ_tangents;
    break;
      }

      /*=== NOT IMPLEMENTED ===*/
    case CO_MOD_MATOUS_ARCINIEGA:
    case CO_MOD_PPR:
    default:
      PGFEM_printerr("Cohesive type not yet implemented! %s:%s\n",
             __func__,__FILE__);
      PGFEM_Abort();
    }
  }

 exit_function:
  free(line);

  if(err){
    PGFEM_printerr("ERROR: problem reading file! %s:%s:%d\n",
           __func__,__FILE__,__LINE__);
    PGFEM_Abort();
  }
  return err;
}

void destroy_cohesive_props(const int n_mat,
                cohesive_props *props)
{
  for(int i=0; i<n_mat; i++){
    free(props[i].props);
  }
  free(props);
}

/*** PPR see "A unified potential-based cohesive model for mixed-mode
     fracture", JMPS 2009 ***/
int PPR_potential(double *p_pot,
          const double *jump,
          const double *normal,
          const double *props,
          const double *vars)
{
  int err = -1;
  /* static const struct PPR_PROP_KEY *pk = &PPR_pkey; */
  /*** NOT IMPLEMENTED ***/
  return err;
}

int PPR_traction(double *traction,
         const double *jump,
         const double *normal,
         const double *props,
         const double *vars)
{
  int err = -1;
  /* static const struct PPR_PROP_KEY *pk = &PPR_pkey; */
  /*** NOT IMPLEMENTED ***/
  return err;
}

int PPR_tangents(double *mat_tan,
         double *geo_tan,
         const double *jump,
         const double *normal,
         const double *props,
         const double *vars)
{
  int err = -1;
  /* static const struct PPR_PROP_KEY *pk = &PPR_pkey; */
  /*** NOT IMPLEMENTED ***/
  return err;
}

/*** Needleman exponential potential ***/
int Needleman_potential(double *p_pot,
            const double *jump,
            const double *normal,
            const double *props,
            const double *vars)
{
  int err = -1;
  /* static const struct NEEDLEMAN_PROP_KEY *pk = &Needleman_pkey; */
  /*** NOT IMPLEMENTED ***/
  return err;
}

static int Needleman_compute_t_hat(double *p_t_hat,
                   const double X_eff,
                   const double X_eff_max,
                   const double t_max,
                   const double *props)
{
  /* t_hat = 1/X_eff (d Psi/d X_eff X_eff) */
  int err = 0;
  static const struct NEEDLEMAN_PROP_KEY *pk = &Needleman_pkey;
  if (X_eff >= X_eff_max || X_eff_max < props[pk->chi_c]){ /* loading */
    *p_t_hat = (props[pk->sig_c]/props[pk->chi_c]
        *exp(1 - X_eff/props[pk->chi_c]));
  } else { /* Unloading */
    *p_t_hat = t_max/X_eff_max;
  }
  return err;
}

int Needleman_traction(double *traction,
               const double *jump,
               const double *normal,
               const double *props,
               const double *vars)
{
  int err = 0;
  static const struct NEEDLEMAN_PROP_KEY *pk = &Needleman_pkey;
  double Xn = 0.0;
  double Xs = 0.0;
  double X_eff = 0.0;
  const double b2 = props[pk->beta]*props[pk->beta];
  const double X_eff_max = vars[0];
  const double t_max = vars[1];

  /* compute opening */
  err += get_jump_nt(&Xn,&Xs,jump,normal);
  err += get_eff_jump(&X_eff,Xn,Xs,props[pk->beta]);

  /* compute t_hat */
  double t_hat = 0.0;
  err += Needleman_compute_t_hat(&t_hat,X_eff,X_eff_max,t_max,props);

  /* determine if we should kill element. Still returns actual value
     in case needed */
  if(Xn > 0
     && X_eff > props[pk->chi_c]
     && t_hat < props[pk->sig_c]/100.){
    err = 100;
  }

  /* compute traction */
  double tn = 0.0;
  if(Xn >= 0.0){
    tn = t_hat*Xn;
  } else { /* Contact */
    double coeff = props[pk->sig_c]/props[pk->alpha];
    double hy_term = props[pk->alpha]*Xn/props[pk->chi_c];
    if(fabs(hy_term) > max_hy_term) {
      hy_term = hy_term/fabs(hy_term)*max_hy_term;
    }

    /** FOR TESTING **/
    double chi = X_eff_max;
    /* double chi = X_eff; */

    tn = coeff*exp(1-chi/props[pk->chi_c])*sinh(hy_term);
  }

  for(int i=0; i<ndim; i++){
    traction[i] = tn*normal[i] + t_hat*b2*(jump[i]-Xn*normal[i]);
  }

  return err;
}

static int Needleman_compute_dt_hat(double *p_dt_hat,
                    const double X_eff,
                    const double t_hat,
                    const double X_eff_max,
                    const double t_max,
                    const double *props)
{
  /* Compute (d t_hat/d Xn) and (d t_hat/d Xs) */
  int err = 0;
  static const struct NEEDLEMAN_PROP_KEY *pk = &Needleman_pkey;

  /* Note: can write the tangents in terms of t_hat */
  if(X_eff >= X_eff_max || X_eff_max < props[pk->chi_c]){ /* loading */
    *p_dt_hat = -t_hat/props[pk->chi_c];
  } else { /* unloading */
    *p_dt_hat = t_max/X_eff_max;
  }
  return err;
}

int Needleman_tangents(double *mat_tan,
               double *geo_tan,
               const double *jump,
               const double *normal,
               const double *props,
               const double *vars)
{
  int err = 0;
  static const struct NEEDLEMAN_PROP_KEY *pk = &Needleman_pkey;
  double Xn = 0.0;
  double Xs = 0.0;
  double X_eff = 0.0;
  const double b2 = props[pk->beta]*props[pk->beta];
  const double X_eff_max = vars[0];
  const double t_max = vars[1];

  /* compute opening */
  err += get_jump_nt(&Xn,&Xs,jump,normal);
  err += get_eff_jump(&X_eff,Xn,Xs,props[pk->beta]);

  /* Compute t_hat and normal/shear derrivatives*/
  double t_hat = 0.0;
  double dt_hat = 0.0;

  err += Needleman_compute_t_hat(&t_hat,X_eff,X_eff_max,t_max,props);

  err += Needleman_compute_dt_hat(&dt_hat,X_eff,t_hat,
                  X_eff_max,t_max,props);

  /* compute derivatives jump */
  double *Chi_s = PGFEM_calloc(double, ndim);
  double *dChi_dX = PGFEM_calloc(double, ndim);
  double *dChi_dN = PGFEM_calloc(double, ndim);
  double *dChi_s_dX = PGFEM_calloc(double, ndim*ndim);
  double *dChi_s_dN = PGFEM_calloc(double, ndim*ndim);

  if(X_eff > 0){
    for(int i=0; i<ndim; i++){
      Chi_s[i] = jump[i] - Xn*normal[i];
      dChi_dX[i] = (b2*jump[i]+(1-b2)*Xn*normal[i])/X_eff;
      dChi_dN[i] = (1-b2)*Xn/X_eff*jump[i];
    }
  }

  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      double dij = (i==j)?1.0:0.0;
      dChi_s_dX[idx_2(i,j)] = dij - normal[i]*normal[j];
      dChi_s_dN[idx_2(i,j)] = -(normal[i]*jump[j] + Xn*dij);
    }
  }

  /* compute shear components which are the same for all cases */
  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      mat_tan[idx_2(i,j)] = (dt_hat*b2*Chi_s[i]*dChi_dX[j]
                 + t_hat*b2*dChi_s_dX[idx_2(i,j)]);

      geo_tan[idx_2(i,j)] = (dt_hat*b2*Chi_s[i]*dChi_dN[j]
                 + t_hat*b2*dChi_s_dN[idx_2(i,j)]);
    }
  }

  /* compute normal component depending on contact condition */
  if(Xn >= 0.0){ /*** OPENING ***/
    for(int i=0; i<ndim; i++){
      for(int j=0; j<ndim; j++){
    double dij = (i==j)?1.0:0.0;
    mat_tan[idx_2(i,j)] = (mat_tan[idx_2(i,j)]
                   + dt_hat*Xn*normal[i]*dChi_dX[j]
                   + t_hat*normal[i]*normal[j]);

    geo_tan[idx_2(i,j)] = (geo_tan[idx_2(i,j)]
                   + dt_hat*Xn*normal[i]*dChi_dN[j]
                   + t_hat*normal[i]*jump[j]
                   + t_hat*Xn*dij);
      }
    }
  } else { /*** CONTACT ***/

    double coeff_1 = props[pk->sig_c]/props[pk->alpha];
    double coeff_2 = props[pk->sig_c]/props[pk->chi_c];

    /** FOR TESTING **/
    double chi = X_eff_max;
    double dchi_dxn = 0.0;
    /* double chi = X_eff; */
    /* double dchi_dxn = (1-b2)*Xn/X_eff; */
    if(X_eff <= 0.0) dchi_dxn = 0.0;

    double coeff_3 = (props[pk->sig_c]*dchi_dxn
              /(props[pk->chi_c]*props[pk->alpha]));

    double exp_term = exp(1-chi/props[pk->chi_c]);
    double hy_term = Xn*props[pk->alpha]/props[pk->chi_c];
    if(fabs(hy_term) > max_hy_term) {
      hy_term = hy_term/fabs(hy_term)*max_hy_term;
    }

    double tn_hat = coeff_1*exp_term*sinh(hy_term);
    double dtn_hat = (coeff_2*exp_term*cosh(hy_term)
              - coeff_3*exp_term*sinh(hy_term));

    for(int i=0; i<ndim; i++){
      for(int j=0; j<ndim; j++){
    double dij = (i==j)?1.0:0.0;
    mat_tan[idx_2(i,j)] = (mat_tan[idx_2(i,j)]
                   + dtn_hat*normal[i]*normal[j]);

    geo_tan[idx_2(i,j)] = (geo_tan[idx_2(i,j)]
                   + dtn_hat*normal[i]*jump[j]
                   + tn_hat*dij );
      }
    }
  }

  free(Chi_s);
  free(dChi_dX);
  free(dChi_dN);
  free(dChi_s_dX);
  free(dChi_s_dN);
  return err;
}

/*** Sc : Xc : b : k -> Matous-Arciniega
     1 < k < 5.04 - parameter shifts the Xc to the right,
     but Gc is preserved ***/
int MA_potential(double *p_pot,
         const double *jump,
         const double *normal,
         const double *props,
         const double *vars)
{
  int err = -1;
  /* static const struct MA_PROP_KEY *pk = &MA_pkey; */
  /*** NOT IMPLEMENTED ***/
  return err;
}

int MA_traction(double *traction,
        const double *jump,
        const double *normal,
        const double *props,
        const double *vars)
{
  int err = 0;
  /* static const struct MA_PROP_KEY *pk = &MA_pkey; */
  return err;
}


int MA_tangents(double *mat_tan,
        double *geo_tan,
        const double *jump,
        const double *normal,
        const double *props,
        const double *vars)
{
  int err = 0;
  /* static const struct MA_PROP_KEY *pk = &MA_pkey; */
  return err;
}

int vdWaals_potential(double *p_pot,
              const double *jump,
              const double *normal,
              const double *props,
              const double *vars)
{
  int err = -1;
  /* static const struct VDW_PROP_KEY *pk = &VDW_pkey; */
  /*** NOT IMPLEMENTED ***/
  return err;
}

static inline int compute_vdWaals_eff_trac(double *eff_trac,
                       const double A,
                       const double X_eff)
{
  int err = 0;
  *eff_trac = A/(6*pi*X_eff*X_eff*X_eff);
  return err;
}

static inline int compute_vdWaals_eff_tan(double *eff_tan,
                      const double A,
                      const double X_eff)
{
  int err = 0;
  *eff_tan = -A/(2*pi*X_eff*X_eff*X_eff*X_eff);
  return err;
}

int vdWaals_traction(double *traction,
             const double *jump,
             const double *normal,
             const double *props,
             const double *vars)
{
  int err = 0;
  static const struct VDW_PROP_KEY *pk = &VDW_pkey;
  double Xn = 0.0;
  double Xs = 0.0;
  double X_eff = 0.0;

  /* NOTE: The normal is hard coded because the attractive force is
     always in the vertical direction */
  double N[3] = {0.0,0.0,1.0};

  /* compute opening */
  err += get_jump_nt(&Xn,&Xs,jump,N);
  err += get_eff_jump(&X_eff,Xn,Xs,0.0);

  /* compute eff_trac */
  double eff_trac = 0.0;
  err += compute_vdWaals_eff_trac(&eff_trac,props[pk->A],X_eff);

  /* compute traction */
  for(int i=0; i<ndim; i++){
    traction[i] = eff_trac*N[i];
  }

  return err;
}

int vdWaals_tangents(double *mat_tan,
             double *geo_tan,
             const double *jump,
             const double *normal,
             const double *props,
             const double *vars)
{
  int err = 0;
  static const struct VDW_PROP_KEY *pk = &VDW_pkey;
  double Xn = 0.0;
  double Xs = 0.0;
  double X_eff = 0.0;

  /* NOTE: The normal is hard coded because the attractive force is
     always in the vertical direction */
  double N[3] = {0.0,0.0,1.0};

  /* compute opening */
  err += get_jump_nt(&Xn,&Xs,jump,N);
  err += get_eff_jump(&X_eff,Xn,Xs,0.0);

  double eff_tan = 0.0;
  err += compute_vdWaals_eff_tan(&eff_tan,props[pk->A],X_eff);

  /* compute tangents */
  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      mat_tan[idx_2(i,j)] = eff_tan*N[i]*N[j];
      geo_tan[idx_2(i,j)] = 0.0;
    }
  }

  return err;
}

int LJ_potential(double *p_pot,
         const double *jump,
         const double *normal,
         const double *props,
         const double *vars)
{
  int err = -1;
  /* static const struct LJ_PROP_KEY *pk = &LJ_pkey; */
  /*** NOT IMPLEMENTED ***/
  return err;
}

static inline int compute_LJ_eff_trac(double *eff_trac,
                      const double A,
                      const double X0,
                      const double X_eff)
{
  int err = 0;
  *eff_trac = A/(6*pi*X0*X0*X0)*(pow(X0/X_eff,3)-pow(X0/X_eff,9));
  return err;
}

static inline int compute_LJ_eff_tan(double *eff_tan,
                     const double A,
                     const double X0,
                     const double X_eff)
{
  int err = 0;
  *eff_tan = A/(2*pi*pow(X_eff,10))*(3*pow(X0,6) - pow(X_eff,6));
  return err;
}

int LJ_traction(double *traction,
        const double *jump,
        const double *normal,
        const double *props,
        const double *vars)
{
  int err = 0;
  static const struct LJ_PROP_KEY *pk = &LJ_pkey;
  double Xn = 0.0;
  double Xs = 0.0;
  double X_eff = 0.0;

  /* NOTE: The normal is hard coded because the attractive force is
     always in the vertical direction */
  double N[3] = {0.0,0.0,1.0};

  /* compute opening */
  err += get_jump_nt(&Xn,&Xs,jump,N);
  err += get_eff_jump(&X_eff,Xn,Xs,0.0);

  /* compute eff_trac */
  double eff_trac = 0.0;
  err += compute_LJ_eff_trac(&eff_trac,props[pk->A],props[pk->X0],X_eff);

  /* compute traction */
  for(int i=0; i<ndim; i++){
    traction[i] = eff_trac*N[i];
  }

  return err;
}

int LJ_tangents(double *mat_tan,
        double *geo_tan,
        const double *jump,
        const double *normal,
        const double *props,
        const double *vars)
{
  int err = 0;
  static const struct LJ_PROP_KEY *pk = &LJ_pkey;
  double Xn = 0.0;
  double Xs = 0.0;
  double X_eff = 0.0;

  /* NOTE: The normal is hard coded because the attractive force is
     always in the vertical direction */
  double N[3] = {0.0,0.0,1.0};

  /* compute opening */
  err += get_jump_nt(&Xn,&Xs,jump,N);
  err += get_eff_jump(&X_eff,Xn,Xs,0.0);

  /* compute eff_trac & eff_tan*/
  double eff_trac = 0.0;
  double eff_tan = 0.0;
  err += compute_LJ_eff_trac(&eff_trac,props[pk->A],props[pk->X0],X_eff);
  err += compute_LJ_eff_tan(&eff_tan,props[pk->A],props[pk->X0],X_eff);

  /* compute tangents */
  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      mat_tan[idx_2(i,j)] = eff_tan*N[i]*N[j];
      geo_tan[idx_2(i,j)] = 0.0;
    }
  }

  return err;
}

int Multiscale_potential(double *p_pot,
             const double *jump,
             const double *normal,
             const double *props,
             const double *vars)
{
  *p_pot = 0.0;
  return 0;
}

int Multiscale_traction(double *traction,
            const double *jump,
            const double *normal,
            const double *props,
            const double *vars)
{
  memset(traction,0,ndim*sizeof(double));
  return 0;
}

int Multiscale_tangents(double *mat_tan,
            double *geo_tan,
            const double *jump,
            const double *normal,
            const double *props,
            const double *vars)
{
  memset(mat_tan,0,ndim*ndim*sizeof(double));
  memset(geo_tan,0,ndim*ndim*sizeof(double));
  return 0;
}
