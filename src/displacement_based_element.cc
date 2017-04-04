/* HEADER */
/**
 * AUTHORS:
 * Matthew Mosby
 */

/** This is the element definition for a TOTAL LAGRANGIAN
    displacement-based formulation. */
#include "displacement_based_element.h"
#include <string.h>
#include <errno.h>
#include <math.h>
#include "mkl_cblas.h"

#include "PGFEM_io.h"
#include "cast_macros.h"
#include "utils.h"
#include "allocation.h"
#include "incl.h"
#include "tensors.h"
#include "def_grad.h"
#include "elem3d.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "new_potentials.h"
#include "cohesive_element_utils.h"
#include "transform_coordinates.h"
#include "index_macros.h"
#include "data_structure_c.h"
#include "constitutive_model.h"
#include <ttl/ttl.h>

#ifndef DISP_DEBUG
#define DISP_DEBUG 0
#endif

#ifndef DISP_DEBUG_NO_LM
#define DISP_DEBUG_NO_LM 1
#endif

Define_Matrix(double);

//ttl declarations
namespace {
  template<int R, int D = 3, class S = double>
  using Tensor = ttl::Tensor<R, D, S>;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
}

static const int ndn = 3;
static const double LAGRANGE_THRESH = 0.1;
static const double DIFF_THRESH = 0.01;
static const damage empty_damage = {0};

static inline double del(const int i, const int j){ return (i==j?1.0:0.0); }

/*==============================================================*/
/*                DECLARE LOCAL FUNCTIONS                       */
/*==============================================================*/

/** get the material stress tensor (effective) */
static void get_material_stress(const double kappa,
				const HOMMAT *hommat, /* pointer to single hommat */
				const double *C,
				const double *C_I,
				const double J,
				double *Sbar);

/** get the material stiffness tensor (TOTAL) */
static void get_material_stiffness(const double kappa,
				   const HOMMAT *hommat, /* p_mat */
				   const damage *dam,
				   const double *C,
				   const double *C_I,
				   const double J,
				   const double *Sbar,
				   double *L);

static void get_Lbar(const double kappa,
		     const HOMMAT *p_mat,
		     const double *C,
		     const double *C_I,
		     const double J,
		     double *Lbar);

/** Get the material sensitivity tensor (TOTAL) */
static void get_material_sensitivity(const double kappa,
				     const HOMMAT *p_mat,
				     const damage *p_dam,
				     const double *C,
				     const double *C_I,
				     const double J,
				     const double *Sbar,
				     const double *Lbar,
				     double *K);


/** Compute the residual for the displacement-based element at an
    integration point. */
static void disp_based_resid_at_ip(double *R,
				   const int nne,
				   const double *ST,
				   const double *F,
				   const double *Sbar,
				   const damage *dam,
				   const double jj,
				   const double wt);

static void disp_based_bnd_resid_at_ip(double *R,
				       const int nne_ve,
				       const int ndofn,
				       const int n_lm_vecs,
				       const double *N_lm,
				       const double *ST,
				       const double *normal,
				       const double *F,
				       const double *Sbar,
				       const double *L,
				       const double *LM,
				       const double ome,
				       const double jj,
				       const double wt);

/** Compute the tangent for the displacement-based element at an
    integration point. */
static void disp_based_tan_at_ip(double *K,
				 const int nne,
				 const double *ST,
				 const double *F,
				 const double *Sbar,
				 const double *L,
				 const damage *dam,
				 const double jj,
				 const double wt);

static void disp_based_bnd_tan_at_ip(double *K,
				     const int nn_ve,
				     const int ndofn,
				     const int n_lm_vecs,
				     const double *N_lm,
				     const double *ST,
				     const double *normal,
				     const double *F,
				     const double *Sbar,
				     const double *L,
				     const double *KK,
				     const double *LM,
				     const double ome,
				     const double jj,
				     const double wt);

static void disp_based_bnd_Kuu_at_ip(double *Kuu,
				     const int nn_ve,
				     const int ndofn,
				     const double *ST,
				     const double *normal,
				     const double *F,
				     const double *L,
				     const double *KK,
				     const double *LM,
				     const double ome,
				     const double jj,
				     const double wt);

static void disp_based_bnd_Kul_Klu_at_ip(double *Kul,
					 double *Klu,
					 const int nn_ve,
					 const int ndofn,
					 const int n_lm_vecs,
					 const double *N_lm,
					 const double *ST,
					 const double *normal,
					 const double *F,
					 const double *Sbar,
					 const double *L,
					 const double ome,
					 const double jj,
					 const double wt);

/** Simple helper function to initialize variables for the integration
    loop(s) */
static int integration_help(const int elem_id,
			    const int nne,
			    const int ndofn,
			    const int i,
			    const int j,
			    const int k,
			    const double *x,
			    const double *y,
			    const double *z,
			    const double *int_pt_ksi,
			    const double *int_pt_eta,
			    const double *int_pt_zet,
			    const double *weights,
			    const double *disp,
			    const SUPP sup,
			    double *wt,
			    double *jj,
			    double *Na,
			    double *N_x,
			    double *N_y,
			    double *N_z,
			    double *ST,
			    double *Fr,
			    double *C,
			    double *C_I,
			    double *J);

/** helper function which handles all of the domain mapping and
    computes function values at the 2D integration points */
static int bnd_integration_help(const BOUNDING_ELEMENT *ptr_be,
				const ELEMENT *ptr_ve,
				const damage *ptr_dam,
				const SUPP sup,
				const HOMMAT *ptr_mat,
				const int ndofn,
				const int ip_2d,
				const double *int_pt_ksi_2d,
				const double *int_pt_eta_2d,
				const double *wts_2d,
				const double *x_ve,
				const double *y_ve,
				const double *z_ve,
				const double *KSI,
				const double *ETA,
				const double *ZET,
				const double *ve_disp,
				const int n_lm_vecs,
				const double *elem_lm,
				const double *x_be,
				const double *y_be,
				const double *z_be,
				const double kappa,
				double *LM,
				double *ome,
				double *Sbar,
				double *F,
				double *J,
				double *C,
				double *C_I,
				double *ST,
				double *normal,
				double *N3d,
				double *N2d,
				double *N_lm,
				double *jj,
				double *wt);

/** compute microscale contribution to K_00_e at an integration
    point */
static int compute_K_00_e_at_ip_2(double *K_00_e,
				  /* macro information */
				  const int macro_nnode,
				  const int macro_ndofn,
				  const double macro_int_wt,
				  const double *macro_shape_func,
				  const double *macro_normal,
				  const double layer_thickness,
				  const double *gNoxN,
				  /* micro information */
				  const int nne,
				  const double micro_volume0,
				  const double jj,
				  const double wt,
				  const double *F,
				  const double *ST,
				  const damage *p_dam,
				  const double *Sbar,
				  const double *L);

/** compute microscale contribution to K_01_e at an integration
    point */
static int compute_K_01_e_at_ip(double *K_01_e,
				/* macro information */
				const int macro_nnode,
				const int macro_ndofn,
				const double macro_int_wt,
				const double *macro_shape_func,
				const double *macro_normal,
				const double *gNoxN,
				/* micro information */
				const int nne,
				const double micro_volume0,
				const double jj,
				const double wt,
				const double *F,
				const double *ST,
				const damage *p_dam,
				const double *Sbar,
				const double *L);

/** compute microscale contribution to K_10_e at an integration
    point */
static int compute_K_10_e_at_ip(double *K_10_e,
				/* macro information */
				const int macro_nnode,
				const int macro_ndofn,
				const double macro_int_wt,
				const double *macro_shape_func,
				const double *macro_normal,
				const double *gNoxN,
				/* micro information */
				const int nne,
				const double micro_volume0,
				const double jj,
				const double wt,
				const double *F,
				const double *ST,
				const damage *p_dam,
				const double *Sbar,
				const double *L);


static int disp_cm_material_response(double *S,
                                     double *L,
                                     const Constitutive_model *m,
                                     const double *F,
                                     const double dt,
                                     const int get_L)
{
  int err = 0;
  Matrix_double MF, ML, MS;
  Matrix_construct(double, MF);
  Matrix_init_w_array(MF, ndn, ndn, F);
  Matrix_construct_init(double, MS, ndn, ndn, 0.0);
  Matrix_construct_init(double, ML, ndn*ndn, ndn*ndn, 0.0);

  err += constitutive_model_default_update_elasticity(m, &MF, &ML, &MS, get_L);
  memcpy(S, MS.m_pdata, 9 * sizeof(*S));
  if(get_L) memcpy(L, ML.m_pdata, 81 * sizeof(*L));

  Matrix_cleanup(MF);
  Matrix_cleanup(MS);
  Matrix_cleanup(ML);
  return err;
}

/*==============================================================*/
/*               DEFINE  API  FUNCTIONS                         */
/*==============================================================*/

int DISP_get_material_potential(const double kappa,
				const HOMMAT *hommat,/* pointer to 1 hommat */
				const double *C,
				const double J,
				double *Wbar)
{
  int err = 0;
  devPotentialFuncPtr DevPot = getDevPotentialFunc(0,hommat);
  UFuncPtr VolPot = getUFunc(0,hommat);
  *Wbar = 0.0;
  double U = 0.0;

  DevPot(C,hommat,Wbar);
  VolPot(J,hommat,&U);

  *Wbar = (*Wbar+kappa*U);

  return err;
}/* DISP_get_material_potential() */

int DISP_stiffmat_el(double *Ks,
		     const int ii,
		     const int ndofn,
		     const int nne,
		     const double *x,
		     const double *y,
		     const double *z,
		     const ELEMENT *elem,
		     const HOMMAT *hommat,
		     const long *nod,
		     const NODE *node,
		     const EPS *eps,
		     const SIG *sig,
		     const SUPP sup,
		     const double *disp,
                     const double dt) 
{
  int err = 0;
  const int mat = elem[ii].mat[2];
  const double kappa = hommat_get_kappa(&(hommat[mat]));
  const HOMMAT *ptrMat = &hommat[mat];

  int ndofe = 0;
  for(int i=0; i<nne; i++){
    ndofe += ndofn;
  }

  /* make sure the stiffenss matrix contains all zeros */
  memset(Ks,0,ndofe*ndofe*sizeof(double));

  double *F, *C, *C_I, *Sbar, *L, J;
  F = aloc1(ndn*ndn);
  C = aloc1(ndn*ndn);
  C_I = aloc1(ndn*ndn);
  Sbar = aloc1(ndn*ndn);
  L = aloc1(ndn*ndn*ndn*ndn);

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *N_x, *N_y, *N_z, *ST;
  Na = aloc1 (nne);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  N_z = aloc1 (nne);
  ST = aloc1(3*3*ndn*nne);

  /*=== INTEGRATION LOOP ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
	    int_pt_ksi,int_pt_eta,int_pt_zet,
	    weights);
  int ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){
	double wt,jj;
	err = integration_help(ii,nne,ndofn,i,j,k,x,y,z,
				int_pt_ksi,int_pt_eta,int_pt_zet,
			       weights,disp,sup,&wt,&jj,Na,N_x,N_y,N_z,ST,
				F,C,C_I,&J);

	/* inverted element detected, exit and return error */
	if(err != 0 ) goto exit_function;

	const damage *ptrDam = NULL;
        if (eps[ii].model == NULL) {
          ptrDam = &(eps[ii].dam[ip]);
          get_material_stress(kappa,&hommat[mat],C,C_I,J,Sbar);
          get_material_stiffness(kappa,ptrMat,ptrDam,C,C_I,J,Sbar,L);
        } else {
          ptrDam = &empty_damage;
          err += disp_cm_material_response(Sbar, L, eps[ii].model + ip,
                                           F, dt, 1);
        }
	disp_based_tan_at_ip(Ks,nne,ST,F,Sbar,L,ptrDam,jj,wt);
	ip++;
      }
    }
  }

 exit_function:
  free(F);
  free(C);
  free(C_I);
  free(Sbar);
  free(L);

  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);

  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  free(ST);

  return err;
}

int DISP_resid_el(double *R,
		  const int ii,
		  const int ndofn,
		  const int nne,
		  const double *x,
		  const double *y,
		  const double *z,
		  const ELEMENT *elem,
		  const HOMMAT *hommat,
		  const long *nod,
		  const NODE *node,
		  const EPS *eps,
		  const SIG *sig,
		  const SUPP sup,
		  const double *disp,
                  const double dt)
{
  int err = 0;
  const int mat = elem[ii].mat[2];
  const double kappa = hommat_get_kappa(&(hommat[mat]));
  const HOMMAT *ptrMat = &hommat[mat];

  int ndofe = 0;
  for(int i=0; i<nne; i++){
    ndofe += ndofn;
  } 

  /* Make sure that the residual vector contains zeros */
  memset(R,0,ndofe*sizeof(double));

  double *F, *C, *C_I, *Sbar, J;
  F = aloc1(ndn*ndn);
  C = aloc1(ndn*ndn);
  C_I = aloc1(ndn*ndn);
  Sbar = aloc1(ndn*ndn);

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *N_x, *N_y, *N_z, *ST;
  Na = aloc1 (nne);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  N_z = aloc1 (nne);
  ST = aloc1(3*3*ndn*nne);

  /*=== INTEGRATION LOOP ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
	    int_pt_ksi,int_pt_eta,int_pt_zet,
	    weights);
  int ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){
	double wt,jj;
	err = integration_help(ii,nne,ndofn,i,j,k,x,y,z,
				int_pt_ksi,int_pt_eta,int_pt_zet,
			       weights,disp,sup,&wt,&jj,Na,N_x,N_y,N_z,ST,
				F,C,C_I,&J);

	/* inverted element detected, exit and return error */
	if(err != 0) goto exit_function;

	const damage *ptrDam = NULL;
        if (eps[ii].model == NULL) {
          ptrDam = &(eps[ii].dam[ip]);
          get_material_stress(kappa,ptrMat,C,C_I,J,Sbar);
        } else {
          ptrDam = &empty_damage;
          void *ctx = NULL;
          construct_model_context(&ctx, eps[ii].model[ip].param->type, F, dt, 0.5, NULL);
          eps[ii].model[ip].param->integration_algorithm(&eps[ii].model[ip], ctx);
          err += disp_cm_material_response(Sbar, NULL, eps[ii].model + ip,
                                           F, dt, 0);
          eps[ii].model[ip].param->destroy_ctx(&ctx);
        }
	disp_based_resid_at_ip(R,nne,ST,F,Sbar,ptrDam,jj,wt);
	ip++;
      }
    }
  }

 exit_function:
  free(F);
  free(C);
  free(C_I);
  free(Sbar);

  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);

  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  free(ST);

  return err;
}

int DISP_resid_bnd_el(double *R,
		      const int b_el_id,
		      const int ndofn,
		      const int ndof_ve,
		      const double *x_ve,
		      const double *y_ve,
		      const double *z_ve,
		      const BOUNDING_ELEMENT *b_elems,
		      const ELEMENT *elem,
		      const HOMMAT *hommat,
		      const NODE *node,
		      const EPS *eps,
		      const SIG *sig,
		      const SUPP sup,
		      const double *vol_disp)
{
  /* NOTE: This function assumes constant damage on the element. More
     work is required to properly account for higher order damage
     functions */

  int err = 0;

  /* Get const pointers */
  const BOUNDING_ELEMENT *ptr_be = &b_elems[b_el_id];
  const ELEMENT *ptr_ve = &elem[ptr_be->vol_elem_id];
  const HOMMAT *ptrMat = &hommat[ptr_ve->mat[2]];
  const damage *ptr_dam = &eps[ptr_be->vol_elem_id].dam[0];
  const long *loc_nod = ptr_be->loc_nodes;

  /* compute constant values */
  const double kappa = hommat_get_kappa(ptrMat);
  const int nne_ve = ptr_ve->toe;
  const int nne_be = ptr_be->nnodes;

  /* clear contents of R */
  memset(R,0,ndof_ve*sizeof(double));

  if(DISP_DEBUG_NO_LM){
    return err;
  }

  /* separate the displacement and Lagrange multiplier degrees of
     freedom. The displacements are stored first in order of nodal
     connectivity, followed by the LMs (see get_dof_ids.h) */
 
  /* how many LM "nodes" are there for interpolation */
  const int n_lm_vecs = ptr_be->n_dofs/ndn; 
  const int node_dof = nne_ve*ndofn;
  const double *ptr_disp = vol_disp;
  const double *ptr_lm = vol_disp + node_dof;

  /* get the volume element parent coordinates */
  double *KSI = aloc1(nne_ve);
  double *ETA = aloc1(nne_ve);
  double *ZET = aloc1(nne_ve);
  get_element_node_parent_coords(nne_ve,KSI,ETA,ZET);

  /* get the boundary element xyz coords */
  double *x_be = aloc1(nne_be);
  double *y_be = aloc1(nne_be);
  double *z_be = aloc1(nne_be);
  for(int i=0; i<nne_be; i++){
    int n_id = loc_nod[i];
    x_be[i] = x_ve[n_id];
    y_be[i] = y_ve[n_id];
    z_be[i] = z_ve[n_id];
  }

  /* compute values from volumetric element */
  double *F = aloc1(9);
  double *C = aloc1(9);
  double *C_I = aloc1(9);
  double *Sbar = aloc1(9);
  double *ST = aloc1(nne_ve*ndn*ndn*ndn);
  double *LM = aloc1(ndn);
  double *L = aloc1(81);

  /* Integrate over the boundary element */
  int n_int_pt = int_pointC(nne_be);
  long n_int_ksi_2d = 0, n_int_eta_2d = 0;
  double *N2d = aloc1(nne_be);
  double *N3d = aloc1(nne_ve);
  double *N_lm = aloc1(n_lm_vecs);
  double *wts_2d = aloc1(n_int_pt);
  double *int_pt_ksi_2d = aloc1(n_int_pt);
  double *int_pt_eta_2d = aloc1(n_int_pt);

  /* get integration points and weights for a 2D element */
  integrateC(ptr_be->nnodes,&n_int_ksi_2d,&n_int_eta_2d,
	     int_pt_ksi_2d,int_pt_eta_2d,wts_2d);

  /* allocate for coordinate transformation */
  double *normal = aloc1(ndn);

  int ip = 0;
  for(int i=0; i<n_int_ksi_2d; i++){
    for(int j=0; j<n_int_eta_2d; j++){
      /* null quantites */
      nulld(F,9); nulld(C,9); nulld(C_I,9); nulld(Sbar,9);
      nulld(ST,ndn*ndn*ndn*nne_ve); nulld(N3d,nne_ve);
      nulld(normal,ndn); nulld(L,81);
      double J=0.0, jj=0.0, wt=0.0, ome=0.0;

      /* compute function values at integration point */
      err += bnd_integration_help(ptr_be,ptr_ve,ptr_dam,sup,ptrMat,ndofn,
				  ip,int_pt_ksi_2d,int_pt_eta_2d,wts_2d,
				  x_ve,y_ve,z_ve,KSI,ETA,ZET,ptr_disp,n_lm_vecs,
				  ptr_lm,x_be,y_be,z_be,kappa,LM,&ome,Sbar,
				  F,&J,C,C_I,ST,normal,N3d,N2d,N_lm,&jj,&wt);

      get_material_stiffness(kappa,ptrMat,ptr_dam,C,C_I,J,Sbar,L);
      disp_based_bnd_resid_at_ip(R,nne_ve,ndofn,n_lm_vecs,N_lm,ST,normal,
				 F,Sbar,L,LM,ome,jj,wt);

      ip++;
    }
  } /* end integration loop */

  free(KSI);
  free(ETA);
  free(ZET);

  free(x_be);
  free(y_be);
  free(z_be);

  free(F);
  free(C);
  free(C_I);
  free(Sbar);
  free(ST);
  free(LM);
  free(L);

  free(N2d);
  free(N3d);
  free(N_lm);
  free(wts_2d);
  free(int_pt_ksi_2d);
  free(int_pt_eta_2d);

  free(normal);

  return err;
} /* DISP_resid_bnd_el() */

int DISP_stiffmat_bnd_el(double *Ks,
			 const int b_el_id,
			 const int ndofn,
			 const int ndof_ve,
			 const double *x_ve,
			 const double *y_ve,
			 const double *z_ve,
			 const BOUNDING_ELEMENT *b_elems,
			 const ELEMENT *elem,
			 const HOMMAT *hommat,
			 const NODE *node,
			 const EPS *eps,
			 const SIG *sig,
			 const SUPP sup,
			 const double *vol_disp)
{
  /* NOTE: This function is coded for constant damage variable! Must
     be extended for general damage function for use with elements
     other than linear tetra. */

  int err = 0;

  /* Get const pointers */
  const BOUNDING_ELEMENT *ptr_be = &b_elems[b_el_id];
  const ELEMENT *ptr_ve = &elem[ptr_be->vol_elem_id];
  const HOMMAT *ptrMat = &hommat[ptr_ve->mat[2]];
  const damage *ptr_dam = &eps[ptr_be->vol_elem_id].dam[0];
  const long *loc_nod = ptr_be->loc_nodes;

  /* compute constant values */
  const double kappa = hommat_get_kappa(ptrMat);
  const int nne_be = ptr_be->nnodes;
  const int nne_ve = ptr_ve->toe;

  /* clear contents of Ks */
  memset(Ks,0,ndof_ve*ndof_ve*sizeof(double));

  if(DISP_DEBUG_NO_LM){
    /* set Kll to identity */
    for(int i=nne_ve*ndofn;i<ndof_ve; i++){
      Ks[idx_2_gen(i,i,ndof_ve,ndof_ve)] = 1.0;
    }
    return err;
  }

  /* get the volume element parent coordinates */
  double *KSI = aloc1(nne_ve);
  double *ETA = aloc1(nne_ve);
  double *ZET = aloc1(nne_ve);
  get_element_node_parent_coords(nne_ve,KSI,ETA,ZET);

  /* get the boundary element xyz coords */
  double *x_be = aloc1(nne_be);
  double *y_be = aloc1(nne_be);
  double *z_be = aloc1(nne_be);
  for(int i=0; i<nne_be; i++){
    int n_id = loc_nod[i];
    x_be[i] = x_ve[n_id];
    y_be[i] = y_ve[n_id];
    z_be[i] = z_ve[n_id];
  }

 /* how many LM "nodes" are there for interpolation */
  const int n_lm_vecs = ptr_be->n_dofs/ndn; 
  const int node_dof = nne_ve*ndofn;
  const double *ptr_disp = vol_disp;
  const double *ptr_lm = vol_disp + node_dof;

 /* compute values from volumetric element */
  double *F = aloc1(9);
  double *C = aloc1(9);
  double *C_I = aloc1(9);
  double *Sbar = aloc1(9);
  double *Lbar = aloc1(81);
  double *L = aloc1(81);
  double *K = aloc1(729);
  double *ST = aloc1(nne_ve*ndn*ndn*ndn);
  double *LM = aloc1(ndn);

  /* Integrate over the boundary element */
  int n_int_pt = int_pointC(nne_be);
  long n_int_ksi_2d = 0, n_int_eta_2d = 0;
  double *N2d = aloc1(nne_be);
  double *N3d = aloc1(nne_ve);
  double *N_lm = aloc1(n_lm_vecs);
  double *wts_2d = aloc1(n_int_pt);
  double *int_pt_ksi_2d = aloc1(n_int_pt);
  double *int_pt_eta_2d = aloc1(n_int_pt);

  /* get integration points and weights for a 2D element */
  integrateC(ptr_be->nnodes,&n_int_ksi_2d,&n_int_eta_2d,
	     int_pt_ksi_2d,int_pt_eta_2d,wts_2d);

  /* allocate for coordinate transformation */
  double *normal = aloc1(ndn);

  int ip = 0;
  for(int i=0; i<n_int_ksi_2d; i++){
    for(int j=0; j<n_int_eta_2d; j++){
      /* null quantites */
      nulld(F,9); nulld(C,9); nulld(C_I,9); nulld(Sbar,9);
      nulld(ST,ndn*ndn*ndn*nne_ve); nulld(N3d,nne_ve);
      nulld(normal,ndn); nulld(L,81); nulld(K,729);
      double J=0.0, jj=0.0, wt=0.0, ome=0.0;

      /* compute function values at integration point */
      err += bnd_integration_help(ptr_be,ptr_ve,ptr_dam,sup,ptrMat,ndofn,
				  ip,int_pt_ksi_2d,int_pt_eta_2d,wts_2d,
				  x_ve,y_ve,z_ve,KSI,ETA,ZET,ptr_disp,n_lm_vecs,
				  ptr_lm,x_be,y_be,z_be,kappa,LM,&ome,Sbar,F,&J,
				  C,C_I,ST,normal,N3d,N2d,N_lm,&jj,&wt);

      get_Lbar(kappa,ptrMat,C,C_I,J,Lbar);
      get_material_stiffness(kappa,ptrMat,ptr_dam,C,C_I,J,Sbar,L);
      get_material_sensitivity(kappa,ptrMat,ptr_dam,C,C_I,J,Sbar,Lbar,K);
      disp_based_bnd_tan_at_ip(Ks,nne_ve,ndofn,n_lm_vecs,N_lm,ST,
			       normal,F,Sbar,L,K,LM,ome,jj,wt);

      ip++;
    }
  } /* end integration loop */

  free(KSI);
  free(ETA);
  free(ZET);

  free(x_be);
  free(y_be);
  free(z_be);

  free(F);
  free(C);
  free(C_I);
  free(Sbar);
  free(Lbar);
  free(L);
  free(K);
  free(ST);
  free(LM);

  free(N2d);
  free(N3d);
  free(N_lm);
  free(wts_2d);
  free(int_pt_ksi_2d);
  free(int_pt_eta_2d);

  free(normal);

   return err;
} /* DISP_stiffmat_bnd_el() */

void DISP_increment_el(const ELEMENT *elem,
		       const int ii, /* id of element working on */
		       const int nne,
		       const NODE *node,
		       const long *nod,
		       const int ndofn,
		       const double *x,
		       const double *y,
		       const double *z,
		       EPS *eps,
		       SIG *sig,
		       const SUPP sup,
		       const HOMMAT *hommat,
		       const double *disp)
{
  const int mat = elem[ii].mat[2];
  const double kappa = hommat_get_kappa(&(hommat[mat]));

  /* compute volume */
  double nc_vol;
  if(nne == 4){ /* linear tet */
    nc_vol = Tetra_V(x,y,z);
  } else if(nne == 10){ /* quadradic tet */
    nc_vol = Tetra_qv_V(nne,ndofn,x,y,z);
  } else if(nne == 8){ /* trilinear hex */
    nc_vol = Hexa_V(x,y,z);
  } else nc_vol = 0.0;

  /* get constant variable */
  const double volume = nc_vol;

  /* Zero the stress & strain on the element */
  memset(sig[ii].el.o,0,6*sizeof(double));
  memset(eps[ii].el.o,0,6*sizeof(double));

  int ndofe = 0;
  for(int i=0; i<nne; i++){
    ndofe += ndofn;
  }

  double *F, *C, *C_I, *Sbar, J;
  F = aloc1(ndn*ndn);
  C = aloc1(ndn*ndn);
  C_I = aloc1(ndn*ndn);
  Sbar = aloc1(ndn*ndn);

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *N_x, *N_y, *N_z, *ST;
  Na = aloc1 (nne);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  N_z = aloc1 (nne);
  ST = aloc1(3*3*ndn*nne);

  double *sigma, *b, *b_I;
  sigma = aloc1(9);
  b = aloc1(9);
  b_I = aloc1(9);

  double *ident;
  ident = aloc1(9);
  ident[0] = ident[4] = ident[8] = 1.0;

  /*=== INTEGRATION LOOP ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
	    int_pt_ksi,int_pt_eta,int_pt_zet,
	    weights);
  int ip = 0;
  for(int i=0; i<npt_x; i++){
    for(int j=0; j<npt_y; j++){
      for(int k=0; k<npt_z; k++){
	double wt,jj;
	integration_help(ii,nne,ndofn,i,j,k,x,y,z,
			 int_pt_ksi,int_pt_eta,int_pt_zet,
			 weights,disp,sup,&wt,&jj,Na,N_x,N_y,N_z,ST,
			 F,C,C_I,&J);

	get_material_stress(kappa,&hommat[mat],C,C_I,J,Sbar);
	/* store deformation gradient */
	memcpy(eps[ii].il[ip].F,F,9*sizeof(double));

	/* store the potential */
	DISP_get_material_potential(kappa,&hommat[mat],C,J,&eps[ii].il[ip].Y);

	/* Compute the Cauchy Stress sigma = (1-w)/J F Sbar F' */
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		    3,3,3,(1.-eps[ii].dam[ip].w)/J,F,3,Sbar,3,0.0,b,3);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
		    3,3,3,1.0,b,3,F,3,0.0,sigma,3);

	/* store S at ip */
	sig[ii].il[ip].o[0] = (1.-eps[ii].dam[ip].w)*Sbar[idx_2(0,0)]; /* XX */
	sig[ii].il[ip].o[1] = (1.-eps[ii].dam[ip].w)*Sbar[idx_2(1,1)]; /* YY */
	sig[ii].il[ip].o[2] = (1.-eps[ii].dam[ip].w)*Sbar[idx_2(2,2)]; /* ZZ */
	sig[ii].il[ip].o[3] = (1.-eps[ii].dam[ip].w)*Sbar[idx_2(1,2)]; /* YZ */
	sig[ii].il[ip].o[4] = (1.-eps[ii].dam[ip].w)*Sbar[idx_2(0,2)]; /* XZ */
	sig[ii].il[ip].o[5] = (1.-eps[ii].dam[ip].w)*Sbar[idx_2(0,1)]; /* XY */

	/* avg(sig) = 1/V int(sig) */
	sig[ii].el.o[0] += (wt*jj/volume*sigma[idx_2(0,0)]); /* XX */
	sig[ii].el.o[1] += (wt*jj/volume*sigma[idx_2(1,1)]); /* YY */
	sig[ii].el.o[2] += (wt*jj/volume*sigma[idx_2(2,2)]); /* ZZ */
	sig[ii].el.o[3] += (wt*jj/volume*sigma[idx_2(1,2)]); /* YZ */
	sig[ii].el.o[4] += (wt*jj/volume*sigma[idx_2(0,2)]); /* XZ */
	sig[ii].el.o[5] += (wt*jj/volume*sigma[idx_2(0,1)]); /* XY */

	/* Compute the logarithmic strain e = 1/2(I - inv(FF'))*/
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
		    3,3,3,1.0,F,3,F,3,0.0,b,3);
	inverse(b,3,b_I);
	/* e <-- b is the Euler strain */
	for(int i=0; i<9; i++){
	  b[i] = 0.5*(ident[i]-b_I[i]);
	}

	/* avg(e) = 1/V int(e) */
	eps[ii].el.o[0] += wt*jj/volume*b[idx_2(0,0)];
	eps[ii].el.o[1] += wt*jj/volume*b[idx_2(1,1)];
	eps[ii].el.o[2] += wt*jj/volume*b[idx_2(2,2)];
	
	eps[ii].el.o[3] += 2.*wt*jj/volume*b[idx_2(1,2)];
	eps[ii].el.o[4] += 2.*wt*jj/volume*b[idx_2(0,2)];
	eps[ii].el.o[5] += 2.*wt*jj/volume*b[idx_2(0,1)];

	/* This is TOTAL LAGRANGIAN formulation. There is no need to
	   store deformation gradients etc. */

	/* update damage variables */
	update_damage(&eps[ii].dam[ip]);

        /* update constitutive model */
        if (eps[ii].model != NULL) {
          eps[ii].model[ip].param->update_state_vars(&eps[ii].model[ip]);
        }

	ip++;
      }
    }
  }

  free(F);
  free(C);
  free(C_I);
  free(Sbar);

  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);

  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  free(ST);
  free(sigma);
  free(b);
  free(b_I);
  free(ident);

}/* DISP_increment_el() */

void DISP_increment(const ELEMENT *elem,
		    const int nelem,
		    NODE *node,
		    const int nnodes,
		    const int ndofn,
		    SUPP sup,
		    EPS *eps,
		    SIG *sig,
		    const HOMMAT *hommat,
		    const double *sol_incr,
		    const double *sol,
		    MPI_Comm mpi_comm,
		    const int mp_id)
{
  /* for each element */
  for (int i=0; i<nelem; i++){
    const ELEMENT *ptr_elem = &elem[i];
    const int nne = ptr_elem->toe;
    const int ndofe = nne*ndofn;

    /* allocate */
    long *cn = aloc1l(ndofe);
    double *x = aloc1(nne);
    double *y = aloc1(nne);
    double *z = aloc1(nne);
    double *disp_incr = aloc1(ndofe);
    double *disp = aloc1(ndofe);

    /* get node ids on element */
    long *nod = ptr_elem->nod;

    /* get local dof ids on elemnt */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* Get the node and bubble coordinates */
    nodecoord_total(nne,nod,node,x,y,z);

    /* TOTAL LAGRANGIAN: get the TOTAL displacement */
    def_elem_total(cn,ndofe,sol,sol_incr,elem,node,sup,disp);

    /* increment the element quantities */
    DISP_increment_el(elem,i,nne,node,nod,ndofn,
		      x,y,z,eps,sig,sup,hommat,disp);

    free(cn);
    free(x);
    free(y);
    free(z);
    free(disp_incr);
    free(disp);
  } /* For each element */

  /*** Update the coordinates FOR COMPUTING CURRENT VOLUME ONLY ***/
  /* Total Lagrangian formaulation takes gradients w.r.t. xi_fd */
  for (int ii=0;ii<nnodes;ii++){
    for (int i=0;i<ndn;i++){
      int II = node[ii].id_map[mp_id].id[i];
      if (II > 0){
	if (i == 0) node[ii].x1 = node[ii].x1_fd + sol[II-1] + sol_incr[II-1];
	else if (i == 1) node[ii].x2 = node[ii].x2_fd + sol[II-1] + sol_incr[II-1];
	else if (i == 2) node[ii].x3 = node[ii].x3_fd + sol[II-1] + sol_incr[II-1];
      }
      else if (II < 0){
	if (i == 0) node[ii].x1 = (node[ii].x1_fd
				   + sup->defl[abs(II)-1]
				   + sup->defl_d[abs(II)-1]);

	else if (i == 1) node[ii].x2 = (node[ii].x2_fd
					+ sup->defl[abs(II)-1]
					+ sup->defl_d[abs(II)-1]);

	else if (i == 2) node[ii].x3 = (node[ii].x3_fd
					+ sup->defl[abs(II)-1]
					+ sup->defl_d[abs(II)-1]);
      }
    }
  }/* end ii < nn */

  /* update the coordinates including macroscale deformations */
  if(sup->multi_scale){
    const double *F = sup->F0;
    double X[ndn];
    double Y[ndn];
    for(int i=0; i<nnodes; i++){
      X[0] = node[i].x1_fd;
      X[1] = node[i].x2_fd;
      X[2] = node[i].x3_fd;
      cblas_dgemv(CblasRowMajor,CblasNoTrans,ndn,ndn,1.0,
		  F,ndn,X,1,0.0,Y,1);
      node[i].x1 += Y[0];
      node[i].x2 += Y[1];
      node[i].x3 += Y[2];
    }
  }

  double PL, GPL;
  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);
  PL = T_VOLUME (nelem,ndn,elem,node);
  /* Gather Volume from all domains */
  MPI_Allreduce(&PL,&GPL,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  if (myrank == 0) {
    PGFEM_printf ("AFTER DEF - VOLUME = %12.12f\n",GPL);
  }
}/* DISP_increment() */

int DISP_cohe_micro_terms_el(double *K_00_e,
			     double *K_01_e,
			     double *K_10_e,
			     double *traction_res_e,
			     double *traction_e,
			     /* macro information */
			     const int macro_nnode,
			     const int macro_ndofn,
			     const double macro_int_wt,
			     const double *macro_shape_func,
			     const double *macro_normal,
			     const double layer_thickness,
			     /* micro information */
			     const double micro_volume0,
			     const int elem_id,
			     const int ndofn,
			     const int nne,
			     const double *x,
			     const double *y,
			     const double *z,
			     const ELEMENT *elem,
			     const HOMMAT *hommat,
			     const long *nod,
			     const NODE *node,
			     const EPS *eps,
			     const SIG *sig,
			     const SUPP sup,
			     const double *disp,
                             const double dt)
{
  int err = 0;
  const ELEMENT *p_elem = elem + elem_id;
  const HOMMAT *p_hommat = &hommat[p_elem->mat[2]];
  const double kappa = hommat_get_kappa(p_hommat);
  const int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);
  const int macro_ndof = macro_nnode*macro_ndofn;

  /* ensure null values of matrices */
  nulld(K_00_e,macro_ndof*macro_ndof);
  nulld(K_01_e,macro_ndof*ndofe);
  nulld(K_10_e,ndofe*macro_ndof);
  nulld(traction_res_e,macro_ndof);

  /* set up integration */
  long int_pt_x = 0;
  long int_pt_y = 0;
  long int_pt_z = 0;
  /* compute number of integration points -> int_pt_x */
  int_point(nne,&int_pt_x);

  /* allocate information for integration */
  double *int_pt_ksi = PGFEM_calloc(int_pt_x,sizeof(double));
  double *int_pt_eta = PGFEM_calloc(int_pt_x,sizeof(double));
  double *int_pt_zet = PGFEM_calloc(int_pt_x,sizeof(double));
  double *weights = PGFEM_calloc(int_pt_x,sizeof(double));
  double *Na = PGFEM_calloc(nne,sizeof(double));
  double *N_x = PGFEM_calloc(nne,sizeof(double));
  double *N_y = PGFEM_calloc(nne,sizeof(double));
  double *N_z = PGFEM_calloc(nne,sizeof(double));

  /* allocate variables */
  double *ST = PGFEM_calloc(nne*ndn*ndn*ndn,sizeof(double));
  double *F = PGFEM_calloc(ndn*ndn,sizeof(double));
  double *C = PGFEM_calloc(ndn*ndn,sizeof(double));
  double *C_I = PGFEM_calloc(ndn*ndn,sizeof(double));
  double *Sbar = PGFEM_calloc(ndn*ndn,sizeof(double));
  double *L = PGFEM_calloc(ndn*ndn*ndn*ndn,sizeof(double));

  /* compute constant terms */
  /*== Indexing of gNoxN ==
    idx_4_gen(a,b,i,j,macro_nnode,macro_ndofn,ndn,ndn)
    ==*/
  double *gNoxN = PGFEM_calloc(macro_ndof*ndn*ndn,sizeof(double));
  {
    for(int a=0; a<macro_nnode; a++){
      for(int b=0; b<macro_ndofn; b++){
	for(int i=0; i<ndn; i++){
	  if(i != b) continue; /* N^a_{i,b} = N^a d_{i,b} */
	  for(int j=0; j<ndn; j++){
	    const int idx = idx_4_gen(a,b,i,j,macro_nnode,
				      macro_ndofn,ndn,ndn);
	    gNoxN[idx] += macro_shape_func[a]*macro_normal[j];
	  }
	}
      }
    }
  }

  /* integration loop */
  {
    integrate(nne,&int_pt_x,&int_pt_y,&int_pt_z,
	      int_pt_ksi,int_pt_eta,int_pt_zet,weights);
    int ip = 0;
    for(int i=0; i<int_pt_x; i++){
      for(int j=0; j<int_pt_y; j++){
	for(int k=0; k<int_pt_z; k++){
	  /* reset variables */
	  double J = 0.0;
	  double jj = 0.0;
	  double wt = 0.0;
	  nulld(Na,nne);
	  nulld(N_x,nne);
	  nulld(N_y,nne);
	  nulld(N_z,nne);
	  nulld(ST,nne*ndn*ndn*ndn);
	  nulld(F,ndn*ndn);
	  nulld(C,ndn*ndn);
	  nulld(C_I,ndn*ndn);
	  nulld(Sbar,ndn*ndn);
	  nulld(L,ndn*ndn*ndn*ndn);

	  err += integration_help(elem_id,nne,ndofn,i,j,k,x,y,z,
				  int_pt_ksi,int_pt_eta,int_pt_zet,
				  weights,disp,sup,&wt,&jj,Na,
				  N_x,N_y,N_z,ST,F,C,C_I,&J);

	  /* inverted element detected, exit and return error */
	  if(err != 0 ) goto exit_function;

	  /* get material stress and stiffness */
          const damage *p_dam = NULL;
          if (eps[elem_id].model == NULL) {
            p_dam = &(eps[elem_id].dam[ip]);
            get_material_stress(kappa,p_hommat,C,C_I,J,Sbar);
            get_material_stiffness(kappa,p_hommat,p_dam,C,C_I,J,Sbar,L);
          } else {
            p_dam = &empty_damage;
            err += disp_cm_material_response(Sbar, L, eps[elem_id].model + ip,
                                             F, dt, 1);
          }

	  /*=== OPTIMIZATION: Combine following function calls into
	    one loop. This was not done originally for simplified
	    testing/proof of concept code */

	  err += compute_K_00_e_at_ip_2(K_00_e,macro_nnode,macro_ndofn,
					macro_int_wt,macro_shape_func,
					macro_normal,layer_thickness,
					gNoxN,nne,micro_volume0,
					jj,wt,F,ST,p_dam,Sbar,L);

	  /* compute K_01_e */
	  err += compute_K_01_e_at_ip(K_01_e,macro_nnode,macro_ndofn,
				      macro_int_wt,macro_shape_func,
				      macro_normal,gNoxN,nne,micro_volume0,
				      jj,wt,F,ST,p_dam,Sbar,L);
	  /* compute K_10_e */
	  err += compute_K_10_e_at_ip(K_10_e,macro_nnode,macro_ndofn,
				      macro_int_wt,macro_shape_func,
				      macro_normal,gNoxN,nne,micro_volume0,
				      jj,wt,F,ST,p_dam,Sbar,L);

	  /* compute traction_res_e */
	  {
	    double *P = PGFEM_calloc(ndn*ndn,sizeof(double));
	    /* compute 1st PK stress */
	    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
			3,3,3,(1-p_dam->w),F,3,Sbar,3,0.0,P,3);

	    /* compute traction */
	    for(int I=0; I<ndn; I++){
	      traction_e[I] = 0.0;
	      for(int J=0; J<ndn; J++){
		traction_e[I] += (1./micro_volume0*jj*wt
				  *P[idx_2(I,J)]*macro_normal[J]);
	      }
	    }

	    for(int a=0; a<macro_nnode; a++){
	      for(int b=0; b<macro_ndofn; b++){
		const int idx = idx_2_gen(a,b,macro_nnode,macro_ndofn);
		for(int I=0; I<ndn; I++){
		  if(I != b) continue;
		  for(int J=0; J<ndn; J++){
		    traction_res_e[idx] +=
		      1./micro_volume0*jj*wt*macro_int_wt
		      *P[idx_2(I,J)]*macro_normal[J]*macro_shape_func[a];
		  }
		}
	      }
	    }

	    free(P);
	  }

	  /* end of integration loop */
	}
      }
    }
    /* end of integration loop */
  }

 exit_function:
  /* clean up */
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(N_x);
  free(N_y);
  free(N_z);
  free(ST);
  free(F);
  free(C);
  free(C_I);
  free(gNoxN);
  free(Sbar);
  free(L);

  return err;
}/* DISP_cohe_micro_terms_el */

/*=========================================================*/
/*             DEFINE  LOCAL FUNCTIONS                     */
/*=========================================================*/

static void get_material_stress(const double kappa,
				const HOMMAT *hommat, /* pointer to 1 hommat */
				const double *C,
				const double *C_I,
				const double J,
				double *Sbar)
{
  devStressFuncPtr Stress = getDevStressFunc(0,hommat);
  dUdJFuncPtr DUDJ = getDUdJFunc(0,hommat);
  double dUdJ = 0.0;

  Stress(C,hommat,Sbar);
  DUDJ(J,hommat,&dUdJ);

  cblas_daxpy(9,kappa*J*dUdJ,C_I,1,Sbar,1);

} /* get_material_stress() */

static void get_material_stiffness(const double kappa,
				   const HOMMAT *hommat, /* pointer to 1 hommat */
				   const damage *dam,
				   const double *C,
				   const double *C_I,
				   const double J,
				   const double *Sbar,
				   double *L)
{
  double *CIoxCI, *CICI, *SoxS;
  CIoxCI = (double*) PGFEM_calloc (81,sizeof(double));
  CICI = (double*) PGFEM_calloc (81,sizeof(double));
  SoxS = (double*) PGFEM_calloc (81,sizeof(double));

  /* Get the potential stuff */
  double H = 0.0;
  double dUdJ = 0.0;
  double d2UdJ2 = 0.0;
  matStiffFuncPtr Stiff = getMatStiffFunc(0,hommat);
  dUdJFuncPtr DUDJ = getDUdJFunc(0,hommat);
  d2UdJ2FuncPtr D2UDJ2 = getD2UdJ2Func(0,hommat);

  Stiff(C,hommat,L);
  DUDJ(J,hommat,&dUdJ); 
  D2UDJ2(J,hommat,&d2UdJ2);
 
  /* if(dam->damaged){ */
  /*   DISP_get_material_potential(kappa,hommat,C,J,&Ybar); */
  /*   H = dam->dmu/(1+dam->dmu)*dam->evolution(Ybar,&(dam->params)); */
  /* } */

  /* use stored damage evolution parameter */
  if(dam->damaged){
    H = dam->dmu/(1+dam->dmu)*dam->H;
  }

  /* compute CIoxCI, CICI, & SoxS */
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	for(int l=0; l<3; l++){
	  CIoxCI[idx_4(i,j,k,l)] = C_I[idx_2(i,j)]*C_I[idx_2(k,l)];
	  SoxS[idx_4(i,j,k,l)] = Sbar[idx_2(i,j)]*Sbar[idx_2(k,l)];
	  CICI[idx_4(i,j,k,l)] = C_I[idx_2(i,k)]*C_I[idx_2(l,j)];
	}
      }
    }
  }
  
  for(int i=0; i<81; i++){
    L[i] = ((1-dam->w)*(
			L[i]
			+ kappa*(J*dUdJ + J*J*d2UdJ2)*CIoxCI[i]
			- 2.0*kappa*J*dUdJ*CICI[i])
	    -H*SoxS[i]);
  }

  free(CIoxCI);
  free(CICI);
  free(SoxS);
}/* get_material_stiffness() */

static void get_Lbar(const double kappa,
		     const HOMMAT *p_mat,
		     const double *C,
		     const double *C_I,
		     const double J,
		     double *Lbar)
{
  double dUdJ = 0.0;
  double d2UdJ2 = 0.0;
  matStiffFuncPtr Stiff = getMatStiffFunc(0,p_mat);
  dUdJFuncPtr DUDJ = getDUdJFunc(0,p_mat);
  d2UdJ2FuncPtr D2UDJ2 = getD2UdJ2Func(0,p_mat);
  Stiff(C,p_mat,Lbar);
  DUDJ(J,p_mat,&dUdJ); 
  D2UDJ2(J,p_mat,&d2UdJ2);
  for(int i=0; i<ndn; i++){
    for(int j=0; j<ndn; j++){
      for(int k=0; k<ndn; k++){
	for(int l=0; l<ndn; l++){
	  const int idx4 = idx_4(i,j,k,l);
	  /* Deviatoric + Volumetric stiffness */
	  Lbar[idx4] += ((kappa*J*(dUdJ+J*d2UdJ2)*C_I[idx_2(i,j)]*C_I[idx_2(k,l)])
			 -(2.*kappa*J*dUdJ*C_I[idx_2(i,k)]*C_I[idx_2(l,j)]));
	}
      }
    }
  }
}/* get_Lbar() */

static void get_material_sensitivity(const double kappa,
				     const HOMMAT *p_mat,
				     const damage *p_dam,
				     const double *C,
				     const double *C_I,
				     const double J,
				     const double *Sbar,
				     const double *Lbar,
				     double *K)
{
  /* Get the potential stuff */
  double H = 0.0;
  double Hp = 0.0;
  double dUdJ = 0.0;
  double d2UdJ2 = 0.0;
  double d3UdJ3 = 0.0;
  matSensFuncPtr Sens = getMatSensFunc(0,p_mat);
  dUdJFuncPtr DUDJ = getDUdJFunc(0,p_mat);
  d2UdJ2FuncPtr D2UDJ2 = getD2UdJ2Func(0,p_mat);
  d3UdJ3FuncPtr D3UDJ3 = getD3UdJ3Func(0,p_mat);

  Sens(C,p_mat,K);
  DUDJ(J,p_mat,&dUdJ); 
  D2UDJ2(J,p_mat,&d2UdJ2);
  D3UDJ3(J,p_mat,&d3UdJ3);

  /* compute tensor products */
  if(p_dam->damaged){
    H = p_dam->dmu/(1+p_dam->dmu)*p_dam->H;
    Hp = p_dam->dmu/(1+p_dam->dmu)*p_dam->Hp;
  }
  for(int i=0; i<ndn; i++){
    for(int j=0; j<ndn; j++){
      const int ij = idx_2(i,j);

      for(int k=0; k<ndn; k++){
	const int ik = idx_2(i,k);

	for(int l=0; l<ndn; l++){
	  const int kl = idx_2(k,l);
	  const int lj = idx_2(l,j);
	  const int ijkl = idx_4(i,j,k,l);

	  for(int r=0; r<ndn; r++){
	    const int ir = idx_2(i,r);
	    const int kr = idx_2(k,r);
	    const int lr = idx_2(l,r);

	    for(int s=0; s<ndn; s++){
	      const int sj = idx_2(s,j);
	      const int sk = idx_2(s,k);
	      const int sl = idx_2(s,l);
	      const int rs = idx_2(r,s);
	      const int ijrs = idx_4(i,j,r,s);
	      const int klrs = idx_4(k,l,r,s);
	      const int ijklrs = idx_6(i,j,k,l,r,s);

	      const double K_vol = 2.*(
				       (0.5*kappa*J*(dUdJ+3*J*d2UdJ2+J*J*d3UdJ3)
					*C_I[ij]*C_I[kl]*C_I[rs])
				       
				       -(kappa*J*(dUdJ+J*d2UdJ2)
					 *(C_I[ir]*C_I[sj]*C_I[kl]
					   + C_I[ij]*C_I[kr]*C_I[sl]
					   + C_I[ik]*C_I[lj]*C_I[rs]))
				       
				       +2.*kappa*J*dUdJ*(C_I[ir]*C_I[sk]*C_I[lj]
							 + C_I[ik]*C_I[lr]*C_I[sj])
				       );
	      const double damage_term = (-H*(Lbar[ijkl]*Sbar[rs]
					      + Lbar[ijrs]*Sbar[kl]
					      + Sbar[ij]*Lbar[klrs])
					  - Hp*(Sbar[ij]*Sbar[kl]*Sbar[rs]));

	      K[ijklrs] = ((1.-p_dam->w)*(K[ijklrs] + K_vol)
			   + damage_term);

	    }
	  }
	}
      }
    }
  }

}/* get_material_sensitivity() */

static void disp_based_resid_at_ip(double *R,
				   const int nne,
				   const double *ST,
				   const double *F,
				   const double *Sbar,
				   const damage *dam,
				   const double jj,
				   const double wt)
{
  double *AA;
  AA = aloc1(9);

  for(int a=0; a<nne; a++){
    for(int b=0; b<ndn; b++){
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,
						   nne,ndn,ndn,ndn)];
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		  3,3,3,1.0,F,3,ptrST_ab,3,0.0,AA,3);

      R[a*ndn+b] += ((1 - dam->w)*cblas_ddot(9,Sbar,1,AA,1)*jj*wt);
    }
  }

  free(AA);
} /* disp_based_resid_at_ip() */

static void disp_based_bnd_resid_at_ip(double *R,
				       const int nne_ve,
				       const int ndofn,
				       const int n_lm_vecs,
				       const double *N_lm,
				       const double *ST,
				       const double *normal,
				       const double *F,
				       const double *Sbar,
				       const double *L,
				       const double *LM,
				       const double ome,
				       const double jj,
				       const double wt)
{
  double *FSTs = aloc1(9);

  /* compute R_lm */
  double *R_lm = R + nne_ve*ndofn;
  for(int a=0; a<n_lm_vecs; a++){
    if(N_lm[a] == 0) continue;
    for(int b=0; b<ndn; b++){
      double result = 0.;
      for(int i=0; i<ndn; i++){
	for(int j=0; j<ndn; j++){
	  for(int k=0; k<ndn; k++){
	    result += ((1.-ome)*(F[idx_2(i,k)]*Sbar[idx_2(k,j)]
				 *normal[j])*N_lm[a]*del(i,b));
	  }
	}
      }
      R_lm[a*ndn+b] += result*jj*wt;
      /* R_lm[a*ndn+b] -= result*jj*wt; */
    }
  }

  /* compute res */
  for(int a=0; a<nne_ve; a++){
    for(int b=0; b<ndofn; b++){
      const double *ptr_ST = &ST[idx_4_gen(a,b,0,0,nne_ve,ndn,ndn,ndn)];
      /* sym(F'ptr_ST) */
      for(int i=0; i<ndn; i++){
	for(int j=0; j<ndn; j++){
	  const int idx_ij = idx_2(i,j);
	  FSTs[idx_ij] = 0.0;
	  for(int k=0; k<ndn; k++){
	    FSTs[idx_ij] += 0.5*(F[idx_2(k,i)]*ptr_ST[idx_2(k,j)]
				 + ptr_ST[idx_2(k,i)]*F[idx_2(k,j)]);
	  }
	}
      }

      double term_I = 0.0;
      double term_II = 0.0;
      for(int i=0; i<ndn; i++){
	for(int j=0; j<ndn; j++){
	  for(int k=0; k<ndn; k++){
	    const int ik = idx_2(i,k);
	    const int kj = idx_2(k,j);
	    term_I += LM[i]*((1-ome)*ptr_ST[ik]*Sbar[kj])*normal[j];
	    for(int l=0; l<ndn; l++){
	      for(int m=0; m<ndn; m++){
		const int kjlm = idx_4(k,j,l,m);
		const int lm = idx_2(l,m);
		term_II += LM[i]*(F[ik]*L[kjlm]*FSTs[lm])*normal[j];
	      }
	    }
	  }
	}
      }

      R[a*ndofn+b] += (term_I+term_II)*jj*wt;
    }
  }

  free(FSTs);
}/* disp_based_bnd_resid_at_ip() */

static void disp_based_tan_at_ip(double *K,
				 const int nne,
				 const double *ST,
				 const double *F,
				 const double *Sbar,
				 const double *L,
				 const damage *dam,
				 const double jj,
				 const double wt)
{
  Tensor<2, 3, double> AA, BB, CC, DD, LA_wg;
  Tensor<2, 3, const double*> F_ttl(F);  // convert F and L to ttl tensors
  Tensor<4, 3, const double*> L_ttl(L);

  /* POSSIBLE OPTIMIZATION: unroll loops and change abwg order */

  for(int a=0; a<nne; a++){
    for(int b=0; b<ndn; b++){
      const Tensor<2, 3, const double*> ptrST_ab(&ST[idx_4_gen(a,b,0,0,
        				 	    nne,ndn,ndn,ndn)]);

      /* AA = sym(F'Grad(del u)) */
      BB(i,j) = F_ttl(k,i).to(i,k) * ptrST_ab(k,j);
      AA(i,j) = .5 * (BB(i,j) + BB(j,i).to(i,j));

      for(int w=0; w<nne; w++){
	for(int g=0; g<ndn; g++){
	  const Tensor<2, 3, const double*> ptrST_wg(&ST[idx_4_gen(w,g,0,0,
	  						nne,ndn,ndn,ndn)]);
	 
	  /* BB = F' * ST_wg */
	  BB(i,j) = F_ttl(k,i).to(i,k) * ptrST_wg(k,j);
	  /* CC = ST_wg' * ST_ab */
	  CC(i,j) = ptrST_wg(k,i).to(i,k) * ptrST_ab(k,j);

	  /* DD = sym(BB) */
	  DD = .5 * (BB(i,j) + BB(j,i).to(i,j));

	  /* compute LA_wg = L:sym(F'ST_wg) */
	  LA_wg(i,j) = L_ttl(i,j,k,l) * DD(k,l);

	  const int K_idx = idx_K(a,b,w,g,nne,ndn);
	  for(int i=0; i<9; i++){
	    K[K_idx] += jj*wt*(AA.data[i]*LA_wg.data[i] + (1-dam->w)*Sbar[i]*CC.data[i]);
	  }
	}
      }
    }
  }
}/* disp_based_tan_at_ip() */

static void disp_based_bnd_tan_at_ip(double *K,
				     const int nn_ve,
				     const int ndofn,
				     const int n_lm_vecs,
				     const double *N_lm,
				     const double *ST,
				     const double *normal,
				     const double *F,
				     const double *Sbar,
				     const double *L,
				     const double *KK,
				     const double *LM,
				     const double ome,
				     const double jj,
				     const double wt)
{
  const int node_dof = nn_ve*ndofn;
  const int lm_dof = n_lm_vecs*ndn;
  const int ndofe = node_dof + lm_dof;

  /* construct tangent by blocks */
  double *Kuu = aloc1(node_dof * node_dof);
  double *Kul = aloc1(node_dof * lm_dof);
  double *Klu = aloc1(lm_dof * node_dof);

  /*=== compute Kuu ===*/
  disp_based_bnd_Kuu_at_ip(Kuu,nn_ve,ndofn,ST,normal,
			   F,L,KK,LM,ome,jj,wt);

  /*=== compute Kul & Klu ===*/
  disp_based_bnd_Kul_Klu_at_ip(Kul,Klu,nn_ve,ndofn,n_lm_vecs,
			       N_lm,ST,normal,F,Sbar,L,ome,jj,wt);

  /*=== Assemble the stiffness matrix from blocks ===*/
  for(int row=0; row<ndofe; row++){
    for(int col=0; col<ndofe; col++){
      const int idx = idx_2_gen(row,col,ndofe,ndofe);
      if(row<node_dof && col<node_dof){ /* Kuu */
  	const int idx_Kuu = idx_2_gen(row,col,node_dof,node_dof);
  	K[idx] += Kuu[idx_Kuu];
      } else if(row<node_dof && col >=node_dof){ /* Kul */
  	const int idx_Kul = idx_2_gen(row,col-node_dof,node_dof,lm_dof);
  	K[idx] += Kul[idx_Kul];
      } else if(row>=node_dof && col<node_dof){ /* Klu */
  	const int idx_Klu = idx_2_gen(row-node_dof,col,lm_dof,node_dof);
  	K[idx] += Klu[idx_Klu];
  	/* K[idx] -= Klu[idx_Klu]; */
      } else if(row>=node_dof && col>=node_dof){ /* Kll = NULL */
	continue;
      } else {
	PGFEM_printerr("Entered undefinded branch in %s\n",__func__);
      }
    }
  }

  free(Kuu);
  free(Kul);
  free(Klu);
} /* disp_based_bnd_tan_at_ip() */

static int integration_help(const int elem_id,
			    const int nne,
			    const int ndofn,
			    const int i,
			    const int j,
			    const int k,
			    const double *x,
			    const double *y,
			    const double *z,
			    const double *int_pt_ksi,
			    const double *int_pt_eta,
			    const double *int_pt_zet,
			    const double *weights,
			    const double *disp,
			    const SUPP sup,
			    double *wt,
			    double *jj,			     
			    double *Na,
			    double *N_x,
			    double *N_y,
			    double *N_z,
			    double *ST,
			    double *Fr,
			    double *C,
			    double *C_I,
			    double *J)
{
  double ksi,eta,zet;
  int err;

  const int ip = k;
  if(nne == 8){/* hexahedron */
    ksi = int_pt_ksi[i];
    eta = int_pt_ksi[j];
    zet = int_pt_ksi[k];
    (*wt) = weights[i]*weights[j]*weights[k];
  } else { /* tetrahedron type */
    ksi = int_pt_ksi[ip];
    eta = int_pt_eta[ip];
    zet = int_pt_zet[ip];
    (*wt) = weights[ip];
  }

  double ****ST_tensor, **Fr_mat;
  ST_tensor = aloc4 (3,3,ndn,nne);
  Fr_mat = aloc2(3,3);

  shape_func(ksi,eta,zet,nne,Na);
  (*jj) = deriv(ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
  shape_tensor (nne,ndn,N_x,N_y,N_z,ST_tensor);
  shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);

  def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,disp,Fr_mat);
  mat2array(Fr,CONST_2(double) Fr_mat,3,3);

  /* add the macroscopic contribution to F */
  if(sup->multi_scale){
    cblas_daxpy(9,1.0,sup->F0,1,Fr,1);
  }

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
	      3,3,3,1.0,Fr,3,Fr,3,0.0,C,3);

  inverse(C,3,C_I);
  *J = getJacobian(Fr,elem_id,&err);

  dealoc4(ST_tensor,3,3,ndn);
  dealoc2(Fr_mat,3);

  return err;
}/* integration_help */

static int bnd_integration_help(const BOUNDING_ELEMENT *ptr_be,
				const ELEMENT *ptr_ve,
				const damage *ptr_dam,
				const SUPP sup,
				const HOMMAT *ptr_mat,
				const int ndofn,
				const int ip_2d,
				const double *int_pt_ksi_2d,
				const double *int_pt_eta_2d,
				const double *wts_2d,
				const double *x_ve,
				const double *y_ve,
				const double *z_ve,
				const double *KSI,
				const double *ETA,
				const double *ZET,
				const double *ve_disp,
				const int n_lm_vecs,
				const double *elem_lm,
				const double *x_be,
				const double *y_be,
				const double *z_be,
				const double kappa,
				double *LM,
				double *ome,
				double *Sbar,
				double *F,
				double *J,
				double *C,
				double *C_I,
				double *ST,
				double *normal,
				double *N3d,
				double *N2d,
				double *N_lm,
				double *jj,
				double *wt)
{
  int err = 0;

  const int nne_2d = ptr_be->nnodes;
  const int nne_3d = ptr_ve->toe;
  const long *loc_nod = ptr_be->loc_nodes;

  /* get simple quantites */
  *ome = ptr_dam->w;
  *wt = wts_2d[ip_2d];

  {
    /* compute 2D transformation and jj */
    double *e1 = aloc1(ndn);
    double *e2h = aloc1(ndn);
    double *e2 = aloc1(ndn);
    double *xl = aloc1(nne_2d);
    double *yl = aloc1(nne_2d);
    double *zl = aloc1(nne_2d);
    double *N_x = aloc1(nne_2d);
    double *N_y = aloc1(nne_2d);
    base_vec(nne_2d,int_pt_ksi_2d[ip_2d],int_pt_eta_2d[ip_2d],
	     x_be,y_be,z_be,e1,e2,e2h,normal,0);
    transform_coordinates(nne_2d,x_be,y_be,z_be,e1,e2,normal,0,xl,yl,zl);
    *jj = dN3_xy(int_pt_ksi_2d[ip_2d],int_pt_eta_2d[ip_2d],nne_2d,xl,yl,zl,N_x,N_y);
    free(e1);
    free(e2h);
    free(e2);
    free(xl);
    free(yl);
    free(zl);
    free(N_x);
    free(N_y);
  }

  /* compute shape functions for 2D element to get ksi eta zet
     coordinates for 3D element. */
  shape_2DC(nne_2d,int_pt_ksi_2d[ip_2d],int_pt_eta_2d[ip_2d],N2d);
  double ksi = 0.0, eta = 0.0, zet = 0.0;
  for(int i=0; i<nne_2d; i++){
    ksi += N2d[i]*KSI[loc_nod[i]];
    eta += N2d[i]*ETA[loc_nod[i]];
    zet += N2d[i]*ZET[loc_nod[i]];
  }

  /* compute the Lagrange multiplier shape function and interpolate
     the Lagrange multiplier value at the integration point */
  shape_2DC(n_lm_vecs,int_pt_ksi_2d[ip_2d],int_pt_eta_2d[ip_2d],N_lm);
  nulld(LM,ndn);
  for(int i=0; i<n_lm_vecs; i++){
    LM[0] += N_lm[i]*elem_lm[i*ndn+0];
    LM[1] += N_lm[i]*elem_lm[i*ndn+1];
    LM[2] += N_lm[i]*elem_lm[i*ndn+2];
  }

  /* begin evaluating functions using 3D shape functions */
  double *N_x = aloc1(nne_3d);
  double *N_y = aloc1(nne_3d);
  double *N_z = aloc1(nne_3d);
  double ****ST_tensor = aloc4(ndn,ndn,ndn,nne_3d);
  double **F_mat = aloc2(ndn,ndn);

  shape_func(ksi,eta,zet,nne_3d,N3d);
  deriv(ksi,eta,zet,nne_3d,x_ve,y_ve,z_ve,N_x,N_y,N_z);
  shape_tensor(nne_3d,ndn,N_x,N_y,N_z,ST_tensor);
  shapeTensor2array(ST,CONST_4(double) ST_tensor,nne_3d);
  def_grad_get(nne_3d,ndofn,CONST_4(double) ST_tensor,ve_disp,F_mat);
  mat2array(F,CONST_2(double) F_mat,ndn,ndn);

  if(sup->multi_scale){
    cblas_daxpy(9,1.0,sup->F0,1,F,1);
  }

  free(N_x);
  free(N_y);
  free(N_z);
  dealoc4(ST_tensor,ndn,ndn,ndn);
  dealoc2(F_mat,ndn);

  *J = getJacobian(F,ptr_be->vol_elem_id,&err);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
	      ndn,ndn,ndn,1.0,F,ndn,F,ndn,0.0,C,ndn);
  err += inverse(C,3,C_I);
  get_material_stress(kappa,ptr_mat,C,C_I,*J,Sbar);

  return err;
}/* bnd_integration_help() */

static void disp_based_bnd_Kuu_at_ip(double *Kuu,
				     const int nn_ve,
				     const int ndofn,
				     const double *ST,
				     const double *normal,
				     const double *F,
				     const double *L,
				     const double *KK,
				     const double *LM,
				     const double ome,
				     const double jj,
				     const double wt)
{
  double *FSTab = aloc1(9);
  double *FSTwg = aloc1(9);
  double *STST = aloc1(9);

  for(int a=0; a<nn_ve; a++){
    for(int b=0; b<ndofn; b++){
      const double *pSTab = &ST[idx_4_gen(a,b,0,0,nn_ve,ndn,ndn,ndn)];

      for(int i=0; i<ndn; i++){
	for(int j=0; j<ndn; j++){
	  const int ij = idx_2(i,j);
	  FSTab[ij] = 0.0;
	  for(int k=0; k<ndn; k++){
	    const int ki = idx_2(k,i);
	    const int kj = idx_2(k,j);
	    FSTab[ij] += 0.5*(F[ki]*pSTab[kj] + pSTab[ki]*F[kj]);
	  }
	}
      }

      for(int w=0; w<nn_ve; w++){
  	for(int g=0; g<ndofn; g++){
	  const double *pSTwg = &ST[idx_4_gen(w,g,0,0,nn_ve,ndn,ndn,ndn)];

	  for(int i=0; i<ndn; i++){
	    for(int j=0; j<ndn; j++){
	      const int ij = idx_2(i,j);
	      FSTwg[ij] = 0.0;
	      STST[ij] = 0.0;
	      for(int k=0; k<ndn; k++){
		const int ki = idx_2(k,i);
		const int kj = idx_2(k,j);
		FSTwg[ij] += 0.5*(F[ki]*pSTwg[kj] + pSTwg[ki]*F[kj]);
		STST[ij] += 0.5*(pSTwg[ki]*pSTab[kj] + pSTab[ki]*pSTwg[kj]);
	      }
	    }
	  }

	  double term_I = 0.0;
	  double term_II = 0.0;
	  double term_III = 0.0;
	  double term_IV = 0.0;
	  for(int i=0; i<ndn; i++){
	    for(int j=0; j<ndn; j++){
	      for(int k=0; k<ndn; k++){
		const int ik = idx_2(i,k);
		for(int l=0; l<ndn; l++){
		  for(int m=0; m<ndn; m++){
		    const int kjlm = idx_4(k,j,l,m);
		    const int lm = idx_2(l,m);
		    term_I += LM[i]*pSTab[ik]*L[kjlm]*FSTwg[lm]*normal[j];
		    term_II += LM[i]*pSTwg[ik]*L[kjlm]*FSTab[lm]*normal[j];
		    term_III += LM[i]*F[ik]*L[kjlm]*STST[lm]*normal[j];

		    for(int r=0; r<ndn; r++){
		      for(int s=0; s<ndn; s++){
		    	const int kjlmrs = idx_6(k,j,l,m,r,s);
		    	const int rs = idx_2(r,s);
		    	term_IV += (LM[i]*normal[j]
		    		    *F[ik]*KK[kjlmrs]*FSTab[lm]*FSTwg[rs]);
		      }
		    }
		  }
		}
	      }
	    }
	  }

	  Kuu[idx_K(a,b,w,g,nn_ve,ndofn)] = 
	    (term_I + term_II + term_III + term_IV)*jj*wt;
  	} /* g */
      } /* w */
    } /* b */
  } /* a */

  free(FSTab);
  free(FSTwg);
  free(STST);

}/* disp_based_bnd_Kuu_at_ip() */

static void disp_based_bnd_Kul_Klu_at_ip(double *Kul,
					 double *Klu,
					 const int nn_ve,
					 const int ndofn,
					 const int n_lm_vecs,
					 const double *N_lm,
					 const double *ST,
					 const double *normal,
					 const double *F,
					 const double *Sbar,
					 const double *L,
					 const double ome,
					 const double jj,
					 const double wt)
{
  double *FSTab = aloc1(9);

  for(int a=0; a<nn_ve; a++){
    for(int b=0; b<ndofn; b++){
      /* pre-compute on nodes dofs */
      const double *pST = &ST[idx_4_gen(a,b,0,0,nn_ve,ndn,ndn,ndn)];

      for(int i=0; i<ndn; i++){
	for(int j=0; j<ndn; j++){
	  const int idx_ij = idx_2(i,j);
	  FSTab[idx_ij] = 0.0;
	  for(int k=0; k<ndn; k++){
	    const int ki = idx_2(k,i);
	    const int kj = idx_2(k,j);
	    FSTab[idx_ij] += 0.5*(F[ki]*pST[kj] + pST[ki]*F[kj]);
	  }
	}
      }

      for(int w=0; w<n_lm_vecs; w++){
	if(N_lm[w] == 0) continue;
	for(int g=0; g<ndn; g++){
	  double term_I = 0.0;
	  double term_II = 0.0;
	  for(int i=0; i<ndn; i++){
	    for(int j=0; j<ndn; j++){
	      for(int k=0; k<ndn; k++){
		const int ik = idx_2(i,k);
		const int kj = idx_2(k,j);
		term_I += N_lm[w]*del(i,g)*((1.-ome)*pST[ik]*Sbar[kj])*normal[j];
		for(int l=0; l<ndn; l++){
		  for(int m=0; m<ndn; m++){
		    const int kjlm = idx_4(k,j,l,m);
		    const int lm = idx_2(l,m);
		    term_II += N_lm[w]*del(i,g)*F[ik]*L[kjlm]*FSTab[lm]*normal[j];
		  }
		}
	      }
	    }
	  }

	  Kul[idx_K_gen(a,b,w,g,nn_ve,ndofn,
			n_lm_vecs,ndn)] = (term_I + term_II)*jj*wt;
	  Klu[idx_K_gen(w,g,a,b,n_lm_vecs,
			ndn,nn_ve,ndofn)] = (term_I + term_II)*jj*wt;
	}
      }
    }
  }

  free(FSTab);
}/* disp_based_bnd_Kul_Klu_at_ip() */

static int compute_K_00_e_at_ip_2(double *K_00_e,
				/* macro information */
				const int macro_nnode,
				const int macro_ndofn,
				const double macro_int_wt,
				const double *macro_shape_func,
				const double *macro_normal,
				const double layer_thickness,
				const double *gNoxN,
				/* micro information */
				const int nne,
				const double micro_volume0,
				const double jj,
				const double wt,
				const double *F,
				const double *ST,
				const damage *p_dam,
				const double *Sbar,
				const double *L)
{
  int err = 0;
  /* compute the terms from the microscale, then integrate at the
     macroscale. Use a bunch of nested loops to ensure correctness,
     optimize later! */
  const double scale = 1.0/(layer_thickness*micro_volume0);
  double val = 0.0;
  for(int k=0; k<ndn; k++){
    if(macro_normal[k] == 0.0) continue;
    for(int j=0; j<ndn; j++){
      val += (macro_normal[k]*(1.0-p_dam->w)*Sbar[idx_2(k,j)]
	      *macro_normal[j]*jj*wt*scale);
    }
  }

  double *term_I = PGFEM_calloc(ndn*ndn,sizeof(double));
  double *term_II = PGFEM_calloc(ndn*ndn,sizeof(double));
  for(int i=0; i<ndn; i++){
    for(int r=0; r<ndn; r++){
      /* get constant pointer to term computing */
      double *const T = term_II + idx_2(i,r);
      if(i==r) term_I[idx_2(i,r)] = val;
      for(int j=0; j<ndn; j++){
	const double Nj = macro_normal[j];
	if(Nj == 0.0) continue;
	for(int k=0; k<ndn; k++){
	  const int ik = idx_2(i,k);
	  for(int ll=0; ll<ndn; ll++){
	    const double Nl =  macro_normal[ll];
	    const int rl = idx_2(r,ll);
	    for(int m=0; m<ndn; m++){
	      const double Nm =  macro_normal[m];
	      const int rm = idx_2(r,m);
	      const int kjlm = idx_4(k,j,ll,m);

	      /* compute term for current index */
	      *T += (scale*0.5*F[ik]*L[kjlm]
		     *(F[rm]*Nl + F[rl]*Nm)*Nj*jj*wt);

	    }
	  }
	}
      }
    }
  }

  /* integrate at macroscale */
  for(int b=0; b<ndn; b++){
    for(int g=0; g<ndn; g++){
      const int bg = idx_2(b,g);

      for(int a=0; a<macro_nnode; a++){
	const double Na = macro_shape_func[a];
	for(int w=0; w<macro_nnode; w++){
	  const double Nw = macro_shape_func[w];
	  const int idx = idx_K(a,b,w,g,macro_nnode,macro_ndofn);
	  K_00_e[idx] += macro_int_wt*Na*Nw*(term_I[bg]+term_II[bg]);
	}
      }
    }
  }

  free(term_I);
  free(term_II);
  return err;
}

/** compute the micro term: d trac_i/ d 1u_p N_p|wg */
static int compute_K_01_micro_term(double *result,
				   const double *F,
				   const double *ST_wg,
				   const double *normal,
				   const double *Sbar,
				   const double *L,
				   const double damage,
				   const double volume,
				   const double wt,
				   const double jj)
{
  int err = 0;

  /* compute symmetric F' ST_wg (3 index) */
  /* compute ST_wg S (3 index) normal */
  double *FSTs = PGFEM_calloc(ndn*ndn,sizeof(double));
  for(int i=0; i<ndn; i++){
   for(int k=0; k<ndn; k++){
     const int ik = idx_2(i,k);
     const int ki = idx_2(k,i);
     for(int j=0; j<ndn; j++){
       const int kj = idx_2(k,j);
       const int ij = idx_2(i,j);
       result[i] += ST_wg[ik]*(1.0 - damage)*Sbar[kj]*normal[j];
       FSTs[ij] += 0.5*(F[ki]*ST_wg[kj] + ST_wg[ki]*F[kj]);
      }
    }
  }

  /* compute F (L:FSTs) normal */
  for(int j=0; j<ndn; j++){
    if(normal[j] == 0.0) continue;
    for(int i=0; i<ndn; i++){
      for(int k=0; k<ndn; k++){
	const int ik = idx_2(i,k);
	if(F[ik] == 0.0) continue;
	const double *L_kj = L + idx_4(k,j,0,0);
	for(int lm=0; lm<ndn*ndn; lm++){
	  result[i] += F[ik]*L_kj[lm]*FSTs[lm]*normal[j];
	}
      }
    }
  }

  /* integrate value and scale */
  cblas_dscal(ndn,jj*wt/volume,result,1);

  free(FSTs);
  return err;
}
				   
				   

static int compute_K_01_e_at_ip(double *K_01_e,
				/* macro information */
				const int macro_nnode,
				const int macro_ndofn,
				const double macro_int_wt,
				const double *macro_shape_func,
				const double *macro_normal,
				const double *gNoxN,
				/* micro information */
				const int nne,
				const double micro_volume0,
				const double jj,
				const double wt,
				const double *F,
				const double *ST,
				const damage *p_dam,
				const double *Sbar,
				const double *L)
{
  int err = 0;
  for(int w=0; w<nne; w++){
    for(int g=0; g<ndn; g++){
      const int wg = idx_4_gen(w,g,0,0,nne,ndn,ndn,ndn);
      double *dt_d1u_wg = PGFEM_calloc(ndn,sizeof(double));
      err += compute_K_01_micro_term(dt_d1u_wg,F,ST+wg,macro_normal,Sbar,L,
				     p_dam->w,micro_volume0,wt,jj);

      /* assemble to element matrix */
      for(int a=0; a<macro_nnode; a++){
	for(int b=0; b<macro_ndofn; b++){
	  const int abwg = idx_K_gen(a,b,w,g,macro_nnode,
				     macro_ndofn,nne,ndn);
	  K_01_e[abwg] = macro_int_wt*macro_shape_func[a]*dt_d1u_wg[b];
	}
      }
      free(dt_d1u_wg);
    }
  }

  return err;

}

static int compute_K_10_micro_term(double *result,
				   const double *F,
				   const double *ST_ab,
				   const double *normal,
				   const double *Sbar,
				   const double *L,
				   const double damage,
				   const double volume,
				   const double wt,
				   const double jj)
{
  int err = 0;
  /* compute symmetric F' ST_ab (3 index) */
  double *FSTs = PGFEM_calloc(ndn*ndn,sizeof(double));
  double *FN = PGFEM_calloc(ndn*ndn*ndn,sizeof(double));
  double *NST = PGFEM_calloc(ndn*ndn*ndn,sizeof(double));
  for(int i=0; i<ndn; i++){
    for(int k=0; k<ndn; k++){
      const int ik = idx_2(i,k);
      const int ki = idx_2(k,i);
      for(int j=0; j<ndn; j++){
	const int kj = idx_2(k,j);
	const int ij = idx_2(i,j);
	const int ijk = idx_3(i,j,k);
	FSTs[ij] += 0.5*(F[ki]*ST_ab[kj] + ST_ab[ki]*F[kj]);
	FN[ijk] += 0.5*(F[ik]*normal[j] + F[ij]*normal[k]);
	NST[ijk] += 0.5*(normal[i]*ST_ab[ik] + ST_ab[ij]*normal[k]);
      }
    }
  }

  /* big loop */
  double int_and_scale = jj*wt/volume;
  for(int p=0; p<ndn; p++){
    result[p] = 0.0;
    const double *pFN = FN + idx_3(p,0,0);
    const double *pNST = NST + idx_3(p,0,0);
    for(int ij=0; ij<ndn*ndn; ij++){
      result[p] += (1.0-damage)*Sbar[ij]*pNST[ij];
      const double *Lij = L + ij*ndn*ndn;
      for(int lm=0; lm<ndn*ndn; lm++){
	result[p] += FSTs[ij]*Lij[lm]*pFN[lm];
      }
    }
    result[p] *= int_and_scale;
  }

  free(FSTs);
  free(FN);
  free(NST);
  return err;
}

static int compute_K_10_e_at_ip(double *K_10_e,
				/* macro information */
				const int macro_nnode,
				const int macro_ndofn,
				const double macro_int_wt,
				const double *macro_shape_func,
				const double *macro_normal,
				const double *gNoxN,
				/* micro information */
				const int nne,
				const double micro_volume0,
				const double jj,
				const double wt,
				const double *F,
				const double *ST,
				const damage *p_dam,
				const double *Sbar,
				const double *L)
{
  int err = 0;
  for(int a=0; a<nne; a++){
    for(int b=0; b<ndn; b++){
      const int ab = idx_4_gen(a,b,0,0,nne,ndn,ndn,ndn);
      double *micro_term_ab = PGFEM_calloc(ndn,sizeof(double));
      err += compute_K_10_micro_term(micro_term_ab,F,ST+ab,macro_normal,
				     Sbar,L,p_dam->w,micro_volume0,wt,jj);

      for(int w=0; w<macro_nnode; w++){
	for(int g=0; g<macro_ndofn; g++){
	  const int abwg = idx_K_gen(a,b,w,g,macro_nnode,
				     macro_ndofn,nne,ndn);
	  K_10_e[abwg] = macro_shape_func[w]*micro_term_ab[g];
	}
      }
      free(micro_term_ab);
    }
  }

  return err;
}
