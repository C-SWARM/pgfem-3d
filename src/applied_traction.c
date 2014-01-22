/* HEADER */
/* This file contains functions to compute the applied tractions to
   model features. Currently, only surface tractions are supported */

#include "applied_traction.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef INDEX_MACROS_H
#include "index_macros.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef INTEGRATE_SURFACE_H
#include "integrate_surface.h"
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

#ifndef GET_NDOF_ON_ELEM_H
#include "get_ndof_on_elem.h"
#endif

#ifndef GET_DOF_IDS_ON_ELEM_H
#include "get_dof_ids_on_elem.h"
#endif

int read_applied_surface_tractions_fname(char *fname,
					 int *n_feats,
					 int **feat_type,
					 int **feat_id,
					 double **loads)
{
  int err = 0;
  FILE *in = PGFEM_fopen(fname,"r");
  err = read_applied_surface_tractions(in,n_feats,feat_type,feat_id,loads);
  PGFEM_fclose(in);
  return err;
}

int read_applied_surface_tractions(FILE *in,
				   int *n_feats,
				   int **feat_type,
				   int **feat_id,
				   double **loads)
{
  int err = 0;
  (*n_feats) = 0;
  (*feat_type) = NULL;
  (*feat_id) = NULL;
  (*loads) = NULL;

  /* read number of features and allocate */
  fscanf(in,"%d",n_feats);
  (*feat_type) = PGFEM_calloc(*n_feats,sizeof(int));
  (*feat_id) = PGFEM_calloc(*n_feats,sizeof(int));
  (*loads) = PGFEM_calloc(3*(*n_feats),sizeof(double));

  int *ft = &(*feat_type)[0];
  int *fi = &(*feat_id)[0];
  double *ld = &(*loads)[0];
  for(int i=0; i<*n_feats; i++){
    int idx = idx_2_gen(i,0,*n_feats,3);
    fscanf(in,"%d %d %lf %lf %lf",
	   ft+i,fi+i,ld+idx,ld+idx+1,ld+idx+2);
  }

  return err;
}

int generate_applied_surface_traction_list(const int ne,
					   const ELEMENT *elem,
					   const int n_feats,
					   const int *feat_type,
					   const int *feat_id,
					   int *n_sur_trac_elem,
					   SUR_TRAC_ELEM **sur_trac_elem)
{
  int err = 0;
  (*n_sur_trac_elem) = 0;
  (*sur_trac_elem) = NULL;

  int *elem_has_trac = PGFEM_calloc(ne,sizeof(int));
  int *elem_n_faces = PGFEM_calloc(ne,sizeof(int));
  int *n_faces = PGFEM_calloc(ne,sizeof(int));

  /* determine how many elements have tractions, which ones, and how
     many faces on the element have an applied traction */
  for(int i=0; i<ne; i++){
    const ELEMENT *el = &elem[i];

    /* get number of element faces */
    switch(el->toe){
      /* tetrahedron */
    case 4: case 5: case 10: n_faces[i] = 4; break;
      /* hexas */
    case 8: n_faces[i] = 6; break;
    default:
      PGFEM_printerr("WARNING: Unsupported elemet type,"
		     " no tractions will be applied! %s:%s:%d\n",
		     __func__,__FILE__,__LINE__);
      break;
    }

    for(int k=0; k<n_feats; k++){
      for(int j=0; j<n_faces[i]; j++){
	if(el->bnd_type[j] == feat_type[k]
	   && el->bnd_id[j] == feat_id[k]){
	  elem_has_trac[i] = 1;
	  elem_n_faces[i] ++;
	}
      }/* faces */
    }/* features */

    /* running total of how many elements */
    if(elem_has_trac[i]){
      (*n_sur_trac_elem)++;
    }
  }

  /* get short list of elements */
  int *elem_ids = NULL;
  if(*n_sur_trac_elem > 0){
    elem_ids = PGFEM_calloc(*n_sur_trac_elem,sizeof(int));
    {
      int idx = 0;
      for(int i=0; i<ne; i++){
	if(elem_has_trac[i]){
	  elem_ids[idx] = i;
	  idx++;
	}      
      }
    }

    /* allocate and populate */
    (*sur_trac_elem) = PGFEM_calloc(*n_sur_trac_elem,sizeof(SUR_TRAC_ELEM));
    for(int i=0; i<(*n_sur_trac_elem); i++){
      int idx = elem_ids[i];
      const ELEMENT *el = &elem[idx];
      SUR_TRAC_ELEM *ste = &(*sur_trac_elem)[i];

      ste->elem_id = idx;
      ste->n_faces = elem_n_faces[idx];
      ste->faces = PGFEM_calloc(ste->n_faces,sizeof(int));
      ste->feat_num = PGFEM_calloc(ste->n_faces,sizeof(int));

      /* loop through faces on element */
      int fidx = 0;
      for(int k=0; k<n_feats; k++){
	for(int j=0; j<n_faces[i]; j++){
	  if(el->bnd_type[j] == feat_type[k]
	     && el->bnd_id[j] == feat_id[k]){
	    ste->faces[fidx] = j;
	    ste->feat_num[fidx] =  k;
	    fidx++;
	  }
	}/* faces */
      }/* features */
    }/* for each loaded element */
  }
  /* free memory */
  free(elem_has_trac);
  free(elem_n_faces);
  free(n_faces);

  /* exit function */
  return err;
}/* generate_applied_surface_traction_list() */

int destroy_applied_surface_traction_list(const int n_sur_trac_elem,
					  SUR_TRAC_ELEM *sur_trac_elem)
{
  int err = 0;
  for(int i=0; i< n_sur_trac_elem; i++){
    free(sur_trac_elem[i].faces);
    free(sur_trac_elem[i].feat_num);
  }
  free(sur_trac_elem);
  return err;
}/* destroy_applied_surface_traction_list() */

/* Distributed dead load is computed in reference configuration */
int compute_applied_traction_res(const int ndofn,
				 const NODE *nodes,
				 const ELEMENT *elem,
				 const int n_ste,
				 const SUR_TRAC_ELEM *ste,
				 const int n_feats,
				 const double *loads,
				 double *res)
{
  int err = 0;
  const int int_order = 1;
  for(int i=0; i<n_ste; i++){
    const SUR_TRAC_ELEM *pste = &ste[i];

    const ELEMENT *el = &elem[pste->elem_id];
    const long *nod_3D = el->nod;
    const int nne_3D = el->toe;

    double *x = PGFEM_calloc(nne_3D,sizeof(double));
    double *y = PGFEM_calloc(nne_3D,sizeof(double));
    double *z = PGFEM_calloc(nne_3D,sizeof(double));
    nodecoord_total(nne_3D,nod_3D,nodes,x,y,z);

    const int ndofe = get_ndof_on_elem_nodes(nne_3D,nod_3D,nodes);
    long *cn = PGFEM_calloc(ndofe,sizeof(long));
    double *res_el = PGFEM_calloc(ndofe,sizeof(double));
    get_dof_ids_on_elem_nodes(0,nne_3D,ndofn,nod_3D,nodes,cn);

    double *N_3D = PGFEM_calloc(nne_3D,sizeof(double));

    for(int j=0; j<pste->n_faces; j++){
      const double *trac = &loads[idx_2_gen(pste->feat_num[j],
					    0,n_feats,3)];
      int n_ip = 0;
      double *ksi_3D = NULL;
      double *eta_3D = NULL;
      double *zet_3D = NULL;
      double *ksi_2D = NULL;
      double *eta_2D = NULL;
      double *wt_2D = NULL;
      int *nod_2D = NULL;
      int nne_2D = 0;

      err += integrate_surface(nne_3D,pste->faces[j],int_order,
			       &n_ip,&ksi_3D,&eta_3D,&zet_3D,
			       &ksi_2D,&eta_2D,&wt_2D,
			       &nne_2D,&nod_2D);

      for(int ip=0; ip<n_ip; ip++){
	double jj = compute_surface_jacobian(nne_2D,nod_2D,x,y,z,
					     ksi_2D[ip],eta_2D[ip]);
	shape_func(ksi_3D[ip],eta_3D[ip],zet_3D[ip],nne_3D,N_3D);

	for(int k=0; k<nne_3D; k++){
	  for(int m=0; m<3; m++){
	    int idx = idx_2_gen(k,m,nne_3D,ndofn);
	    res_el[idx] += N_3D[k]*trac[m]*jj*wt_2D[ip];
	  }/* dims */
	}/* nodes */

      }/* for each integration point */

      free(ksi_3D);
      free(eta_3D);
      free(zet_3D);
      free(ksi_2D);
      free(eta_2D);
      free(wt_2D);
      free(nod_2D);
    }/* for each face w/ app. trac. */

    /* assembly */
    for(int j=0; j<ndofe; j++){
      if(cn[j] - 1 >= 0){
	res[cn[j] - 1] += res_el[j];
      }
    }

    free(x);
    free(y);
    free(z);
    free(N_3D);
    free(cn);
    free(res_el);
  }/* for elements */

  return err;
}/* compute_applied_traction_res() */
