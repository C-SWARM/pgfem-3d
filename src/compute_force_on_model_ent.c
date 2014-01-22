#include "compute_force_on_model_ent.h"
#include "mkl_cblas.h"

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef DISP_BASED_ELEM_H
#include "displacement_based_element.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef INCL_H
#include "incl.h"
#endif

#ifndef GET_NDOF_ON_ELEM_H
#include "get_ndof_on_elem.h"
#endif

#ifndef GET_DOF_IDS_ON_ELEM_H
#include "get_dof_ids_on_elem.h"
#endif

int read_model_entity_list(const char *filename,
			   MODEL_ENTITY **me,
			   const int nnodes,
			   const int nelems)
{
  int err = 0;
  FILE *in;
  if((in = fopen(filename,"r")) != NULL){
   err += read_model_entity_listf(in,me,nnodes,nelems);
  } else {
    PGFEM_printerr("ERROR opening %s!\n",filename);
    err = 1;
  }
  return err;
} /* read_model_entity_list() */

int read_model_entity_listf(FILE *in,
			    MODEL_ENTITY **me,
			    const int nnodes,
			    const int nelems)
{
  int err = 0;
  int n_entities;
  fscanf(in,"%d",&n_entities);
  err += create_model_entity(n_entities,3,nnodes,nelems,me);
  if(!err){/* successfully created, read types and ids */
    int i=0;
    while(i<n_entities && !feof(in)){
      fscanf(in,"%d %d",(*me)->ent_type+i,(*me)->ent_id+i);
      i++;
    }
  }
  return err;
}/*  */

int create_model_entity(const int n_entities,
			const int n_dim,
			const int nnode,
			const int nelem,
			MODEL_ENTITY **me)
{
  int err = 0;
  (*me) = PGFEM_calloc(1,sizeof(MODEL_ENTITY));
  (*me)->n_entities = n_entities;
  (*me)->n_dim = n_dim;
  (*me)->ent_type = PGFEM_calloc(n_entities,sizeof(int));
  (*me)->ent_id = PGFEM_calloc(n_entities,sizeof(int));
  (*me)->node_ent = PGFEM_calloc(nnode,sizeof(int));
  (*me)->elem_related = PGFEM_calloc(nelem,sizeof(int));
  if((*me)->elem_related == NULL){/* out of memory */
    err = 1;
  }
  return err;
} /* create_model_entity() */

int destroy_model_entity(MODEL_ENTITY *me)
{
  if(me != NULL){
    free(me->ent_type);
    free(me->ent_id);
    free(me->node_ent);
    free(me->elem_related);
  }
  free(me);
  return 0;
}

int model_entity_mark_nodes_elems(const int nnodes,
				  const NODE *nodes,
				  const int nelems,
				  const ELEMENT *elems,
				  MODEL_ENTITY *me)
{
  int err = 0;
  /* mark nodes first */
  for(int i=0; i<nnodes; i++){
    int T = nodes[i].model_type;
    int I = nodes[i].model_id;
    int j = 0;
    for(j=0; j<me->n_entities; j++){
      if(T == me->ent_type[j] && I == me->ent_id[j]){
	break;
      }
    }
    me->node_ent[i] = j;
    /* note j == n_entities --> skip this index */
  }

  /* mark elements */
  for(int i=0; i<nelems; i++){
    const int nne = elems[i].toe;
    const long *nod = elems[i].nod;
    int related = 0;
    for(int j=0; j<nne; j++){
      if(me->node_ent[nod[j]] < me->n_entities) related = 1;
    }
    me->elem_related[i] = related;
  }

  return err;
}

/* this function must be called before element update functions */
int compute_force_increment_on_model_entity(const int nelems,
					    const int ndofn,
					    const double *total_disp_n,
					    const double *disp_k,
					    const double *del_disp_k,
					    const MODEL_ENTITY *me,
					    const ELEMENT *elems,
					    const NODE *nodes,
					    const SUPP sup,
					    const HOMMAT *hommat,
					    const EPS *eps,
					    const SIG *sig,
					    const MPI_Comm mpi_comm,
					    const int analysis,
					    double *forces)
{
  int err;
  /* note forces must be of length (n_entities+1)*me->n_dim */
  for(int i=0; i<nelems; i++){
    if(!me->elem_related[i]) continue;
    const int nne = elems[i].toe;
    const long *nod = elems[i].nod;
    const int ndofe = get_ndof_on_elem_nodes(nne,nod,nodes);
    double *K = PGFEM_calloc(ndofe*ndofe,sizeof(double));
    double *elem_f = PGFEM_calloc(ndofe,sizeof(double));


    long *dof_ids = PGFEM_calloc(ndofe,sizeof(long));
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,nodes,dof_ids);

    double *x = PGFEM_calloc(nne,sizeof(double));
    double *y = PGFEM_calloc(nne,sizeof(double));
    double *z = PGFEM_calloc(nne,sizeof(double));

    switch(analysis){
    case DISP:
      nodecoord_total(nne,nod,nodes,x,y,z);
      break;
    default:
      nodecoord_updated(nne,nod,nodes,x,y,z);
      break;
    }

    /* displacement vectors u_n, delta u, u_n+1 */
    double *n_disp = PGFEM_calloc(ndofe,sizeof(double));
    double *km1_disp = PGFEM_calloc(ndofe,sizeof(double));
    double *k_disp = PGFEM_calloc(ndofe,sizeof(double));
    double *del_disp = PGFEM_calloc(ndofe,sizeof(double));

    /* get displacements at N+1(k-1), and N+1(k) */
    if(analysis == DISP){ /* total Lagrangian */
      def_elem(dof_ids,ndofe,total_disp_n,NULL,NULL,n_disp,sup,1);
      def_elem(dof_ids,ndofe,disp_k,NULL,NULL,k_disp,sup,0);
      def_elem(dof_ids,ndofe,del_disp_k,NULL,NULL,del_disp,sup,0);
    } else { /* updated Lagrangian */
      nulld(n_disp,ndofe); /* u_n = 0 */
      def_elem(dof_ids,ndofe,disp_k,NULL,NULL,k_disp,sup,0);
      def_elem(dof_ids,ndofe,del_disp_k,NULL,NULL,del_disp,sup,0);
    }
    for(int j=0; j<ndofe; j++){
      k_disp[j] += n_disp[j];
      km1_disp[j] = k_disp[j] - del_disp[j];
    }

    /* compute tangent at (k-1) */
    switch(analysis){
    case DISP:
      err += DISP_stiffmat_el(K,i,ndofn,nne,x,y,z,elems,hommat,
			      nod,nodes,eps,sig,sup,km1_disp);
      break;
    default:
      /* do nothing for now */
      break;
    }

    /* compute incremental force on element */
    cblas_dgemv(CblasRowMajor,CblasNoTrans,ndofe,ndofe,1.0,
		K,ndofe,del_disp,1,0.0,elem_f,1);

    /* assemeble the force */
    for(int j=0; j<nne; j++){
      const int ent_id = me->node_ent[nod[j]];
      if(ent_id >= me->n_entities) continue; /* skip non-marked node */
      for(int k=0; k<me->n_dim; k++){
	forces[ent_id*me->n_dim + k] += elem_f[j*ndofn + k];
      }
    }

    /* deallocate */
    free(K);
    free(elem_f);
    free(x);
    free(y);
    free(z);
    free(n_disp);
    free(km1_disp);
    free(del_disp);
    free(k_disp);
  }

  /* the forces on each domain must be accumulated only after
     convergence is reached */

  return err;
}
