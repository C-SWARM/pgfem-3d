/*******************************************
 * Program FEM3d ver. 2.0                  *
 * FEM - 3D analysis                       *
 * Karel Matous & Jaroslav Kruis           *
 *******************************************/

/*****************/
/* November 2000 */
/*****************/

#include "load.h"

#ifndef ENUMERATIONS_H
#include "enumerations.h"
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

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

/* #ifndef INCL_H */
/* #include "incl.h" */
/* #endif */

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

#ifndef STABILIZED_H
#include "stabilized.h"
#endif

#ifndef STIFFMATEL_FD_H
#include "stiffmatel_fd.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef MINI_ELEMENT_H
#include "MINI_element.h"
#endif

#ifndef MINI_3F_ELEMENT_H
#include "MINI_3f_element.h"
#endif

#ifndef DISP_BASED_ELEM_H
#include "displacement_based_element.h"
#endif

long* compute_times_load (FILE *in1,
			  long nt,
			  long nlod_tim)
     /*
       
     */
{
  long i,j,*tim_load;
  long *help;
  
  if(nt == 0){
    tim_load = aloc1l (1);
  } else {
    tim_load = aloc1l (nt);
  }
  
  if(nlod_tim == 0){
    help = aloc1l (1);
  } else {
    help = aloc1l (nlod_tim);
  }
  
  for (i=0;i<nlod_tim;i++)
    fscanf (in1,"%ld",&help[i]);
  
  for (i=0;i<nlod_tim;i++){
    for (j=0;j<nt;j++){
      if (j == help[i]){
	tim_load[j] = 1;
	continue;
      }
    }
  }
  
  dealoc1l (help);
  
  return (tim_load);
}

void load_vec_node (double *f,
		    long nln,
		    long ndofn,
		    ZATNODE *znode,
		    NODE *node)
	  /*
	    
	  */
{
  long i,j,ii;
  
  for (i=0;i<nln;i++){
    for (j=0;j<ndofn;j++){
      ii = node[znode[i].nod].id[j]-1;
      if (ii < 0)  continue;
      f[ii] += znode[i].load[j];
    }
  }
}

int load_vec_node_defl (double *f,
			long ne,
			long ndofn,
			ELEMENT *elem,
			BOUNDING_ELEMENT *b_elems,
			NODE *node,
			HOMMAT *hommat,
			MATGEOM matgeom,
			SUPP sup,
			long npres,
			double nor_min,
			SIG *sig,
			EPS *eps,
			double dt,
			CRPL *crpl,
			double stab,
			double *r,
			const PGFem3D_opt *opts)
     /*
       
     */
{
  int err = 0;
  double *fe = NULL;

  /* 3D */
  const int ndn = 3;

  for (int i=0;i<sup->nde;i++){
    
    /* Number of element nodes */
    int nne = elem[sup->lepd[i]].toe;
    int nne_t = nne + elem[sup->lepd[i]].n_bub;
    /* Nodes on element */
    long *nod = aloc1l (nne);
    elemnodes (sup->lepd[i],nne,nod,elem);
    
    /* Element Dof */
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node);

    /* Allocation */
    long *cn = aloc1l (ndofe);
    double *lk= aloc1 (ndofe*ndofe);

    double *x,*y,*z;
    if (opts->analysis_type == MINI
	|| opts->analysis_type == MINI_3F){ /* P1+B/P1 */
      x = aloc1 (nne_t);
      y = aloc1 (nne_t);
      z = aloc1 (nne_t);
    } else {
      x = aloc1 (nne);
      y = aloc1 (nne);
      z = aloc1 (nne);
    }

    double *floc = aloc1 (ndofe);
    double *rloc = aloc1 (ndofe); 
    
    /* Coordinates of nodes */
    switch(opts->analysis_type){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);    
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);    
      break;
    }
    if (opts->analysis_type == MINI
	|| opts->analysis_type == MINI_3F){ /* P1+B/P1 */
      element_center(nne,x,y,z);
    }
   
    nulld (lk,ndofe*ndofe);
    double *r_e, *r_r;

    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    r_e = aloc1 (ndofe); 
    r_r = aloc1 (ndofe); /* for TOTAL LAGRANGIAN */
    switch(opts->analysis_type){
    case FS_CRPL:
    case FINITE_STRAIN:
      stiffmatel_fd (sup->lepd[i],ndofn,nne,nod,x,y,z,elem,matgeom,
		     hommat,node,sig,eps,r_e,npres,nor_min,lk,dt,
		     crpl,0,0.0,fe,opts->analysis_type);
      break;
    case STABILIZED:
      stiffmatel_st (sup->lepd[i],ndofn,nne,x,y,z,elem,
		     hommat,nod,node,sig,eps,r_e,npres,
		     nor_min,lk,dt,stab,0,0.0,fe);
      break;
    case MINI:
      MINI_stiffmat_el(lk,sup->lepd[i],ndofn,nne,x,y,z,elem,
		       hommat,nod,node,eps,sig,r_e);
      break;
    case MINI_3F:
      MINI_3f_stiffmat_el(lk,sup->lepd[i],ndofn,nne,x,y,z,elem,
			  hommat,nod,node,eps,sig,r_e);
      break;
    case DISP:
	/* Total a-vector */
	def_elem (cn,ndofe,r,elem,node,r_r,sup,1);
	err = DISP_stiffmat_el(lk,sup->lepd[i],ndofn,nne,x,y,z,elem,
			       hommat,nod,node,eps,sig,sup,r_r);
      break;
    default:
      stiffmatel (sup->lepd[i],x,y,z,nne,ndofn,elem,hommat,node,lk,opts);
      break;
    } /* switch analysis */
    dealoc1 (r_e);
    dealoc1 (r_r);
  
    
    /* get the disp increment from BC */
    {
      int k = 0;
      int jj = 0;
      for (int ii=0;ii<nne;ii++){
	for (jj=0;jj<node[nod[ii]].ndofn;jj++){
	  if (jj  < ndn && cn[k+jj] <= -1)
	    rloc[k+jj] = sup->defl_d[abs(cn[k+jj])-1];
	  else
	    rloc[k+jj] = 0.0;
	}
	k += jj;
      }
    }
    
    /* Matrix vector multiplication */
    mv (lk,rloc,floc,ndofe,ndofe);
    
    /* Localization */
    {
      int k = 0;
      int jj = 0;
      int II = 0;
      for (jj=0;jj<nne;jj++){
	for (int ii=0;ii<node[nod[jj]].ndofn;ii++){
	  II = node[nod[jj]].id[ii]-1;
	  if (II < 0)  continue;
	  f[II] += floc[k+ii];
	}/*end ii*/
	k += node[nod[jj]].ndofn;
      }/*end jj*/
    }

    /*  dealocation  */
    dealoc1l (cn);
    dealoc1l (nod);
    dealoc1 (lk);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1(floc);
    dealoc1(rloc);
    
    if(err != 0) return err;
  }/* end i (each volume element) */

  /* ADDED 1/7/2013 MM */
  for(int i=0; i<sup->nd_be; i++){
    /* get pointers and constant vaules for the bounding element we
       are working on. */
    const int be_id = sup->lbepd[i];
    const BOUNDING_ELEMENT *ptr_be = &b_elems[be_id];
    const int ve_id = ptr_be->vol_elem_id;
    const ELEMENT *ptr_ve = &elem[ve_id];
    const long *ve_nod = ptr_ve->nod;
    const int nne_ve = ptr_ve->toe;
    const int nne_ve_t = nne_ve + ptr_ve->n_bub;

    /* get ndofs on element */
    int ndof_ve = get_ndof_on_bnd_elem(node,ptr_be,elem);

    /* get coordinates of volume element nodes */
    double *x = aloc1(nne_ve);
    double *y = aloc1(nne_ve);
    double *z = aloc1(nne_ve);

    switch(opts->analysis_type){
    case DISP:
      nodecoord_total (nne_ve,ve_nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne_ve,ve_nod,node,x,y,z);    
      break;
    }

    /* allocate space for stiffness and localization. */
    double *lk = aloc1(ndof_ve*ndof_ve);
    double *floc = aloc1(ndof_ve);
    double *rloc = aloc1(ndof_ve);
    long *cn_ve = aloc1l(ndof_ve);
    get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn_ve);

    /* get displacements */
    double *ve_disp = aloc1(ndof_ve);
    if(opts->analysis_type == DISP){ /* TOTAL LAGRANGIAN formulation */
      def_elem(cn_ve,ndof_ve,r,NULL,NULL,ve_disp,sup,1);
    }

    /* compute element stiffness */
    if(opts->analysis_type == DISP){
      err += DISP_stiffmat_bnd_el(lk,be_id,ndofn,ndof_ve,x,y,z,b_elems,
				  elem,hommat,node,eps,sig,sup,ve_disp);
    } else {
      /* not implemented, do nothing */
    }

    /* get the local disp increment from BC */
    {
      int k = 0;
      int jj = 0;
      for (int ii=0;ii<nne_ve;ii++){
	for (jj=0;jj<node[ve_nod[ii]].ndofn;jj++){
	  if (jj  < ndn && cn_ve[k+jj] <= -1)
	    rloc[k+jj] = sup->defl_d[abs(cn_ve[k+jj])-1];
	  else
	    rloc[k+jj] = 0.0;
	}
	k += jj;
      }
    }

    /* matvec mult */
    mv(lk,rloc,floc,ndof_ve,ndof_ve);

    /* Localization */
    {
      int j = 0;
      int II = 0;
      for (int jj=0;jj<nne_ve;jj++){
	for (int ii=0;ii<node[ve_nod[jj]].ndofn;ii++){
	  II = node[ve_nod[jj]].id[ii]-1;
	  if (II < 0)  continue;
	  f[II] += floc[j+ii];
	}/*end ii*/
	j += node[ve_nod[jj]].ndofn;
      }/*end jj*/
    }

    free(x);
    free(y);
    free(z);

    free(lk);
    free(floc);
    free(rloc);
    free(cn_ve);

    free(ve_disp);

    if(err != 0) return err;
  }/* end for each bounding element (i) in the list */
  return err;
}

void load_vec_elem_sur (double *f,
			long nle_s,
			long ndofn,
			ELEMENT *elem,
			ZATELEM *zele_s)
     /*
       NOT IMPLIMENTED
     */
{
  
  
}
