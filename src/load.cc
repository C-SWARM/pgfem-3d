/* HEADER */


#include "PGFem3D_data_structure.h"
#include "load.h"
#include "enumerations.h"
#include "incl.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "elem3d.h"
#include "allocation.h"
#include "matice.h"
#include "stabilized.h"
#include "stiffmatel_fd.h"
#include "utils.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "three_field_element.h"
#include "stiffmat_fd.h"
#include "dynamics.h"
#include "constitutive_model.h"
#include "energy_equation.h"

long* compute_times_load (FILE *in1,
			  const long nt,
			  const long nlod_tim)
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
		    const long nln,
		    const long ndofn,
		    const ZATNODE *znode,
		    const NODE *node,
		    const int mp_id)
{
  long i,j,ii;
  
  for (i=0;i<nln;i++){
    for (j=0;j<ndofn;j++){
      ii = node[znode[i].nod].id_map[mp_id].id[j]-1;
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
                        double *r_n,
                        const PGFem3D_opt *opts,
                        double alpha,
                        const int mp_id)
{
  int err = 0;
  double *fe = NULL;

  /* 3D */
  const int ndn = 3;

  for (int i=0;i<sup->nde;i++){
    
    const int mat = elem[sup->lepd[i]].mat[2];
    double rho = hommat[mat].density;
    long include_inertia = 1;
    int nsd = 3;
  
    if(fabs(rho)<MIN_DENSITY)
    {
      include_inertia = 0;
    }    
    
    /* Number of element nodes */
    int nne = elem[sup->lepd[i]].toe;
    int nne_t = nne + elem[sup->lepd[i]].n_bub;
    /* Nodes on element */
    long *nod = elem[sup->lepd[i]].nod;
    
    /* Element Dof */
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);
    long *cn = aloc1l (ndofe);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* element tangent */
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
    double *r_e = aloc1 (ndofe); 

    /* Coordinates of nodes */
    if(sup->multi_scale){
      nodecoord_total (nne,nod,node,x,y,z);
      def_elem (cn,ndofe,r,elem,node,r_e,sup,1);
    } else {
      switch(opts->analysis_type){
      case DISP: case TF: /* total Lagrangian formulations */
        nodecoord_total(nne,nod,node,x,y,z);
        def_elem(cn,ndofe,r,elem,node,r_e,sup,1);
        break;

      case CM:
        {
          switch(opts->cm) {
            case UPDATED_LAGRANGIAN:
              nodecoord_updated(nne,nod,node,x,y,z);
              break;
            case TOTAL_LAGRANGIAN: /* total Lagrangian */
              nodecoord_total(nne,nod,node,x,y,z);
              def_elem(cn,ndofe,r,elem,node,r_e,sup,1);
              break;
            case MIXED_ANALYSIS_MODE: /* total Lagrangian */
              nodecoord_total(nne,nod,node,x,y,z);
              def_elem(cn,ndofe,r,elem,node,r_e,sup,1);
              break;
            default: assert(0 && "should never reach this case"); break;
          }
        }
        break;

      default: /* updated Lagrangian */
        nodecoord_updated(nne,nod,node,x,y,z);
        break;
      }
    }

    if (opts->analysis_type == MINI
	|| opts->analysis_type == MINI_3F){ /* P1+B/P1 */
      element_center(nne,x,y,z);
    }

    int nVol = N_VOL_TREE_FIELD;
    long FNR = 0;
    double lm = 0.0;
    err += el_compute_stiffmat(sup->lepd[i],lk,ndofn,nne,npres,nVol,nsd,
                               elem,node,hommat,matgeom,sig,eps,sup,
                               dt,nor_min,stab,crpl,FNR,lm,
                               x,y,z,fe,nod,r_n,r_e,
                               alpha,include_inertia,opts->analysis_type, opts->cm);

    /* get the disp increment from BC */
    {
      int k = 0;
      int jj = 0;
      for (int ii=0;ii<nne;ii++){
	for (jj=0;jj<ndofn;jj++){
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
	for (int ii=0;ii<ndofn;ii++){
	  II = node[nod[jj]].id_map[mp_id].id[ii]-1;
	  if (II < 0)  continue;
	  f[II] += floc[k+ii];
	}/*end ii*/
	k += ndofn;
      }/*end jj*/
    }

    /*  dealocation  */
    dealoc1l (cn);
    dealoc1 (lk);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1(floc);
    dealoc1(rloc);
    dealoc1 (r_e);

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

    /* get ndofs on element */
    int ndof_ve = get_ndof_on_bnd_elem(node,ptr_be,elem,ndofn);

    /* get coordinates of volume element nodes */
    double *x = aloc1(nne_ve);
    double *y = aloc1(nne_ve);
    double *z = aloc1(nne_ve);

    switch(opts->analysis_type){
    case DISP: case TF:
      nodecoord_total (nne_ve,ve_nod,node,x,y,z);
      break;
	  case CM:
	  {
      switch(opts->cm)
      {
        case UPDATED_LAGRANGIAN:
          nodecoord_updated(nne_ve,ve_nod,node,x,y,z);
          break;        
        case TOTAL_LAGRANGIAN:
          nodecoord_total(nne_ve,ve_nod,node,x,y,z);
          break;
        case MIXED_ANALYSIS_MODE:
        default:
          nodecoord_total(nne_ve,ve_nod,node,x,y,z);
          break;
      }
      break;
    }        
    default:
      nodecoord_updated (nne_ve,ve_nod,node,x,y,z);    
      break;
    }

    /* allocate space for stiffness and localization. */
    double *lk = aloc1(ndof_ve*ndof_ve);
    double *floc = aloc1(ndof_ve);
    double *rloc = aloc1(ndof_ve);
    long *cn_ve = aloc1l(ndof_ve);
    get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn_ve,mp_id);
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
	for (jj=0;jj<ndofn;jj++){
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
	for (int ii=0;ii<ndofn;ii++){
	  II = node[ve_nod[jj]].id_map[mp_id].id[ii]-1;
	  if (II < 0)  continue;
	  f[II] += floc[j+ii];
	}/*end ii*/
	j += ndofn;
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


/// Compute load vector for prescribed BCs(Dirichlet)
/// This compute load, f, as below:
/// [Kii Kio]<ui>   <bi>
/// [Koi Koo]<uo> = <bo>
/// [Kii][ui] = <bi> - [Kio]<uo>
/// where f = [Kio]<uo>, uo is Drichlet BCs
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] sol solution scheme object
/// \param[in] load object for loading
/// \param[in] dt time step
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return non-zero on internal erro
int momentum_equation_load4pBCs(GRID *grid,
                                MATERIAL_PROPERTY *mat,
                                FIELD_VARIABLES *fv,
                                SOLVER_OPTIONS *sol,
                                LOADING_STEPS *load,
                                double dt,
                                CRPL *crpl,
                                const PGFem3D_opt *opts,
                                MULTIPHYSICS *mp,
                                int mp_id,
                                int myrank)
{
  int err = 0;
  double *be = NULL;
  int intg_order = 0;
  
  ELEMENT *elem = grid->element;
  NODE *node = grid->node;
  SUPP sup = load->sups[mp_id];
  
  
  // multi-scale simulation, displacement base and three-fied-mixed (elastic only)
  // are total Lagrangian based
  int total_Lagrangian = 0;
  switch(opts->analysis_type)
  {
    case DISP: // intented to flow
    case TF:
      total_Lagrangian = 1;
      break;
    case CM:
      if(opts->cm != UPDATED_LAGRANGIAN)
        total_Lagrangian = 1;
      
      break;
  }
  
  if(sup->multi_scale)
    total_Lagrangian = 1;
  
  for (int i=0;i<sup->nde;i++)
  {
    int eid = sup->lepd[i];
    // set FEMLIB
    FEMLIB fe;
    if (opts->analysis_type == MINI || opts->analysis_type == MINI_3F)
      FEMLIB_initialization_by_elem_w_bubble(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);
    else
      FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);
    
    long *nod = (fe.node_id).m_pdata; // list of node ids in this element
    
    /* Element Dof */
    long ndofe = get_ndof_on_elem_nodes(fe.nne,nod,grid->node,fv->ndofn);
    
    long *cn = aloc1l (ndofe);
    get_dof_ids_on_elem_nodes(0,fe.nne,fv->ndofn,nod,node,cn,mp_id);
    
    // element tangent
    double *floc = aloc1 (ndofe);
    double *rloc = aloc1 (ndofe);
    double *r_e = aloc1 (ndofe);
    
    /* Coordinates of nodes */
    if(total_Lagrangian)
      def_elem (cn,ndofe,fv->u_np1,elem,node,r_e,sup,1);
    
    int nVol = N_VOL_TREE_FIELD;
    double lm = 0.0;
    Matrix(double) lk;
    Matrix_construct_redim(double,lk,ndofe,ndofe);
    Matrix_init(lk, 0.0);

    err += el_compute_stiffmat_MP(&fe,lk.m_pdata,grid,mat,fv,sol,load,
                                  crpl,opts,mp,mp_id,dt,lm,be,r_e);    
    
    /* get the disp increment from BC */
    {
      int k = 0;
      int jj = 0;
      for (int ii=0;ii<fe.nne;ii++)
      {
        for(jj=0;jj<fv->ndofn;jj++){
          if (jj  < fe.nsd && cn[k+jj] <= -1)
            rloc[k+jj] = sup->defl_d[abs(cn[k+jj])-1];
          else
            rloc[k+jj] = 0.0;
        }
        k += jj;
      }
    }
            
    /* Matrix vector multiplication */
    mv(lk.m_pdata,rloc,floc,ndofe,ndofe);
            
    /* Localization */
    {
      int k = 0;
      int jj = 0;
      int II = 0;
      for (jj=0;jj<fe.nne;jj++){
        for (int ii=0;ii<fv->ndofn;ii++){
          II = node[nod[jj]].id_map[mp_id].id[ii]-1;
          if (II < 0)  continue;
          fv->f_defl[II] += floc[k+ii];
        }/*end ii*/
        k += fv->ndofn;
      }/*end jj*/
    }
    
    /*  dealocation  */
    dealoc1l (cn);
    Matrix_cleanup(lk);
    dealoc1(floc);
    dealoc1(rloc);
    dealoc1 (r_e);
    FEMLIB_destruct(&fe);
    if(err != 0) return err;
  }/* end i (each volume element) */
  
  int ndn = 3;
  /* ADDED 1/7/2013 MM */
  for(int i=0; i<sup->nd_be; i++)
  {
    /* get pointers and constant vaules for the bounding element we
     * are working on. */
    const int be_id = sup->lbepd[i];
    const BOUNDING_ELEMENT *ptr_be = &grid->b_elems[be_id];
    const int ve_id = ptr_be->vol_elem_id;
    const ELEMENT *ptr_ve = &elem[ve_id];
    const long *ve_nod = ptr_ve->nod;
    const int nne_ve = ptr_ve->toe;
    
    /* get ndofs on element */
    int ndof_ve = get_ndof_on_bnd_elem(node,ptr_be,elem,fv->ndofn);

    /* get coordinates of volume element nodes */
    double *x = aloc1(nne_ve);
    double *y = aloc1(nne_ve);
    double *z = aloc1(nne_ve);
    
    switch(opts->analysis_type){
      case DISP: case TF:
        nodecoord_total (nne_ve,ve_nod,node,x,y,z);
        break;
      case CM:
      {
        switch(opts->cm)
        {
          case UPDATED_LAGRANGIAN:
            nodecoord_updated(nne_ve,ve_nod,node,x,y,z);
            break;
          case TOTAL_LAGRANGIAN:
            nodecoord_total(nne_ve,ve_nod,node,x,y,z);
            break;
          case MIXED_ANALYSIS_MODE:
          default:
            nodecoord_total(nne_ve,ve_nod,node,x,y,z);
            break;
        }
        break;
      }
      default:
        nodecoord_updated (nne_ve,ve_nod,node,x,y,z);
        break;
    }
    
    /* allocate space for stiffness and localization. */
    double *lk = aloc1(ndof_ve*ndof_ve);
    double *floc = aloc1(ndof_ve);
    double *rloc = aloc1(ndof_ve);
    long *cn_ve = aloc1l(ndof_ve);
    get_dof_ids_on_bnd_elem(0,fv->ndofn,node,ptr_be,elem,cn_ve,mp_id);
    /* get displacements */
    double *ve_disp = aloc1(ndof_ve);
    if(opts->analysis_type == DISP){ /* TOTAL LAGRANGIAN formulation */
      def_elem(cn_ve,ndof_ve,fv->u_np1,NULL,NULL,ve_disp,sup,1);
    }
    
    /* compute element stiffness */
    if(opts->analysis_type == DISP){
      err += DISP_stiffmat_bnd_el(lk,be_id,fv->ndofn,ndof_ve,x,y,z,grid->b_elems,
              elem,mat->hommat,node,fv->eps,fv->sig,sup,ve_disp);
    } else {
      /* not implemented, do nothing */
    }
    
    /* get the local disp increment from BC */
    {
      int k = 0;
      int jj = 0;
      for (int ii=0;ii<nne_ve;ii++){
        for (jj=0;jj<fv->ndofn;jj++){
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
        for (int ii=0;ii<fv->ndofn;ii++){
          II = node[ve_nod[jj]].id_map[mp_id].id[ii]-1;
          if (II < 0)  continue;
          fv->f_defl[II] += floc[j+ii];
        }/*end ii*/
        j += fv->ndofn;
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

/// Compute load vector for prescribed BCs(Dirichlet)
/// This compute load, f, as below:
/// [Kii Kio]<ui>   <bi>
/// [Koi Koo]<uo> = <bo>
/// [Kii][ui] = <bi> - [Kio]<uo>
/// where f = [Kio]<uo>, uo is Drichlet BCs
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] sol solution scheme object
/// \param[in] load object for loading
/// \param[in] dt time step
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int compute_load_vector_for_prescribed_BC(GRID *grid,
                                          MATERIAL_PROPERTY *mat,
                                          FIELD_VARIABLES *fv,
                                          SOLVER_OPTIONS *sol,
                                          LOADING_STEPS *load,
                                          double dt,                                          
                                          CRPL *crpl,
                                          const PGFem3D_opt *opts,
                                          MULTIPHYSICS *mp,
                                          int mp_id,
                                          int myrank)
{
  int err = 0;
  switch(mp->physics_ids[mp_id])
  {
    case MULTIPHYSICS_MECHANICAL:
      err + momentum_equation_load4pBCs(grid,
                                        mat,
                                        fv,
                                        sol,
                                        load,
                                        dt,                                          
                                        crpl,
                                        opts,
                                        mp,
                                        mp_id,
                                        myrank);
      break;                                  
    case MULTIPHYSICS_THERMAL:
      err += energy_equation_compute_load4pBCs(grid,
                                               mat,
                                               fv,
                                               sol,
                                               load,
                                               myrank,
                                               opts,
                                               mp_id,
                                               dt);
      break;                                         
    default:
      printf("%s is not supported\n", mp->physicsname[mp_id]);                                               
  
  }
  return err;
}

void load_vec_elem_sur (double *f,
			const long nle_s,
			const long ndofn,
			const ELEMENT *elem,
			const ZATELEM *zele_s)
{
}