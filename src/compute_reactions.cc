/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "compute_reactions.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "allocation.h"
#include "displacement_based_element.h"
#include "elem3d.h"
#include "enumerations.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "matice.h"
#include "resid_on_elem.h"
#include "stabilized.h"
#include "utils.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

int compute_reactions(long ne,
                      long ndofn,
                      long npres,
                      double *r,
                      Node *node,
                      Element *elem,
                      MATGEOM matgeom,
                      HOMMAT *hommat,
                      SUPP sup,
                      EPS *eps,
                      SIG *sig,
                      double nor_min,
                      CRPL *crpl,
                      double dt,
                      double stab,
                      CommunicationStructure *com,
                      const int analysis,
                      const int mp_id)
{
  int err = 0;

  if(sup->npd > 0){
    double *rxn = aloc1(sup->npd);

    int myrank;
    myrank = com->rank;

    for (int i=0;i<sup->nde;i++){
      const int elem_id = sup->lepd[i];
      const Element *ptr_elem = &elem[elem_id];

      /* Number of element nodes */
      const int nne = ptr_elem->toe;
      const int nne_t = nne + ptr_elem->n_bub;
      /* Nodes on element */
      long *nod = ptr_elem->nod;

      /* Element Dof */
      const int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

      /* allocation */
      double *r_e = aloc1 (ndofe);
      double *fe = aloc1 (ndofe);
      long *cn = aloc1l (ndofe);

      double *x,*y,*z;
      switch(analysis){
       case MINI:
       case MINI_3F:
        x = aloc1 (nne_t);
        y = aloc1 (nne_t);
        z = aloc1 (nne_t);
        break;
       default:
        x = aloc1 (nne);
        y = aloc1 (nne);
        z = aloc1 (nne);
        break;
      }

      /* coordinates */
      if(analysis == DISP){
        nodecoord_total (nne,nod,node,x,y,z);
      } else {
        nodecoord_updated (nne,nod,node,x,y,z);
      }
      if(analysis == MINI
         || analysis == MINI_3F){/* P1+B/P1 */
        element_center(nne,x,y,z);
      }

      get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

      /* After convergence, no increment of deformation, r_e = 0 */

      /* Residuals on element */
      switch(analysis){
       case STABILIZED:
        err = resid_st_elem (elem_id,ndofn,nne,elem,nod,node,hommat,
                             x,y,z,eps,sig,sup,r_e,nor_min,fe,dt,stab);
        break;
       case MINI:
        MINI_resid_el(fe,elem_id,ndofn,nne,x,y,z,elem,
                      nod,node,hommat,eps,sig,r_e);
        break;
       case MINI_3F:
        MINI_3f_resid_el(fe,elem_id,ndofn,nne,x,y,z,elem,
                         nod,node,hommat,eps,sig,r_e);
        break;
       case DISP:
         {
           /* Get TOTAL deformation on element; r_e already contains
              INCREMENT of deformation, add the deformation from previous. */
           double *r_en;
           r_en = aloc1(ndofe);
           def_elem (cn,ndofe,r,elem,node,r_en,sup,1);
           vvplus(r_e,r_en,ndofe);
           err =  DISP_resid_el(fe,elem_id,ndofn,nne,x,y,z,elem,
                                hommat,nod,node,eps,sig,sup,r_e,dt);
           free(r_en);
         }
         break;
       default:
        resid_on_elem (elem_id,ndofn,nne,nod,elem,node,matgeom,
                       hommat,x,y,z,eps,sig,r_e,npres,
                       nor_min,fe,crpl,dt,analysis);
        break;
      }

      /* fe contains the local residual vector on the element. */
      int j = 0;
      for (int k=0; k<nne; k++){
        for(int kk=0; kk<ndofn; kk++){
          long II = node[nod[k]].id_map[mp_id].id[kk];
          if (II <= -1){ /* dof is prescribed displacement */
            rxn[abs(II + 1)] -= fe[j+kk]; /* add value to appropriate pre. disp. */
          }
        }
        j += ndofn;
      }

      /* deallocate */
      /* free(nod); */
      free(r_e);
      free(x);
      free(y);
      free(z);
      free(fe);
      free(cn);

    }/* for each element in list */

    double *Grxn = aloc1(sup->npd);
    com->net->reduce(rxn,Grxn,sup->npd,NET_DT_DOUBLE,NET_OP_SUM,0,com->comm);
		     
    if(myrank == 0){
      printf("Reactions: ");
      print_array_d(stdout,Grxn,sup->npd,1,sup->npd);
    }

    free(rxn);
    free(Grxn);
  }

  return err;
}

