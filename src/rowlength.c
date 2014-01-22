#include "rowlength.h"
#include <stdlib.h>

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef GET_DOF_IDS_ON_ELEM_H
#include "get_dof_ids_on_elem.h"
#endif

#ifndef GET_NDOF_ON_ELEM_H
#include "get_ndof_on_elem.h"
#endif

void rowlength (long *adr,
		long ne,
		long ndofn,
		long ndof,
		NODE *node,
		ELEMENT *elem,
		long gr4)
/** Computes lengths of the rows. */
{
  long i,j,k,nne,ndofe,min;
  long *nod,*cn;
  
  nod = (long*) PGFEM_calloc (10,sizeof(long));
  cn = (long*) PGFEM_calloc (30,sizeof(long));
  
  for (i=0;i<ne;i++){
    
    /* Number of element nodes */
    nne = elem[i].toe;
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    ndofe = get_ndof_on_elem_nodes(nne,nod,node);
    
    if (gr4 == 0) { for (k=0;k<nne;k++) cn[k] = nod[k]+1; }
    else get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
	 
    min = ndof+4;
    
    for (k=0;k<ndofe;k++){
      if (cn[k] < min && cn[k] != 0 && cn[k] > -1)  min = cn[k];
    }
    for (k=0;k<ndofe;k++){
      if (cn[k]-1 <= -1)  continue;
      if (adr[cn[k]-1] < cn[k]-min+1)  adr[cn[k]-1] = cn[k]-min+1;
    }
  }
  free (nod); free (cn);
}
