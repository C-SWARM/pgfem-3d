#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "rowlength.h"
#include "allocation.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "utils.h"
#include <cstdlib>

void rowlength (long *adr,
                long ne,
                long ndofn,
                long ndof,
                Node *node,
                Element *elem,
                long gr4,
                const int mp_id)
{
  long i,k,nne,ndofe,min;
  long *nod,*cn;

  nod = PGFEM_calloc (long, 10);
  cn = PGFEM_calloc (long, 30);

  for (i=0;i<ne;i++){

    /* Number of element nodes */
    nne = elem[i].toe;
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

    if (gr4 == 0) {
      for (k=0;k<nne;k++) {
        cn[k] = nod[k]+1;
      }
    }
    else {
      get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);
    }

    min = ndof+4;

    for (k=0;k<ndofe;k++){
      if (cn[k] < min && cn[k] != 0 && cn[k] > -1)  {
        min = cn[k];
      }
    }
    for (k=0;k<ndofe;k++){
      if (cn[k]-1 <= -1)  {
        continue;
      }
      if (adr[cn[k]-1] < cn[k]-min+1) {
        adr[cn[k]-1] = cn[k]-min+1;
      }
    }
  }
  PGFEM_free (nod);
  PGFEM_free (cn);
}
