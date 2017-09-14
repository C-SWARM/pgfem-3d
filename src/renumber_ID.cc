#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "renumber_ID.h"

void renumber_ID (int ndofn, int nn, Node *node, int *g_order, MPI_Comm comm,
                  const int mp_id)
{
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  for (int i = 0, e = nn; i < e; ++i) {
    if (node[i].Dom != myrank) continue;
    for (int k = 0, e = ndofn; k < e; ++k) {
      if (node[i].id_map[mp_id].Gid[k] > 0) {
        node[i].id_map[mp_id].Gid[k] = g_order[node[i].id_map[mp_id].Gid[k]-1] +1;
      }
    }
  }
}
