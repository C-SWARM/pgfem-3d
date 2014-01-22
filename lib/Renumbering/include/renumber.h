#ifndef RENUMBER_H
#define RENUMBER_H

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef _MOPT2
#define _MOPT2 0
#undef _MOPT1
#endif /* #ifndef _MOPT2 */

#ifndef _MOPT3
#define _MOPT3 0
#undef _MOPT1
#endif /* #ifndef _MOPT3 */

#ifndef _MOPT1
#define _MOPT1 0
#endif /* #ifndef _MOPT1 */

  /** Runs the ParMETIS fill reducing function. */
void renumber(int *dist, int *Ap, int *Ai, int *l_order,
	      int *sizes, int *order, MPI_Comm comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef RENUMBER_H */
