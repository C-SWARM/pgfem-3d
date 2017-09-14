#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "incl.h"
#include "elem3d.h"
#include "enumerations.h"
#include <cstring>

void build_elem_inelas (long ne,Element *elem)
{
  long i,j,II,nne;

  for (i=0;i<ne;i++){

    nne = elem[i].toe;

    /* Integration */
    int_point (nne,&II);

    elem[i].L  = PGFEM_calloc (double*, II);
    elem[i].LO = PGFEM_calloc (long, II);

    for (j=0;j<II;j++) {
      elem[i].L[j] = PGFEM_calloc (double, 21);
    }
  }

}

void build_pressure_nodes (long ne,
                           long npres,
                           Element *elem,
                           SIG *sig,
                           EPS *eps,
                           const int analysis)
{
  long i;

  /* allocate 1 if npres = 0 */
  if(npres == 0) npres = 1;

  for (i=0;i<ne;i++){
    switch(analysis){
     case FS_CRPL:
     case FINITE_STRAIN:
      eps[i].T   = PGFEM_calloc (double, npres);
      eps[i].d_T = PGFEM_calloc (double, npres);
      break;
     case MINI_3F:
      eps[i].T   = PGFEM_calloc (double, 4);
      eps[i].d_T = PGFEM_calloc (double, 4);
      break;
     default:
      eps[i].T   = NULL;
      eps[i].d_T = NULL;
      break;
    }

    sig[i].p    = PGFEM_calloc (double, npres);
    sig[i].d_p  = PGFEM_calloc (double, npres);
    sig[i].pn_1 = PGFEM_calloc (double, npres);
  }
}

void build_crystal_plast (long ne,
                          Element *elem,
                          SIG *sig,
                          EPS *eps,
                          CRPL *crpl,
                          const int analysis,
                          const int plc)
/*

 */
{
  long i,II,nne,k;
  if (analysis == FS_CRPL) {
    for (i=0;i<ne;i++){

      nne = elem[i].toe;

      /* Integration */
      int_point (nne,&II);

      for (k=0;k<II;k++) {
        sig[i].il[k].Tau = PGFEM_calloc (double, crpl[elem[i].mat[2]].nss);
        eps[i].il[k].GA  = PGFEM_calloc (double, crpl[elem[i].mat[2]].nss);
        eps[i].il[k].GA1 = PGFEM_calloc (double, crpl[elem[i].mat[2]].nss);
        if (plc == 1) {
          eps[i].il[k].PLC_B = PGFEM_calloc (long, crpl[elem[i].mat[2]].nss);
        }
      }
    }
  }
}

void nulld (double *a,long n)
{
  memset(a,0,n*sizeof(double));
}

void nulld2 (double **a,long m,long n)
{
  long i;

  for (i=0;i<m;i++){
    memset(a[i],0,n*sizeof(double));
  }
}

// #ifndef PGFEM_MACRO_ALLOCATION
// int* aloc1i (int m)
//      /*
//        m - Size of alocated array
//      */
// {
//   int *pom;

//   pom = (int*) calloc (m,sizeof(int));

//   if (pom == NULL){
//     PGFEM_printf ("\n Memory is full.\n"); fflush(stdout);
//     abort ();
//   }

//   return (pom);
// }

// void dealoc1i (int *a)
//      /*
//        Array for dealocation
//      */
// {
//   free (a);
// }

// long* aloc1l (long m)
//      /*
//        m - Size of alocated array
//      */
// {
//   long *pom;

//   if (m == 0) m = 1;

//   pom = (long*) calloc (m,sizeof(long));

//   if (pom == NULL){
//     PGFEM_printf ("\n Memory is full.\n"); fflush(stdout);
//     abort ();
//   }

//   return (pom);
// }

// void dealoc1l (long *a)
//      /*
//        Array for dealocation
//      */
// {
//   free (a);
// }

// long** aloc2l (long m,long n)
//      /*
//        m - Size of alocated array in one direction
//        n - Size of alocated array in two direction
//      */
// {
//   long i;
//   long **pom;

//   pom = (long**) calloc (m,sizeof(long*));

//   for (i=0;i<m;i++){
//     pom[i]= (long*) calloc (n,sizeof(long));
//   }

//   if (pom == NULL){
//     PGFEM_printf ("\n Memory is full.\n"); fflush(stdout);
//     abort ();
//   }

//   return (pom);
// }

// void dealoc2l (long **a,long m)
// {
//   long i;

//   for (i=0;i<m;i++){
//     free (a[i]);
//   }
//   free (a);
// }

// long*** aloc3l (long m,long n,long p)
//      /*

//       */
// {
//   long i,j;
//   long ***pom;

//   pom = (long***) calloc (m,sizeof(long**));

//   for (i=0;i<m;i++){
//     pom[i]= (long**) calloc (n,sizeof(long*));
//     for (j=0;j<n;j++){
//       pom[i][j] = (long*) calloc (p,sizeof(long));
//     }
//   }

//   if (pom == NULL){
//     PGFEM_printf ("\n Memory is full.\n");  fflush(stdout);
//     abort ();
//   }

//   return (pom);
// }

// void dealoc3l (long ***a,long m,long n)
// {
//   long i,j;

//   for (i=0;i<m;i++){
//     for (j=0;j<n;j++){
//       free (a[i][j]);
//     }
//     free (a[i]);
//   }
//   free (a);
// }

// double* aloc1 (long m)
//      /*

//       */
// {
//   double *pom;

//   if (m == 0) m = 1;

//   pom = (double*) calloc (m,sizeof(double));

//   if (pom == NULL){
//     PGFEM_printf ("\n Memory is full. \n");  fflush(stdout);
//     abort ();
//   }

//   return (pom);
// }

// void dealoc1 (double *a)
// {
//   free (a);
// }

// double** aloc2 (long m,long n)
//      /*

//       */
// {
//   long i;
//   double **pom;

//   pom = (double**) calloc (m,sizeof(double*));

//   for (i=0;i<m;i++){
//     pom[i] = (double*) calloc (n,sizeof(double));
//   }

//   if (pom == NULL){
//     PGFEM_printf ("\n Memory is full.\n");  fflush(stdout);
//     abort ();
//   }

//   return (pom);
// }

// void dealoc2 (double **a,long m)
// {
//   long i;

//   for (i=0;i<m;i++){
//     free (a[i]);
//   }
//   free (a);
// }

// double*** aloc3 (long m,long n,long p)
//      /*

//       */
// {
//   long i,j;
//   double ***pom;

//   pom = (double***) calloc (m,sizeof(double**));

//   for (i=0;i<m;i++){
//     pom[i]= (double**) calloc (n,sizeof(double*));
//     for (j=0;j<n;j++){
//       pom[i][j] = (double*) calloc (p,sizeof(double));
//     }
//   }

//   if (pom == NULL){
//     PGFEM_printf ("\n Memory is full.\n");  fflush(stdout);
//     abort ();
//   }

//   return (pom);
// }

// void dealoc3 (double ***a,long m,long n)
// {
//   long i,j;

//   for (i=0;i<m;i++){
//     for (j=0;j<n;j++){
//       free (a[i][j]);
//     }
//     free (a[i]);
//   }
//   free (a);
// }

// double**** aloc4 (long m,long n,long p,long q)
//      /*

//      */
// {
//   long i,j,k;
//   double ****pom;

//   pom = (double****) calloc (m,sizeof(double***));

//   for (i=0;i<m;i++){
//     pom[i] = (double***) calloc (n,sizeof(double**));
//     for (j=0;j<n;j++){
//       pom[i][j] = (double**) calloc (p,sizeof(double*));
//       for (k=0;k<p;k++){
//     pom[i][j][k] = (double*) calloc (q,sizeof(double));
//       }
//     }
//   }

//   if (pom==NULL){
//     PGFEM_printf ("\n Memory is full.\n");  fflush(stdout);
//     abort ();
//   }

//   return pom;
// }

// void dealoc4 (double ****a,long m,long n,long p)
// {
//   long i,j,k;

//   for (i=0;i<m;i++){
//     for (j=0;j<n;j++){
//       for (k=0;k<p;k++){
//     free (a[i][j][k]);
//       }
//       free (a[i][j]);
//     }
//     free (a[i]);
//   }
//   free (a);
// }
// #endif /* #ifndef PGFEM_MACRO_ALLOCATION */

