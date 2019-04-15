/**********************************
 *  Matrix and Vectors operations  *
 *  Karel Matous & Jaroslav Kruis  *
 **********************************/

#include "matice.h"
#include <math.h>

#include "pgfem3d/Communication.hpp"
#include "incl.h"
#include "resice.h"

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

using namespace pgfem3d;
using namespace multiscale::net;

int compare_val_w_key(const void *a,
              const void *b)
{
  return ( ((val_key*)a)->val - ((val_key*)b)->val);
}


int compare_long (const void * a, const void * b)
{
  return (*(long*)a - *(long*)b);
}

int compare_int (const void *a,const void *b)
{
  return (*(int*)a - *(int*)b);
}

int compare_double (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

long round1 (double a)
{
  if ( a >= 0 ) return (long)(a + 0.5);
  return (long)(a - 0.5);
}

void nas_AB (double **A,double **B,double **C,long m,long n,long p)
/*
   Procedura na nasobeni matic typu A(m,n)*B(n,p)=C(m,p), plne matice
*/
{

  long i,j,o;

  for(i=0;i<m;i++){
    for(j=0;j<p;j++){
      C[i][j]=0;
    }
  }

  for(i=0;i<m;i++){
    for(j=0;j<p;j++){
      for(o=0;o<n;o++){
    C[i][j] += A[i][o]*B[o][j];
      }
    }
  }
}

void nas_ATB (double **AT,double **B,double **C,long m,long n,long p)
/*
   Procedura na nasobeni matic typu
   A_T(m,n)*B(n,p)=C(m,p), plne matice
*/
{

  long i,j,o;

  for(i=0;i<m;i++){
    for(j=0;j<p;j++){
      C[i][j]=0;
    }
  }

  for(i=0;i<m;i++){
    for(j=0;j<p;j++){
      for(o=0;o<n;o++){
    C[i][j] += AT[o][i]*B[o][j];
      }
    }
  }
}

void nas_ABT (double **AT,double **B,double **C,long m,long n,long p)
/*
   Procedura na nasobeni matic typu A(m,n)*B_T(n,p)=C(m,p), plne
   matice
*/
{

  long i,j,o;

  for(i=0;i<m;i++){
    for(j=0;j<p;j++){
      C[i][j]=0;
    }
  }

  for(i=0;i<m;i++){
    for(j=0;j<p;j++){
      for(o=0;o<n;o++){
    C[i][j] += AT[i][o]*B[j][o];
      }
    }
  }
}

void inv_I (double **A,double **I,long m)
/*
   Procedura na sestaveni Inverzni matice typu A(m,n), plna matice
*/
{

  long i,j,pom;
  double *a,*e,*c;

  a = PGFEM_calloc (double, m*m);
  e = PGFEM_calloc (double, m*m);
  c = PGFEM_calloc (double, m*m);

  pom=0;
  for (i=0;i<m;i++){
    for (j=0;j<m;j++){
      *(a+pom)=A[i][j];
      pom++;
    }
  }

  pom=0;
  for (i=0;i<m*m;i++){
    if ((m+1)*pom==i){
      *(e+i)=1;
      pom++;
    }
    else
      *(e+i)=0;
  }

  reseni_rovnic (a,c,e,m,m,2);

  pom=0;
  for (i=0;i<m;i++){
    for (j=0;j<m;j++){
      I[i][j]=*(c+pom);
      pom++;
    }
  }

  PGFEM_free(a);PGFEM_free(e);PGFEM_free(c);

}

void vvplus (double *a,double *b,long n)
/*
  a = a + b
*/
{
  long i;
  for (i=0;i<n;i++) a[i] += b[i];
}

void vvminus (double *a,double *b,long n)
/*
  a = a - b
*/
{
  long i;
  for (i=0;i<n;i++) a[i] -= b[i];
}

void mv (double *a,double *b,double *c,long m,long n)
/*
  funkce nasobi matici A(m,n) s vektorem b(n,1)
  A(m,n).b(n,1)=c(m,1)

  vystup
  c - vysledny vektor

  vstupy
  a - matice A
  b - vektor b
  m,n - rozmer matice A; A stored by rows

  19.2.1997
*/
{
  long i,j,aca;
  double s;

  aca=0;
  for (i=0;i<m;i++){
    s=0.0;
    for (j=0;j<n;j++){
      s+=a[aca]*b[j];
      aca++;
    }
    c[i]=s;
  }
}

void mtv (double *a,double *b,double *c,long m,long n)
/*
  A^T . b = c
  A(m,n), b(m,1), c(n,1)

  A is stored by rows
*/
{
  long i,j,ii;
  double s;

  for (i=0;i<n;i++){
    s=0.0;  ii=i;
    for (j=0;j<m;j++){
      s+=a[ii]*b[j];
      ii+=n;
    }
    c[i]=s;
  }
}

void mvc (double *a,double *b,double *c,long m,long n)
/*
  funkce nasobi matici A vektorem b, vysledek je vektor c
  matice A je ulozena po sloupcich

  vstupy
  a - pole obsahujici matici A(m,n)
  b - pole obsahujici vektor b(n,1)
  m,n - rozmery matice a vektoru

  vystup
  c - pole obsahujici vysledny vektor c(m,1)

  15.12.1998
*/
{
  long i,j,k;
  double s;

  for (i=0;i<m;i++){
    s=0.0;  k=i;
    for (j=0;j<n;j++){
      s+=a[k]*b[j];  k+=m;
    }
    c[i]=s;
  }
}

void mtvc (double *a,double *b,double *c,long m,long n)
/*
  funkce nasobi transponovanou matici a vektor A^T b = c

  matice A je ulozena po sloupcich

  vstupy
  a - pole obsahujici matici A(m,n)
  b - pole obsahujici vektor b(m,1)
  m,n - rozmery matice A

  vystup
  c - pole obsahujici vektor c(n,1)

  7.12.1998
*/
{
  long i,j,k;
  double s;

  for (i=0;i<n;i++){
    s=0.0;  k=i*m;
    for (j=0;j<m;j++){
      s+=a[k]*b[j];  k++;
    }
    c[i]=s;
  }
}

void mm (double *a,double *b,double *c,long l,long m,long n)
/*
  funkce pocita soucin A.B=C

  A, B stored by rows

  19.2.1997
*/
{
  long i,j,k,ac,acb,acu,acl;
  double s;

  acl=0;  ac=0;
  for (i=0;i<l;i++){
    acu=acl+m;
    for (j=0;j<n;j++){
      s=0.0;
      acb=j;
      for (k=acl;k<acu;k++){
    s+=a[k]*b[acb];
    acb+=n;
      }
      c[ac]=s;
      ac++;
    }
    acl+=m;
  }
}

void mmt (double *a,double *b,double *c,long ra,long ca,long rb)
/*
  funkce provadi soucin A . B^T = C

  A, B stored by rows

  a(ra,ca) - matice A
  b(rb,ca) - matice B
  c(ra,rb) - matice C

  9.11.1999
*/
{
  long i,j,k,ii,kk,lk,uk;
  double s;

  ii=0;
  for (i=0;i<ra;i++){
    lk=i*ca;  uk=lk+ca;
    for (j=0;j<rb;j++){
      s=0.0;  kk=j*ca;
      for (k=lk;k<uk;k++){
    s+=a[k]*b[kk];  kk++;
      }
      c[ii]=s;  ii++;
    }
  }
}

void mtm (double *a,double *b,double *c,long ra,long ca,long cb)
/*
  funkce provadi soucin A^T . B = C

  A, B stored by rows

  a(ra,ca), b(ra,cb), c(ca,cb)
*/
{
  long i,j,k,ii;
  double s;

  ii=0;
  for (i=0;i<ca;i++){
    for (j=0;j<cb;j++){
      s=0.0;
      for (k=0;k<ra;k++){
    s+=a[k*ca+i]*b[k*cb+j];
      }
      c[ii]=s;  ii++;
    }
  }
}

void mtmccr (double *a,double *b,double *c,long m,long n,long p)
/*
  funkce nasobi matice A^T . B = C

  matice A a B jsou ulozeny po sloupcich
  vysledna matice C je ulozena po radcich

  vstupy
  a - pole obsahujici matici A(m,n)
  b - pole obsahujici matici B(m,p)
  m,n,p - rozmery matic A, B, C

  vystup
  c - pole obsahujici matici C(n,p)

  7.12.1998
*/
{
  long i,j,k,ii,lk,uk,ac;
  double s;

  ac=0;
  for (i=0;i<n;i++){
    for (j=0;j<p;j++){
      s=0.0;  lk=i*m;  uk=lk+m;  ii=j*m;
      for (k=lk;k<uk;k++){
    s+=a[k]*b[ii];  ii++;
      }
      c[ac]=s;  ac++;
    }
  }

}

void copyi (long *a,long *b,long n)
/*
  funkce kopiruje celociselne pole b do pole a
  obe pole maji n slozek

  31.5.1998
*/
{
  long i;
  for (i=0;i<n;i++){
    a[i]=b[i];
  }
}

void copyd (double *a,double *b,long n)
/*
  funkce kopiruje celociselne pole b do pole a
  obe pole maji n slozek

  31.5.1998
*/
{
  long i;
  for (i=0;i<n;i++){
    a[i]=b[i];
  }
}

double ss (const double *a, const double *b,const long n)
/*
  funkce pocita skalarni soucin vektoru a(n) a b(n)

  // TRANSLATION: computes scalar product of vectors

  5.4.1996
*/
{
  long i;
  double s;

  s = 0.0; for (i=0;i<n;i++) s += a[i]*b[i];

  return (s);
}

void nor_vec (double *BS_a,long n,CommunicationStructure *com)
/*
  Return normalized global vector BS_a
*/
{
  double nor,tmp;
  long i;
  int myrank = com->rank;
  
  tmp = ss (BS_a,BS_a,n);
  com->net->allreduce(&tmp,&nor,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  if (nor<1.0e-20){
    if (myrank == 0) PGFEM_printf ("Zero norm in routine normal\n");
    PGFEM_Comm_code_abort (com, 0);
  }

  tmp = sqrt(nor);
  for (i=0;i<n;i++) BS_a[i] /= tmp;
}

void nor_vec_serial (double *a,long n,int myrank)
/*
  Return normalized local vector a
*/
{
  double nor;
  long i;

  nor = ss (a,a,n);
  if (nor<1.0e-20){
    if (myrank == 0) PGFEM_printf ("Zero norm in routine %s\n",__func__);
    PGFEM_Abort();
  }

  nor = sqrt(nor);
  for (i=0;i<n;i++) a[i] /= nor;
}

void normalization (double *a,double nor,long n)
/*

 */
{
  long i;

  if (nor<1.0e-20){
    PGFEM_printf ("\n\n do funkce normovani byla zaslana nekladna norma");
    PGFEM_printf ("\n program konci.\n");
    abort ();
  }

  for (i=0;i<n;i++){
    a[i]/=nor;
  }
}

void mv_sky (double *a,double *b,double *c,long *adr,long n)
/*
  funkce nasobi matici A(n,n) s vektorem B(n,1)
  matice A je ulozena ve skylinu

  a - matice ulozena ve skylinu
  b - vektor
  c - vysledny vektor
  adr - pole adres diagonalnich prvku
  n - rozmer matice

  funkce byla testovana 24.7.1996 programem test_fin_res.c v /u/jk/TESTY/METODY

  **** otestovano ****
  */
{
  long i,j,acb,aci,aci1;
  double s,g;

  for (i=0;i<n;i++){
    aci=adr[i];  aci1=adr[i+1]-1;
    g=b[i];  s=0.0;  acb=i-aci1+aci;
    for (j=aci1;j>aci;j--){
      s+=a[j]*b[acb];
      c[acb]+=a[j]*g;
      acb++;
    }
    c[i]=s+a[aci]*g;
  }
}

void utv (double *a,double *b,double *c,long *adr,long n)
/*
  a - upper triangular matrix in skyline storage
  b - vector
  c - vector
  adr - addresses of diagonal elements
  n - dimension

  | 1 x x x x x x |
  | 0 1 x x x x x |
  | 0 0 1 x x x x |
  A = | 0 0 0 1 x x x |
  | 0 0 0 0 1 x x |
  | 0 0 0 0 0 1 x |
  | 0 0 0 0 0 0 1 |
*/
{
  long i,j,k,lj,uj;
  double s;

  nulld (c,n);

  for (i=n-1;i>-1;i--){
    s=b[i];  k=i-1;
    c[i]+=b[i];
    lj=adr[i]+1;  uj=adr[i+1];
    for (j=lj;j<uj;j++){
      c[k]+=a[j]*s;  k--;
    }
  }
}

void ltv (double *a,double *b,double *c,long *adr,long n)
/*
  a - lower triagular matrix in skyline storage
  b - vector
  c - vector
  adr - addresses of diagonal elements
  n - dimension

  | 1 0 0 0 0 0 0 |
  | x 1 0 0 0 0 0 |
  A = | x x 1 0 0 0 0 |
  | x x x 1 0 0 0 |
  | x x x x 1 0 0 |
  | x x x x x 1 0 |
  | x x x x x x 1 |
*/
{
  long i,j,k,lj,uj;
  double s;

  nulld (c,n);

  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1]-1;
    k=i-(uj-lj);
    s=0.0;
    for (j=uj;j>lj;j--){
      s+=a[j]*b[k];  k++;
    }
    c[i]=s+b[i];
  }
}

void per_sym (double ***e)
/*
  Function returns permutation symbol Eijk
  http://mathworld.wolfram.com/PermutationSymbol.html
*/
{

  e[0][1][2] = e[1][2][0] = e[2][0][1] = +1.0;
  e[0][2][1] = e[2][1][0] = e[1][0][2] = -1.0;
}

void cross_product (double *a,double *b,double *c)
/*
  Cross product
  c = a x b
*/
{
  long i,j,k;
  double ***e;

  e = aloc3 (3,3,3);

  per_sym (e);

  for (i=0;i<3;i++){
    c[i] = 0.0;
    for (j=0;j<3;j++){
      for (k=0;k<3;k++){
    c[i] += e[i][j][k]*a[j]*b[k];
      }
    }
  }

  dealoc3 (e,3,3);
}

void tred2 (double **a,int n,double *d,double *e)
/*
  Householder reduction of a real, symmetric matrix
  a[1..n][1..n]. On output, a is replaced by the orthogonal
  matrix Q effecting the transformation. d[1..n] returns the
  diagonal elements of the tridiagonal matrix, and e[1..n] the
  off-diagonal elements, with e[1]=0. Several statements, as
  noted in comments, can be omitted if only eigenvalues are to be
  found, in which case a contains no useful information on
  output. Otherwise they are to be included.
*/
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
    scale += fabs(a[i][k]);
      if (scale == 0.0)
    e[i]=a[i][l];
      else {
    for (k=1;k<=l;k++) {
      a[i][k] /= scale;
      h += a[i][k]*a[i][k];
    }
    f=a[i][l];
    g = f>0 ? -sqrt(h) : sqrt(h);
    e[i]=scale*g;
    h -= f*g;
    a[i][l]=f-g;
    f=0.0;
    for (j=1;j<=l;j++) {
      /* Next statement can be omitted if eigenvectors not wanted */
      a[j][i]=a[i][j]/h;
      g=0.0;
      for (k=1;k<=j;k++)
        g += a[j][k]*a[i][k];
      for (k=j+1;k<=l;k++)
        g += a[k][j]*a[i][k];
      e[j]=g/h;
      f += e[j]*a[i][j];
    }
    hh=f/(h+h);
    for (j=1;j<=l;j++) {
      f=a[i][j];
      e[j]=g=e[j]-hh*f;
      for (k=1;k<=j;k++)
        a[j][k] -= (f*e[k]+g*a[i][k]);
    }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  /* Next statement can be omitted if eigenvectors not wanted */
  d[1]=0.0;
  e[1]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not wanted
     except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
    g=0.0;
    for (k=1;k<=l;k++)
      g += a[i][k]*a[k][j];
    for (k=1;k<=l;k++)
      a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
}

void tqli (double *d,double *e,int n,double **z)
/*
  QL algorithm with implicit shifts, to determine the eigenvalues
  and eigenvectors of a real, symmetric, tridiagonal matrix, or
  of a real, symmetric matrix previously reduced by tred2
  ¡×11.2. On input, d[1..n] contains the diagonal elements of the
  tridiagonal matrix. On output, it returns the eigenvalues. The
  vector e[1..n] inputs the subdiagonal elements of the
  tridiagonal matrix, with e[1] arbitrary. On output e is
  destroyed. When finding only the eigenvalues, several lines may
  be omitted, as noted in the comments. If the eigenvectors of a
  tridiagonal matrix are desired, the matrix z[1..n][1..n] is
  input as the identity matrix. If the eigenvectors of a matrix
  that has been reduced by tred2 are required, then z is input as
  the matrix output by tred2.  In either case, the kth column of
  z returns the normalized eigenvector corresponding to d[k].
*/
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
    dd=fabs(d[m])+fabs(d[m+1]);
    if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
    if (iter++ == 30) {
      PGFEM_printf("Too many iterations in TQLI");
      abort();
    }
    g=(d[l+1]-d[l])/(2.0*e[l]);
    r=sqrt((g*g)+1.0);
    g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
    s=c=1.0;
    p=0.0;
    for (i=m-1;i>=l;i--) {
      f=s*e[i];
      b=c*e[i];
      if (fabs(f) >= fabs(g)) {
        c=g/f;
        r=sqrt((c*c)+1.0);
        e[i+1]=f*r;
        c *= (s=1.0/r);
      } else {
        s=f/g;
        r=sqrt((s*s)+1.0);
        e[i+1]=g*r;
        s *= (c=1.0/r);
      }
      g=d[i+1]-p;
      r=(d[i]-g)*s+2.0*c*b;
      p=s*r;
      d[i+1]=g+p;
      g=c*r-b;
      /* Next loop can be omitted if eigenvectors not wanted */
      for (k=1;k<=n;k++) {
        f=z[k][i+1];
        z[k][i+1]=s*z[k][i]+c*f;
        z[k][i]=c*z[k][i]-s*f;
      }
    }
    d[l]=d[l]-p;
    e[l]=g;
    e[m]=0.0;
      }
    } while (m != l);
  }
}
