/********************
 *  Solvers         *
 *  Jaroslav Kruis  *
 ********************/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "resice.h"

#include "allocation.h"
#include "enumerations.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "matice.h"
#include "utils.h"
#include <cmath>

void reseni_rovnic (double *a,double *x,double *y,long n,long m,long as)
/* funkce provadi reseni maticove rovnice A.X=Y
   lze ji tedy bez problemu pouzit na vypocet inverzni matice
   matice A je plna

   A(n,n),X(n,m),Y(n,m)
   as - rozhodovaci konstanta
   as=1 - pivot se hleda jen v pripade, ze A(i,i)=0
   as=2 - pivot se hleda pokazde

   testovano 24.7.1996 pomoci programu test_fin_res.c v /u/jk/TESTY/METODY
   procedura dava stejne vysledky jako ldl,ldlkon,ldlblok,congrad_sky,congrad_comp

   **** otestovano ****
   */
{
  long    i,j,k,ac,acr,acc{},aca,aca1,acx,acy,acy1,aci,acj;
  long    *av;
  double  s,g;

  av = PGFEM_calloc (long, n);

  /*************************************************************************/
  /*  nastaveni hodnot vektoru, ktery udava poradi jednotlivych neznamych  */
  /*************************************************************************/
  for (i=0;i<n;i++){
    av[i]=i;
  }

  for (i=0;i<n-1;i++){
    acr=i;  acc=i;
    if (as==1){
      /******************************************/
      /*  pivot se hleda jen pokud je A(i,i)=0  */
      /******************************************/
      if (fabs(a[i*n+i])<1.0e-5){
        /*  vyber pivota  */
        s=0.0;
        /*  smycka pres radky  */
        for (j=i;j<n;j++){
          aca=j*n+i;
          /*  smycka pres sloupce  */
          for (k=i;k<n;k++){
            if (s<fabs(a[aca])){
              s=fabs(a[aca]);  acr=j;  acc=k;
            }
            aca++;
          }
        }
        if (s==0.0){
          PGFEM_printf ("\n Singularni matice v procedure gemp (krok %ld).\n",i);
          abort ();
        }
      }
    }
    if (as==2){
      /****************************/
      /*  pivot se hleda pokazde  */
      /****************************/
      s=0.0;
      /*  smycka pres radky  */
      for (j=i;j<n;j++){
        aca=j*n+i;
        /*  smycka pres sloupce  */
        for (k=i;k<n;k++){
          if (s<fabs(a[aca])){
            s=fabs(a[aca]);  acr=j;  acc=k;
          }
          aca++;
        }
      }
      if (s<1.0e-15){
        PGFEM_printf ("\n Singularni matice v procedure reseni_rovnic(krok %ld).\n",i);
        abort ();
      }
    }

    /******************/
    /*  vymena radku  */
    /******************/
    if (acr!=i){
      aca=i*n+i;  aca1=acr*n+i;
      for (j=i;j<n;j++){
        s=a[aca];
        a[aca]=a[aca1];
        a[aca1]=s;
        aca++;  aca1++;
      }
      acy=i*m;  acy1=acr*m;
      for (j=0;j<m;j++){
        s=y[acy];
        y[acy]=y[acy1];
        y[acy1]=s;
        acy++;  acy1++;
      }
    }
    /********************/
    /*  vymena sloupcu  */
    /********************/
    if (acc!=i){
      ac=av[i];
      av[i]=av[acc];
      av[acc]=ac;

      aca=i;  aca1=acc;
      for (j=0;j<n;j++){
        s=a[aca];
        a[aca]=a[aca1];
        a[aca1]=s;
        aca+=n;  aca1+=n;
      }
    }
    /***************/
    /*  eliminace  */
    /***************/

    for (j=i+1;j<n;j++){
      acj=j*n+i;  aci=i*n+i;
      s=a[acj]/a[aci];
      /*  modifikace matice A  */
      for (k=i;k<n;k++){
        a[acj]-=s*a[aci];
        acj++;  aci++;
      }
      acj=j*m;  aci=i*m;
      /*  modifikace matice pravych stran Y  */
      for (k=0;k<m;k++){
        y[acj]-=s*y[aci];
        acj++;  aci++;
      }
    }
  }

  /*****************/
  /*  zpetny chod  */
  /*****************/

  for (i=n-1;i>-1;i--){
    g=a[i*n+i];  acx=i*m;
    for (j=0;j<m;j++){
      s=0.0;  aca=i*n+i+1;  acy=(i+1)*m+j;
      for (k=i+1;k<n;k++){
        s+=a[aca]*x[acy];
        aca++;  acy+=m;
      }
      x[acx]=(y[acx]-s)/g;
      acx++;
    }
  }

  /***********************************/
  /*  prerovnani do puvodniho stavu  */
  /***********************************/
  for (i=0;i<n;i++){
    if (av[i]!=i){
      for (j=i;j<n;j++){
        if (av[j]==i){
          acc=j;  break;
        }
      }

      ac=av[i];
      av[i]=av[acc];
      av[acc]=ac;

      aca=i*m;  aca1=acc*m;
      for (j=0;j<m;j++){
        s=x[aca];
        x[aca]=x[aca1];
        x[aca1]=s;
        aca++;  aca1++;
      }
    }
  }

  PGFEM_free (av);
}


int resic_sky (double *a,double *x,double *y,long *adr,long n,long tc,const int analysis)

/*
  funkce resi soustavu linearnich algebraickych rovnic
  reseni se provadi rozkladem LDL
  matice soustavy je ulozena ve skylinu

  a - matice soustavy
  x - vektor reseni
  y - vektor prave strany
  adr - pole adres diagonalnich prvku
  n - pocet neznamych
  tc - typ vypoctu  tc=1 - provede se eliminace i zpetny chod
  tc=2 - provede se pouze eliminace
  tc=3 - provede se pouze zpetny chod

  10.7.1996
  funkce je totozna s procedurou ve fortranu, ktera je stejne
  rychla jako colsol od Batheho
  funkce byla testovana s vysledky od Batheho
  ldl ve fortranu
  stare ldl v c
*/
{
  int SINGULAR;
  long i,j,k,ac,ac1,ac2,acs,ack,ack1,acrk,aci,aci1,acri,acj,acj1;
  double s,g;

  SINGULAR = 0;

  if (tc==1 || tc==2){
    /*****************************/
    /*  rozklad matice soustavy  */
    /*****************************/
    for (k=1;k<n;k++){
      /*  smycka pres vsechny radky matice  */
      ack=adr[k];  ack1=adr[k+1];
      ac1=k+ack;   acrk=ac1-ack1+1;
      acj=ack1-2;
      /*  uprava mimodiagonalnich prvku k-teho sloupce  */
      for (i=acrk+1;i<k;i++){
        /*  smycka pres prvky k-teho sloupce  */
        aci=adr[i];  aci1=adr[i+1];
        ac2=i+aci;   acri=ac2-aci1+1;
        if (acri<acrk)  ac=acrk;
        else            ac=acri;
        acj1=ac1-ac;  acs=ac2-ac;
        s=0.0;
        for (j=acj1;j>acj;j--){
          s+=a[j]*a[acs];
          acs--;
        }
        a[acj]-=s;  acj--;
      }
      /*  uprava diagonalniho prvku v k-tem sloupci  */
      s=0.0;
      for (i=ack1-1;i>ack;i--){
        /*  smycka pres mimodiagonalni prvky k-teho sloupce  */
        ac1=adr[acrk];  acrk++;
        g=a[i];
        a[i]/=a[ac1];
        s+=a[i]*g;
      }
      a[ack]-=s;
      if (fabs(a[ack]) < 1.0e-10){
        if (SINGULAR == 0) PGFEM_printf ("\n\n Singular matrix in the LDL factorization. Progam Aborted. \n\n");
        SINGULAR = 1;
        if (analysis == ELASTIC) abort ();
      }
    }
  }
  if (tc==1 || tc==3){
    /********************/
    /*  vypocet reseni  */
    /********************/
    /*  vypocet  Lz=y => z (y se prepisuji na z) */
    for (k=1;k<n;k++){
      /*  smycka pres nezname  */
      ack=adr[k]+1;  ack1=adr[k+1];
      ac1=k-1;  s=0.0;
      for (i=ack;i<ack1;i++){
        s+=a[i]*y[ac1];
        ac1--;
      }
      y[k]-=s;
    }
    /*  deleni zbyle soustavy diagonalnimi prvky (DLx=z => Lx=1/Dz) */
    for (k=0;k<n;k++){
      ac=adr[k];  y[k]/=a[ac];
    }
    /*  vypocet Lx=1/Dz => x  */
    for (k=n-1;k>-1;k--){
      /*  smycka pres nezname  */
      ack=adr[k]+1;  ack1=adr[k+1];
      x[k]=y[k];  g=x[k];
      ac=k-1;
      for (i=ack;i<ack1;i++){
        y[ac]-=a[i]*g;
        ac--;
      }
    }
  }

  return (SINGULAR);
}

void lokalizace_scr (double *a,double *b,long *lcn,long *adr,long *ci,long n)
/*
  funkce lokalizuje matici B do matice A
  uklada se do modifikovanych kompresovanych radku

  a - matice, do ktere se lokalizuje
  b - matice, ktera se lokalizuje
  ci - pole sloupcovych indexu
  adr - pole adres zacatku radku
  lcn - pole kodovych cisel jednoho prvku
  n - rozmer matice b

  11.6.1997
*/
{
  long i,j,k,aci,acj,ack,acp,lk,uk;

  for (i=0;i<n;i++){
    aci=lcn[i]-1;
    if (aci<0) continue;
    acp=0;  lk=adr[aci];  uk=adr[aci+1];  ack=lk;
    for (j=0;j<n;j++){
      acj=lcn[j]-1;
      if (acj<0) continue;
      if (aci<acj)  continue;
      if (acp<acj){
        acp=acj;
        for (k=ack;k<uk;k++){
          if (ci[k]!=acj)  continue;
          else{
            a[k]+=b[i*n+j];
            ack=k;
            break;
          }
        }
      }
      else{
        acp=acj;
        for (k=ack;k>=lk;k--){
          if (ci[k]!=acj)  continue;
          else{
            a[k]+=b[i*n+j];
            ack=k;
            break;
          }
        }
      }
    }
  }
}


void strci_scr (long *adr,long ne,long ndofn,Element *elem,Node *node, const int mp_id)
/*
  funkce zjistuje prispevek k velikosti pomocneho pole
  sloupcovych indexu v ukladani symmetric compressed rows
  od jednoho typu prvku

  vystup
  adr - pole poctu prispevku v pomocnem poli sloupcovych indexu

  vstupy
  cn - pole kodovych cisel prvku
  ne - pocet prvku zpracovavaneho typu
  ndofe - pocet stupnu volnosti jednoho prvku

  29.5.1998
*/
{
  long i,j,k,lcnj,lcnk,*cn,nne,ndofe,*nod;

  cn = aloc1l (30); nod = aloc1l (10);

  for (i=0;i<ne;i++){

    /* Number of element nodes */
    nne = elem[i].toe;
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);
    /* Id numbers */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    for (j=0;j<ndofe;j++){
      lcnj=cn[j]-1;
      if (lcnj<0)  continue;
      for (k=0;k<ndofe;k++){
        lcnk=cn[k]-1;
        if (lcnk<0)  continue;
        if (lcnj<lcnk)  continue;
        adr[lcnj]++;
      }
    }
  }
  dealoc1l (cn);  dealoc1l (nod);
}

void aci_scr (long *ci,long *adr,long ne,long ndofn,Element *elem,Node *node, const int mp_id)
/*
  funkce sestavuje prispevky do pomocneho pole sloupcovych
  indexu v ukladani symmetric compressed rows od jednoho typu prvku

  vystup
  ci - pomocne pole sloupcovych indexu
  adr - pole aktualnich pozic v pomocnem poli sloupcovych indexu

  vstupy
  node - pole cisel uzlu na jednom typu prvku
  cn - pole kodovych cisel uzlu
  ne - pocet prvku zpracovavaneho typu
  ndofe - pocet stupnu volnosti jednoho prvku
  ndofn - pocet stupnu volnosti jednoho uzlu
  nne - pocet uzlu na jednom prvku

  29.5.1998
*/
{
  long i,j,k,lcnj,lcnk,*cn,nne,ndofe,*nod;

  cn = aloc1l (30); nod = aloc1l (10);

  for (i=0;i<ne;i++){

    /* Number of element nodes */
    nne = elem[i].toe;
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);
    /* Id numbers */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    for (j=0;j<ndofe;j++){
      lcnj=cn[j]-1;
      if (lcnj<0)  continue;
      for (k=0;k<ndofe;k++){
        lcnk=cn[k]-1;
        if (lcnk<0)  continue;
        if (lcnj<lcnk)  continue;
        ci[adr[lcnj]]=lcnk;
        adr[lcnj]++;
      }
    }
  }
  dealoc1l (cn);  dealoc1l (nod);
}

void sort_cr (long *ci,long *adrb,long *adre,long ndof)
/*
  funkce tridi a usporadava podle velikosti pomocne pole
  sloupcovych indexu v ukladani compressed row

  ci - pomocne pole sloupcovych indexu
  adrb - pole adres zacatku v poli ci
  adre - pole adres koncu v poli ci
  ndof - pocet stupnu volnosti problemu

  29.5.1998
*/
{
  long i,j,k,ii{},jj,lj,uj,min,prev;

  for (i=0;i<ndof;i++){
    lj=adrb[i];  uj=adre[i];  prev=-1;
    for (j=lj;j<uj;j++){
      min=100000000;
      for (k=j;k<uj;k++){
        if (ci[k]<min){
          min=ci[k];  ii=k;
        }
      }
      if (min==prev){
        uj--;  j--;
        ci[ii]=ci[uj];
      }
      else{
        jj=ci[j];
        ci[j]=min;
        ci[ii]=jj;
        prev=min;
      }
    }
    adre[i]=uj;
  }
}

long nonzero_sky (double *a,long *adr,double limit,long n)
/*
  funkce zjistuje pocet prvku v matici a, ktera je ulozena ve skyline,
  vetsich nez je limit

  a - matice ve skyline
  adr - pole adres diagonalnich prvku
  limit - velikost pro testovani
  n - rozmer matice

  8.11.1998
*/
{
  long i,j,lj,uj,m;

  m=0;
  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1];
    for (j=lj;j<uj;j++){
      if (fabs(a[j])>limit)  m++;
    }
  }
  return m;
}

void sky_scr (double *sky,long *adrs,double *scr,long *adrc,long *ci,
              double limit,long n)
/*
  funkce ulozi matici ve skyline do symmetric compressed row

  sky - matice ve skyline
  adrs - pole adres diagonalnich prvku
  scr - matice v symmetric compressed rows
  adrc - pole adres zacatku
  ci - pole sloupcovych indexu
  limit - hodnota pro testovani
  n - rozmer matice

  8.11.1998
*/
{
  long i,j,k,lj,uj,m;

  m=0;  adrc[0]=0;
  for (i=0;i<n;i++){
    lj=adrs[i];  uj=adrs[i+1]-1;  k=i-(uj-lj);
    for (j=uj;j>=lj;j--){
      if (fabs(sky[j])>limit){
        scr[m]=sky[j];  ci[m]=k;  m++;
      }
      k++;
    }
    adrc[i+1]=m;
  }
}

long size (long *adrb,long *adre,long n)
/*
  funkce pocita velikost pole z poli adres zacatku a koncu

  adrb - pole zacatku
  adre - pole koncu
  n - pocet slozek v polich adrb a adre

  31.5.1998
*/
{
  long i,j;

  j=0;
  for (i=0;i<n;i++){
    j+=adre[i]-adrb[i];
  }
  return j;
}

void colindex_cr (long *ci,long *adr,long *aci,long *adre,long ndof)
/*
  funkce sestavuje vysledne pole sloupcovych indexu

  vystupy
  ci - pole sloupcovych indexu
  adr - pole adres zacatku v poli ci

  vstupy
  adr - pole adres zacatku v poli aci
  aci - pomocne pole sloupcovych indexu
  adre - pole adres koncu v poli aci
  ndof - pocet stupnu volnosti celeho problemu

  29.5.1998
*/
{
  long i,j,ii,lj,uj;

  ii=0;  lj=0;  adr[0]=0;
  for (i=0;i<ndof;i++){
    uj=adre[i];
    for (j=lj;j<uj;j++){
      ci[ii]=aci[j];  ii++;
    }
    lj=adr[i+1];
    adr[i+1]=ii;
  }
}

void mv_scr (double *a,double *b,double *c,long *adr,long *ci,long n)
/*
  funkce nasobi matici a s vektorem b, vysledek je vektor c
  matice a je ulozena v modifikovanem compress rows

  vystup
  c - vysledny vektor

  vstupy
  a - matice v modifikovanem compress rows
  b - vektor
  adr - pole adres prvnich prvku v radku
  ci - pole sloupcovych indexu
  n - rozmer matice a

  28.7.1997
*/
{
  long i,j,ii,lj,uj;
  double s,d;

  for (i=0;i<n;i++){
    lj=adr[i];  uj=adr[i+1];
    s=0.0;  d=b[i];
    for (j=lj;j<uj;j++){
      ii=ci[j];
      s+=a[j]*b[ii];
      c[ii]+=a[j]*d;
    }
    c[i]=s;
  }
}

void minimize_cr (double *a,long *b,long *adrn,long *adro,long n,double limit)
/*
  funkce minimalizuje velikost poli a a b

  a - pole realnych cisel (prvky matice)
  b - pole celych cisel (sloupcove indexy)
  adrn - pole adres prvnich prvku v poli a po minimalizaci
  adro - pole adres prvnich prvku v poli a pred minimalizaci
  n - pocet radku matice a
  limit - konstanta vyberu

  30.7.1997
*/
{
  long i,j,n1,cor;

  cor=-1;  n1=0;
  for (i=0;i<n;i++){
    adrn[i]=adro[i]-n1;
    for (j=adro[i];j<adro[i+1];j++){
      cor++;
      if (fabs(a[j])>limit){
        a[cor]=a[j];
        b[cor]=b[j];
      }
      else{
        cor--;  n1++;
      }
    }
  }
  adrn[n]=adro[n]-n1;
  /*  return n1; */
}

double energie_scr (double *a,double *x,double *b,long *adr,long *ci,long n)
/*
  funkce pocita hodnotu energetickeho funkcionalu
  E = x A x - x b

  a - matice soustavy ulozena jako symmetric compressed rows
  x - vektor posunu
  b - vektor prave strany
  adr - pole adres prvnich prvku v radcich
  ci - pole sloupcovych indexu
  n - pocet neznamych

  7.4.1998
*/
{
  double e,*p;

  p = aloc1 (n);

  mv_scr (a,x,p,adr,ci,n);
  e = 0.5*ss(x,p,n) - ss(x,b,n);

  PGFEM_free (p);
  return e;
}


void cg_scr (double *a,double *x,double *y,long *adr,long *ci,long n,long ni,double err,long *ani,double *ares,double limit,long iv)
/*
  funkce resi soustavu linearnich algebraickych rovnic
  metodou sdruzenych gradientu (podle Axelssona)

  matice soustavy je ulozena v compresovane forme
  ukladaji se pouze nenulove prvky dolni trojuhelniku vcetne
  diagonalnich, ukladani je provedeno po radcich

  a - matice soustavy
  x - vektor reseni
  y - vektor prave strany
  adr - pole adres prvnich prvku v poli ci
  ci - pole sloupcovych indexu
  n - pocet neznamych
  ni - maximalni pocet iteraci
  err - maximalni pripustna chyba
  ani - pocet skutecne provedenych operaci
  ares - norma vektoru rezidui
  limit - konstanta pro testovani na nulu
  iv - prepinac pocatecni aproximace
  iv=0 - pocatecni aproximace je nulovy vektor
  iv=1 - pocatecni aproximace se do funkce posila zvnejsku

  30.4.1997
*/
{
  long i,j;
  double nom,denom,nory,alpha,beta,ener,res{};
  double *d,*r,*p;

  /* double  en; */
  /*  FILE *str; */
  /*  str = fopen ("../out/kuk","w"); */

  d = aloc1 (n);
  r = aloc1 (n);
  p = aloc1 (n);

  /**********************************/
  /*  nastaveni pocatecnich hodnot  */
  /**********************************/
  if (iv==0){
    for (i=0;i<n;i++)
      x[i]=0.0;
  }
  mv_scr (a,x,p,adr,ci,n);
  nory=0.0;  nom=0.0;
  for (i=0;i<n;i++){
    nory+=y[i]*y[i];
    r[i]=p[i]-y[i];
    nom+=r[i]*r[i];
    d[i]=-1.0*r[i];
  }

  if (nory<limit){
    PGFEM_printf ("\n\n norma vektoru prave strany v metode sdruzenych gradientu");
    PGFEM_printf ("\n je mensi nez %e",limit);
    *ares=nory;  *ani=0;
    return;
  }

  /*************/
  /*  iterace  */
  /*************/
  for (i=0;i<ni;i++){

    /*  vypocet nove alphy  */
    mv_scr (a,d,p,adr,ci,n);

    denom = ss (d,p,n);
    if (fabs(denom)<limit){
      PGFEM_printf ("\n V metode sdruzenych gradientu je nulovy jmenovatel u alfa");
      abort ();
    }

    alpha = nom/denom;

    /*  nova aproximace x a r */
    ener=0.0;  res=0.0;
    for (j=0;j<n;j++){
      ener+=alpha*d[j]*alpha*d[j];
      x[j]+=alpha*d[j];
      r[j]+=alpha*p[j];
      res+=r[j]*r[j];
    }

    /* en = energie_scr (a,x,y,adr,ci,n); */

    PGFEM_printf ("Iteration [%ld], Relative residual error |r|/|b| %e\n",i,res/nory);
    /*    PGFEM_printf ("%ld   %e   %e\n",i,res/nory,ener/nory); */

    /* PGFEM_printf ("iterace   %ld   %e   %e   %e\n",i,res/nory,ener/nory,en); */
    /* PGFEM_fprintf (str,"iterace   %ld   %le   %le   %le\n",i,res/nory,ener/nory,en); */

    denom=nom;

    /*  vypocet beta  */
    nom = ss (r,r,n);

    if (res/nory < err)  break;
    /*   if (ener/nory<err)  break; */
    /*    if (fabs(nom) < limit)  break; */

    beta = nom/denom;

    /*  vypocet noveho d */
    for (j=0;j<n;j++){
      d[j]=beta*d[j]-r[j];
    }
  }/* end i<ni */

  *ani=i;  *ares=res/nory;

  /*  fclose (str); */

  PGFEM_free (p);  PGFEM_free (r);  PGFEM_free (d);
}
