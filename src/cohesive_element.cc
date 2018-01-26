#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "cohesive_element.h"
#include "PGFEM_mpi.h"
#include "allocation.h"
#include "cohesive_element_utils.h"
#include "enumerations.h"
#include "get_dof_ids_on_elem.h"
#include "index_macros.h"
#include "matice.h"
#include "utils.h"
#include <cmath>

static void update_state_variables_co (COEL *cel, /* ptr to single el */
                                       double *x,
                                       double *y,
                                       double *z,
                                       double *r_e);

/*** MAIN ROUTINES ***/
void destroy_coel(COEL* coel, long nce)
{
  for(long i=0; i<nce; i++){
    free(coel[i].nod);
    free(coel[i].e1);
    free(coel[i].e2);
    free(coel[i].n);
    free(coel[i].x);
    free(coel[i].y);
    free(coel[i].z);
    free(coel[i].Xmax);
    free(coel[i].tmax);
    free(coel[i].Xi);
    free(coel[i].ti);

    int ip = int_pointC(coel[i].toe/2);
    for(int j=0; j<ip; j++){
      free(coel[i].vars[j]);
    }
    free(coel[i].vars);
  }
  free(coel);
}

void reset_coel_props(const cohesive_props *co_props,
                      COEL *p_coel)
{
  p_coel->props = &co_props[p_coel->mat];
}

COEL* read_cohe_elem (FILE *in1,
                      long ncom,
                      long ndofn,
                      long nn,
                      Node *node,
                      long *NCE,
                      double **comat,
                      Ensight *ensight,
                      long gr2,
                      int myrank,
                      const cohesive_props *co_props)
/*

 */
{
  long i,j,k,l,nce,nne,ndofe,*nod,II,*Enodes;
  double *r_e,*x,*y,*z,*e2h;
  COEL *coel;

  /*
    We have currently two cohesive laws:

    Sc : Xc : b :   -> Needleman
    Sc : Xc : b : k -> Matous-Arciniega
    1 < k < 5.04 - parameter shifts the Xc to the right,
    but Gc is preserved
  */

  for (i=0;i<ncom;i++) {
    CHECK_SCANF (in1,"%lf %lf %lf %lf\n",&comat[i][0],&comat[i][1],
                 &comat[i][2],&comat[i][3]); /* Sc : Xc : b : k */

    if (comat[i][3] < 1.0 || comat[i][3] >= 3.00) {
      /* Why is this different from the explination of the law? */
      PGFEM_printf ("Incorrect input for cohesive element TYPE 1\n");
      /* PGFEM_Abort(); */
    }

  }

  CHECK_SCANF (in1,"%ld\n",&nce);

  if (nce == 0) coel = PGFEM_calloc (COEL, 1);
  else coel = PGFEM_calloc (COEL, nce);

  for (i=0;i<nce;i++) CHECK_SCANF (in1,"%ld\n",&coel[i].toe);

  for (i=0;i<nce;i++){
    nne = coel[i].toe/2;

    /* Number of intergration points */
    II = int_pointC (nne);

    coel[i].nod = PGFEM_calloc (long, coel[i].toe);
    coel[i].e1  = PGFEM_calloc (double, 3);
    coel[i].e2  = PGFEM_calloc (double, 3);
    coel[i].n   = PGFEM_calloc (double, 3);
    coel[i].ti  = PGFEM_calloc (double, 3);
    coel[i].Xi  = PGFEM_calloc (double, 3);

    coel[i].x = PGFEM_calloc (double, nne);
    coel[i].y = PGFEM_calloc (double, nne);
    coel[i].z = PGFEM_calloc (double, nne);

    coel[i].Xmax = PGFEM_calloc (double, II);
    coel[i].tmax = PGFEM_calloc (double, II);

    /* state variables */
    coel[i].vars = PGFEM_calloc(double*, II);

  }/* end i<nce */

  /*
    0   30 37 33 31 36 32      0      0       0
    id  connectivity table   TYPE MATERIAL PROPERTY
  */

  for (i=0;i<nce;i++){
    CHECK_SCANF (in1,"%ld",&j);
    for (l=0;l<coel[j].toe;l++){
      CHECK_SCANF (in1,"%ld",&coel[j].nod[l]);
    }
    CHECK_SCANF (in1,"%ld %ld %ld",&coel[j].typ,&coel[j].mat,&coel[j].pr);

    /* Cohesive elements are always integrated in reference configuration */
    coel[j].Jjn = 1.0;

    /* Set the element properties */
    reset_coel_props(co_props,coel + j);

    /* set internal state variables */
    switch(coel[j].props->type){
     case CO_MOD_MS:
      coel[j].nvar = 5;
      /* state vars are: max_traction, max_jump, t_1,t_2,t_3 */
      break;
     case CO_MOD_NEEDLEMAN:
       {
         coel[j].nvar = 2;
         break;
       }
     case CO_MOD_VDW: case CO_MOD_LJ:
       {
         coel[j].nvar = 0;
         break;
       }

       /*=== NOT IMPLEMENTED ===*/
     case CO_MOD_MATOUS_ARCINIEGA:
     case CO_MOD_PPR:
     default:
      PGFEM_printerr("Cohesive type not yet implemented! %s:%s\n",
                     __FILE__,__func__);
      PGFEM_Abort();
    }

    II = int_pointC(coel[j].toe/2);
    for(int ip=0; ip<II; ip++){
      if(coel[j].nvar > 0){
        coel[j].vars[ip] = PGFEM_calloc(double, coel[j].nvar);
      } else {
        coel[j].vars[ip] = NULL;
      }
    }
  }

  /* ENSIGHT NODES */
  if (gr2 == VIS_ENSIGHT || gr2 == VIS_VTK) {

    Enodes = aloc1l (nn);
    for (i=0;i<nn;i++) Enodes[i] = -1;

    k = 0;
    for (i=0;i<nce;i++){
      for (j=0;j<coel[i].toe;j++) {
        /* mark nodes on - surface for Ensight (-) */
        if (j < coel[i].toe/2 && Enodes[coel[i].nod[j]] == -1) {
          Enodes[coel[i].nod[j]] = k;
          k++;
        }
      }
    }

    ensight->ncn = 0;
    for (i=0;i<nn;i++) {
      if (Enodes[i] != -1)
        ensight->ncn++;
      Enodes[i] = -1;
    }
    if(ensight->ncn > 0){
      ensight->Sm = new long[ensight->ncn]{};
      ensight->Sp = new long[ensight->ncn]{};
    } else {
      ensight->Sm = ensight->Sp = nullptr;
    }
    ensight->No = new long[nn];

    k = 0;
    for (i=0;i<nce;i++){
      for (j=0;j<coel[i].toe;j++) {
        /* mark nodes on - surface for Ensight (-) */
        if (j < coel[i].toe/2 && Enodes[coel[i].nod[j]] == -1) {
          Enodes[coel[i].nod[j]] = k;
          ensight->Sm[k] = coel[i].nod[j];
          ensight->Sp[k] = coel[i].nod[j+coel[i].toe/2];
          k++;
        }
      }
    }

    for (i=0;i<nn;i++) ensight->No[i] = Enodes[i];

    dealoc1l (Enodes);
  }/* end ensight */

  e2h = aloc1 (3);

  /* Compute base vectors for cohesive elements */
  for (i=0;i<nce;i++){

    nne = coel[i].toe/2;
    ndofe = coel[i].toe*ndofn;

    nod = aloc1l (coel[i].toe);
    r_e = aloc1 (ndofe);
    x = aloc1 (coel[i].toe);
    y = aloc1 (coel[i].toe);
    z = aloc1 (coel[i].toe);

    for (j=0;j<coel[i].toe;j++) nod[j] = coel[i].nod[j];
    nodecoord_updated (coel[i].toe,nod,node,x,y,z);

    /* Mean mapping */
    mean_map (nne,x,y,z,r_e,coel[i].x,coel[i].y,coel[i].z);

    /* Transformation base vectors :: normal vector */
    base_vec (nne,0.0,0.0,coel[i].x,coel[i].y,coel[i].z,
              coel[i].e1,coel[i].e2,e2h,coel[i].n,myrank);

    dealoc1l (nod);
    dealoc1 (r_e);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
  }/* end i < nce */

  *NCE = nce; dealoc1 (e2h);

  if (coel == NULL){PGFEM_printf ("\n Memory is full.\n");abort ();}

  return (coel);
}

void stiff_mat_coh (long ii,
                    long ndofn,
                    long nne,
                    long *nod,
                    double *x,
                    double *y,
                    double *z,
                    COEL *coel,
                    double *r_e,
                    double *Kch,
                    double nor_min,
                    EPS *eps,
                    long FNR,
                    double lm,
                    double *fe,
                    int myrank)
/*
*** NOTE: integration is on A0, even for updated Lagrangian
*** formulation!
*/
{
  long ip,i,j,II,JJ,M,N,P,R,U,ndofe;
  double ****Kt,****Kn,*gk,*ge,*w,ksi{},eta{},ai{},aj{},*Nf,*N_x,*N_y,
                                                        *xl,*yl,*zl,**TX,***AA,**TN,
                                                        ***Nxb,*xb,*yb,*zb,*Xi,*e1,
                                                        *e2,*e2h,*n;
  double **STIFF,J=0.0,Jjn;

  static const int ndn = 3;
  COEL *const cel = &coel[ii];
  const cohesive_props *props = cel->props;
  ndofe = ndofn*coel[ii].toe;
  // @todo Commented as dead code. @cp should review. LD
  // double bet = coel[ii].b;
  Jjn = coel[ii].Jjn;

  Kt = aloc4 (ndn,ndn,nne,nne);
  Kn = aloc4 (ndn,ndn,nne,nne);
  Nf = aloc1 (nne);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  xl = aloc1 (nne);
  yl = aloc1 (nne);

  zl = aloc1 (nne);
  TX = aloc2 (ndn,ndn);
  AA = aloc3 (ndn,ndn,nne);
  TN = aloc2 (ndn,ndn);
  Nxb = aloc3 (ndn,ndn,nne);
  xb = aloc1 (nne);
  yb = aloc1 (nne);
  zb = aloc1 (nne);
  Xi = aloc1 (ndn);
  e1 = aloc1 (3);
  e2 = aloc1 (3);
  e2h = aloc1 (3);
  n = aloc1 (3);
  STIFF = aloc2 (ndofe,ndofe);

  /* Number of intergration points */
  II = int_pointC (nne);

  /* int. point and weights */
  gk = aloc1 (II);
  ge = aloc1 (II);
  w = aloc1 (II);

  /* Transform coordinates (at time n) from G->L in E1;E2 */
  tran_coord (nne,coel[ii].x,coel[ii].y,coel[ii].z,coel[ii].e1,
              coel[ii].e2,coel[ii].n,xl,yl,zl,1);

  /* Mean mapping */
  mean_map (nne,x,y,z,r_e,xb,yb,zb);

  /* Numerical integration */
  integrateC (nne,&II,&JJ,gk,ge,w);


  /* allocate arrays for tractions and tangents */
  double *traction = PGFEM_calloc(double, ndn);
  double *mat_tan =  PGFEM_calloc(double, ndn*ndn);
  double *geo_tan =  PGFEM_calloc(double, ndn*ndn);

  ip = 0;
  for(i=0;i<II;i++){
    for(j=0;j<JJ;j++){

      /* Coordinates ksi, eta */
      if (nne == 3){
        ksi = gk[j];
        eta = ge[j];
        ai = w[j];
        aj = 1.00;
      }
      if (nne == 4){
        ksi = gk[i];
        eta = gk[j];
        ai = w[i];
        aj = w[j];
      }

      /* Transformation base vectors :: normal vector */
      base_vec (nne,ksi,eta,xb,yb,zb,e1,e2,e2h,n,myrank);

      /* Derivatives of the shape functions */
      J = dN3_xy (ksi,eta,nne,xl,yl,zl,N_x,N_y);

      /* Shape functions */
      shape_2DC (nne,ksi,eta,Nf);

      /* Get diplacement jump */
      get_jump (nne,x,y,z,r_e,Nf,Xi);

      /* compute the traction and tangents */
      int err = 0;
      err += props->get_traction(traction,Xi,n,props->props,cel->vars[ip]);
      err += props->get_tangents(mat_tan,geo_tan,Xi,n,props->props,
                                 cel->vars[ip]);

      /* Remove damaged elements */
      if(err >= 100){ /* Killed element */
        ip++;
        continue;
      } else if(err != 0){
        PGFEM_printerr("Error in %s:%s(%d)\n",
                       __FILE__,__func__,__LINE__);
      }


      /* Get dN/dxb */
      dN_dxb (nne,ksi,eta,e1,e2h,n,Nxb);
      for (M=0;M<ndn;M++){
        for (N=0;N<ndn;N++){
          for (P=0;P<nne;P++){

            AA[M][N][P] = 0.0;
            for (U=0;U<ndn;U++){
              AA[M][N][P] += geo_tan[idx_2(M,U)]*Nxb[U][N][P];
            }

            for (R=0;R<nne;R++){
              Kt[M][N][P][R] += (ai*aj*J*Jjn*mat_tan[idx_2(M,N)]
                                 *Nf[P]*Nf[R]);
            }
          }
        }
      }/* end M < ndn */

      for (M=0;M<ndn;M++){
        for (N=0;N<ndn;N++){
          for (P=0;P<nne;P++){
            for (R=0;R<nne;R++){
              Kn[M][N][P][R] += (ai*aj*J * 1./2.*Jjn*AA[M][N][R]*Nf[P]);
            }
          }
        }
      }/* M < ndn */

      ip++;
    } /*end j < JJ */
  } /*end i < II */

  /* Composition I */
  for (M=0;M<ndn;M++){
    for (N=0;N<ndn;N++){
      for (P=0;P<nne;P++){
        for (R=0;R<nne;R++){/* aebr */
          STIFF[P*ndofn+M][R*ndofn+N] =
          Kt[M][N][P][R] - Kn[M][N][P][R]; /* -- */

          STIFF[nne*ndofn + P*ndofn+M][nne*ndofn + R*ndofn+N] =
          Kt[M][N][P][R] + Kn[M][N][P][R]; /* ++ */

          STIFF[P*ndofn+M][nne*ndofn + R*ndofn+N] =
          -Kt[M][N][P][R] - Kn[M][N][P][R]; /* -+ */

          STIFF[nne*ndofn + P*ndofn+M][R*ndofn+N] =
          -Kt[M][N][P][R] + Kn[M][N][P][R]; /* +- */
        }
      }
    }
  }

  /* Composition II */
  for (M=0;M<ndofe;M++)
    for (N=0;N<ndofe;N++)
      Kch[M*ndofe+N] = STIFF[M][N];

  dealoc4 (Kt,ndn,ndn,nne);
  dealoc4 (Kn,ndn,ndn,nne);
  dealoc1 (Nf);
  dealoc1 (N_x);
  dealoc1 (N_y);
  dealoc1 (xl);
  dealoc1 (yl);
  dealoc1 (zl);
  dealoc2 (TX,ndn);
  dealoc3 (AA,ndn,ndn);
  dealoc2 (TN,ndn);
  dealoc3 (Nxb,ndn,ndn);
  dealoc1 (xb);
  dealoc1 (yb);
  dealoc1 (zb);
  dealoc1 (Xi);
  dealoc1 (e1);
  dealoc1 (e2);
  dealoc1 (e2h);
  dealoc1 (n);
  dealoc2 (STIFF,ndofe);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1(w);

  free(traction);
  free(mat_tan);
  free(geo_tan);
}

void resid_co_elem (long ii,
                    long ndofn,
                    long nne,
                    long *nod,
                    double *x,
                    double *y,
                    double *z,
                    COEL *coel,
                    double *r_e,
                    double *fe,
                    double nor_min,
                    int myrank)
/*
*** NOTE: integration is on A0, even for updated Lagrangian
*** formulation!
*/
{
  long ip,i,j,II,JJ,M,N;
  double *gk,*ge,*w,ksi{},eta{},ai{},aj{},*T,**Rc,*Nf,*N_x,*N_y,
                                          *xl,*yl,*zl,*xb,*yb,*zb,*Xi,*e1,*e2,*e2h,
                                          *n,J,Jjn;

  static const int ndn = 3;
  COEL *const cel = &coel[ii];
  const cohesive_props *props = cel->props;
  // @todo Commented as dead code. @cp should review. LD
  // long ndofe = ndofn*coel[ii].toe;
  // double bet = coel[ii].b;
  Jjn = coel[ii].Jjn;

  T = aloc1 (ndn);
  Rc = aloc2 (ndn,nne);
  Nf = aloc1 (nne);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  xl = aloc1 (nne);
  yl = aloc1 (nne);
  zl = aloc1 (nne);
  xb = aloc1 (nne);
  yb = aloc1 (nne);
  zb = aloc1 (nne);
  Xi = aloc1 (ndn);
  e1 = aloc1 (3);
  e2 = aloc1 (3);
  e2h = aloc1 (3);
  n = aloc1 (3);

  /* Number of intergration points */
  II = int_pointC (nne);

  /* int. point and weights */
  gk = aloc1 (II);
  ge = aloc1 (II);
  w = aloc1 (II);

  /* Transform coordinates (at time n) from G->L in E1;E2 */
  tran_coord (nne,coel[ii].x,coel[ii].y,coel[ii].z,coel[ii].e1,
              coel[ii].e2,coel[ii].n,xl,yl,zl,1);

  /* element area */
  /* AA = area (nne,xl,yl); */

  /* Mean mapping */
  mean_map (nne,x,y,z,r_e,xb,yb,zb);

  /* Numerical integration */
  integrateC (nne,&II,&JJ,gk,ge,w);

  /* allocate arrays for tractions and tangents */
  double *traction = PGFEM_calloc(double, ndn);

  ip = 0;
  for(i=0;i<II;i++){
    for(j=0;j<JJ;j++){

      /* Coordinates ksi, eta */
      if (nne == 3){
        ksi = gk[j];
        eta = ge[j];
        ai = w[j];
        aj = 1.00;
      }
      if (nne == 4){
        ksi = gk[i];
        eta = gk[j];
        ai = w[i];
        aj = w[j];
      }

      /* Transformation base vectors :: normal vector */
      base_vec (nne,ksi,eta,xb,yb,zb,e1,e2,e2h,n,myrank);

      /* Derivatives of the shape functions */
      J = dN3_xy (ksi,eta,nne,xl,yl,zl,N_x,N_y);

      /* Shape functions */
      shape_2DC (nne,ksi,eta,Nf);

      /* Get diplacement jump */
      get_jump (nne,x,y,z,r_e,Nf,Xi);

      /* compute the traction and tangents */
      int err = 0;
      err += props->get_traction(traction,Xi,n,props->props,cel->vars[ip]);

      /* Remove damaged elements */
      if(err >= 100){ /* Killed element */
        ip++;
        continue;
      } else if(err != 0){
        PGFEM_printerr("Error in %s:%s(%d)\n",
                       __FILE__,__func__,__LINE__);
      }

      for (M=0;M<ndn;M++){
        for (N=0;N<nne;N++){
          Rc[M][N] += ai*aj*J * Jjn*traction[M]*Nf[N];
        }
      }

      ip++;
    }/* end j < JJ */
  }/* end i < II */

  /* Residuals composition */
  for (N=0;N<nne;N++){
    for (M=0;M<ndn;M++){
      fe[N*ndofn+M] = -Rc[M][N]; /* - */
      fe[ndofn*nne + N*ndofn + M] = Rc[M][N]; /* + */
    }
  }

  free(traction);
  dealoc1 (T);
  dealoc2 (Rc,ndn);
  dealoc1 (Nf);
  dealoc1 (N_x);
  dealoc1 (N_y);
  dealoc1 (xl);
  dealoc1 (yl);
  dealoc1 (zl);
  dealoc1 (xb);
  dealoc1 (yb);
  dealoc1 (zb);
  dealoc1 (Xi);
  dealoc1 (e1);
  dealoc1 (e2);
  dealoc1 (e2h);
  dealoc1 (n);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (w);
}

int increment_cohesive_elements(const int nce,
                                COEL *coel,
                                double *pores,
                                const Node *node,
                                const SUPP sup,
                                const double *d_r,
                                const int mp_id)
{
  int err = 0;
  const int ndofc = 3;
  *pores = 0.0;

  for (int i=0;i<nce;i++){
    COEL *cel = &coel[i];
    const int nne_t = cel->toe;
    const int ndofe = nne_t*ndofc;

    double *r_e = aloc1 (ndofe);
    double *x = aloc1 (nne_t);
    double *y = aloc1 (nne_t);
    double *z = aloc1 (nne_t);
    long *cn = aloc1l (ndofe);

    /* coordinates */
    const long *nod = cel->nod;
    nodecoord_updated (nne_t,nod,node,x,y,z);

    /* Id numbers */
    get_dof_ids_on_elem_nodes(0,nne_t,ndofc,nod,node,cn,mp_id);

    /* deformation on element */
    def_elem (cn,ndofe,d_r,NULL,node,r_e,sup,0);

    /* Update Cohesive elements */
    update_state_variables_co (cel,x,y,z,r_e);

    /* Voids volume */
    *pores += cel->vo;

    dealoc1 (r_e);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1l (cn);
  }/* end i < nce */

  return err;
}/* increment_cohesive_elements */

static void update_state_variables_co (COEL *cel, /* ptr to single el */
                                       double *x,
                                       double *y,
                                       double *z,
                                       double *r_e)
/*
*** NOTE: integration is on A0, even for updated Lagrangian
*** formulation!
*/
{
  static const int ndn = 3;
  const int nne = cel->toe/2;

  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);

  double bet = cel->b;

  /* Allocate */
  double *Nf = aloc1 (nne);
  double *xb = aloc1 (nne);
  double *yb = aloc1 (nne);
  double *zb = aloc1 (nne);
  double *e1 = aloc1 (3);
  double *e2 = aloc1 (3);
  double *e2h = aloc1 (3);
  double *n = aloc1 (3);

  double *Xi = aloc1 (ndn);
  double *xl = aloc1 (nne);
  double *yl = aloc1 (nne);
  double *zl = aloc1 (nne);
  double *N_x = aloc1 (nne);
  double *N_y = aloc1 (nne);
  double *vec1 = aloc1 (3);
  double *vec2 = aloc1 (3);
  double *vec3 = aloc1 (3);

  /* allocate traction */
  double *traction = PGFEM_calloc(double, ndn);

  /* Number of intergration points */
  long II = int_pointC (nne);
  long JJ = 0;

  double *gk = aloc1 (II);
  double *ge = aloc1 (II);
  double *w = aloc1 (II);

  /* Transform coordinates (at time n) from G->L in E1;E2 */
  tran_coord (nne,cel->x,cel->y,cel->z,cel->e1,
              cel->e2,cel->n,xl,yl,zl,1);

  /* Area of cohesive element */
  double aa = area (nne,xl,yl);

  /* Mean mapping */
  mean_map (nne,x,y,z,r_e,xb,yb,zb);

  /* Null fields */
  cel->txi = cel->Xxi = cel->ts = 0.0;
  cel->tn = cel->Xn = cel->vo = 0.0;
  cel->Xn = cel->Xs = 0.0;
  for (int M=0;M<ndn;M++){
    cel->ti[M] = cel->Xi[M] = 0.0;
  }

  /* Numerical integration */
  integrateC (nne,&II,&JJ,gk,ge,w);

  int ip = 0;
  for(int i=0;i<II;i++){
    for(int j=0;j<JJ;j++){
      double ksi, eta, ai, aj;

      /* Coordinates ksi, eta */
      switch(nne){
       case 3:
        ksi = gk[j];
        eta = ge[j];
        ai = w[j];
        aj = 1.00;
        break;
       case 4:
        ksi = gk[i];
        eta = gk[j];
        ai = w[i];
        aj = w[j];
        break;
       default:
        PGFEM_printerr("[%d]ERROR: incorrect number of nodes! %s:%s:%d\n",
                       err_rank,__func__,__FILE__,__LINE__);
        PGFEM_Abort();
      }

      /* Transformation base vectors :: normal vector */
      base_vec (nne,ksi,eta,xb,yb,zb,e1,e2,e2h,n,err_rank);

      /* Derivatives of the shape functions */
      double J = dN3_xy (ksi,eta,nne,xl,yl,zl,N_x,N_y);

      /* Shape functions */
      shape_2DC (nne,ksi,eta,Nf);

      /* Get diplacement jump */
      get_jump (nne,x,y,z,r_e,Nf,Xi);

      /* Normal opening displacement */
      double Xn = ss (Xi,n,3);

      /* Get effective opening displacement */
      double Xxi = get_Xxi (bet,Xi,n);

      /* get traction */
      cel->props->get_traction(traction,Xi,n,cel->props->props,
                               cel->vars[ip]);

      /* Average cohesive element fields */
      double txi = 0.0;
      for (int M=0;M<ndn;M++){
        cel->ti[M] += ai*aj*J/aa *traction[M];
        cel->Xi[M] += ai*aj*J/aa * Xi[M];
        txi += traction[M]*traction[M];
      }

      /* Average traction and effective opening displacement */
      txi = sqrt(txi);
      cel->txi += ai*aj*J/aa*txi;
      cel->Xxi += ai*aj*J/aa*Xxi;

      /* Normal and Shear tractions || Normal opening */
      for (int M=0;M<ndn;M++){
        /* vec1[M] = traction[M]*n[M]; */
        /* vec2[M] = EE1*bet*bet*(Xi[M] - Xn*n[M]); */
        vec3[M] = Xi[M] - Xn*n[M]; /* shear direction */
      }
      double Xs = sqrt (ss(vec3,vec3,3)); /* shear magnitude */
      double tn = sqrt (pow(ss(traction,n,3),2));
      double ts = 0.0;
      if(Xs > 0){
        ts = sqrt (pow(ss(traction,vec3,3)/Xs,2));
      }
      cel->tn += ai*aj*J/aa*tn;
      cel->ts += ai*aj*J/aa*ts;
      cel->Xn += ai*aj*J/aa*Xn;
      cel->Xs += ai*aj*J/aa*Xs;

      /* For loading update state variables */

      if (Xxi >= cel->Xmax[ip] && Xn > 0.0){
        cel->Xmax[ip] = Xxi;
        cel->tmax[ip] = txi;
      }

      switch(cel->props->type){
       case CO_MOD_MS:
        /* do nothing, multiscale state variables are updated
           elsewhere */
        break;
       default:
        if(cel->nvar > 0){
          if (Xxi >= cel->Xmax[ip] && Xn > 0.0){
            cel->vars[ip][0] = Xxi;
            cel->vars[ip][1] = txi;
          }
        }
        break;
      }

      /* Void concentration */
      if (Xn > 0.0 && cel->Xmax[ip] > cel->Xc)
        cel->vo += ai*aj*J* Xn;

      ip++;
    }/* end j < JJ */
  }/* end i < II */

  /* Positive or negative normal traction */
  cel->tn *= cel->Xn/fabs(cel->Xn);

  free(traction);
  dealoc1 (Nf);
  dealoc1 (xb);
  dealoc1 (yb);
  dealoc1 (zb);
  dealoc1 (e1);
  dealoc1 (e2);
  dealoc1 (e2h);
  dealoc1 (n);
  dealoc1 (Xi);
  dealoc1 (xl);

  dealoc1 (yl);
  dealoc1 (zl);
  dealoc1 (N_x);
  dealoc1 (N_y);
  dealoc1 (vec1);
  dealoc1 (vec2);
  dealoc1 (vec3);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (w);
}

size_t coel_list_get_state_length_bytes(const int nce,
                                        const COEL *coel)
{
  size_t len_b = 0;
  for(int i = 0; i < nce; i++){
    int nip = int_pointC(coel[i].toe/2);
    /* length in bytes/elem = n_int_points*nvars*sizeof(vars) */
    len_b += nip*coel[i].nvar*sizeof(**(coel[i].vars));
  }
  return len_b;
}

void coel_list_pack_state(const int nce,
                          const COEL *coel,
                          char *buffer,
                          size_t *buf_pos)
{
  if(nce <= 0 ) return;
  const size_t len_vars = sizeof(**(coel[0].vars));
  for(int i = 0; i < nce; i++){
    int nip = int_pointC(coel[i].toe/2);
    for(int j = 0; j < nip; j++){
      pack_data(coel[i].vars[j],buffer,buf_pos,coel[i].nvar,len_vars);
    }
  }
}

void coel_list_unpack_state(const int nce,
                            COEL *coel,
                            const cohesive_props *co_props,
                            const char *buffer,
                            size_t *buf_pos)
{
  if(nce <= 0 ) return;
  const size_t len_vars = sizeof(**(coel[0].vars));
  for(int i = 0; i < nce; i++){
    int nip = int_pointC(coel[i].toe/2);
    for(int j = 0; j < nip; j++){
      unpack_data(buffer,coel[i].vars[j],buf_pos,coel[i].nvar,len_vars);
    }
    reset_coel_props(co_props,coel+i);
  }
}
