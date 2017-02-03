/* HEADER */

#include "elem3d.h"
#include "enumerations.h"
#include "incl.h"
#include "quadrature_rules.h"
#include "get_dof_ids_on_elem.h"
#include "homogen.h"
#include "localizat.h"
#include "matice.h"
#include "utils.h"
#include "resice.h"

#ifndef NO_BUBBLE
#define NO_BUBBLE 0
#endif

void int_point (const long nne,
		long *II)
     /*
       
     */
{
  if (nne == 4)  *II = 1; /* *II = 1; */ /* Linear tetra */
  if (nne == 5)  *II = 4; /* P1+B/P1 */
  if (nne == 10) *II = 4; /* Quadratic tetra */
  if (nne == 8)  *II = 8; /* Linear Hexa */
}

void integrate (const long nne,
		long *II,
		long *JJ,
		long *KK,
		double *gk,
		double *ge,
		double *gz,
		double *w)
     /*
       Determine the integration points (gk, ge, gz)
       Determine the integration weights (w)
       Number of summations in each direction (II, JJ, KK)
     */
{
  if (nne == 4)  {
    int_tetra_1 (gk,ge,gz,w);
    /* int_tetra_4 (gk,ge,gz,w); */
    *II = 1;
    *JJ = 1;
    *KK = 1;
    /* *KK = 4; */
  } /* Linear tetra */
  
  if (nne == 5){
    //int_tetra_24 (gk,ge,gz,w);
    int_tetra_4 (gk,ge,gz,w);
    *II = 1;
    *JJ = 1;
    //*KK = 24;
    *KK = 4; 
  } /* P1 + B/P1 || bubble function is 4th degree polynomial!!! */

  if (nne == 10) {
    int_tetra_4 (gk,ge,gz,w);
    *II = 1;
    *JJ = 1;
    *KK = 4;
  } /* Quadratic tetra */

  if (nne == 8)  {
    intpoints_2 (gk,w);
    *II = 2;
    *JJ = 2;
    *KK = 2;
  } /* Linear Hexa */
}



void int_tria_1 (double *gs,double *w)
/*
  Function return integration points for degenerated triangle
*/
{
  *(w+0) = 3.0;
  *(gs+0)  = -0.333333333333333;
}

void intpoints_2 (double *gs,double *w)
     /*
       Function return Gauss integration points
     */
{
  *(w+0) = 1.0;
  *(w+1) = 1.0;

  *(gs+0) = -0.577350269189626;
  *(gs+1) =  0.577350269189626;
}

void intpoints_3 (double *gs,double *w)
     /*
       Function return Gauss integration points
     */
{
  *(w+0) = 0.555555555555555;
  *(w+1) = 0.888888888888888;
  *(w+2) = 0.555555555555555;
  
  *(gs+0) = -0.774596669241483;
  *(gs+1) =  0.0;
  *(gs+2) =  0.774596669241483;
}

void int_tetra_1 (double *gk,double *ge,double *gz,double *w)
     /*
       Returns integration points and integration weights
     */
{
  w[0] = 1./6.;  gk[0] = 1./4.;  ge[0] = 1./4.;  gz[0] = 1./4.;
}

void int_tetra_4 (double *gk,double *ge,double *gz,double *w)
     /*
       
     */
{
  w[0] = (1./4.)/6.;
  gk[0] = 0.58541020;
  ge[0] = 0.13819660;
  gz[0] = 0.13819660;

  w[1] = (1./4.)/6.;
  gk[1] = 0.13819660;
  ge[1] = 0.58541020;
  gz[1] = 0.13819660;

  w[2] = (1./4.)/6.;
  gk[2] = 0.13819660;
  ge[2] = 0.13819660;
  gz[2] = 0.58541020;

  w[3] = (1./4.)/6.;
  gk[3] = 0.13819660;
  ge[3] = 0.13819660;
  gz[3] = 0.13819660;

  
}

void int_tetra_5 (double *gk,double *ge,double *gz,double *w)
     /*
       
     */
{
  w[0] = (-4./5.)/6.;  gk[0] = 1./4.;  ge[0] = 1./4.;  gz[0] = 1./4.;
  w[1] = (9./20.)/6.;  gk[1] = 1./2.;  ge[1] = 1./6.;  gz[1] = 1./6.;
  w[2] = (9./20.)/6.;  gk[2] = 1./6.;  ge[2] = 1./2.;  gz[2] = 1./6.;
  w[3] = (9./20.)/6.;  gk[3] = 1./6.;  ge[3] = 1./6.;  gz[3] = 1./2.;
  w[4] = (9./20.)/6.;  gk[4] = 1./6.;  ge[4] = 1./6.;  gz[4] = 1./6.;
}

void int_tetra_11 (double *gk,
		   double *ge,
		   double *gz,
		   double *w)
{

  /* 11 point formula for exact integration of 4th degree polunomial
     on tetrahedron */

  /* Formula from Keast 1985. Note that the factor of 1/6 is already
     incorperated into the weights */

  w[0] = -0.01315555556;  gk[0] = 1./4.;  ge[0] = 1./4.;  gz[0] = 1./4.;

  w[1] = 0.007622222222;  gk[1] = ge[1] = gz[1] = 0.07142857143;
  w[2] = 0.007622222222; gk[2] = 0.7857142857; ge[2] = gz[2] = 0.07142857143;
  w[3] = 0.007622222222; ge[3] = 0.7857142857; gk[3] = gz[3] = 0.07142857143;
  w[4] = 0.007622222222; gz[4] = 0.7857142857; gk[4] = ge[4] = 0.07142857143;

  w[5] = 0.02488888889; gk[5] = ge[5] = 0.3994035762; gz[5] = 0.1005964238;
  w[6] = 0.02488888889; gk[6] = 0.3994035762; ge[6] = gz[6] = 0.1005964238;
  w[7] = 0.02488888889; gz[7] = 0.3994035762; gk[7] = ge[7] = 0.1005964238;
  w[8] = 0.02488888889; ge[8] = gz[8] = 0.3994035762; gk[8] = 0.1005964238;
  w[9] = 0.02488888889; gk[9] = gz[9] = 0.3994035762; ge[9] = 0.1005964238;
  w[10] = 0.02488888889; ge[10] = 0.3994035762; gz[10]=gk[10] = 0.1005964238;

}

void int_tetra_24(double *gk,
		  double *ge,
		  double *gz,
		  double *w)
{
  /* 24 point formula for exact integration of 6th degree polunomial
     on tetrahedron */

  /* Formula from Keast 1985. Note that the factor of 1/6 is already
     incorperated into the weights */

  double a1, a2, a3, a4;

  /* block 1 */
  w[0] = w[1] = w[2] = w[3] = 0.665379170969464506e-2;
  a1 = a2 = a3 = 0.214602871259151684;
  a4 = 0.356191386222544953;

  gk[0] = a1;
  ge[0] = a2;
  gz[0] = a3;

  gk[1] = a4;
  ge[1] = a2;
  gz[1] = a3;

  gk[2] = a1;
  ge[2] = a4;
  gz[2] = a3;

  gk[3] = a1;
  ge[3] = a2;
  gz[3] = a4;

  /* block 2 */
  w[4] = w[5] = w[6] = w[7] = 0.167953517588677620e-2;
  a1 = a2 = a3 = 0.406739585346113397e-1;
  a4 =  0.877978124396165982;

  gk[4] = a1;
  ge[4] = a2;
  gz[4] = a3;
  
  gk[5] = a4;
  ge[5] = a2;
  gz[5] = a3;
  
  gk[6] = a1;
  ge[6] = a4;
  gz[6] = a3;

  gk[7] = a1;
  ge[7] = a2;
  gz[7] = a4;

  /* block 3*/
  w[8] = w[9] = w[10] = w[11] = 0.922619692394239843e-2;
  a1 = a2 = a3 = 0.322337890142275646;
  a4 = 0.329863295731730594e-1;

  gk[8] = a1;
  ge[8] = a2;
  gz[8] = a3;

  gk[9] = a4;
  ge[9] = a2;
  gz[9] = a3;

  gk[10] = a1;
  ge[10] = a4;
  gz[10] = a3;

  gk[11] = a1;
  ge[11] = a2;
  gz[11] = a4;

  /* block 4 */
  w[12] = w[13] = w[14] = w[15] = w[16]
    = w[17] = w[18] = w[19] = w[20] = w[21]
    = w[22] = w[23] = 0.803571428571428248e-2;
  a1 = a2 = 0.636610018750175299e-1;
  a3 = 0.269672331458315867;
  a4 = 0.603005664791649076;

  /* rotate through a1 a2 a3 */
  gk[12] = a1;
  ge[12] = a2;
  gz[12] = a3;

  gk[13] = a3;
  ge[13] = a1;
  gz[13] = a2;

  gk[14] = a2;
  ge[14] = a3;
  gz[14] = a1;

  /* replace a1 with a4 */
  gk[15] = a4;
  ge[15] = a2;
  gz[15] = a3;

  gk[16] = a3;
  ge[16] = a4;
  gz[16] = a2;

  gk[17] = a2;
  ge[17] = a3;
  gz[17] = a4;

  /* replace a2 with a4 */
  gk[18] = a1;
  ge[18] = a4;
  gz[18] = a3;

  gk[19] = a3;
  ge[19] = a1;
  gz[19] = a4;

  gk[20] = a4;
  ge[20] = a3;
  gz[20] = a1;

  /* replace a3 with a4 */
  gk[21] = a1;
  ge[21] = a2;
  gz[21] = a4;

  gk[22] = a4;
  ge[22] = a1;
  gz[22] = a2;

  gk[23] = a2;
  ge[23] = a4;
  gz[23] = a1;

}

void shape_tria (const double ksi,
		 const double eta,
		 const double zet,
		 double *N)
/*
       
     */
{
  N[0] = (2.*(1.-ksi-eta)-1.)*(1.-ksi-eta);
  N[1] = (2.*ksi-1.)*ksi;
  N[2] = (2.*eta-1.)*eta;
  N[3] = 4.*(1.-ksi-eta)*ksi;
  N[4] = 4.*ksi*eta;
  N[5] = 4.*(1.-ksi-eta)*eta;
}

void shape_func (const double ksi,
		 const double eta,
		 const double zet,
		 const long nne,
		 double *N)
     /*
       ksi -souradnice integracniho bodu (integration points)
       eta -souradnice integracniho bodu
       N - bazove funkce (basis functions)
     */
{
  if (nne == 1){ /* Constant */
    N[0] = 1.0;
  }
  if (nne == 4){/* Linear tetrahedron */
    N[0] = (1.-ksi-eta-zet);
    N[1] = ksi;
    N[2] = eta;
    N[3] = zet;
  }
  if(nne == 5){ /* linear tetra + bubble */
    N[0] = (1.-ksi-eta-zet);
    N[1] = ksi;
    N[2] = eta;
    N[3] = zet;
    if(NO_BUBBLE){
      N[4] = 0;
    } else {
      N[4] = 256.*ksi*eta*zet*(1.-ksi-eta-zet);
    }
    
    /* N[0] = 1 - eta - ksi - zet - 64.0*eta*ksi*(1 - eta - ksi - zet)*zet; */
    /* N[1] = ksi - 64.0*eta*ksi*(1 - eta - ksi - zet)*zet; */
    /* N[2] = eta - 64.0*eta*ksi*(1 - eta - ksi - zet)*zet; */
    /* N[3] = zet - 64.0*eta*ksi*(1 - eta - ksi - zet)*zet; */
    /* N[4] = 256.0*eta*ksi*(1 - eta - ksi - zet)*zet; */

  }
  if (nne == 8){/* Linear brick */
    N[0] = 1./8.*(1.-ksi)*(1.-eta)*(1.-zet);
    N[1] = 1./8.*(1.+ksi)*(1.-eta)*(1.-zet);
    N[2] = 1./8.*(1.+ksi)*(1.+eta)*(1.-zet);
    N[3] = 1./8.*(1.-ksi)*(1.+eta)*(1.-zet);

    N[4] = 1./8.*(1.-ksi)*(1.-eta)*(1.+zet);
    N[5] = 1./8.*(1.+ksi)*(1.-eta)*(1.+zet);
    N[6] = 1./8.*(1.+ksi)*(1.+eta)*(1.+zet);
    N[7] = 1./8.*(1.-ksi)*(1.+eta)*(1.+zet);
  }
  if (nne == 10){/* Quadratic tetrahedron */
    N[0] = (2.*(1.-ksi-eta-zet)-1.)*(1.-ksi-eta-zet);
    N[1] = (2.*ksi-1.)*ksi;
    N[2] = (2.*eta-1.)*eta;
    N[3] = (2.*zet-1.)*zet;
    N[4] = 4.*(1.-ksi-eta-zet)*ksi;
    N[5] = 4.*ksi*eta;
    N[6] = 4.*(1.-ksi-eta-zet)*eta;
    N[7] = 4.*(1.-ksi-eta-zet)*zet;
    N[8] = 4.*ksi*zet;
    N[9] = 4.*eta*zet;
  }

  if(nne == 6){ /* linear wdge */
    N[0] = 0.5*(1-zet)*ksi;
    N[1] = 0.5*(1-zet)*eta;
    N[2] = 0.5*(1-zet)*(1-ksi-eta);
    N[3] = 0.5*(1+zet)*ksi;
    N[4] = 0.5*(1+zet)*eta;
    N[5] = 0.5*(1+zet)*(1-ksi-eta);
  }
}

void shape_2D (const long nne,
	       const double ksi,
	       const double eta,
	       double *N)
     /*
       ksi - integration poin
       eta - integration poin
       N - shape functions
     */
{
  if (nne == 3){/* Degenerated Triangle */
    *(N+0) = 0.25*(1.-ksi)*(1.-eta);
    *(N+1) = 0.25*(1.+ksi)*(1.-eta);		  
    *(N+2) = 0.50*(1.+eta);
  }
  if (nne == 4){/*Quadrilateral */
    *(N+0) = 0.25*(1.-ksi)*(1.-eta);
    *(N+1) = 0.25*(1.+ksi)*(1.-eta);
    *(N+2) = 0.25*(1.+ksi)*(1.+eta);
    *(N+3) = 0.25*(1.-ksi)*(1.+eta);
  }
}

void dN_kez (const double ksi,
	     const double eta,
	     const double zet,
	     const long nne,
	     double *N_ksi,
	     double *N_eta,
	     double *N_zet)
     /*
       
     */
{
  if (nne == 4){/* Linear tetrahedron */
    N_ksi[0] = -1.;
    N_eta[0] = -1.;
    N_zet[0] = -1.;

    N_ksi[1] = +1.;
    N_eta[1] =  0.;
    N_zet[1] =  0.;

    N_ksi[2] =  0.;
    N_eta[2] = +1.;
    N_zet[2] =  0.;

    N_ksi[3] =  0.;
    N_eta[3] =  0.;
    N_zet[3] = +1.;

  }
  if (nne == 5){ /* Linear tetra + bubble */
    N_ksi[0] = -1.;
    N_eta[0] = -1.;
    N_zet[0] = -1.;

    N_ksi[1] = +1.;
    N_eta[1] =  0.;
    N_zet[1] =  0.;

    N_ksi[2] =  0.;
    N_eta[2] = +1.;
    N_zet[2] =  0.;

    N_ksi[3] =  0.;
    N_eta[3] =  0.;
    N_zet[3] = +1.;

    if(NO_BUBBLE){
      N_ksi[4] = N_eta[4] = N_zet[4] = 0;
    } else {
      N_ksi[4] = 256.*eta*zet*(1-2*ksi-eta-zet);
      N_eta[4] = 256.*ksi*zet*(1-ksi-2*eta-zet);
      N_zet[4] = 256.*eta*ksi*(1-ksi-eta-2*zet);
    }

    /* N_ksi[0] = (-1 + 64*eta*zet*(-1 + eta + 2*ksi + zet)); */
    /* N_ksi[1] = (1 + 64*eta*zet*(-1 + eta + 2*ksi + zet)); */
    /* N_ksi[2] = 64*eta*zet*(-1 + eta + 2*ksi + zet); */
    /* N_ksi[3] = 64*eta*zet*(-1 + eta + 2*ksi + zet); */
    /* N_ksi[4] = -256*eta*zet*(-1 + eta + 2*ksi + zet); */

    /* N_eta[0] = -1 + 64*ksi*zet*(-1 + 2*eta + ksi + zet); */
    /* N_eta[1] = 64*ksi*zet*(-1 + 2*eta + ksi + zet); */
    /* N_eta[2] = 1 + 64*ksi*zet*(-1 + 2*eta + ksi + zet); */
    /* N_eta[3] = 64*ksi*zet*(-1 + 2*eta + ksi + zet); */
    /* N_eta[4] = -256*ksi*zet*(-1 + 2*eta + ksi + zet); */

    /* N_zet[0] = -1 + 64*eta*ksi*(-1 + eta + ksi + 2*zet); */
    /* N_zet[1] = 64*eta*ksi*(-1 + eta + ksi + 2*zet); */
    /* N_zet[2] = 64*eta*ksi*(-1 + eta + ksi + 2*zet); */
    /* N_zet[3] = 1 + 64*eta*ksi*(-1 + eta + ksi + 2*zet); */
    /* N_zet[4] = -256*eta*ksi*(-1 + eta + ksi + 2*zet); */

  }
  if (nne == 8){/* Linear brick */
    N_ksi[0] = -1./8.*(1.-eta)*(1.-zet);
    N_eta[0] = -1./8.*(1.-ksi)*(1.-zet);
    N_zet[0] = -1./8.*(1.-ksi)*(1.-eta);

    N_ksi[1] = +1./8.*(1.-eta)*(1.-zet);
    N_eta[1] = -1./8.*(1.+ksi)*(1.-zet);
    N_zet[1] = -1./8.*(1.+ksi)*(1.-eta);

    N_ksi[2] = +1./8.*(1.+eta)*(1.-zet);
    N_eta[2] = +1./8.*(1.+ksi)*(1.-zet);
    N_zet[2] = -1./8.*(1.+ksi)*(1.+eta);

    N_ksi[3] = -1./8.*(1.+eta)*(1.-zet);
    N_eta[3] = +1./8.*(1.-ksi)*(1.-zet);
    N_zet[3] = -1./8.*(1.-ksi)*(1.+eta);

    N_ksi[4] = -1./8.*(1.-eta)*(1.+zet);
    N_eta[4] = -1./8.*(1.-ksi)*(1.+zet);
    N_zet[4] = +1./8.*(1.-ksi)*(1.-eta);

    N_ksi[5] = +1./8.*(1.-eta)*(1.+zet);
    N_eta[5] = -1./8.*(1.+ksi)*(1.+zet);
    N_zet[5] = +1./8.*(1.+ksi)*(1.-eta);

    N_ksi[6] = +1./8.*(1.+eta)*(1.+zet);
    N_eta[6] = +1./8.*(1.+ksi)*(1.+zet);
    N_zet[6] = +1./8.*(1.+ksi)*(1.+eta);

    N_ksi[7] = -1./8.*(1.+eta)*(1.+zet);
    N_eta[7] = +1./8.*(1.-ksi)*(1.+zet);
    N_zet[7] = +1./8.*(1.-ksi)*(1.+eta);

  }
  if (nne == 10){/* Quadratic tetrahedron */
    N_ksi[0] = 4.*(ksi+eta+zet)-3.;
    N_eta[0] = 4.*(eta+ksi+zet)-3.;
    N_zet[0] = 4.*(zet+ksi+eta)-3.;

    N_ksi[1] = 4.*ksi-1.;
    N_eta[1] = 0.0;
    N_zet[1] = 0.0;

    N_ksi[2] = 0.0;
    N_eta[2] = 4.*eta-1.;
    N_zet[2] = 0.0;

    N_ksi[3] = 0.0;
    N_eta[3] = 0.0;
    N_zet[3] = 4.*zet-1.;

    N_ksi[4] = 4.*(-2.*ksi-eta-zet+1.);
    N_eta[4] = -4.*ksi;
    N_zet[4] = -4.*ksi;

    N_ksi[5] = 4.*eta;
    N_eta[5] = 4.*ksi;
    N_zet[5] = 0.0;

    N_ksi[6] = -4.*eta;
    N_eta[6] = 4.*(-ksi-2.*eta-zet+1.);
    N_zet[6] = -4.*eta;

    N_ksi[7] = -4.*zet;
    N_eta[7] = -4.*zet;
    N_zet[7] = 4.*(-ksi-eta-2.*zet+1.);

    N_ksi[8] = 4.*zet;
    N_eta[8] = 0.0;
    N_zet[8] = 4.*ksi;

    N_ksi[9] = 0.0;
    N_eta[9] = 4.*zet;
    N_zet[9] = 4.*eta;

  }

  if(nne == 6){
    N_ksi[0] = 0.5*(1-zet);
    N_eta[0] = 0.0;
    N_zet[0] = -0.5*ksi;

    N_ksi[1] = 0.0;
    N_eta[1] = 0.5*(1-zet);
    N_zet[1] = -0.5*eta;

    N_ksi[2] = -0.5*(1-zet);
    N_eta[2] = -0.5*(1-zet);
    N_zet[2] = -0.5*(1-ksi-eta);

    N_ksi[3] = 0.5*(1+zet);
    N_eta[3] = 0.0;
    N_zet[3] = 0.5*ksi;

    N_ksi[4] = 0.0;
    N_eta[4] = 0.5*(1+zet);
    N_zet[4] = 0.5*eta;

    N_ksi[5] = -0.5*(1+zet);
    N_eta[5] = -0.5*(1+zet);
    N_zet[5] = 0.5*(1-ksi-eta);
  }
}

void dxyz_kez (const double ksi,
	       const double eta,
	       const double zet,
	       const long nne,
	       const double *x,
	       const double *y,
	       const double *z,
	       const double *N_ksi,
	       const double *N_eta,
	       const double *N_zet,
	       double *dx,
	       double *dy,
	       double *dz)
     /*
       dx - derivace x podle ksi, eta, zet
       dy - derivace y podle ksi, eta, zet
       dz - derivace z podle ksi, eta, zet
     */
{
  long i;
  
  dx[0] = 0.0;
  dx[1] = 0.0;
  dx[2] = 0.0;

  dy[0] = 0.0;
  dy[1] = 0.0;
  dy[2] = 0.0;

  dz[0] = 0.0;
  dz[1] = 0.0;
  dz[2] = 0.0;
  
  for (i=0;i<nne;i++){
    dx[0] += N_ksi[i]*x[i];
    dx[1] += N_eta[i]*x[i];
    dx[2] += N_zet[i]*x[i];

    dy[0] += N_ksi[i]*y[i];
    dy[1] += N_eta[i]*y[i];
    dy[2] += N_zet[i]*y[i];

    dz[0] += N_ksi[i]*z[i];
    dz[1] += N_eta[i]*z[i];
    dz[2] += N_zet[i]*z[i];
  }
}

double Jacobi (const double ksi,
	       const double eta,
	       const double zet,
	       const double *x,
	       const double *y,
	       const double *z,
	       const double *dx,
	       const double *dy,
	       const double *dz)
{
  double J;
  
  J = ((dx[0]*dy[1]*dz[2]) +
       (dy[0]*dz[1]*dx[2]) + 
       (dz[0]*dx[1]*dy[2]) - 
       (dz[0]*dy[1]*dx[2]) - 
       (dx[0]*dz[1]*dy[2]) - 
       (dy[0]*dx[1]*dz[2]));
  
  if (J <= 0.0) {
    PGFEM_printf ("Negative determinant is isoparametric"
		  " mapping : J = %12.12e || Bye Bye\n",J);
    PGFEM_Abort();
  }
  
  return (J);
}

double deriv (const double ksi,
	      const double eta,
	      const double zet,
	      const long nne,
	      const double *x,
	      const double *y,
	      const double *z,
	      double *N_x,
	      double *N_y,
	      double *N_z)
     /*
       
     */
{
  double *N_ksi,*N_eta,*N_zet,*dx,*dy,*dz,J,**Jac;
  long i;
  
  N_ksi = aloc1(nne); 
  N_eta = aloc1(nne);
  N_zet = aloc1(nne);
  dx = aloc1(3);
  dy = aloc1(3);
  dz = aloc1(3);
  Jac = aloc2(3,3);
  
  dN_kez (ksi,eta,zet,nne,N_ksi,N_eta,N_zet);
  
  dxyz_kez (ksi,eta,zet,nne,x,y,z,N_ksi,N_eta,N_zet,dx,dy,dz);
  
  J = Jacobi (ksi,eta,zet,x,y,z,dx,dy,dz);
  
  Jac[0][0] = dy[1]*dz[2] - dz[1]*dy[2];
  Jac[0][1] = dz[0]*dy[2] - dy[0]*dz[2];
  Jac[0][2] = dy[0]*dz[1] - dz[0]*dy[1];

  Jac[1][0] = dz[1]*dx[2] - dx[1]*dz[2];
  Jac[1][1] = dx[0]*dz[2] - dz[0]*dx[2];
  Jac[1][2] = dz[0]*dx[1] - dx[0]*dz[1];

  Jac[2][0] = dx[1]*dy[2] - dy[1]*dx[2];
  Jac[2][1] = dy[0]*dx[2] - dx[0]*dy[2];
  Jac[2][2] = dx[0]*dy[1] - dy[0]*dx[1];
  
  for (i=0;i<nne;i++){
    N_x[i] = 1./J * (Jac[0][0]*N_ksi[i] 
		     + Jac[0][1]*N_eta[i] 
		     + Jac[0][2]*N_zet[i]);

    N_y[i] = 1./J * (Jac[1][0]*N_ksi[i] 
		     + Jac[1][1]*N_eta[i] 
		     + Jac[1][2]*N_zet[i]);

    N_z[i] = 1./J * (Jac[2][0]*N_ksi[i] 
		     + Jac[2][1]*N_eta[i] 
		     + Jac[2][2]*N_zet[i]);
  }
  
  dealoc1 (N_ksi);
  dealoc1 (N_eta);
  dealoc1 (N_zet);
  dealoc1 (dx);
  dealoc1 (dy);
  dealoc1 (dz);
  dealoc2 (Jac,3);
  
  return (J);
}

void get_element_node_parent_coords(const int nne,
				    double *ksi,
				    double *eta,
				    double *zet)
{
  /* NOTE: When adding new elements, be sure that the coordinates are
     consistend with the shape function order!!! */
  switch(nne){

    /* Linear tetrahedrons */
  case 5:
    /* bubble does not change the element coordinate
       transformation. fall through to case 4 */
  case 4:
    ksi[0] = 0.0; eta[0] = 0.0; zet[0] = 0.0;
    ksi[1] = 1.0; eta[1] = 0.0; zet[1] = 0.0;
    ksi[2] = 0.0; eta[2] = 1.0; zet[2] = 0.0;
    ksi[3] = 0.0; eta[3] = 0.0; zet[3] = 1.0;
    break;

    /* Tri-linear brick */
  case 8:
    ksi[0] = -1.0; eta[0] = -1.0; zet[0] = -1.0;
    ksi[1] =  1.0; eta[1] = -1.0; zet[1] = -1.0;
    ksi[2] =  1.0; eta[2] =  1.0; zet[2] = -1.0;
    ksi[3] = -1.0; eta[3] =  1.0; zet[3] = -1.0;

    ksi[4] = -1.0; eta[4] = -1.0; zet[4] =  1.0;
    ksi[5] =  1.0; eta[5] = -1.0; zet[5] =  1.0;
    ksi[6] =  1.0; eta[6] =  1.0; zet[6] =  1.0;
    ksi[7] = -1.0; eta[7] =  1.0; zet[7] =  1.0;
    break;

  default:
    PGFEM_printerr("%s does not currently support element type %d.\n",
	    __func__,nne);
    PGFEM_Abort();
  }
}

void element_center(const int nne,
		    double *x,
		    double *y,
		    double *z)
{
  /* The centroid of an general polyhedron consisting of a finite set
     of points is simply the average of its nodal coordinates */
  x[nne] = y[nne] = z[nne] = 0;
  for(int i=0; i<nne; i++){
    x[nne] += x[i];
    y[nne] += y[i];
    z[nne] += z[i];
  }
  x[nne] /= nne;
  y[nne] /= nne;
  z[nne] /= nne;
}

int element_center_kez(const int nne,
		       double *x,
		       double *y,
		       double *z)
{
  /*** NOTE: length(x_i) >= nne + 1 ***/
  /* compute the element center in (x,y,z) from the center in
     (k,e,z). x_i = sum(Na(k,e,z)*x_i) */
  int err = 0;
  double ksi,eta,zet;
  double *Na;

  /* get (k,e,z) coordinates of element center */
  switch(nne){
  case 4: /* linear tetra */
    ksi = eta = zet = 0.25;
    break;

  default:
    PGFEM_printf("ERROR: element type not supported for computing center\n");
    return 1;
  }

  /* Compute the element center in (x,y,z) */
  Na = aloc1(nne);
  shape_func(ksi,eta,zet,nne,Na);
  x[nne] = y[nne] = z[nne] = 0;
  for(int i=0; i<nne; i++){
    x[nne] += Na[i]*x[i];
    y[nne] += Na[i]*y[i];
    z[nne] += Na[i]*z[i];
  }
  free(Na);
  return err;
}

void get_bubble_grad(const int nne_t, /* The bubble is the last node */
		     const double ksi,
		     const double eta,
		     const double zet,
		     const double *x,
		     const double *y,
		     const double *z,
		     double *N_x,
		     double *N_y,
		     double *N_z)
{
  double *N_ksi, *N_eta, *N_zet,*dx,*dy,*dz,**Jac;

  N_ksi = aloc1(nne_t); 
  N_eta = aloc1(nne_t);
  N_zet = aloc1(nne_t);
  dx = aloc1(3);
  dy = aloc1(3);
  dz = aloc1(3);
  Jac = aloc2(3,3);

  dN_kez (ksi,eta,zet,nne_t,N_ksi,N_eta,N_zet);

  /* only get dx, dy, dz for linear element */
  dxyz_kez (ksi,eta,zet,nne_t-1,
	    x,y,z,N_ksi,N_eta,N_zet,
	    dx,dy,dz);

  double J = Jacobi (ksi,eta,zet,x,y,z,dx,dy,dz);

  Jac[0][0] = dy[1]*dz[2] - dz[1]*dy[2];
  Jac[0][1] = dz[0]*dy[2] - dy[0]*dz[2];
  Jac[0][2] = dy[0]*dz[1] - dz[0]*dy[1];

  Jac[1][0] = dz[1]*dx[2] - dx[1]*dz[2];
  Jac[1][1] = dx[0]*dz[2] - dz[0]*dx[2];
  Jac[1][2] = dz[0]*dx[1] - dx[0]*dz[1];

  Jac[2][0] = dx[1]*dy[2] - dy[1]*dx[2];
  Jac[2][1] = dy[0]*dx[2] - dx[0]*dy[2];
  Jac[2][2] = dx[0]*dy[1] - dy[0]*dx[1];


  N_x[nne_t-1] = 1./J * (Jac[0][0]*N_ksi[nne_t-1] 
			 + Jac[0][1]*N_eta[nne_t-1] 
			 + Jac[0][2]*N_zet[nne_t-1]);

  N_y[nne_t-1] = 1./J * (Jac[1][0]*N_ksi[nne_t-1] 
			 + Jac[1][1]*N_eta[nne_t-1] 
			 + Jac[1][2]*N_zet[nne_t-1]);

  N_z[nne_t-1] = 1./J * (Jac[2][0]*N_ksi[nne_t-1] 
			 + Jac[2][1]*N_eta[nne_t-1] 
			 + Jac[2][2]*N_zet[nne_t-1]);

  
  dealoc1 (N_ksi);
  dealoc1 (N_eta);
  dealoc1 (N_zet);
  dealoc1 (dx);
  dealoc1 (dy);
  dealoc1 (dz);
  dealoc2 (Jac,3);
}

double Bmat (const double ksi,
	     const double eta,
	     const double zet,
	     const long nne,
	     const double *x,
	     const double *y,
	     const double *z,
	     double **B)
     /*
       
     */
{
  long i;
  double *N_x,*N_y,*N_z,J;
  
  N_x = aloc1(nne);  N_y = aloc1(nne);  N_z = aloc1(nne);
  
  J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
  
  for (i=0;i<nne;i++){
    B[i*3+0][0] = N_x[i];
    B[i*3+0][1] = 0.0;
    B[i*3+0][2] = 0.0;
    B[i*3+0][3] = 0.0;
    B[i*3+0][4] = N_z[i];
    B[i*3+0][5] = N_y[i];
 
    B[i*3+1][0] = 0.0;
    B[i*3+1][1] = N_y[i];
    B[i*3+1][2] = 0.0;
    B[i*3+1][3] = N_z[i];
    B[i*3+1][4] = 0.0;
    B[i*3+1][5] = N_x[i];
 
    B[i*3+2][0] = 0.0;
    B[i*3+2][1] = 0.0;
    B[i*3+2][2] = N_z[i];
    B[i*3+2][3] = N_y[i];
    B[i*3+2][4] = N_x[i];
    B[i*3+2][5] = 0.0;
 
  }
  
  dealoc1 (N_x);  dealoc1 (N_y);  dealoc1 (N_z);
  
  return (J);
}

void B_BAR (double **B,
	    const long nne,
	    const double *x,
	    const double *y,
	    const double *z)
     /*
       
     */
{
  long i;
  double *N_x,*N_y,*N_z,J,ksi,eta,zet,B1,B2,B3;
  
  N_x = aloc1(nne);  N_y = aloc1(nne);  N_z = aloc1(nne);
  
  if (nne == 4)  {ksi = 1./4.; eta = 1./4.; zet = 1./4.;}
  if (nne == 8)  {ksi = 0.0;   eta = 0.0;   zet = 0.0;}
  if (nne == 10) {ksi = 1./4.; eta = 1./4.; zet = 1./4.;}
  
  J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
  
  for (i=0;i<nne;i++){
    B1 = B[i*3+0][0];
    B2 = B[i*3+1][1];
    B3 = B[i*3+2][2];

    
    B[i*3+0][0] = B1 + (N_x[i] - B1)/3.;
    B[i*3+0][1] =      (N_x[i] - B1)/3.;
    B[i*3+0][2] =      (N_x[i] - B1)/3.;

    B[i*3+1][0] =      (N_y[i] - B2)/3.;
    B[i*3+1][1] = B2 + (N_y[i] - B2)/3.;
    B[i*3+1][2] =      (N_y[i] - B2)/3.;

    B[i*3+2][0] =      (N_z[i] - B3)/3.;
    B[i*3+2][1] =      (N_z[i] - B3)/3.;
    B[i*3+2][2] = B3 + (N_z[i] - B3)/3.;

  }
  
  dealoc1 (N_x);  dealoc1 (N_y);  dealoc1 (N_z);
}

void stiffmatel (long ii,
		 double *x,
		 double *y,
		 double *z,
		 long nne,
		 long ndofn,
		 ELEMENT *elem,
		 HOMMAT *hommat,
		 NODE *node,
		 double *K,
		 const PGFem3D_opt *opts)
     /*
       
     */
{
  long i,j,k,jj,kk,II,JJ,KK,ndofe,ip;
  double ksi,eta,zet,ai,aj,ak,J,**A,**B_T,**BDB;
  double *gk,*ge,*gz,*w,**D;
  
  ndofe = ndofn*nne;
  
  gk = aloc1(5);
  ge = aloc1(5);
  gz = aloc1(5);
  w = aloc1(5);
  A = aloc2(ndofe,6);
  B_T = aloc2(ndofe,6);
  BDB = aloc2(ndofe,ndofe);
  D = aloc2 (6,6);

  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  nulld (K,ndofe*ndofe);  ip = 0;
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){
	
	if (nne == 4)  {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 10) {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 8)  {
	  ksi = *(gk+i);
	  eta = *(gk+j);
	  zet = *(gk+k);
	  ai = *(w+i);
	  aj = *(w+j);
	  ak = *(w+k);
}
	
	J = Bmat (ksi,eta,zet,nne,x,y,z,B_T);
	
	/*  Matice B_BAR */
	B_BAR (B_T,nne,x,y,z);
	
	/**************************/
	/* Vypocet matice tuhosti */
	/**************************/
	
	Stiffness_Matrix_3D (ii,ip,elem,hommat,D,opts->analysis_type);
	
	nas_AB  (B_T,D,A,ndofe,6,6);
	nas_ABT (A,B_T,BDB,ndofe,6,ndofe);
	
	for (jj=0;jj<ndofe;jj++){
	  for (kk=0;kk<ndofe;kk++){
	    K[jj*ndofe+kk] += ai*aj*ak*BDB[jj][kk]*J;
	  }
	}
	ip++;
      }/*end of k*/
    }/*end of j*/
  }/*end of i*/
  
  dealoc2 (A,ndofe);
  dealoc2 (B_T,ndofe);
  dealoc2 (BDB,ndofe);
  dealoc2 (D,6);
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  
}

void stiffmat (long *adr,
	       long ne,
	       long ndofn,
	       ELEMENT *elem,
	       NODE *node,
	       HOMMAT *hommat,
	       long *ci,
	       long typsolveru,
	       double *k,
	       const PGFem3D_opt *opts,
	       const int mp_id)
     /*
       
     k - stiffness matrix of the subdomain
     adr - array of diagonal addresses
     toe - array of element types
     ne - number of elements
     ndofn - number of DOF of one node
     
     */
{
  long i,nne,*cn,*nod,ndofe;
  double *lk,*x,*y,*z;
  
 
  /*  pro prvky s vice nez 10 uzly a 30 stupni volnosti je treba
      predelat nasledujici alokaci */
  
  cn = aloc1l (30);
  nod = aloc1l (10);
  lk= aloc1 (30*30);
  x = aloc1 (10);
  y = aloc1 (10);
  z = aloc1 (10);

  
  for (i=0;i<ne;i++){
    
    nne = elem[i].toe;  
    ndofe = ndofn*nne;
    elemnodes (i,nne,nod,elem);
    switch(opts->analysis_type){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
    
    /* Stiffness matrix of one element */
    stiffmatel (i,x,y,z,nne,ndofn,elem,hommat,node,lk,opts);
    
    /* ID numbers */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);
    
    /* Assembly */
    if (typsolveru == 0) localizat (k,lk,adr,cn,ndofe);
    if (typsolveru == 1) lokalizace_scr (k,lk,cn,adr,ci,ndofe);
  }
  dealoc1l (cn);
  dealoc1l (nod);
  dealoc1 (lk);
  dealoc1 (x);
  dealoc1(y);
  dealoc1 (z);

}
