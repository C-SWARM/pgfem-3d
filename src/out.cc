#include "out.h"
#include <sys/time.h> 
#include <sys/resource.h>
#include <math.h>
#include "enumerations.h"
#include "def_grad.h"
#include "incl.h"
#include "allocation.h"
#include "utils.h"
#include "elem3d.h"
#include "gen_path.h"
#include "cast_macros.h"

static const int periodic = 0;
static const int ndim = 3;

static const char *PGFEM_LOGO = 
  " _______    ______   ________                        ______   _______  \n"
  "/       \\  /      \\ /        |                      /      \\ /       \\ \n"
  "$$$$$$$  |/$$$$$$  |$$$$$$$$/______   _____  ____  /$$$$$$  |$$$$$$$  |\n"
  "$$ |__$$ |$$ | _$$/ $$ |__  /      \\ /     \\/    \\ $$ ___$$ |$$ |  $$ |\n"
  "$$    $$/ $$ |/    |$$    |/$$$$$$  |$$$$$$ $$$$  |  /   $$< $$ |  $$ |\n"
  "$$$$$$$/  $$ |$$$$ |$$$$$/ $$    $$ |$$ | $$ | $$ | _$$$$$  |$$ |  $$ |\n"
  "$$ |      $$ \\__$$ |$$ |   $$$$$$$$/ $$ | $$ | $$ |/  \\__$$ |$$ |__$$ |\n"
  "$$ |      $$    $$/ $$ |   $$       |$$ | $$ | $$ |$$    $$/ $$    $$/ \n"
  "$$/        $$$$$$/  $$/     $$$$$$$/ $$/  $$/  $$/  $$$$$$/  $$$$$$$/  \n";

void logo (FILE *out)
     /*
       
     */
{
  PGFEM_fprintf (out,"\n%s\n",PGFEM_LOGO);
}

void coordinates (FILE *out, NODE *node, long nn)
/* This function prints out the nodal coordinates.  Created to aid in
   comparing parallel output to single processor output */
{
  long i;

  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"*****************\n");
  PGFEM_fprintf (out,"*  COORDINATES  *\n");
  PGFEM_fprintf (out,"*****************\n");
  PGFEM_fprintf (out,"\n");

  PGFEM_fprintf(out,"*NODE [Gnn][Dom][Lnn]*\n");

  for(i=0;i<nn;i++)
    {
      PGFEM_fprintf (out,"*NODE [%ld][%ld][%ld]*  ", node[i].Gnn, node[i].Dom, i);
      PGFEM_fprintf(out,"%12.12f %12.12f  %12.12f\n", node[i].x1, node[i].x2, node[i].x3);
    }
}

void deform (FILE *out,NODE *node,ELEMENT *elem,long nn,long ne,long ndofn,SUPP sup,double *r)
     /*
       
     */
     
{
  long i,j;
  double *rl;
  
  rl = aloc1 (ndofn);
  
  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"****************\n");
  PGFEM_fprintf (out,"* DISPLACEMENT *\n");
  PGFEM_fprintf (out,"****************\n");
  PGFEM_fprintf (out,"\n");
  
  for (i=0;i<nn;i++){
    PGFEM_fprintf (out,"*NODE [%ld][%ld][%ld]*  ", node[i].Gnn, node[i].Dom, i);
    
    nulld (rl,ndofn);
    
    for (j=0;j<ndofn;j++){
      
      if (node[i].id[j] > 0)  
	rl[j] = r[node[i].id[j]-1];
      if (node[i].id[j] < 0) 
	rl[j] = sup->defl[abs(node[i].id[j])-1];
      
      PGFEM_fprintf (out,"%12.12f  ",rl[j]);
    }
    PGFEM_fprintf (out,"\n");
  }
  dealoc1 (rl);
}

void stress_out (FILE *out,long ne,long nn,ELEMENT *elem,SIG *sig_e,SIG *sig_n,long gr4)
     /*
      */
{
  long k,jj;
  
  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"*****************************************************\n");
  PGFEM_fprintf (out,"* STRESS on ELEMENT sig={s11,s22,s33,s23,s13,s12}^T *\n");
  PGFEM_fprintf (out,"*****************************************************\n");
  PGFEM_fprintf (out,"\n");
  
  for (jj=0;jj<ne;jj++){
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"ELEMENT [%ld], Pr.[%ld]\n",jj,elem[jj].pr);
    PGFEM_fprintf (out,"Mises  %12.8f\n",sig_e[jj].el.eq);
    for (k=0;k<6;k++)
      PGFEM_fprintf (out,"%12.8f  ",sig_e[jj].el.o[k]);
    PGFEM_fprintf (out,"\n");
  }
  
  if (gr4 == 0){
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"**************************************************\n");
    PGFEM_fprintf (out,"* SIGMA on NODES sig={s11,s22,s33,s23,s13,s12}^T *\n");
    PGFEM_fprintf (out,"**************************************************\n");
    PGFEM_fprintf (out,"\n");
    
    for (jj=0;jj<nn;jj++){
      PGFEM_fprintf (out,"\n");
      PGFEM_fprintf (out,"NODE [%ld]\n",jj);
      PGFEM_fprintf (out,"Mises  %12.8f\n",sig_n[jj].el.eq);
      for (k=0;k<6;k++)
	PGFEM_fprintf (out,"%12.8f  ",sig_n[jj].el.o[k]);
      PGFEM_fprintf (out,"\n");
    }
  }
}

void strain_out (FILE *out,long ne,ELEMENT *elem,EPS *eps,const PGFem3D_opt *opts)
     /*
      */
{
  long k,jj;
  
  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"*****************************************************\n");
  PGFEM_fprintf (out,"* STRAIN on ELEMENT sig={e11,e22,e33,e23,e13,e12}^T *\n");
  PGFEM_fprintf (out,"*****************************************************\n");
  PGFEM_fprintf (out,"\n");
  
  if (opts->analysis_type == FS_CRPL && periodic == 1) PGFEM_fprintf (out,"Macro-scale Effective Plastic Strain || Eeq = %12.8f\n",eps[0].eff);
  
  for (jj=0;jj<ne;jj++){
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"ELEMENT [%ld], Pr.[%ld]\n",jj,elem[jj].pr);
    PGFEM_fprintf (out,"Mises  %12.8f\n",eps[jj].el.eq);
    for (k=0;k<6;k++) PGFEM_fprintf (out,"%12.8f  ",eps[jj].el.o[k]);
    PGFEM_fprintf (out,"\n");
    if (opts->analysis_type == FS_CRPL){
      PGFEM_fprintf (out,"Plastic strain\n");
      PGFEM_fprintf (out,"Mises  %12.8f || %12.8f\n",eps[jj].pl.eq[0],eps[jj].pl.eq[1]);
      for (k=0;k<6;k++) PGFEM_fprintf (out,"%12.8f  ",eps[jj].pl.o[k]);
      PGFEM_fprintf (out,"\n");
    }
  }/* end jj < ne */
}

void deform_grad_out (FILE *out,long ne,ELEMENT *elem,EPS *eps)
     /*
       
     */
{
  /** output the deformation gradient at the integration points **/

}

void macro_fields_out (FILE *out,EPS *eps,const PGFem3D_opt *opts)
     /*
       
     */
{
  long M,N,P,Q;
  double **s,**e,**F,**FB,**F_I,J,**S,**E,dij,*Sij,S_eq,**FeFp,E_eq;
  
  s = aloc2 (3,3); e = aloc2 (3,3); F_I = aloc2 (3,3); F = aloc2 (3,3); FB = aloc2 (3,3); S = aloc2 (3,3); E = aloc2 (3,3); Sij = aloc1 (6);
  FeFp = aloc2(3,3);
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      FB[M][N] = eps[0].FB[M][N];
      F[M][N] = eps[0].F[M][N];
    }
  }
  
  J = def_grad_det (CCONST_2(double) F);
  def_grad_inv (CCONST_2(double) F,F_I);
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      if (M == N) dij = 1.; else dij = 0.;
      E[M][N] = -1./2.*dij;
      for (P=0;P<3;P++){
	E[M][N] += 1./2.*F[P][M]*F[P][N];
	S[M][N] += F_I[M][P]*eps[0].P[P][N];
	if (opts->analysis_type == FS_CRPL) FeFp[M][N] += eps[0].Fe[M][P]*eps[0].Fp[P][N];
      }
    }
  }
  
  /* Compute Logarithmic strain */
  Logarithmic_strain (F,e);
  
  /* Solve for Cauchy stress */
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      s[M][N] = 0.0;
      for (P=0;P<3;P++){
	for (Q=0;Q<3;Q++){
	  s[M][N] += 1./J*F[M][P]*S[P][Q]*F[N][Q];
	}
      }
    }
  }
  
  Sij[0] = s[0][0] - (s[0][0] + s[1][1] + s[2][2])/3.;
  Sij[1] = s[1][1] - (s[0][0] + s[1][1] + s[2][2])/3.;
  Sij[2] = s[2][2] - (s[0][0] + s[1][1] + s[2][2])/3.;
  Sij[3] = s[1][2];
  Sij[4] = s[0][2];
  Sij[5] = s[0][1];
  
  S_eq = sqrt (3./2.*(Sij[0]*Sij[0] + Sij[1]*Sij[1] + Sij[2]*Sij[2] + 2.*(Sij[3]*Sij[3] + Sij[4]*Sij[4]+ Sij[5]*Sij[5])));
  
  Sij[0] = e[0][0] - (e[0][0] + e[1][1] + e[2][2])/3.;
  Sij[1] = e[1][1] - (e[0][0] + e[1][1] + e[2][2])/3.;
  Sij[2] = e[2][2] - (e[0][0] + e[1][1] + e[2][2])/3.;
  Sij[3] = e[1][2];
  Sij[4] = e[0][2];
  Sij[5] = e[0][1];
  
  E_eq = sqrt (2./3.*(Sij[0]*Sij[0] + Sij[1]*Sij[1] + Sij[2]*Sij[2] + 2.*(Sij[3]*Sij[3] + Sij[4]*Sij[4]+ Sij[5]*Sij[5])));
  
  PGFEM_fprintf (out,"\n\n");
  PGFEM_fprintf (out,"****************\n");
  PGFEM_fprintf (out,"* Macro fields *\n");
  PGFEM_fprintf (out,"****************\n");
  PGFEM_fprintf (out,"\n");
  
  PGFEM_fprintf (out,"**                    F                                               FB                    **\n");
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",eps[0].F[0][0],eps[0].F[0][1],eps[0].F[0][2],FB[0][0],FB[0][1],FB[0][2]);
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",eps[0].F[1][0],eps[0].F[1][1],eps[0].F[1][2],FB[1][0],FB[1][1],FB[1][2]);
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",eps[0].F[2][0],eps[0].F[2][1],eps[0].F[2][2],FB[2][0],FB[2][1],FB[2][2]);
  PGFEM_fprintf (out,"\n");
  if (opts->analysis_type == FS_CRPL){
    PGFEM_fprintf (out,"**                    Fe                                              Fp                                         Fe*Fp               **\n");
    PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",
	     eps[0].Fe[0][0],eps[0].Fe[0][1],eps[0].Fe[0][2],eps[0].Fp[0][0],
	     eps[0].Fp[0][1],eps[0].Fp[0][2],FeFp[0][0],FeFp[0][1],FeFp[0][2]);
    PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",
	     eps[0].Fe[1][0],eps[0].Fe[1][1],eps[0].Fe[1][2],eps[0].Fp[1][0],
	     eps[0].Fp[1][1],eps[0].Fp[1][2],FeFp[1][0],FeFp[1][1],FeFp[1][2]);
    PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",
	     eps[0].Fe[2][0],eps[0].Fe[2][1],eps[0].Fe[2][2],eps[0].Fp[2][0],
	     eps[0].Fp[2][1],eps[0].Fp[2][2],FeFp[2][0],FeFp[2][1],FeFp[2][2]);
    PGFEM_fprintf (out,"\n");
  }/* end opts->analysis_type == FS_CRPL */
  PGFEM_fprintf (out,"**                    E0                                              S0                                                     P         **\n");
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",E[0][0],E[0][1],E[0][2],S[0][0],S[0][1],S[0][2],
	   eps[0].P[0][0],eps[0].P[0][1],eps[0].P[0][2]);
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",E[1][0],E[1][1],E[1][2],S[1][0],S[1][1],S[1][2],
	   eps[0].P[1][0],eps[0].P[1][1],eps[0].P[1][2]);
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f \t %12.12f %12.12f %12.12f\n",E[2][0],E[2][1],E[2][2],S[2][0],S[2][1],S[2][2],
	   eps[0].P[2][0],eps[0].P[2][1],eps[0].P[2][2]);
  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"**       Macro eps (Logarithmic strain)                        sig (Cauchy stress)                   **\n");
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f\t%12.12f %12.12f %12.12f\n",e[0][0],e[0][1],e[0][2],s[0][0],s[0][1],s[0][2]);
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f\t%12.12f %12.12f %12.12f\n",e[1][0],e[1][1],e[1][2],s[1][0],s[1][1],s[1][2]);
  PGFEM_fprintf (out,"%12.12f %12.12f %12.12f\t%12.12f %12.12f %12.12f\n",e[2][0],e[2][1],e[2][2],s[2][0],s[2][1],s[2][2]);
  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"Mises stress = %12.12f | Macro Effective strain = %12.12f || Effective strain = %12.12f\n",S_eq,E_eq,eps[0].eff);
  PGFEM_fprintf (out,"\n");
  
  dealoc2 (s,3); dealoc2 (e,3); dealoc2 (F_I,3); dealoc2 (F,3); dealoc2 (S,3); dealoc2 (E,3); dealoc1 (Sij); dealoc2 (FeFp,3); dealoc2 (FB,3);
}

void cohesive_out (FILE *out,long nce,COEL *coel)
     /*
       
     */
{
  long i,k;
  
  PGFEM_fprintf (out,"\n\n");
  PGFEM_fprintf (out,"*********************\n");
  PGFEM_fprintf (out,"* COHESIVE ELEMENTS *\n");
  PGFEM_fprintf (out,"*********************\n");
  PGFEM_fprintf (out,"\n");
  
  for (i=0;i<nce;i++){
    PGFEM_fprintf (out,"\n");
    if (coel[i].Xn < 0.0) PGFEM_fprintf (out,"COHESIVE ELEMENT [%ld], Pr.[%ld] || COMPRESSION\n",i,coel[i].pr);
    else                  PGFEM_fprintf (out,"COHESIVE ELEMENT [%ld], Pr.[%ld]\n",i,coel[i].pr);
    PGFEM_fprintf (out,"Effective traction  %12.12f : Effective opening %12.12f || Normal traction %12.12f : Shear traction %12.12f\n",
	     coel[i].txi,coel[i].Xxi,coel[i].tn,coel[i].ts);
    PGFEM_fprintf (out,"Tractions\n");
    for (k=0;k<3;k++) PGFEM_fprintf (out,"%12.12f  ",coel[i].ti[k]);
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Effective opening\n");
    for (k=0;k<3;k++) PGFEM_fprintf (out,"%12.12f  ",coel[i].Xi[k]);
    PGFEM_fprintf (out,"\n");
  }
  
}
void damage_out(FILE *out,
		const long ne,
		const ELEMENT *elem,
		const EPS *eps)
{
  /* print out the damage variable for each integration point of the
     element */
  PGFEM_fprintf (out,
	   "\n\n"
	   "*********************\n"
	   "*  DAMAGE VARIABLES  *\n"
	   "*********************\n"
	   "\n");

  /* for each element */
  for(long i=0; i<ne; i++){
    /* get the number of integration points */
    long n_ip = 0;
    int_point(elem[i].toe,&n_ip);

    PGFEM_fprintf(out,
	    "\nELEMENT [%ld], Pr.[%ld] \n"
	    "Damage variable at int points (w || X)\n",
	    i,elem[i].pr);

    /* loop over integration points and print damage variables */
    for(long j=0; j<n_ip; j++){
      /* Print values at n since those are the updated values after
	 the solve */
      PGFEM_fprintf(out,"(%12.12f || %12.12f) ",
	      eps[i].dam[j].wn,eps[i].dam[j].Xn);
    }
    PGFEM_fprintf(out,"\n");
  }/* end each element */
}/* damage_out() */

void elixir (char jmeno[50],
	     long nn,
	     long ne,
	     long ndofn,
	     NODE *node,
	     ELEMENT *elem,
	     SUPP sup,
	     double *r,
	     SIG *sig_e,
	     SIG *sig_n,
	     EPS *eps,
	     long gr4,
	     long nce,
	     COEL *coel,
	     const PGFem3D_opt *opts)
     /*
       
     */
	  
{
  long i,j,jj,k,*nod,nne,M,N;
  double *rl,*X,*xx;
  
  FILE *out;
  
  nod = aloc1l (10); rl = aloc1 (ndofn); X = aloc1 (3); xx = aloc1 (3);
  
  if ((out = fopen(jmeno,"w")) == NULL ){
    PGFEM_printf("Vystupni soubor pro ELIXIR nelze otevrit\n");
    PGFEM_Abort();
  } 
  
  PGFEM_fprintf (out,"%ld %ld %ld\n",nn,ne,ndofn);
  
  PGFEM_fprintf (out,"\n");
  
  for (i=0;i<ne;i++){
    nne = elem[i].toe;
    PGFEM_fprintf (out,"%ld ",nne);
  }
  
  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"\n");
  
  /* Souradnice uzlu prvku ==> Coordinates of the node element*/
  
  for (i=0;i<nn;i++){
    PGFEM_fprintf (out,"%ld  %ld  %ld  ",node[i].Gnn,node[i].Dom,i);
    
    switch(opts->analysis_type){
    case FS_CRPL:
    case FINITE_STRAIN:
    case STABILIZED:
      if (periodic == 1){/* Homogeneous deformation */
	X[0] = node[i].x1_fd; X[1] = node[i].x2_fd; X[2] = node[i].x3_fd;
	
	for (M=0;M<3;M++){
	  xx[M] = 0.0;
	  for (N=0;N<3;N++){
	    xx[M] += eps[0].F[M][N]*X[N];
	  }
	}
	PGFEM_fprintf (out,"%12.12f %12.12f %12.12f  %ld\n",xx[0],xx[1],xx[2],node[i].pr);
      } else { /* Nondeformed configuration */
	PGFEM_fprintf (out,"%12.12f %12.12f %12.12f  %ld\n",
		 node[i].x1_fd,node[i].x2_fd,node[i].x3_fd,node[i].pr);
      }
      break;
    default:
      PGFEM_fprintf (out,"%12.8f %12.8f %12.8f  %ld\n",
	       node[i].x1,node[i].x2,node[i].x3,node[i].pr);
      break;
    }
  }
  
  PGFEM_fprintf (out,"\n");

  /* Uzly na prvku ==> The nodes on the element */
  
  for (i=0;i<ne;i++){
    nne = elem[i].toe;
    elemnodes (i,nne,nod,elem);
    for (j=0;j<nne;j++)
      PGFEM_fprintf (out,"%ld  ",nod[j]);
    PGFEM_fprintf (out,"%ld",elem[i].pr);
    PGFEM_fprintf (out,"\n");
  }
  
  PGFEM_fprintf (out,"\n");
  
  /* Vektor deformace */
  
  for (i=0;i<nn;i++){
    
    nulld (rl,ndofn);
    
    for (j=0;j<ndofn;j++){
      
      if (node[i].id[j] > 0)  
	rl[j] = r[node[i].id[j]-1];
      if (node[i].id[j] < 0) 
	rl[j] = sup->defl[abs(node[i].id[j])-1];
      
      PGFEM_fprintf (out,"%12.12f  ",rl[j]);
    }
    PGFEM_fprintf (out,"\n");
  }
  
  PGFEM_fprintf (out,"\n");
  
  /* Stress */
  
  for (jj=0;jj<ne;jj++){
    PGFEM_fprintf (out,"%12.8f\n",sig_e[jj].el.eq);
    for (k=0;k<6;k++)
      PGFEM_fprintf (out,"%12.8f  ",sig_e[jj].el.o[k]);
    PGFEM_fprintf (out,"\n");
  }
  
  PGFEM_fprintf (out,"\n");
  
  if (gr4 == 0){
    for (jj=0;jj<nn;jj++){
      PGFEM_fprintf (out,"%12.8f\n",sig_n[jj].el.eq);
      for (k=0;k<6;k++)
	PGFEM_fprintf (out,"%12.8f  ",sig_n[jj].el.o[k]);
      PGFEM_fprintf (out,"\n");
    }
  }
  
  /* Strain */
  if (opts->analysis_type == FS_CRPL){
    PGFEM_fprintf (out,"%12.8f\n\n",eps[0].eff);
    for (jj=0;jj<ne;jj++){
      PGFEM_fprintf (out,"%12.8f\n",eps[jj].pl.eq[0]);
      for (k=0;k<6;k++)
	PGFEM_fprintf (out,"%12.8f  ",eps[jj].pl.o[k]);
      PGFEM_fprintf (out,"\n");
    }
  }/* end opts->analysis_type == FS_CRPL */
  else{
    PGFEM_fprintf (out,"%12.8f\n\n",eps[0].eff);
    for (jj=0;jj<ne;jj++){
      PGFEM_fprintf (out,"%12.8f\n",eps[jj].el.eq); /* difference */
      for (k=0;k<6;k++)
	PGFEM_fprintf (out,"%12.8f  ",eps[jj].el.o[k]);
      PGFEM_fprintf (out,"\n");
    }
  }
  
  /****************************************/
  /* DATA OUTPUT FOR THE COHESIVE SURFACE */
  /****************************************/
  
  PGFEM_fprintf (out,"\n");
  
  if (opts->cohesive == 1){
    
    /* Number of cohesive elements */
    PGFEM_fprintf (out,"%ld\n\n",nce);
    
    /* Type of cohesive element */
    for (i=0;i<nce;i++) PGFEM_fprintf (out,"%ld ",coel[i].toe);
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"\n");
    
    /* Nodes of cohesive element */
    for (i=0;i<nce;i++){
      for (j=0;j<coel[i].toe;j++)
	PGFEM_fprintf (out,"%ld  ",coel[i].nod[j]);
      PGFEM_fprintf (out,"%ld",coel[i].pr);
      PGFEM_fprintf (out,"\n");
    }
    
    PGFEM_fprintf (out,"\n");
    
    /* Effective opening and effective traction */
    for (i=0;i<nce;i++){
      PGFEM_fprintf (out,"%12.12f %12.12f %12.12f  %12.12f %12.12f %12.12f\n",
	       coel[i].Xxi,coel[i].Xn,coel[i].Xs,coel[i].txi,coel[i].tn,coel[i].ts);
    }
    
    PGFEM_fprintf (out,"\n");
    
  }/* end cohesive elements */

  fclose (out);
  
  dealoc1l (nod); dealoc1 (rl); dealoc1 (X); dealoc1 (xx);
}

void EnSight (char jmeno[500],
	      long tim,
	      long nt,
	      long nn,
	      long ne,
	      long ndofn,
	      NODE *node,
	      ELEMENT *elem,
	      SUPP sup,
	      double *r,
	      SIG *sig_e,
	      SIG *sig_n,
	      EPS *eps,
	      long gr4,
	      long nce,
	      COEL *coel,
	      /*long nge,
		GEEL *geel,
		long ngn,
		GNOD *gnod,*/
	      long FNR,
	      double lm,
	      ENSIGHT ensight,
	      MPI_Comm mpi_comm,
	      const PGFem3D_opt *opts)
{ 
  int part1,part2,*GNVpint,*shift;
  long *GNVp,i,j,k,tmp,*property,M,N,tetra4,tetra10,hexa8,*nod,tria3,quad4,*GVp,GNV;
  double *X,*xx,*u;
  char name[500];
  FILE *out,*out1,*out2;

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);
  
  X = aloc1 (3); xx = aloc1 (3); u = aloc1 (3);
  
  if (opts->cohesive == 0)  part1 = 1;
  if (opts->cohesive == 1) {part1 = 1; part2 = 2;}
  
  /* CASE */
  if (tim < nt) sprintf (name,"%s_%d.case%ld",jmeno,myrank,tim);
  else          sprintf (name,"%s_%d.case",jmeno,myrank);
  
  if ((out = fopen(name,"w")) == NULL ){
    PGFEM_printf("File is not possible to open, EnSight\n");
    PGFEM_Abort();
  }
  
  PGFEM_fprintf(out,"FORMAT\n\n");
  PGFEM_fprintf(out,"type:  ensight gold\n\n");
  PGFEM_fprintf(out,"GEOMETRY\n\n");
  if (tim < nt)
    PGFEM_fprintf(out,"model:                      %s_%d.geo%ld\n\n",
	    jmeno,myrank,tim);
  else 
    PGFEM_fprintf(out,"model:                      %s_%d.geo\n\n",
	    jmeno,myrank);
  PGFEM_fprintf(out,"VARIABLE\n\n");
 
 /* Load multiplier */
  if (FNR == 2 || FNR  == 3)
    PGFEM_fprintf(out,"constant per case:          Lambda %12.12f\n",lm);
  if (periodic == 1){

    /* Displacement */
    if (tim < nt)
      PGFEM_fprintf(out,"vector per node:          Fluct_Displ %s_%d.dis%ld\n",
	      jmeno,myrank,tim);
    else 
      PGFEM_fprintf(out,"vector per node:          Fluct_Displ %s_%d.dis\n",
	      jmeno,myrank);

    /* Displacement */
    if (tim < nt)
      PGFEM_fprintf(out,"vector per node:          Total_Displ %s_%d.dit%ld\n",
	      jmeno,myrank,tim);
    else
      PGFEM_fprintf(out,"vector per node:          Total_Displ %s_%d.dit\n",
	      jmeno,myrank);
  }
  else{

    /* Displacement */
    if (tim < nt)
      PGFEM_fprintf(out,"vector per node:          Displacement %s_%d.dis%ld\n",
	      jmeno,myrank,tim);
    else 
      PGFEM_fprintf(out,"vector per node:          Displacement %s_%d.dis\n",
	      jmeno,myrank);
  }

  /* Pressure */
  if (opts->analysis_type == STABILIZED
      || opts->analysis_type == MINI
      || opts->analysis_type == MINI_3F){
    if (tim < nt)
      PGFEM_fprintf(out,"scalar per node:          Pressure %s_%d.pre%ld\n",
	      jmeno,myrank,tim);
    else
      PGFEM_fprintf(out,"scalar per node:          Pressure %s_%d.pre\n",
	      jmeno,myrank);
  }

  /* Stress ELEMENT */
  if (tim < nt)
    PGFEM_fprintf(out,"tensor per element:       Stress %s_%d.ste%ld\n",
	    jmeno,myrank,tim);
  else      
    PGFEM_fprintf(out,"tensor per element:       Stress %s_%d.ste\n",jmeno,myrank);

  /* Strain */
  if (tim < nt)
    PGFEM_fprintf(out,"tensor per element:       Strain %s_%d.sta%ld\n",
	    jmeno,myrank,tim);
  else
    PGFEM_fprintf(out,"tensor per element:       Strain %s_%d.sta\n",
	    jmeno,myrank);

  /* Stress NODES */
  if (gr4 == 0){
    if (tim < nt)
      PGFEM_fprintf(out,"tensor per node:          Nodal_Stress %s_%d.stn%ld\n",
	      jmeno,myrank,tim);
    else 
      PGFEM_fprintf(out,"tensor per node:          Nodal_Stress %s_%d.stn\n",
	      jmeno,myrank);
  }

  /* Effective Strain */
  if (tim < nt)
    PGFEM_fprintf(out,"scalar per element:       Eff_Strain %s_%d.est%ld\n",
	    jmeno,myrank,tim);
  else      
    PGFEM_fprintf(out,"scalar per element:       Eff_Strain %s_%d.est\n",
	    jmeno,myrank);
  if (opts->analysis_type == FS_CRPL){/* Plastic strains */

    /* Effective Plastic Almansi Strain */
    if (tim < nt)
      PGFEM_fprintf(out,"scalar per element:       Eff_Pl_Al_Strain %s_%d.pla%ld\n",
	      jmeno,myrank,tim);
    else
      PGFEM_fprintf(out,"scalar per element:       Eff_Pl_Al_Strain %s_%d.pla\n",
	      jmeno,myrank);

    /* Effective Plastic Strain int_0^t sqrt(2/3*Dp.Dp) */
    if (tim < nt) PGFEM_fprintf(out,"scalar per element:       Eff_Pl_Strain %s_%d.ple%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"scalar per element:       Eff_Pl_Strain %s_%d.ple\n",jmeno,myrank);
  }
  if (opts->cohesive == 1){
    /* Opening Displacement */
    if (tim < nt) PGFEM_fprintf(out,"vector per element:       Opening %s_%d.odi%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"vector per element:       Opening %s_%d.odi\n",jmeno,myrank);
    /* Tractions */
    if (tim < nt) PGFEM_fprintf(out,"vector per element:       Tractions %s_%d.tra%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"vector per element:       Tractions %s_%d.tra\n",jmeno,myrank);
    /* Effective opening */
    if (tim < nt) PGFEM_fprintf(out,"scalar per element:       Eff_Opening %s_%d.ode%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"scalar per element:       Eff_Opening %s_%d.ode\n",jmeno,myrank);
    /* Normal opening */
    if (tim < nt) PGFEM_fprintf(out,"scalar per element:       No_Opening %s_%d.odn%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"scalar per element:       No_Opening %s_%d.odn\n",jmeno,myrank);
    /* Shear opening */
    if (tim < nt) PGFEM_fprintf(out,"scalar per element:       Sh_Opening %s_%d.ods%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"scalar per element:       Sh_Opening %s_%d.ods\n",jmeno,myrank);
    /* Effective tractions */
    if (tim < nt) PGFEM_fprintf(out,"scalar per element:       Eff_Traction %s_%d.tre%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"scalar per element:       Eff_Traction %s_%d.tre\n",jmeno,myrank);
    /* Normal tractions */
    if (tim < nt) PGFEM_fprintf(out,"scalar per element:       No_Traction %s_%d.trn%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"scalar per element:       No_Traction %s_%d.trn\n",jmeno,myrank);
    /* Shear tractions */
    if (tim < nt) PGFEM_fprintf(out,"scalar per element:       Sh_Traction %s_%d.trs%ld\n",jmeno,myrank,tim);
    else          PGFEM_fprintf(out,"scalar per element:       Sh_Traction %s_%d.trs\n",jmeno,myrank);
  }

  PGFEM_fprintf(out,"\n");

  /* Volumetric elements */
  property = aloc1l (ne);
  
  for (i=0;i<ne;i++) property[i] = elem[i].pr;
  
  /* Sort property */
  for(i=0;i<ne-1;i++){
    for(j=i+1;j<ne;j++){
      if(property[i] > property[j]){
	tmp = property[j]; 
	property[j] = property[i]; 
	property[i] = tmp; 
      }
    }
  }
  
  ensight->NVp = 1; for (i=0;i<ne-1;i++) if (property[i] < property[i+1]) ensight->NVp++;
  ensight->Vp = (long *) PGFEM_calloc (ensight->NVp,sizeof(long));
  k = 1; ensight->Vp[0] = property[0]; for (i=0;i<ne-1;i++) if (property[i] < property[i+1]) {ensight->Vp[k] = property[i+1]; k++;}
  
  dealoc1l (property);
  
  /* Gather the number of volume properties from domains */
  GNVp = aloc1l (nproc);
  GNVpint = aloc1i (nproc);
  MPI_Allgather (&ensight->NVp,1,MPI_LONG,GNVp,1,MPI_LONG,mpi_comm);
  GNV = 0;
  for (i=0;i<nproc;i++){ GNV += GNVp[i]; GNVpint[i] = GNVp[i];}
  
  shift = aloc1i (nproc); for (i=1;i<nproc;i++) shift[i] = shift[i-1] + GNVp[i-1];
  
  /* Gather the volumetric property types */
  GVp = aloc1l (GNV); MPI_Allgatherv (ensight->Vp,ensight->NVp,MPI_LONG,GVp,GNVpint,shift,MPI_LONG,mpi_comm);
  
  property = aloc1l (GNV); for (i=0;i<GNV;i++) property[i] = GVp[i];;
  
  /* Sort property */
  for(i=0;i<GNV-1;i++){
    for(j=i+1;j<GNV;j++){
      if(property[i] > property[j]){
	tmp = property[j]; 
	property[j] = property[i]; 
	property[i] = tmp; 
      }
    }
  }
  
  free (ensight->Vp);
  ensight->NVp = 1; for (i=0;i<GNV-1;i++) if (property[i] < property[i+1]) ensight->NVp++;
  ensight->Vp = (long *) PGFEM_calloc (ensight->NVp,sizeof(long));
  k = 1; ensight->Vp[0] = property[0]; for (i=0;i<GNV-1;i++) if (property[i] < property[i+1]) {ensight->Vp[k] = property[i+1]; k++;}
  
  dealoc1l (GNVp); dealoc1i (GNVpint); dealoc1i (shift); dealoc1l (GVp); dealoc1l (property);

  /* Cohesive elements */
  if (opts->cohesive == 1){
    
    /* Alloc cohesive property */
    if (nce == 0) i = 1; else i = nce; property = aloc1l (i);
    
    for (i=0;i<nce;i++) property[i] = coel[i].pr;
    
    /* Sort property */
    for(i=0;i<nce-1;i++){
      for(j=i+1;j<nce;j++){
	if(property[i] > property[j]){
	  tmp = property[j]; 
	  property[j] = property[i]; 
	  property[i] = tmp; 
	}
      }
    }
    
    ensight->NCp = 1; for (i=0;i<nce-1;i++) if (property[i] < property[i+1]) ensight->NCp++;
    ensight->Cp = (long *) PGFEM_calloc (ensight->NCp,sizeof(long));
    k = 1; ensight->Cp[0] = property[0]; for (i=0;i<nce-1;i++) if (property[i] < property[i+1]) {ensight->Cp[k] = property[i+1]; k++;}

    dealoc1l (property);
    
    /* Gather the number of cohesive properties from domains */
    GNVp = aloc1l (nproc);
    GNVpint = aloc1i (nproc);
    MPI_Allgather (&ensight->NCp,1,MPI_LONG,GNVp,1,MPI_LONG,mpi_comm);
    GNV = 0; for (i=0;i<nproc;i++) {GNV += GNVp[i]; GNVpint[i] = GNVp[i];}
    
    shift = aloc1i (nproc); for (i=1;i<nproc;i++) shift[i] = shift[i-1] + GNVp[i-1];
    
    /* Gather cohesive property types */
    GVp = aloc1l (GNV); MPI_Allgatherv (ensight->Cp,ensight->NCp,MPI_LONG,GVp,GNVpint,shift,MPI_LONG,mpi_comm);
    
    property = aloc1l (GNV); for (i=0;i<GNV;i++) property[i] = GVp[i];;
    
    /* Sort property */
    for(i=0;i<GNV-1;i++){
      for(j=i+1;j<GNV;j++){
	if(property[i] > property[j]){
	  tmp = property[j]; 
	  property[j] = property[i]; 
	  property[i] = tmp; 
	}
      }
    }
    
    free (ensight->Cp);
    ensight->NCp = 1; for (i=0;i<GNV-1;i++) if (property[i] < property[i+1]) ensight->NCp++;
    ensight->Cp = (long *) PGFEM_calloc (ensight->NCp,sizeof(long));
    k = 1; ensight->Cp[0] = property[0]; for (i=0;i<GNV-1;i++) if (property[i] < property[i+1]) {ensight->Cp[k] = property[i+1]; k++;}
    
    dealoc1l (GNVp); dealoc1i (GNVpint); dealoc1i (shift); dealoc1l (GVp); dealoc1l (property);
  }/* end coh */
  
  /* The following section is commented out because our version of Paraview
     does not support material assignment */
  /*
  PGFEM_fprintf(out,"MATERIAL\n\n");
  PGFEM_fprintf(out,"material set number:         1 Material\n");
  PGFEM_fprintf(out,"material id count:           %ld\n",ensight->NVp+ensight->NCp);
  PGFEM_fprintf(out,"material id numbers:         ");
  for (i=0;i<ensight->NVp;i++) PGFEM_fprintf(out,"%ld ",ensight->Vp[i]);
  for (i=0;i<ensight->NCp;i++) PGFEM_fprintf(out,"%ld ",ensight->Cp[i]);
  PGFEM_fprintf(out,"\n");
  PGFEM_fprintf(out,"material id names:           ");
  for (i=0;i<ensight->NVp;i++) PGFEM_fprintf(out,"Vol%ld ",ensight->Vp[i]);
  for (i=0;i<ensight->NCp;i++) PGFEM_fprintf(out,"Coh%ld ",ensight->Cp[i]);
  PGFEM_fprintf(out,"\n");
  if (tim < nt) PGFEM_fprintf(out,"material id per element:     %s_%d.mat%ld\n\n",jmeno,myrank,tim);
  else          PGFEM_fprintf(out,"material id per element:     %s_%d.mat\n\n",jmeno,myrank);
  
  */
  
  /* Matt Mosby 03/29/2010 */

  fclose (out);
  /* CLOSE CASE FILE */
  
  /* GEO */
  if (tim < nt) sprintf (name,"%s_%d.geo%ld",jmeno,myrank,tim);
  else          sprintf (name,"%s_%d.geo",jmeno,myrank);
  
  if ((out = fopen(name,"w")) == NULL ){
    PGFEM_printf("File is not possible to open, EnSight\n");
    PGFEM_Abort();
  }

  PGFEM_fprintf(out,"Geometry file\n");
  PGFEM_fprintf(out,"Non-deformed configuration and/or macroscopic deformation for multi-scale analysis\n");
  PGFEM_fprintf(out,"node id given\n");
  PGFEM_fprintf(out,"element id given\n");
  PGFEM_fprintf(out,"part\n");
  PGFEM_fprintf(out,"%10d\n",part1);
  PGFEM_fprintf(out,"Domain %d\n",myrank);
  PGFEM_fprintf(out,"coordinates\n");
  PGFEM_fprintf(out,"%10ld\n",nn);
  for (i=0;i<nn;i++) PGFEM_fprintf(out,"%10ld\n",i);
  
  /* X coordinate */
  for (i=0;i<nn;i++){
    switch(opts->analysis_type){
    default:
      PGFEM_fprintf(out,"%12.5e\n",node[i].x1_fd);
      break;
    case ELASTIC:
    case TP_ELASTO_PLASTIC:
      PGFEM_fprintf(out,"%12.5e\n",node[i].x1);
      break;
    }
  }/* end X */
  /* Y coordinate */
  for (i=0;i<nn;i++){
    switch(opts->analysis_type){
    default:
      PGFEM_fprintf(out,"%12.5e\n",node[i].x2_fd);
      break;
    case ELASTIC:
    case TP_ELASTO_PLASTIC:
      PGFEM_fprintf(out,"%12.5e\n",node[i].x2);
      break;
    }
  }/* end Y */
  /* Z coordinate */
  for (i=0;i<nn;i++){
    switch(opts->analysis_type){
    default:
      PGFEM_fprintf(out,"%12.5e\n",node[i].x3_fd);
      break;
    case ELASTIC:
    case TP_ELASTO_PLASTIC:
      PGFEM_fprintf(out,"%12.5e\n",node[i].x3);
      break;
    }
  }/* end Z */
  
  tetra4 = tetra10 = hexa8 = 0;
  for (i=0;i<ne;i++){
    if (elem[i].toe == 4) tetra4++;
    if (elem[i].toe == 8) hexa8++;
    if (elem[i].toe == 10) tetra10++;
  }
  
  if (tetra4 > 0){/* tetra4 */
    PGFEM_fprintf(out,"tetra4\n");
    PGFEM_fprintf(out,"%10ld\n",tetra4);
    for (i=0;i<ne;i++) if (elem[i].toe == 4) PGFEM_fprintf(out,"%10ld\n",i);
    for (i=0;i<ne;i++) {
      nod = aloc1l (elem[i].toe);
      /* conectivity table */
      elemnodes (i,elem[i].toe,nod,elem);
      for (j=0;j<elem[i].toe;j++) if (elem[i].toe == 4) PGFEM_fprintf(out,"%10ld",nod[j]+1);
      PGFEM_fprintf(out,"\n");
      dealoc1l (nod);
    }/* end i */
  }/* end tetra4 */
  
  if (hexa8 > 0){/* hexa8 */
    PGFEM_fprintf(out,"hexa8\n");
    PGFEM_fprintf(out,"%10ld\n",hexa8);
    for (i=0;i<ne;i++) if (elem[i].toe == 8) PGFEM_fprintf(out,"%10ld\n",i);
    for (i=0;i<ne;i++) {
      nod = aloc1l (elem[i].toe);
      /* conectivity table */
      elemnodes (i,elem[i].toe,nod,elem);
      for (j=0;j<elem[i].toe;j++){
	if (elem[i].toe == 8) PGFEM_fprintf(out,"%10ld",nod[j]+1);
      }/*end j */
      PGFEM_fprintf(out,"\n");
      dealoc1l (nod);
    }/* end i */
  }/* end hexa8 */
  if (tetra10 > 0){/* tetra10 */
    PGFEM_fprintf(out,"tetra10\n");
    PGFEM_fprintf(out,"%10ld\n",tetra10);
    for (i=0;i<ne;i++) if (elem[i].toe == 10) PGFEM_fprintf(out,"%10ld\n",i);
    for (i=0;i<ne;i++) {
      nod = aloc1l (elem[i].toe);
      /* conectivity table */
      elemnodes (i,elem[i].toe,nod,elem);
      for (j=0;j<elem[i].toe;j++){
	if (elem[i].toe == 10) PGFEM_fprintf(out,"%10ld",nod[j]+1);
      }/*end j */
      PGFEM_fprintf(out,"\n");
      dealoc1l (nod);
    }/* end i */
  }/* end tetra10 */
  
  /* COHESIVE ELEMENTS */
  if (opts->cohesive == 1){
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);
    PGFEM_fprintf(out,"Cohesive elements - D%d\n",myrank);
    PGFEM_fprintf(out,"coordinates\n");
    PGFEM_fprintf(out,"%10ld\n",ensight->ncn);
    for (i=0;i<ensight->ncn;i++) PGFEM_fprintf(out,"%10ld\n",ensight->Sm[i]);
    /* X coordinate */
    for (i=0;i<ensight->ncn;i++){
      switch(opts->analysis_type){
      default:
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(node[ensight->Sm[i]].x1_fd
				      + node[ensight->Sp[i]].x1_fd));
	break;
      case ELASTIC:
      case TP_ELASTO_PLASTIC:
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(node[ensight->Sm[i]].x1
				      + node[ensight->Sp[i]].x1));
	break;
      }
    }/* end X */
    /* Y coordinate */
    for (i=0;i<ensight->ncn;i++){
      switch(opts->analysis_type){
      default:
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(node[ensight->Sm[i]].x2_fd
				      + node[ensight->Sp[i]].x2_fd));
	break;
      case ELASTIC:
      case TP_ELASTO_PLASTIC:
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(node[ensight->Sm[i]].x2
				      + node[ensight->Sp[i]].x2));
	break;
      }
    }/* end Y */
    /* Z coordinate */
    for (i=0;i<ensight->ncn;i++){
      switch(opts->analysis_type){
      default:
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(node[ensight->Sm[i]].x3_fd
				      + node[ensight->Sp[i]].x3_fd));
	break;
      case ELASTIC:
      case TP_ELASTO_PLASTIC:
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(node[ensight->Sm[i]].x3
				      + node[ensight->Sp[i]].x3));
	break;
      }
    }/* end Z */
    
    tria3 = quad4 = 0;
    for (i=0;i<nce;i++){
      if (coel[i].toe == 6) tria3++;
      if (coel[i].toe == 8) quad4++;
    }
    
    if (tria3 > 0){/* tria3 */
      PGFEM_fprintf(out,"tria3\n");
      PGFEM_fprintf(out,"%10ld\n",tria3);
      for (i=0;i<nce;i++) if (coel[i].toe == 6) PGFEM_fprintf(out,"%10ld\n",i);
       for (i=0;i<nce;i++) {
	 /* conectivity table */
	 for (j=0;j<coel[i].toe/2;j++){
	   if (coel[i].toe == 6){
	     PGFEM_fprintf(out,"%10ld",ensight->No[coel[i].nod[j]]+1);
	   }
	 }
	 PGFEM_fprintf(out,"\n");
       }
    }
  }/* opts->cohesive == 1*/
  
  fclose (out);
  /* CLOSE GEO FILE */

  /* MAT */
  if (tim < nt) sprintf (name,"%s_%d.mat%ld",jmeno,myrank,tim);
  else          sprintf (name,"%s_%d.mat",jmeno,myrank);
  
  if ((out = fopen(name,"w")) == NULL ){
    PGFEM_printf("File is not possible to open, EnSight\n");
    PGFEM_Abort();
  }

  PGFEM_fprintf(out,"Material Number/Id File\n");
  PGFEM_fprintf(out,"part\n");
  PGFEM_fprintf(out,"%10d\n",part1);
  if (tetra4 > 0){/* tetra4 */
    PGFEM_fprintf(out,"tetra4\n");
    for (i=0;i<ne;i++) if (elem[i].toe == 4) PGFEM_fprintf(out,"%10ld\n",elem[i].pr);
  }
  if (hexa8 > 0){/* hexa8 */
    PGFEM_fprintf(out,"hexa8\n");
    for (i=0;i<ne;i++) if (elem[i].toe == 8) PGFEM_fprintf(out,"%10ld\n",elem[i].pr);
  }
  if (tetra10 > 0){/* tetra10 */
    PGFEM_fprintf(out,"tetra10\n");
    for (i=0;i<ne;i++) if (elem[i].toe == 10) PGFEM_fprintf(out,"%10ld\n",elem[i].pr);
  }
  if (opts->cohesive == 1){
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);
    if (tria3 > 0){/* tria3 */
      PGFEM_fprintf(out,"tria3\n");
      for (i=0;i<nce;i++) if (coel[i].toe == 6) PGFEM_fprintf(out,"%10ld\n",coel[i].pr);
    }
    if (quad4 > 0){/* quad4 */
      PGFEM_fprintf(out,"quad4\n");
      for (i=0;i<nce;i++) if (coel[i].toe == 8) PGFEM_fprintf(out,"%10ld\n",coel[i].pr);
    }
  }/* end coh */
  
  fclose (out);  
  /* CLOSE MATERIAL FILE */

  /* DISPLACEMENT */
  if (tim < nt) sprintf (name,"%s_%d.dis%ld",jmeno,myrank,tim);
  else          sprintf (name,"%s_%d.dis",jmeno,myrank);
  
  if ((out = fopen(name,"w")) == NULL ){
    PGFEM_printf("File is not possible to open, EnSight\n");
    PGFEM_Abort();
  }
  
  if (periodic == 1) PGFEM_fprintf(out,"Fluctuating Displacement\n"); else PGFEM_fprintf(out,"Displacement\n");
  PGFEM_fprintf(out,"part\n");
  PGFEM_fprintf(out,"%10d\n",part1);
  PGFEM_fprintf(out,"coordinates\n");
  /* X */
  for (i=0;i<nn;i++){
    if (node[i].id[0] == 0) PGFEM_fprintf(out,"%12.5e\n",0.0);
    if (node[i].id[0] >  0) PGFEM_fprintf(out,"%12.5e\n",r[node[i].id[0]-1]);
    if (node[i].id[0] <  0) PGFEM_fprintf(out,"%12.5e\n",sup->defl[abs(node[i].id[0])-1]);
  }/* end X */
  /* Y */
  for (i=0;i<nn;i++){
    if (node[i].id[1] == 0) PGFEM_fprintf(out,"%12.5e\n",0.0);
    if (node[i].id[1] >  0) PGFEM_fprintf(out,"%12.5e\n",r[node[i].id[1]-1]);
    if (node[i].id[1] <  0) PGFEM_fprintf(out,"%12.5e\n",sup->defl[abs(node[i].id[1])-1]);
  }/* end Y */
  /* Z */
  for (i=0;i<nn;i++){
    if (node[i].id[2] == 0) PGFEM_fprintf(out,"%12.5e\n",0.0);
    if (node[i].id[2] >  0) PGFEM_fprintf(out,"%12.5e\n",r[node[i].id[2]-1]);
    if (node[i].id[2] <  0) PGFEM_fprintf(out,"%12.5e\n",sup->defl[abs(node[i].id[2])-1]);
  }/* end Z */
  if (opts->cohesive == 1){
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);
    PGFEM_fprintf(out,"coordinates\n");
    /* X */
    for (i=0;i<ensight->ncn;i++){
      if (node[ensight->Sm[i]].id[0] == 0 && node[ensight->Sp[i]].id[0] == 0)
	PGFEM_fprintf(out,"%12.5e\n",0.0);
      if (node[ensight->Sm[i]].id[0] == 0 && node[ensight->Sp[i]].id[0] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(0.0 + r[node[ensight->Sp[i]].id[0]-1]));
      if (node[ensight->Sm[i]].id[0] == 0 && node[ensight->Sp[i]].id[0] < 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(0.0 + sup->defl[abs(node[ensight->Sp[i]].id[0])-1]));
      if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] == 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[0]-1] + 0.0));
      if (node[ensight->Sm[i]].id[0] <  0 && node[ensight->Sp[i]].id[0] == 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[0])-1] + 0.0));
      if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] <  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[0]-1] + sup->defl[abs(node[ensight->Sp[i]].id[0])-1]));
      if (node[ensight->Sm[i]].id[0] <  0 && node[ensight->Sp[i]].id[0] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[0])-1] + r[node[ensight->Sp[i]].id[0]-1]));
      if (node[ensight->Sm[i]].id[0] <  0 && node[ensight->Sp[i]].id[0] <  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[0])-1] + sup->defl[abs(node[ensight->Sp[i]].id[0])-1]));
      if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[0]-1]+ r[node[ensight->Sp[i]].id[0]-1]));
    }/* end X */
    /* Y */
    for (i=0;i<ensight->ncn;i++){
      if (node[ensight->Sm[i]].id[1] == 0 && node[ensight->Sp[i]].id[1] == 0)
	PGFEM_fprintf(out,"%12.5e\n",0.0);
      if (node[ensight->Sm[i]].id[1] == 0 && node[ensight->Sp[i]].id[1] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(0.0 + r[node[ensight->Sp[i]].id[1]-1]));
      if (node[ensight->Sm[i]].id[1] == 0 && node[ensight->Sp[i]].id[1] < 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(0.0 + sup->defl[abs(node[ensight->Sp[i]].id[1])-1]));
      if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] == 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[1]-1] + 0.0));
      if (node[ensight->Sm[i]].id[1] <  0 && node[ensight->Sp[i]].id[1] == 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[1])-1] + 0.0));
      if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] <  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[1]-1] + sup->defl[abs(node[ensight->Sp[i]].id[1])-1]));
      if (node[ensight->Sm[i]].id[1] <  0 && node[ensight->Sp[i]].id[1] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[1])-1] + r[node[ensight->Sp[i]].id[1]-1]));
      if (node[ensight->Sm[i]].id[1] <  0 && node[ensight->Sp[i]].id[1] <  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[1])-1] + sup->defl[abs(node[ensight->Sp[i]].id[1])-1]));
      if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[1]-1]+ r[node[ensight->Sp[i]].id[1]-1]));
    }/* end Y */
    /* Z */
    for (i=0;i<ensight->ncn;i++){
      if (node[ensight->Sm[i]].id[2] == 0 && node[ensight->Sp[i]].id[2] == 0)
	PGFEM_fprintf(out,"%12.5e\n",0.0);
      if (node[ensight->Sm[i]].id[2] == 0 && node[ensight->Sp[i]].id[2] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(0.0 + r[node[ensight->Sp[i]].id[2]-1]));
      if (node[ensight->Sm[i]].id[2] == 0 && node[ensight->Sp[i]].id[2] < 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(0.0 + sup->defl[abs(node[ensight->Sp[i]].id[2])-1]));
      if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] == 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[2]-1] + 0.0));
      if (node[ensight->Sm[i]].id[2] <  0 && node[ensight->Sp[i]].id[2] == 0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[2])-1] + 0.0));
      if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] <  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[2]-1] + sup->defl[abs(node[ensight->Sp[i]].id[2])-1]));
      if (node[ensight->Sm[i]].id[2] <  0 && node[ensight->Sp[i]].id[2] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[2])-1] + r[node[ensight->Sp[i]].id[2]-1]));
      if (node[ensight->Sm[i]].id[2] <  0 && node[ensight->Sp[i]].id[2] <  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(sup->defl[abs(node[ensight->Sm[i]].id[2])-1] + sup->defl[abs(node[ensight->Sp[i]].id[2])-1]));
      if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] >  0)
	PGFEM_fprintf(out,"%12.5e\n",1./2.*(r[node[ensight->Sm[i]].id[2]-1]+ r[node[ensight->Sp[i]].id[2]-1]));
    }/* end Z */
  }
  
  fclose (out);
  /* CLOSE DISPLACEMENT FILE */
  
  /* TOTAL DISPLACEMENT */
  if (periodic == 1){
    
    /* DISPLACEMENT */
    if (tim < nt) sprintf (name,"%s_%d.dit%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.dit",jmeno,myrank);
    
    if ((out = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out,"Total Displacement\n");
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part1);
    PGFEM_fprintf(out,"coordinates\n");
    /* X */
    for (i=0;i<nn;i++){
      /* Homegeneous displacement */
      X[0] = node[i].x1_fd; X[1] = node[i].x2_fd; X[2] = node[i].x3_fd;
      
      for (M=0;M<3;M++){
	xx[M] = 0.0; u[M] = 0.0;
	for (N=0;N<3;N++){
	  xx[M] += eps[0].F[M][N]*X[N];
	}
      }
      /* Fluctuating displacement */
      if (node[i].id[0] >  0) u[0] = r[node[i].id[0]-1]; 
      if (node[i].id[1] >  0) u[1] = r[node[i].id[1]-1]; 
      if (node[i].id[2] >  0) u[2] = r[node[i].id[2]-1]; 
      
      u[0] += xx[0] - X[0];
      u[1] += xx[1] - X[1];
      u[2] += xx[2] - X[2];
      
      PGFEM_fprintf(out,"%12.5e\n",u[0]);
    }/* end X */
    /* Y */
    for (i=0;i<nn;i++){
      /* Homegeneous displacement */
      X[0] = node[i].x1_fd; X[1] = node[i].x2_fd; X[2] = node[i].x3_fd;

      for (M=0;M<3;M++){
	xx[M] = 0.0; u[M] = 0.0;
	for (N=0;N<3;N++){
	  xx[M] += eps[0].F[M][N]*X[N];
	}
      }
      
      /* Fluctuating displacement */
      if (node[i].id[0] >  0) u[0] = r[node[i].id[0]-1]; 
      if (node[i].id[1] >  0) u[1] = r[node[i].id[1]-1]; 
      if (node[i].id[2] >  0) u[2] = r[node[i].id[2]-1]; 
      
      u[0] += xx[0] - X[0];
      u[1] += xx[1] - X[1];
      u[2] += xx[2] - X[2];
      
      PGFEM_fprintf(out,"%12.5e\n",u[1]);
    }/* end Y */
    /* Z */
    for (i=0;i<nn;i++){
      /* Homegeneous displacement */
      X[0] = node[i].x1_fd; X[1] = node[i].x2_fd; X[2] = node[i].x3_fd;
      
      for (M=0;M<3;M++){
	xx[M] = 0.0; u[M] = 0.0;
	for (N=0;N<3;N++){
	  xx[M] += eps[0].F[M][N]*X[N];
	}
      }

      /* Fluctuating displacement */
      if (node[i].id[0] >  0) u[0] = r[node[i].id[0]-1]; 
      if (node[i].id[1] >  0) u[1] = r[node[i].id[1]-1]; 
      if (node[i].id[2] >  0) u[2] = r[node[i].id[2]-1]; 
      
      u[0] += xx[0] - X[0];
      u[1] += xx[1] - X[1];
      u[2] += xx[2] - X[2];
      
      PGFEM_fprintf(out,"%12.5e\n",u[2]);
    }/* end Z */
    if (opts->cohesive == 1){
      PGFEM_fprintf(out,"part\n");
      PGFEM_fprintf(out,"%10d\n",part2);
      PGFEM_fprintf(out,"coordinates\n");
      /* X */
      for (i=0;i<ensight->ncn;i++){
	/* Homegeneous displacement */
	X[0] = 1./2.*(node[ensight->Sm[i]].x1_fd + node[ensight->Sp[i]].x1_fd);
	X[1] = 1./2.*(node[ensight->Sm[i]].x2_fd + node[ensight->Sp[i]].x2_fd);
	X[2] = 1./2.*(node[ensight->Sm[i]].x3_fd + node[ensight->Sp[i]].x3_fd);
	
	for (M=0;M<3;M++){
	  xx[M] = 0.0; u[M] = 0.0;
	  for (N=0;N<3;N++){
	    xx[M] += eps[0].F[M][N]*X[N];
	  }
	}
	
	/* Fluctuating displacement */
	if (node[ensight->Sm[i]].id[0] == 0 && node[ensight->Sp[i]].id[0] >  0) u[0] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[0]-1]);
	if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] == 0) u[0] = 1./2.*(r[node[ensight->Sm[i]].id[0]-1] + 0.0);
	if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] >  0) u[0] = 1./2.*(r[node[ensight->Sm[i]].id[0]-1]+ r[node[ensight->Sp[i]].id[0]-1]);
	
	if (node[ensight->Sm[i]].id[1] == 0 && node[ensight->Sp[i]].id[1] >  0) u[1] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[1]-1]);
	if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] == 0) u[1] = 1./2.*(r[node[ensight->Sm[i]].id[1]-1] + 0.0);
	if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] >  0) u[1] = 1./2.*(r[node[ensight->Sm[i]].id[1]-1]+ r[node[ensight->Sp[i]].id[1]-1]);
	
	if (node[ensight->Sm[i]].id[2] == 0 && node[ensight->Sp[i]].id[2] >  0) u[2] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[2]-1]);
	if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] == 0) u[2] = 1./2.*(r[node[ensight->Sm[i]].id[2]-1] + 0.0);
	if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] >  0) u[2] = 1./2.*(r[node[ensight->Sm[i]].id[2]-1]+ r[node[ensight->Sp[i]].id[2]-1]);
	
	u[0] += xx[0] - X[0];
	u[1] += xx[1] - X[1];
	u[2] += xx[2] - X[2];
	
	PGFEM_fprintf(out,"%12.5e\n",u[0]);
      }/* end X */
      /* Y */
      for (i=0;i<ensight->ncn;i++){
	/* Homegeneous displacement */
	X[0] = 1./2.*(node[ensight->Sm[i]].x1_fd + node[ensight->Sp[i]].x1_fd);
	X[1] = 1./2.*(node[ensight->Sm[i]].x2_fd + node[ensight->Sp[i]].x2_fd);
	X[2] = 1./2.*(node[ensight->Sm[i]].x3_fd + node[ensight->Sp[i]].x3_fd);

	for (M=0;M<3;M++){
	  xx[M] = 0.0; u[M] = 0.0;
	  for (N=0;N<3;N++){
	    xx[M] += eps[0].F[M][N]*X[N];
	  }
	}

	/* Fluctuating displacement */
	if (node[ensight->Sm[i]].id[0] == 0 && node[ensight->Sp[i]].id[0] >  0) u[0] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[0]-1]);
	if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] == 0) u[0] = 1./2.*(r[node[ensight->Sm[i]].id[0]-1] + 0.0);
	if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] >  0) u[0] = 1./2.*(r[node[ensight->Sm[i]].id[0]-1]+ r[node[ensight->Sp[i]].id[0]-1]);
	
	if (node[ensight->Sm[i]].id[1] == 0 && node[ensight->Sp[i]].id[1] >  0) u[1] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[1]-1]);
	if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] == 0) u[1] = 1./2.*(r[node[ensight->Sm[i]].id[1]-1] + 0.0);
	if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] >  0) u[1] = 1./2.*(r[node[ensight->Sm[i]].id[1]-1]+ r[node[ensight->Sp[i]].id[1]-1]);
	
	if (node[ensight->Sm[i]].id[2] == 0 && node[ensight->Sp[i]].id[2] >  0) u[2] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[2]-1]);
	if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] == 0) u[2] = 1./2.*(r[node[ensight->Sm[i]].id[2]-1] + 0.0);
	if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] >  0) u[2] = 1./2.*(r[node[ensight->Sm[i]].id[2]-1]+ r[node[ensight->Sp[i]].id[2]-1]);
	
	u[0] += xx[0] - X[0];
	u[1] += xx[1] - X[1];
	u[2] += xx[2] - X[2];
	
	PGFEM_fprintf(out,"%12.5e\n",u[1]);
      }/* end Y */
      /* Z */
      for (i=0;i<ensight->ncn;i++){
	/* Homegeneous displacement */
	X[0] = 1./2.*(node[ensight->Sm[i]].x1_fd + node[ensight->Sp[i]].x1_fd);
	X[1] = 1./2.*(node[ensight->Sm[i]].x2_fd + node[ensight->Sp[i]].x2_fd);
	X[2] = 1./2.*(node[ensight->Sm[i]].x3_fd + node[ensight->Sp[i]].x3_fd);
	
	for (M=0;M<3;M++){
	  xx[M] = 0.0; u[M] = 0.0;
	  for (N=0;N<3;N++){
	    xx[M] += eps[0].F[M][N]*X[N];
	  }
	}
	
	/* Fluctuating displacement */
	if (node[ensight->Sm[i]].id[0] == 0 && node[ensight->Sp[i]].id[0] >  0) u[0] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[0]-1]);
	if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] == 0) u[0] = 1./2.*(r[node[ensight->Sm[i]].id[0]-1] + 0.0);
	if (node[ensight->Sm[i]].id[0] >  0 && node[ensight->Sp[i]].id[0] >  0) u[0] = 1./2.*(r[node[ensight->Sm[i]].id[0]-1]+ r[node[ensight->Sp[i]].id[0]-1]);
	
	if (node[ensight->Sm[i]].id[1] == 0 && node[ensight->Sp[i]].id[1] >  0) u[1] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[1]-1]);
	if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] == 0) u[1] = 1./2.*(r[node[ensight->Sm[i]].id[1]-1] + 0.0);
	if (node[ensight->Sm[i]].id[1] >  0 && node[ensight->Sp[i]].id[1] >  0) u[1] = 1./2.*(r[node[ensight->Sm[i]].id[1]-1]+ r[node[ensight->Sp[i]].id[1]-1]);
	
	if (node[ensight->Sm[i]].id[2] == 0 && node[ensight->Sp[i]].id[2] >  0) u[2] = 1./2.*(0.0 + r[node[ensight->Sp[i]].id[2]-1]);
	if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] == 0) u[2] = 1./2.*(r[node[ensight->Sm[i]].id[2]-1] + 0.0);
	if (node[ensight->Sm[i]].id[2] >  0 && node[ensight->Sp[i]].id[2] >  0) u[2] = 1./2.*(r[node[ensight->Sm[i]].id[2]-1]+ r[node[ensight->Sp[i]].id[2]-1]);
	
	u[0] += xx[0] - X[0];
	u[1] += xx[1] - X[1];
	u[2] += xx[2] - X[2];
	
	PGFEM_fprintf(out,"%12.5e\n",u[2]);
      }/* end Z */
    }/* end cohesive */
    
    fclose (out);
  }/* end periodic */
  
  if (opts->analysis_type == STABILIZED
      || opts->analysis_type == MINI
      || opts->analysis_type == MINI_3F){/* Pressure field */
    
    if (tim < nt) sprintf (name,"%s_%d.pre%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.pre",jmeno,myrank);
    
    if ((out = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out,"Pressure\n");
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part1);
    PGFEM_fprintf(out,"coordinates\n");
    for (i=0;i<nn;i++) PGFEM_fprintf(out,"%12.5e\n",r[node[i].id[3]-1]);
  
    if (opts->cohesive == 1){
      PGFEM_fprintf(out,"part\n");
      PGFEM_fprintf(out,"%10d\n",part2);
      PGFEM_fprintf(out,"coordinates partial\n");
      PGFEM_fprintf(out,"%10d\n",0);
    }

    fclose (out);  
  }/* end pressure */
  
  /* STRESS + STRAIN ELEMENT */
  if (tim < nt) sprintf (name,"%s_%d.ste%ld",jmeno,myrank,tim);
  else          sprintf (name,"%s_%d.ste",jmeno,myrank);
  
  if ((out = fopen(name,"w")) == NULL ){
    PGFEM_printf("File is not possible to open, EnSight\n");
    PGFEM_Abort();
  }
  
  PGFEM_fprintf(out,"Stress per element\n");
  PGFEM_fprintf(out,"part\n");
  PGFEM_fprintf(out,"%10d\n",part1);
  
  if (tim < nt) sprintf (name,"%s_%d.sta%ld",jmeno,myrank,tim);
  else          sprintf (name,"%s_%d.sta",jmeno,myrank);
  
  if ((out1 = fopen(name,"w")) == NULL ){
    PGFEM_printf("File is not possible to open, EnSight\n");
    PGFEM_Abort();
  }
  
  PGFEM_fprintf(out1,"Strain per element\n");
  PGFEM_fprintf(out1,"part\n");
  PGFEM_fprintf(out1,"%10d\n",part1);

  if (tetra4 > 0){/* tetra4 */
    PGFEM_fprintf(out,"tetra4\n");
    PGFEM_fprintf(out1,"tetra4\n");
    for (i=0;i<ne;i++)
      if (elem[i].toe == 4) {
	/* S11 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[0]);
	/* E11 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[0]);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 4) {
	/* S22 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[1]);
	/* E22 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[1]);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 4) {
	/* S33 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[2]);
	/* E33 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[2]);
      }
    for (i=0;i<ne;i++)
      if (elem[i].toe == 4) {
	/* S12 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[5]);
	/* E12 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[5]/2.);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 4) {
	/* S13 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[4]);
	/* E13 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[4]/2.);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 4) {
	/* S23 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[3]);
	/* E23 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[3]/2.);
      }
  }
  if (hexa8 > 0){/* hexa8 */
    PGFEM_fprintf(out,"hexa8\n");
    PGFEM_fprintf(out1,"hexa8\n");
    for (i=0;i<ne;i++)
      if (elem[i].toe == 8) {
	/* S11 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[0]);
	/* E11 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[0]);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 8) {
	/* S22 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[1]);
	/* E22 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[1]);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 8) {
	/* S33 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[2]);
	/* E33 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[2]);
      }
    for (i=0;i<ne;i++)
      if (elem[i].toe == 8) {
	/* S12 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[5]);
	/* E12 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[5]/2.);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 8) {
	/* S13 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[4]);
	/* E13 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[4]/2.);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 8) {
	/* S23 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[3]);
	/* E23 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[3]/2.);
      }
  }
  if (tetra10 > 0){/* tetra10 */
    PGFEM_fprintf(out,"tetra10\n");
    PGFEM_fprintf(out1,"tetra10\n");
    for (i=0;i<ne;i++)
      if (elem[i].toe == 10) {
	/* S11 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[0]);
	/* E11 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[0]);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 10) {
	/* S22 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[1]);
	/* E22 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[1]);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 10) {
	/* S33 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[2]);
	/* E33 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[2]);
      }
    for (i=0;i<ne;i++)
      if (elem[i].toe == 10) {
	/* S12 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[5]);
	/* E12 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[5]/2.);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 10) {
	/* S13 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[4]);
	/* E13 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[4]/2.);
      }
    for (i=0;i<ne;i++) 
      if (elem[i].toe == 10) {
	/* S23 */ PGFEM_fprintf(out,"%12.5e\n",sig_e[i].el.o[3]);
	/* E23 */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].el.o[3]/2.);
      }
  }
  
  if (opts->cohesive == 1){
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part2);
  }
  
  fclose (out); fclose (out1);
  /* CLOSE STRESS + STRAIN ELEMENT FILES */
  
  /* STRESS NODES */
  if (gr4 == 0){
    if (tim < nt) sprintf (name,"%s_%d.stn%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.stn",jmeno,myrank);
    
    if ((out = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out,"Stress per nodes (Z-Z smoothing)\n");
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part1);
    PGFEM_fprintf(out,"coordinates\n");
    /* S11 */ for (i=0;i<nn;i++) PGFEM_fprintf(out,"%12.5e\n",sig_n[i].el.o[0]);
    /* S22 */ for (i=0;i<nn;i++) PGFEM_fprintf(out,"%12.5e\n",sig_n[i].el.o[1]);
    /* S33 */ for (i=0;i<nn;i++) PGFEM_fprintf(out,"%12.5e\n",sig_n[i].el.o[2]);
    /* S12 */ for (i=0;i<nn;i++) PGFEM_fprintf(out,"%12.5e\n",sig_n[i].el.o[5]);
    /* S13 */ for (i=0;i<nn;i++) PGFEM_fprintf(out,"%12.5e\n",sig_n[i].el.o[4]);
    /* S23 */ for (i=0;i<nn;i++) PGFEM_fprintf(out,"%12.5e\n",sig_n[i].el.o[3]);
    
    if (opts->cohesive == 1){
      PGFEM_fprintf(out,"part\n");
      PGFEM_fprintf(out,"%10d\n",part2);
      PGFEM_fprintf(out,"coordinates partial\n");
      PGFEM_fprintf(out,"%10d\n",0);
    }
    
    fclose (out);  
    /* CLOSE STRESS NODES FILE */
  }
  
  /* EFFECTIVE STRAIN */
  if (tim < nt) sprintf (name,"%s_%d.est%ld",jmeno,myrank,tim);
  else          sprintf (name,"%s_%d.est",jmeno,myrank);
  
  if ((out = fopen(name,"w")) == NULL ){
    PGFEM_printf("File is not possible to open, EnSight\n");
    PGFEM_Abort();
  }
  
  PGFEM_fprintf(out,"Effective Strain per element\n");
  PGFEM_fprintf(out,"part\n");
  PGFEM_fprintf(out,"%10d\n",part1);
  if (tetra4 > 0){/* tetra4 */
    PGFEM_fprintf(out,"tetra4\n");
    /* Eeq */ for (i=0;i<ne;i++) if (elem[i].toe == 4) PGFEM_fprintf(out,"%12.5e\n",eps[i].el.eq);
  }
  if (hexa8 > 0){/* hexa8 */
    PGFEM_fprintf(out,"hexa8\n");
    /* Eeq */ for (i=0;i<ne;i++) if (elem[i].toe == 8) PGFEM_fprintf(out,"%12.5e\n",eps[i].el.eq);
  }
  if (tetra10 > 0){/* tetra10 */
    PGFEM_fprintf(out,"tetra10\n");
    /* Eeq */ for (i=0;i<ne;i++) if (elem[i].toe == 10) PGFEM_fprintf(out,"%12.5e\n",eps[i].el.eq);
  }

  if (opts->cohesive == 1){
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);
  }

  fclose (out);  
  /* CLOSE EFFECTIVE STRAIN FILE */
  
  if (opts->analysis_type == FS_CRPL){/* Plastic strains */
    /* EFFECTIVE PLASTIC ALMANSI + EFFECTIVE STRAIN */
    if (tim < nt) sprintf (name,"%s_%d.pla%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.pla",jmeno,myrank);
    
    if ((out = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out,"Effective Plastic Almansi Strain\n");
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part1);
    
    /* EFFECTIVE PLASTIC STRAIN */
    if (tim < nt) sprintf (name,"%s_%d.ple%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.ple",jmeno,myrank);
    
    if ((out1 = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out1,"Effective Plastic Strain || Int_0^t (sqrt(2/3*Dp.Dp)) \n");
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part1);

    if (tetra4 > 0){/* tetra4 */
      PGFEM_fprintf(out,"tetra4\n");
      PGFEM_fprintf(out1,"tetra4\n");
      for (i=0;i<ne;i++) 
	if (elem[i].toe == 4) {
	  /* Eeq */ PGFEM_fprintf(out,"%12.5e\n",eps[i].pl.eq[1]);
	  /* Eeq */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].pl.eq[0]);
	}
    }
    if (hexa8 > 0){/* hexa8 */
      PGFEM_fprintf(out,"hexa8\n");
      PGFEM_fprintf(out1,"hexa8\n");
      for (i=0;i<ne;i++) 
	if (elem[i].toe == 8) {
	  /* Eeq */ PGFEM_fprintf(out,"%12.5e\n",eps[i].pl.eq[1]);
	  /* Eeq */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].pl.eq[0]);
	}
    }
    if (tetra10 > 0){/* tetra10 */
      PGFEM_fprintf(out,"tetra10\n");
      PGFEM_fprintf(out1,"tetra10\n");
      for (i=0;i<ne;i++) 
	if (elem[i].toe == 10) {
	  /* Eeq */ PGFEM_fprintf(out,"%12.5e\n",eps[i].pl.eq[1]);
	  /* Eeq */ PGFEM_fprintf(out1,"%12.5e\n",eps[i].pl.eq[0]);
	}
    }

    if (opts->cohesive == 1){
      PGFEM_fprintf(out,"part\n");
      PGFEM_fprintf(out,"%10d\n",part2);
      PGFEM_fprintf(out1,"part\n");
      PGFEM_fprintf(out1,"%10d\n",part2);
    }
    
    fclose (out); fclose (out1);
    /* CLOSE PLASTIC ALMANSI + EFFECTIVE STRAIN FILE */
  }/* end analysis == FS_CRPL */
  
  if (opts->cohesive == 1){
    /* OPENING DISPLACEMENT + TRACTIONS */
    if (tim < nt) sprintf (name,"%s_%d.odi%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.odi",jmeno,myrank);
    
    if ((out = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out,"Cohesive Opening Displacement\n");
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part1);
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);
    
    if (tim < nt) sprintf (name,"%s_%d.tra%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.tra",jmeno,myrank);
    
    if ((out1 = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out1,"Cohesive Tractions\n");
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part1);
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part2);
    
    if (tria3 > 0){/* tria3 */
      PGFEM_fprintf(out,"tria3\n");
      PGFEM_fprintf(out1,"tria3\n");
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 6) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xi[0]);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].ti[0]);
	}
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 6) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xi[1]);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].ti[1]);
	}
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 6) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xi[2]);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].ti[2]);
	}
    }
    if (quad4 > 0){/* quad4 */
      PGFEM_fprintf(out,"quad4\n");
      PGFEM_fprintf(out1,"quad4\n");
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 8) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xi[0]);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].ti[0]);
	}
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 8) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xi[1]);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].ti[1]);
	}
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 8) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xi[2]);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].ti[2]);
	}
    }
    
    fclose (out); fclose (out1);
    /* CLOSE OPENING DISPLACEMENT + TRACTION FILES */
  
    /* OPENING DISPLACEMENT COMPONENT FILES */
    if (tim < nt) sprintf (name,"%s_%d.ode%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.ode",jmeno,myrank);
    
    if ((out = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out,"Effective opening displacement\n");
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part1);
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);

    if (tim < nt) sprintf (name,"%s_%d.odn%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.odn",jmeno,myrank);
    
    if ((out1 = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out1,"Normal opening displacement\n");
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part1);
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part2);
    
    if (tim < nt) sprintf (name,"%s_%d.ods%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.ods",jmeno,myrank);
    
    if ((out2 = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out2,"Shear opening displacement\n");
    PGFEM_fprintf(out2,"part\n");
    PGFEM_fprintf(out2,"%10d\n",part1);
    PGFEM_fprintf(out2,"part\n");
    PGFEM_fprintf(out2,"%10d\n",part2);
    
    if (tria3 > 0){/* tria3 */
      PGFEM_fprintf(out,"tria3\n");
      PGFEM_fprintf(out1,"tria3\n");
      PGFEM_fprintf(out2,"tria3\n");
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 6) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xxi);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].Xn);
	  PGFEM_fprintf(out2,"%12.5e\n",coel[i].Xs);
	}
    }
    if (quad4 > 0){/* quad4 */
      PGFEM_fprintf(out,"quad4\n");
      PGFEM_fprintf(out1,"quad4\n");
      PGFEM_fprintf(out2,"quad4\n");
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 8) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].Xxi);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].Xn);
	  PGFEM_fprintf(out2,"%12.5e\n",coel[i].Xs);
	}
    }
    
    fclose (out); fclose (out1); fclose (out2);
    /* CLOSE OPENING DISPLACEMENT COMPONENT FILES */

    /* OPENING TRACTION COMPONENT FILES */
    if (tim < nt) sprintf (name,"%s_%d.tre%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.tre",jmeno,myrank);
    
    if ((out = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out,"Effective traction\n");
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part1);
    PGFEM_fprintf(out,"part\n");
    PGFEM_fprintf(out,"%10d\n",part2);
    
    if (tim < nt) sprintf (name,"%s_%d.trn%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.trn",jmeno,myrank);
    
    if ((out1 = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out1,"Normal traction\n");
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part1);
    PGFEM_fprintf(out1,"part\n");
    PGFEM_fprintf(out1,"%10d\n",part2);
    
    if (tim < nt) sprintf (name,"%s_%d.trs%ld",jmeno,myrank,tim);
    else          sprintf (name,"%s_%d.trs",jmeno,myrank);
    
    if ((out2 = fopen(name,"w")) == NULL ){
      PGFEM_printf("File is not possible to open, EnSight\n");
      PGFEM_Abort();
    }
    
    PGFEM_fprintf(out2,"Shear traction\n");
    PGFEM_fprintf(out2,"part\n");
    PGFEM_fprintf(out2,"%10d\n",part1);
    PGFEM_fprintf(out2,"part\n");
    PGFEM_fprintf(out2,"%10d\n",part2);
    
    if (tria3 > 0){/* tria3 */
      PGFEM_fprintf(out,"tria3\n");
      PGFEM_fprintf(out1,"tria3\n");
      PGFEM_fprintf(out2,"tria3\n");
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 6) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].txi);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].tn);
	  PGFEM_fprintf(out2,"%12.5e\n",coel[i].ts);
	}
    }
    if (quad4 > 0){/* quad4 */
      PGFEM_fprintf(out,"quad4\n");
      PGFEM_fprintf(out1,"quad4\n");
      PGFEM_fprintf(out2,"quad4\n");
      for (i=0;i<nce;i++) 
	if (coel[i].toe == 8) {
	  PGFEM_fprintf(out,"%12.5e\n",coel[i].txi);
	  PGFEM_fprintf(out1,"%12.5e\n",coel[i].tn);
	  PGFEM_fprintf(out2,"%12.5e\n",coel[i].ts);
	}
    }
    
    fclose (out); fclose (out1); fclose (out2);
    /* CLOSE OPENING TRACTION COMPONENT FILES */
    
  }/* end opts->cohesive == 1 */
  
  dealoc1 (X); dealoc1 (xx); dealoc1 (u);
}

void ASCII_output(const PGFem3D_opt *opts,
		  MPI_Comm comm,
		  long tim,
		  double *times,
		  long  Gnn,
		  long nn,
		  long ne,
		  long nce,
		  long ndofd,
		  long *DomDof,
		  int *Ap,
		  long FNR,
		  double lm,
		  double pores,
		  double VVolume,
		  NODE *node,
		  ELEMENT *elem,
		  SUPP sup,
		  double *r,
		  EPS *eps,
		  SIG *sig_e,
		  SIG *sig_n,
		  COEL *coel)
{
  int myrank = 0;
  MPI_Comm_rank(comm,&myrank);
  struct rusage usage;

  /* compute filename and open file */
  char *filename = PGFEM_calloc(500,sizeof(char));
  sprintf (filename,"%s/ASCII/STEP_%.5ld",opts->opath,tim);
  if(make_path(filename,DIR_MODE) != 0){
    /* could not create, default dir name current directory */
    sprintf(filename,".");
  }
  sprintf (filename,"%s/%s_%d.out%ld",
	   filename,opts->ofname,myrank,tim);
  FILE *out = PGFEM_fopen(filename,"w");
  free(filename);

  /*=== START WRITING ===*/	
  logo (out);
  PGFEM_fprintf (out,"\n");
  PGFEM_fprintf (out,"%s Step - %ld | Time - %f || Stab = %12.12f\n",
	   "FINITE STRAINS: STABILIZED FORMULATION         :",
	   tim,times[tim+1],opts->stab);
  if (opts->analysis_type == FS_CRPL) {
    PGFEM_fprintf (out,"FINITE DEFORMA. + CRYSTAL PLASTICITY ANALYSIS : ");
  } else if (opts->analysis_type == FINITE_STRAIN){
    PGFEM_fprintf (out,"FINITE DEFORMATION                            : ");
  }
  PGFEM_fprintf (out,"Number of total nodes                          : %ld\n",
	   Gnn);
  PGFEM_fprintf (out,"Number of nodes on domain                      : %ld\n",
	   nn);
  PGFEM_fprintf (out,"Number of elements on domain                   : %ld\n",
	   ne);
  if (opts->cohesive == 1){
    PGFEM_fprintf (out,"%s %ld\n",
	     "Number of cohesive elements on domain          :",
	     nce);
  }
  PGFEM_fprintf (out,"Number of equations                            : %ld\n",
	   ndofd);
  PGFEM_fprintf (out,"Number of global degrees of freedom on domain  : %ld\n",
	   DomDof[myrank]);
  PGFEM_fprintf (out,"%s %d\n",
	   "Number of elements in the matrix - SPARSE      :",
	   Ap[DomDof[myrank]]);
  if (FNR == 2 || FNR == 3){
    PGFEM_fprintf (out,
	     "Load multiplier level                          : %12.12f\n",
	     lm);
  }
  if (opts->cohesive == 1){
    PGFEM_fprintf (out,"%s %12.12f || %12.12f\n",
	     "Volume of voids                                :",
	     pores,pores/VVolume);
  }
  PGFEM_fprintf (out,"\n");
  getrusage (RUSAGE_SELF,&usage);
  PGFEM_fprintf (out,"%s %ld.%ld, User %ld.%ld\n",
	   "Time of solution of the system of equations  -  System",
	   usage.ru_stime.tv_sec,usage.ru_stime.tv_usec,
	   usage.ru_utime.tv_sec,usage.ru_utime.tv_usec);
  deform (out,node,elem,nn,ne,ndim,sup,r);

  if(opts->analysis_type == DISP){
    damage_out(out,ne,elem,eps);
  }
  /* Print stress to output file */
  stress_out (out,ne,nn,elem,sig_e,sig_n,opts->smoothing);
  /* Print strain to output file */
  strain_out (out,ne,elem,eps,opts);
  /* Print Fe */
  deform_grad_out (out,ne,elem,eps);
  /* Print cohesive elements */
  if (opts->cohesive == 1){
    cohesive_out (out,nce,coel);
  }
	
  fclose (out); 
}
