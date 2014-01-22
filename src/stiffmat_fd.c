#include "stiffmat_fd.h"

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef GET_NDOF_ON_ELEM_H
#include "get_ndof_on_elem.h"
#endif

#ifndef GET_DOF_IDS_ON_ELEM_H
#include "get_dof_ids_on_elem.h"
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef PLOC_SPARSE_H
#include "PLoc_Sparse.h"
#endif

#ifndef STABILIZED_H
#include "stabilized.h"
#endif

#ifndef STIFFMATEL_FD_H
#include "stiffmatel_fd.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef MINI_ELEMENT_H
#include "MINI_element.h"
#endif

#ifndef MINI_3F_ELEMENT_H
#include "MINI_3f_element.h"
#endif

#ifndef DISP_BASED_ELEM_H
#include "displacement_based_element.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

static const int periodic = 0;

/* This function may not be used outside this file */
static int el_stiffmat (int i, /* Element ID */
			double **Lk,
			BSspmat *K,
			int *Ap,
			int *Ai,
			long ndofn,
			ELEMENT *elem,
			NODE *node,
			HOMMAT *hommat,
			MATGEOM matgeom,
			SIG *sig,
			EPS *eps,
			double *d_r,
			double *r,
			long npres,
			SUPP sup,
			long iter,
			double nor_min,
			double dt,
			CRPL *crpl,
			double stab,
			long FNR,
			double lm,
			double *f_u,
			int myrank,
			int nproc,
			long GDof,
			COMMUN comm,
			int *Ddof,
			int interior,
			const int analysis,
			PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  int err = 0;
  long j,l,nne,ndofe,*cnL,*cnG,*nod,II;
  double *lk,*x,*y,*z,*r_e,*sup_def,*fe;
  long kk;
  
  /* Number of element nodes */
  nne = elem[i].toe;
    
  /* Nodes on element */
  nod = aloc1l (nne);
  elemnodes (i,nne,nod,elem);
    
  /* Element Dof */
  ndofe = get_ndof_on_elem_nodes(nne,nod,node);

  /* allocation */
  cnL = aloc1l (ndofe);
  cnG = aloc1l (ndofe);
  lk = aloc1 ((ndofe*ndofe)); 

  const int nne_t = nne + elem[i].n_bub;

  if (analysis == MINI
      || analysis == MINI_3F){ /* P1+B/P1 element */
    x = aloc1 (nne_t);
    y = aloc1 (nne_t);
    z = aloc1 (nne_t);
  } else {
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);
  }
  r_e = aloc1 (ndofe);
  fe = aloc1 (ndofe);
  if(sup->npd>0){
    sup_def = aloc1(sup->npd);
  } else {
    sup_def = NULL;
  }
    
  /* Coordinates of nodes */
  switch(analysis){
  case DISP:
    nodecoord_total (nne,nod,node,x,y,z);
    break;
  default:
    nodecoord_updated (nne,nod,node,x,y,z);
    break;
  }

  /* if P1+B/P1, get element centroid coords */
  if (analysis == MINI || analysis == MINI_3F){
    element_center(nne,x,y,z);
  }
    
  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cnL);
  get_dof_ids_on_elem_nodes(1,nne,ndofn,nod,node,cnG);
    
  /* deformation on element */
  if (iter == 0) {
    for (j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }
  }

  def_elem (cnL,ndofe,d_r,elem,node,r_e,sup,0);
  if (iter == 0){
    for (j=0;j<sup->npd;j++)
      sup->defl_d[j] = sup_def[j];
   }
    
  nulld (lk,ndofe*ndofe);  
  switch(analysis){
  case STABILIZED:
    err = stiffmatel_st (i,ndofn,nne,x,y,z,elem,hommat,nod,node,sig,eps,
			 r_e,npres,nor_min,lk,dt,stab,FNR,lm,fe);
    break;
  case MINI:
    err = MINI_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			   hommat,nod,node,eps,sig,r_e);
    break;
  case MINI_3F:
    err = MINI_3f_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			      hommat,nod,node,eps,sig,r_e);
    break;
  case DISP:
    {
      /* Get TOTAL deformation on element; r_e already contains
	 INCREMENT of deformation, add the deformation from previous. */
      double *r_en;
      r_en = aloc1(ndofe);
      def_elem (cnL,ndofe,r,elem,node,r_en,sup,1);
      vvplus(r_e,r_en,ndofe);
      err = DISP_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			     hommat,nod,node,eps,sig,sup,r_e);
      free(r_en);
    } 
    break;
  default:
    err = stiffmatel_fd (i,ndofn,nne,nod,x,y,z,elem,matgeom,
			 hommat,node,sig,eps,r_e,npres,
			 nor_min,lk,dt,crpl,FNR,lm,fe,analysis);
    break;
  } /* switch (analysis) */

  if (PFEM_DEBUG){
    char filename[50];
    switch(analysis){
    case STABILIZED:
      sprintf(filename,"stab_stiff_%d.log",myrank);
      break;
    case MINI:
      sprintf(filename,"MINI_stiff_%d.log",myrank);
      break;
    case MINI_3F:
      sprintf(filename,"MINI_3f_stiff_%d.log",myrank);
      break;
    default:
      sprintf(filename,"stiff_%d.log",myrank);
      break;
    }

    FILE *output;
    output = fopen(filename,"a");
    print_array_d(output,lk,ndofe*ndofe,ndofe,ndofe);
    fclose(output);
  }

  /* Localization of TANGENTIAL LOAD VECTOR */
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    for (l=0;l<nne;l++){
      for (kk=0;kk<node[nod[l]].ndofn;kk++){
	II = node[nod[l]].id[kk]-1;
	if (II < 0)  continue;
	f_u[II] += fe[l*node[nod[l]].ndofn+kk];
      }/*end l */
    }/*end kk */
  }/* end periodic */

  /* Assembly */
  PLoc_Sparse (K,Lk,lk,Ai,Ap,cnL,cnG,ndofe,Ddof,GDof,
	       myrank,nproc,comm,interior,PGFEM_hypre,analysis);

  /*  dealocation  */
  free (cnL);
  free (cnG);
  free (nod);
  free (lk);
  free (x);
  free (y);
  free (z);
  free (fe);
  free (r_e);
  free (sup_def);

  return err;

} /* ELEMENT STIFFNESS */

/* This function may not be used outside of this file */
static void coel_stiffmat(int i, /* coel ID */
			  double **Lk,
			  BSspmat *K,
			  int *Ap,
			  int *Ai,
			  long ndofc,
			  ELEMENT *elem,
			  NODE *node,
			  EPS *eps,
			  double *d_r,
			  double *r,
			  long npres,
			  SUPP sup,
			  long iter,
			  double nor_min,
			  double dt,
			  CRPL *crpl,
			  double stab,
			  COEL *coel,
			  long FNR,
			  double lm,
			  double *f_u,
			  int myrank,
			  int nproc,
			  long *DomDof,
			  long GDof,
			  COMMUN comm,
			  int *Ddof,
			  int interior,
			  const int analysis,
			  PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  long j,l,nne,ndofe,*cnL,*cnG,*nod,P,R,II;
  double *lk,*x,*y,*z,*r_e,*sup_def,*fe, *X, *Y;
  long kk;

  nne = coel[i].toe/2;
  ndofe = coel[i].toe*ndofc;
      
  /* Alocation */
  nod = aloc1l (coel[i].toe);
  lk = aloc1 (ndofe*ndofe);
  x = aloc1 (coel[i].toe); 
  y = aloc1 (coel[i].toe);
  z = aloc1 (coel[i].toe);
  r_e = aloc1 (ndofe);
  fe = aloc1 (ndofe);
  cnL = aloc1l (ndofe);
  cnG = aloc1l (ndofe);

  if(sup->npd > 0){
    sup_def = aloc1(sup->npd);
  } else {
    sup_def = NULL;
  }

  /* Element Node */
  for (j=0;j<coel[i].toe;j++)
    nod[j] = coel[i].nod[j];

  /* Coordinates */
  nodecoord_updated (coel[i].toe,nod,node,x,y,z);
      
  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nod,node,cnL);
  get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nod,node,cnG);

  /* deformation on element */
  if (iter == 0){
    for (j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }
  }

  def_elem (cnL,ndofe,d_r,elem,node,r_e,sup,0);
  if (iter == 0)
    for (j=0;j<sup->npd;j++)
      sup->defl_d[j] = sup_def[j];
      
  if (periodic == 1){/* Periodic */

    X = aloc1 (3);
    Y = aloc1 (3);
	
    for (j=0;j<coel[i].toe;j++){
      X[0] = x[j];
      X[1] = y[j];
      X[2] = z[j];
      for (P=0;P<3;P++){
	Y[P] = 0.0;
	for (R=0;R<3;R++){
	  Y[P] += eps[0].F[P][R]*X[R];
	}
      }
      for (P=0;P<3;P++){
	II = node[nod[j]].id[P];
	    
	if (P == 0)
	  x[j] = Y[P];
	if (P == 1)
	  y[j] = Y[P];
	if (P == 2)
	  z[j] = Y[P];
	    
	if (II > 0){
	  if (P == 0)
	    x[j] += r[II-1];
	  if (P == 1)
	    y[j] += r[II-1];
	  if (P == 2)
	    z[j] += r[II-1];
	}
      }/* P < 3 */
    }/* j < coel[i].toe */

    dealoc1 (X);
    dealoc1 (Y);
	
  }/* end periodic */

  /**************************/
  /*     FOR PERIODIC       */
  /*                        */
  /* xn+1 = Fn+1*X + un + u */
  /*    x = Fn+1*X + un     */
  /**************************/

  nulld (lk,ndofe*ndofe); 
  stiff_mat_coh (i,ndofc,nne,nod,x,y,z,coel,r_e,lk,
		 nor_min,eps,FNR,lm,fe,myrank);
      
  /* Assembly */
  PLoc_Sparse (K,Lk,lk,Ai,Ap,cnL,cnG,ndofe,Ddof,
	       GDof,myrank,nproc,comm,interior,PGFEM_hypre,analysis);

  /* Localization of TANGENTIAL LOAD VECTOR */
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    for (l=0;l<coel[i].toe;l++){
      for (kk=0;kk<ndofc;kk++){
	II = node[nod[l]].id[kk]-1;
	if (II < 0)  continue;
	f_u[II] += fe[l*ndofc+kk];
      }/*end l */
    }/*end kk */
  }/* end periodic */

  /*  dealocation  */
  dealoc1l (cnL);
  dealoc1l (cnG);
  dealoc1l (nod);
  dealoc1 (lk); 
  dealoc1 (x);
  dealoc1(y); 
  dealoc1 (z);
  dealoc1 (r_e);
  dealoc1 (fe);
  free(sup_def);

} /* COHESIVE ELEMENT STIFFNESS */

static int bnd_el_stiffmat(int belem_id,
			   double **Lk,
			   int *Ap,
			   int *Ai,
			   long ndofn,
			   ELEMENT *elem,
			   BOUNDING_ELEMENT *b_elems,
			   NODE *node,
			   HOMMAT *hommat,
			   MATGEOM matgeom,
			   SIG *sig,
			   EPS *eps,
			   double *d_r,
			   double *r,
			   long npres,
			   SUPP sup,
			   long iter,
			   double nor_min,
			   double dt,
			   CRPL *crpl,
			   double stab,
			   double FNR,
			   double lm,
			   double *f_u,
			   int myrank,
			   int nproc,
			   long GDof,
			   COMMUN comm,
			   int *Ddof,
			   int interior,
			   const int analysis,
			   PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  int err = 0;
  const BOUNDING_ELEMENT *ptr_be = &b_elems[belem_id];
  const ELEMENT *ptr_ve = &elem[ptr_be->vol_elem_id];
  const long *ptr_vnodes = ptr_ve->nod;
  const int nn_ve = ptr_ve->toe;

  /* get coordinated for BOUNDING ELEMENT */
  double *x = aloc1(nn_ve);
  double *y = aloc1(nn_ve);
  double *z = aloc1(nn_ve);
  switch(analysis){
  case DISP:
    nodecoord_total(nn_ve,ptr_vnodes,node,x,y,z);
    break;
  default:
    nodecoord_updated(nn_ve,ptr_vnodes,node,x,y,z);
    break;
  }

  /* get the local and global dof id's */
  double ndof_ve = get_ndof_on_bnd_elem(node,ptr_be,elem);

  long *cn_ve = aloc1l(ndof_ve);
  long *Gcn_ve = aloc1l(ndof_ve);

  get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn_ve);
  get_dof_ids_on_bnd_elem(1,ndofn,node,ptr_be,elem,Gcn_ve);

  /* compute the deformation on the element */
  double *v_disp = aloc1(ndof_ve);

  if(iter == 0){
    /* on iter == 0, null increment of deflection */
    double *sup_def = aloc1(sup->npd);
    for (int j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }

    def_elem(cn_ve,ndof_ve,d_r,NULL,NULL,v_disp,sup,0);

    for (int j=0;j<sup->npd;j++){
      sup->defl_d[j] = sup_def[j];
    }
    free(sup_def);
  } else {
    def_elem(cn_ve,ndof_ve,d_r,NULL,NULL,v_disp,sup,0);
  }

  if(analysis == DISP){ /* TOTAL LAGRANGIAN formulation */
    double *ve_n = aloc1(ndof_ve);
    def_elem(cn_ve,ndof_ve,r,NULL,NULL,ve_n,sup,1);
    vvplus(v_disp,ve_n,ndof_ve);
    free(ve_n);
  }

  /* compute the local stiffness matrix */
  double *lk = aloc1(ndof_ve*ndof_ve);
  if(analysis == DISP){
    err += DISP_stiffmat_bnd_el(lk,belem_id,ndofn,ndof_ve,
				x,y,z,b_elems,elem,hommat,node,eps,
				sig,sup,v_disp);
  } else {
    /* Not implemented, do nothing */
  }

  /* only assemble to global stiffness if no error */
  if(err == 0){
    /* PLoc_Sparse_rec(Lk,lk,Ai,Ap,Gcn_be,Gcn_ve,ndof_be,ndof_ve,Ddof, */
    /* 		   GDof,myrank,nproc,comm,interior); */
    PLoc_Sparse(NULL,Lk,lk,Ai,Ap,cn_ve,Gcn_ve,ndof_ve,Ddof,
		GDof,myrank,nproc,comm,interior,PGFEM_hypre,analysis);
  }


  free(x);
  free(y);
  free(z);

  free(cn_ve);
  free(Gcn_ve);

  free(v_disp);
  free(lk);

  return err;
} /* Bounding element stiffnes matrix */

/* This is the re-written function which computes elem stiffness on
   boundaries first, then interior elem stiffnesses before
   assembly. */
int stiffmat_fd (BSspmat *K,
		 int *Ap,
		 int *Ai,
		 long ne,
		 int n_be,
		 long ndofn,
		 ELEMENT *elem,
		 BOUNDING_ELEMENT *b_elems,
		 long nbndel,
		 long *bndel,
		 NODE *node,
		 HOMMAT *hommat,
		 MATGEOM matgeom,
		 SIG *sig,
		 EPS *eps,
		 double *d_r,
		 double *r,
		 long npres,
		 SUPP sup,
		 long iter,
		 double nor_min,
		 double dt,
		 CRPL *crpl,
		 double stab,
		 long nce,
		 COEL *coel,
		 long FNR,
		 double lm,
		 double *f_u,
		 int myrank,
		 int nproc,
		 long *DomDof,
		 long GDof,
		 COMMUN comm,
		 MPI_Comm mpi_comm,
		 PGFEM_HYPRE_solve_info *PGFEM_hypre,
		 const PGFem3D_opt *opts)
{
  int err = 0;
  long i,j,k,ndofc;
  int *Ddof;
  long KK;
  double **Lk,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;

  /* interior element counters */
  int idx = 0;
  int skip = 0;

  if(opts->solverpackage == 0){
    PGFEM_printf("BlockSolve no longer supported\n");
    abort();
  }

  err += init_and_post_stiffmat_comm(&Lk,&recieve,&req_r,&sta_r,
				     mpi_comm,comm);

  /* Lk = (double**) PGFEM_calloc (nproc,sizeof(double*)); */
  /* for (i=0;i<nproc;i++) { */
  /*   if (myrank == i || comm->S[i] == 0) */
  /*     k = 1; */
  /*   else  */
  /*     k = comm->AS[i]; */
  /*   Lk[i] = (double*) PGFEM_calloc (k,sizeof(double)); */
  /* } */
  /* if (Lk == NULL){ */
  /*   PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__); */
  /*   fflush(stdout);  */
  /*   PGFEM_Comm_code_abort (mpi_comm,i); */
  /* } */
  
  /* /\* Allocate recieve *\/ */
  /* recieve = (double**) PGFEM_calloc (nproc,sizeof(double*)); */
  /* for (i=0;i<nproc;i++) { */
  /*   if (comm->AR[i] == 0) */
  /*     KK = 1; */
  /*   else */
  /*     KK = comm->AR[i]; */
  /*   recieve[i] = (double*) PGFEM_calloc (KK,sizeof(double)); */
  /* } */
  
  /* /\* Allocate request fields *\/ */
  /* if (comm->Nr == 0) */
  /*   KK = 1; */
  /* else */
  /*   KK = comm->Nr; */
  /* sta_r = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status)); */
  /* req_r = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request)); */
  
  /* /\* Receive data *\/ */
  /* for (i=0;i<comm->Nr;i++){ */
  /*   KK = comm->Nrr[i]; */
  /*   MPI_Irecv (recieve[KK],comm->AR[KK],MPI_DOUBLE,KK, */
  /* 	       MPI_ANY_TAG,mpi_comm,&req_r[i]); */
  /* }/\* end i < nproc *\/ */
  
  /* Allocate */
  Ddof = aloc1i (nproc);
  
  /* Set Ddof */
  Ddof[0] = DomDof[0];
  for (i=1;i<nproc;i++)
    Ddof[i] = Ddof[i-1] + DomDof[i];
  
  /***** COMM BOUNDARY ELEMENTS *****/
  for(i=0; i<nbndel; i++){
    err += el_stiffmat(bndel[i],Lk,K,Ap,Ai,ndofn,elem,node,hommat,
		       matgeom,sig,eps,d_r,r,npres,sup,iter,nor_min,
		       dt,crpl,stab,FNR,lm,f_u,myrank,nproc,GDof,comm,
		       Ddof,0,opts->analysis_type,PGFEM_hypre);

    /* If there is an error, complete communication and exit */
    if(err != 0) goto send;
  }

  /***** COHESIVE ELEMENTS *****/

  /* Need to split into boundary and interior parts as with the
     regular elements */
  if (opts->cohesive == 1){/* COHESIVE STIFFNESS */
    ndofc = 3; 
    
    /* WHY IS THIS HERE */
    if (iter == 0) FNR = 0;
    
    if (nor_min < 1.e-10)
      nor_min = 1.e-10;
    
    for (i=0;i<nce;i++){
      coel_stiffmat(i,Lk,K,Ap,Ai,ndofc,elem,node,eps,
		    d_r,r,npres,sup,iter,nor_min,dt,crpl,
		    stab,coel,FNR,lm,f_u,myrank,nproc,DomDof,
		    GDof,comm,Ddof,0,opts->analysis_type,PGFEM_hypre);
    }
  }


  /**** BOUNDING ELEMENTS ******/
  /* In the future, these elements will be listed as with the
     volumetric elements to properly overlay computation and
     communication. For now, this is the best place for them as they
     are a proportinally smaller group than the interior volume
     elements for typical problems and are guarenteed to be on the
     communication boundary for periodic domains. */

  /* temporary for compile testing */
  /* int n_be = 0; */
  /* BOUNDING_ELEMENT *b_elems = NULL; */
  for(i=0; i<n_be; i++){
    err += bnd_el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,b_elems,node,hommat,
			   matgeom,sig,eps,d_r,r,npres,sup,iter,nor_min,
			   dt,crpl,stab,FNR,lm,f_u,myrank,nproc,GDof,
			   comm,Ddof,0,opts->analysis_type,PGFEM_hypre);

    /* If there is an error, complete communication and exit */
    if(err != 0) goto send;
  }

  if (PFEM_DEBUG){
    char ofile[50];
    switch(opts->analysis_type){
    case STABILIZED:
      sprintf(ofile,"stab_el_stiff_send_%d.log",myrank);
      break;
    case MINI:
      sprintf(ofile,"MINI_el_stiff_send_%d.log",myrank);
      break;
    case MINI_3F:
      sprintf(ofile,"MINI_3f_el_stiff_send_%d.log",myrank);
      break;
    default:
      sprintf(ofile,"el_stiff_send_%d.log",myrank);
      break;
    }

    FILE *out;
    out = fopen(ofile,"a");
    for(int send_proc=0; send_proc<nproc; send_proc++){
      if(send_proc==myrank || comm->S[send_proc] ==0) continue;
      for(int n_data=0; n_data<comm->AS[send_proc]; n_data++){
	PGFEM_fprintf(out,"%12.12e ",Lk[send_proc][n_data]);
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n");
    fclose(out);
  }

  /**********************************/
  /***** SEND BOUNDARY AND COEL *****/
  /**********************************/
 send:
  err += send_stiffmat_comm(&sta_s,&req_s,Lk,mpi_comm,comm);

  /* /\* Allocate status fields *\/ */
  /* if (comm->Ns == 0) */
  /*   KK = 1; */
  /* else */
  /*   KK = comm->Ns; */
  /* sta_s = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status)); */
  /* req_s = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request)); */

  /* /\* Send data *\/ */
  /* for (i=0;i<comm->Ns;i++){ */
  /*   KK = comm->Nss[i]; */
  /*   MPI_Isend (Lk[KK],comm->AS[KK],MPI_DOUBLE,KK,myrank,mpi_comm,&req_s[i]); */
  /* }/\* end i < nproc *\/ */

  /* If error, complete communication and exit */
  if(err != 0) goto wait;

   /**********************************/
  /*****    CONTINUE WORKING    *****/
  /**********************************/

  /***** INTERIOR ELEMENTS *****/

  if(nbndel > 0){/* this is 99% of the time */
    for(i=0; i<ne; i++){
      if(idx < nbndel-1){
	if(i == 0 && idx == 0 && bndel[idx] == 0){
	  idx++;
	  skip++;
	  continue;
	} else if(i == bndel[idx]){
	  idx++;
	  skip++;
	  continue;
	} else if (idx == 0 && i < bndel[idx]){
	  err = el_stiffmat(i,Lk,K,Ap,Ai,ndofn,elem,node,hommat,matgeom,
			    sig,eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,
			    stab,FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre);
	} else if (idx > 0 && bndel[idx-1] < i && i < bndel[idx]){
	  err = el_stiffmat(i,Lk,K,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			    eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,
			    FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre);
	} else {
	  PGFEM_printf("[%d]ERROR: problem in determining if element %ld"
		 " is on interior.\n", myrank, i);
	}
      } else {
	if(i != bndel[nbndel-1]){
	  err = el_stiffmat(i,Lk,K,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			    eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,
			    FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre);
	}
      }

      /* If there is an error, complete communication and exit */
      if(err != 0) goto wait;
    }

    /* Check to make sure I got all of them */
    if(skip != nbndel - 1){
      PGFEM_printf("[%d]WARNING: number skipped elem != nbndel, check code\n",myrank);
    }
  } else { /* communication by coheisve elements only, nbndel = 0 */
    for(i=0; i<ne; i++){
      err = el_stiffmat(i,Lk,K,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,FNR,
			lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			opts->analysis_type,PGFEM_hypre);
      /* If there is an error, complete communication and exit */
      if(err != 0) goto wait;
    }
  }

  /**********************************/
  /*****    RECEIVE AND ADD     *****/
  /**********************************/
  
 wait:
  err += assemble_nonlocal_stiffmat(comm,sta_r,req_r,PGFEM_hypre,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,comm);

  /* /\* Wait to complete the comunications *\/ */
  /* MPI_Waitall (comm->Ns,req_s,sta_s); */
  /* MPI_Waitall (comm->Nr,req_r,sta_r); */

  /* /\* If there is an error, return after completed send *\/ */
  /* if(err != 0) goto exit_function; */

  /* /\* Add Received data *\/ */

  /* /\* Currently, the comm datatype has longs, so everything must be */
  /*    copied.  In the future, all int/long will be Index_t and we can */
  /*    just use pointers.  For now, intelligently allocate/copy/assemble */
  /*    for each processor received from. *\/ */

  /* int proc; */
  /* int nrows, *row_idx, *ncols, *col_idx; */
  /* for (i=0;i<comm->Nr;i++){ */
  /*   proc = comm->Nrr[i]; */
  /*   nrows = comm->R[proc]; */
  /*   row_idx = aloc1i(nrows); */
  /*   ncols = aloc1i(nrows); */
  /*   col_idx = aloc1i(comm->AR[proc]); */

  /*   idx = 0; /\* counter for col_idx index *\/ */
  /*   for(j=0; j<comm->R[proc]; j++){ */
  /*     row_idx[j] = comm->RGID[proc][j]; */
  /*     ncols[j] = comm->RAp[proc][j]; */
  /*     for(k=0; k<comm->RAp[proc][j]; k++){ */
  /* 	col_idx[idx] = comm->RGRId[proc][idx]; */
  /* 	++idx; */
  /*     } */
  /*   } */

  /*   if(PFEM_DEBUG){ */
  /*     char ofile[50]; */
  /*     switch(opts->analysis_type){ */
  /*     case STABILIZED: */
  /* 	sprintf(ofile,"stab_assem_rec_%d.log",myrank); */
  /* 	break; */
  /*     case MINI: */
  /* 	sprintf(ofile,"MINI_assem_rec_%d.log",myrank); */
  /* 	break; */
  /*     case MINI_3F: */
  /* 	sprintf(ofile,"MINI_3f_assem_rec_%d.log",myrank); */
  /* 	break; */
  /*     default: */
  /* 	sprintf(ofile,"el_assem_rec_%d.log",myrank); */
  /* 	break; */
  /*     } */

  /*     FILE *out; */
  /*     out = fopen(ofile,"a"); */
  /*     PGFEM_fprintf(out,"********************************************\n"); */
  /*     print_array_i(out,ncols,nrows,1,nrows); */
  /*     print_array_i(out,row_idx,nrows,1,nrows); */
  /*     print_array_i(out,col_idx,comm->AR[proc],1,comm->AR[proc]); */
  /*     print_array_d(out,recieve[proc],comm->AR[proc],1,comm->AR[proc]); */
  /*     fclose(out); */
  /*   } */

  /*   /\* Add this processor's info to the matrix *\/ */
  /*   int hy_err = HYPRE_IJMatrixAddToValues(PGFEM_hypre->hypre_k, */
  /* 					   nrows,ncols, */
  /* 					   row_idx,col_idx, */
  /* 					   recieve[proc]); */
  /*   if(hy_err != 0){ */
  /*     PGFEM_printerr("[%d]WARNING: HYPRE_IJMatrixAddToValues returned error. %s:%s:%d\n", */
  /* 	      myrank,__func__,__FILE__,__LINE__); */
  /*   } */

  /*   free(row_idx); */
  /*   free(ncols); */
  /*   free(col_idx); */
  /* } */
    
 /* exit_function: */
  /* Deallocate recieve */
  for (i=0;i<nproc;i++)
    free (recieve[i]);
  free (recieve);
  
  /*  dealocation  */
  for (i=0;i<nproc;i++)
    free(Lk[i]);
  
  free (Lk);
  free (Ddof);
  free (sta_s);
  free (sta_r);
  free (req_s);
  free (req_r);

  return err;
}


/** Assemble non-local parts as they arrive */
int assemble_nonlocal_stiffmat(const COMMUN pgfem_comm,
			       MPI_Status *sta_r,
			       MPI_Request *req_r,
			       PGFEM_HYPRE_solve_info *PGFEM_hypre,
			       double **recv)
{
  int err = 0;
  int comm_idx = 0;
  int n_received = 0;
  while (n_received < pgfem_comm->Nr){
    /* get the communication index */
    err += MPI_Waitany(pgfem_comm->Nr,req_r,&comm_idx,sta_r);

    /* convert communication index to proc_id */
    const int proc = pgfem_comm->Nrr[comm_idx];

    /* get number of rows */
    const int nrows = pgfem_comm->R[proc];

    /* allocate rows and cols to receive */
    int *row_idx = PGFEM_calloc(nrows,sizeof(int));
    int *ncols = PGFEM_calloc(nrows,sizeof(int));
    int *col_idx = PGFEM_calloc(pgfem_comm->AR[proc],sizeof(int));

    /* get row and column ids */
    int idx = 0;
    for(int j=0; j<pgfem_comm->R[proc]; j++){
      row_idx[j] = pgfem_comm->RGID[proc][j];
      ncols[j] = pgfem_comm->RAp[proc][j];
      for(int k=0; k<ncols[j]; k++){
	col_idx[idx] = pgfem_comm->RGRId[proc][idx];
	++idx;
      }
    }

    /* assemble to local part of global stiffness */
    err += HYPRE_IJMatrixAddToValues(PGFEM_hypre->hypre_k,
				     nrows,ncols,row_idx,col_idx,
				     recv[proc]);

    /* free memory */
    free(row_idx);
    free(ncols);
    free(col_idx);

    /* increment counter */
    n_received++;
  }
  return err;
}
