/* HEADER */
#include "RNPsparse_ApAi.h"
#include "PGFEM_io.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "matice.h"
#include "utils.h"
#include "allocation.h"
#include "incl.h"

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#ifndef PFEM_DEBUG_ALL
#define PFEM_DEBUG_ALL 0
#endif

#ifndef PFEM_PRINT
#define PFEM_PRINT 0
#endif

int* RNPsparse_ApAi (int nproc,
		     int myrank,
		     long ne,
		     long n_be,
		     long nn,
		     long ndofn,
		     long *ndofd,
		     ELEMENT *elem,
		     BOUNDING_ELEMENT *b_elems,
		     NODE *node,
		     int *Ap,
		     long nce,
		     COEL *coel,
		     long *DomDof,
		     int *GDof,
		     COMMUN comm,
		     MPI_Comm Comm_RN,
		     const int cohesive)
{
  char jmeno[200];
  FILE *out;
  long i,j,k,II,JJ,*cnL,*cnG,nne,ndofe,*nod,*ap,**AA,*ap1,**ID,
       *LG,Nrs=0,*GL,LI1,LI2,KK,*cncL,*cncG,ndofc,*nodc;
  long *send,NRr=0,*GNRr,*ApRr,**GIDRr,**RECI,**SEND;
  int *Ddof,*Ai;
  MPI_Status stat,*sta_s,*sta_r;
  MPI_Request *req_s,*req_r;

  /* my variables */
  int num_missing, *idx_missing;

  if (PFEM_DEBUG_ALL || PFEM_PRINT){
  
    sprintf (jmeno,"%s%d.report","RNApAi_",myrank);
    if ((out = fopen(jmeno,"w")) == NULL ){
      PGFEM_printf("Output file is not possible to open on processor [%d]\n",myrank);
      PGFEM_printf("Check the output file and run program again\n");
      return (0);
    }
  }
  
  /* Alocation has to be changed for elements with more than 10 nodes and 30 degreess of freedom */
  cnL = aloc1l (30); cnG = aloc1l (30); nod = aloc1l (10); cncL = aloc1l (30); cncG = aloc1l (30); nodc = aloc1l (10);
  LG = aloc1l (*ndofd); 

  comm->S = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->R = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->AS = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->AR = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->GL = (long*) PGFEM_calloc (DomDof[myrank],sizeof(long));

  /*****************************************/
  /*          Check GL mapping             */
  /*****************************************/

  /* First we determine the local to global mapping.  Then we
     determine the global to local mapping.  If any mappings are
     missing, we take note and modify things next. */

  if (PFEM_DEBUG || PFEM_DEBUG_ALL){
    PGFEM_printf("[%d]::Checking Global to Local mapping for regular elements...\n",myrank);
  }

  for (i=0;i<ne;i++){/* List of ID number for Aij */
  
    /* Number of element nodes */  
    nne = elem[i].toe;  
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    ndofe = get_total_ndof_on_elem(nne,nod,node,b_elems,&elem[i]);
    /* Id numbers */
    get_all_dof_ids_on_elem(0,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnL);
    get_all_dof_ids_on_elem(1,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnG);
    
    for (j=0;j<ndofe;j++){/* row */
      II = cnL[j]-1;
      if (II < 0)  continue;
      LG[II] = cnG[j]-1;
    }
  }/* end i < ne */

  /* COHESIVE ELEMENTS */
  if (cohesive == 1){

    if (PFEM_DEBUG || PFEM_DEBUG_ALL){
      PGFEM_printf("[%d]::Checking Global to Local mapping for cohesive elements...\n",myrank);
    }

    ndofc = 4;
    for (i=0;i<nce;i++){
      nne = coel[i].toe/2;
      ndofe = coel[i].toe*ndofc;
      for (j=0;j<coel[i].toe;j++) nodc[j] = coel[i].nod[j];
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nodc,node,cncL);
      get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nodc,node,cncG);

      for (j=0;j<ndofe;j++){/* row */
	II = cncL[j]-1;
	if (II < 0)  continue;
	LG[II] = cncG[j]-1;
	}
    }/* end i < nce */
  }/* end coh == 1 */
  
  /* Create GL mapping noting any missing maps */
  Ddof = aloc1i (nproc); Ddof[0] = DomDof[0]; for (i=1;i<nproc;i++) Ddof[i] = Ddof[i-1] + DomDof[i];
  if (myrank == 0) *GDof = 0; else *GDof = Ddof[myrank-1];
  GL = aloc1l (DomDof[myrank]);

  num_missing = 0;
  for (i=0;i<DomDof[myrank];i++){
    for (j=0;j<*ndofd;j++){
      if (i == LG[j] - *GDof) {GL[i] = comm->GL[i] = j; break;}
      if (i != LG[j] - *GDof && j == *ndofd - 1 && GL[i] == 0) num_missing++;
    }
  }

  if (PFEM_DEBUG || PFEM_DEBUG_ALL){
    /* Print total missing dofs (for testing only). */
    int tot_missing = 0;
    MPI_Reduce(&num_missing,&tot_missing,1,MPI_INT,MPI_SUM,0,Comm_RN);
    if(myrank == 0) PGFEM_printf("Total number of maps missing on all domains:\t %d\n",tot_missing);
  }

  if(num_missing != 0)
    {
      idx_missing = aloc1i(num_missing);
      k=0;
      for (i=0;i<DomDof[myrank];i++){
	for (j=0;j<*ndofd;j++){
	  if (i == LG[j] - *GDof) break;
	  if (i != LG[j] - *GDof && j == *ndofd - 1 && GL[i] == 0) {idx_missing[k] = i; k++;}
	}
      }

      /* Change length of LG to include missing maps */
      LG = change_length(LG,*ndofd,*ndofd+num_missing);

      /* Append mappings */
      k = 0;
      for(i=*ndofd;i<*ndofd+num_missing;i++)
	{
	  LG[i] = idx_missing[k] + *GDof;
	  GL[idx_missing[k]] = i;
	  comm->GL[idx_missing[k]] = i;
	  k++;
	}

      *ndofd += num_missing;

      /* Sanity checks */
      for(i=0;i<*ndofd;i++){
	if(LG[i] < 0) {PGFEM_printf("ERROR[%d]::LG mapping has a map that is out of bounds!\n",myrank);break;}
      }

      if (PFEM_DEBUG || PFEM_DEBUG_ALL){
	for(i=0;i<DomDof[myrank];i++){
	  if(GL[i] == 0) PGFEM_printf("WARNING[%d]::GL[%ld] = 0 (GID %ld)\n",myrank,i,i+*GDof);
	}
	PGFEM_printf("[%d]::ndofd = %ld\nidx_missing[num_missing-1] = %d\n",myrank,*ndofd,idx_missing[num_missing-1]);
      }

      free(idx_missing);

    }/*if(num_missing > 0)*/

  /*****************************************/
  /*        Create ApAi as normal          */
  /*    but with LG and GL already known   */
  /*****************************************/
  comm->LG = (long*) PGFEM_calloc (*ndofd,sizeof(long));
  ap = aloc1l (*ndofd); ap1 = aloc1l (*ndofd); 

  for (i=0;i<ne;i++){/* Number of contributions from elements to Aij */
    
    /* Number of element nodes */
    nne = elem[i].toe;
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    ndofe = get_total_ndof_on_elem(nne,nod,node,b_elems,&elem[i]);
    /* Id numbers */
    get_all_dof_ids_on_elem(0,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnL);
    
    for (j=0;j<ndofe;j++){/* row */
      II = cnL[j]-1;
      if (II < 0)  continue;
      for (k=0;k<ndofe;k++){/* column */
	JJ = cnL[k]-1;
	if (JJ < 0)  continue;
	ap[II]++;
      }
    }
  }/* end i < ne */
  
  /* COHESIVE ELEMENTS */
  if (cohesive == 1){
    /* ndofc =  4
       I have two groups of elements and pressure dof is not accounted for
       in the element loop, when cohesive element lives in the other domain.
    */
    ndofc = 4; 
    for (i=0;i<nce;i++){
      
      nne = coel[i].toe/2;
      ndofe = coel[i].toe*ndofc;
      for (j=0;j<coel[i].toe;j++) nodc[j] = coel[i].nod[j];
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nodc,node,cncL);
      
      for (j=0;j<ndofe;j++){/* row */
	II = cncL[j]-1;
	if (II < 0)  continue;
	for (k=0;k<ndofe;k++){/* column */
	  JJ = cncL[k]-1;
	  if (JJ < 0)  continue;
	  ap[II]++;
	}
      }
    }/* end i < nce */
  }/* end coh == 1 */

  AA = (long**) PGFEM_calloc (*ndofd,sizeof(long*));
  for (i=0;i<*ndofd;i++)
    {
      if(ap[i]>0)
	{
	  AA[i]= (long*) PGFEM_calloc (ap[i],sizeof(long));
	  for(j=0;j<ap[i];j++) AA[i][j] = 0;
	}
      else
	{
	  AA[i]= (long*) PGFEM_calloc (1,sizeof(long));
	  AA[i][0] = 0;
	}
    }

  if (AA == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  for (i=0;i<ne;i++){/* List of ID number for Aij */
  
    /* Number of element nodes */  
    nne = elem[i].toe;  
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    ndofe = get_total_ndof_on_elem(nne,nod,node,b_elems,&elem[i]);
    /* Id numbers */
    get_all_dof_ids_on_elem(0,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnL);
    get_all_dof_ids_on_elem(1,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnG);
    
    for (j=0;j<ndofe;j++){/* row */
      II = cnL[j]-1;
      if (II < 0)  continue;
      for (k=0;k<ndofe;k++){/* column */
	JJ = cnG[k]-1;
	if (JJ < 0)  continue;
	AA[II][ap1[II]] = JJ;
	ap1[II]++;
      }
    }
  }/* end i < ne */

  /* COHESIVE ELEMENTS */
  if (cohesive == 1){
    for (i=0;i<nce;i++){

      nne = coel[i].toe/2;
      ndofe = coel[i].toe*ndofc;
      for (j=0;j<coel[i].toe;j++) nodc[j] = coel[i].nod[j];
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nodc,node,cncL);
      get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nodc,node,cncG);

      for (j=0;j<ndofe;j++){/* row */
	II = cncL[j]-1;
	if (II < 0)  continue;
	for (k=0;k<ndofe;k++){/* column */
	  JJ = cncG[k]-1;
	  if (JJ < 0)  continue;
	  AA[II][ap1[II]] = JJ;
	  ap1[II]++;
	}
      }
    }/* end i < nce */
  }/* end coh == 1 */
  

  if (PFEM_DEBUG || PFEM_DEBUG_ALL){
    if(myrank == 0){
      for(i=0;i<*ndofd;i++)
	PGFEM_printf("LG[%ld] :: %ld\n",i,LG[i]);
    }
  }

  for (k=0;k<*ndofd;k++){/* Sort list of IDs for Aij */
    ap[k] = 0; /* null ap */
    qsort (AA[k], ap1[k], sizeof(long*), compare_long);
  }
  
  for (k=0;k<*ndofd;k++){ /* Number of non-zeros in rows */
    comm->LG[k] = LG[k];
    for (j=0;j<ap1[k]-1;j++) {if (AA[k][j] < AA[k][j+1]) ap[k]++;}
    ap[k]++;
  }

  ID = (long**) PGFEM_calloc (*ndofd,sizeof(long*));
  for (i=0;i<*ndofd;i++)
    {
      if(ap[i]>0)
	ID[i]= (long*) PGFEM_calloc (ap[i],sizeof(long));
      else
	ID[i]= (long*) PGFEM_calloc (1,sizeof(long));
      ap[i] = 0; /* null */ 
    }
  if (ID == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  for (k=0;k<*ndofd;k++){ /* Global row indexes */
    for (j=0;j<ap1[k]-1;j++)
      if (AA[k][j] < AA[k][j+1]){
	ID[k][ap[k]] = AA[k][j];
	ap[k]++;
      }
    ID[k][ap[k]] = AA[k][ap1[k]-1];
    ap[k]++;
  }

  /* If we have added global mappings, force the only entry in the ID
     matrix to be the diagonal */
  if(num_missing > 0)
    {
      for(i=*ndofd-num_missing;i<*ndofd;i++)
	ID[i][0] = LG[i];
    }
  
  /* Free unnecessary constants */
  for (i=0;i<*ndofd;i++) free (AA[i]); free (AA); dealoc1l (ap1);

  /*
  Ddof = aloc1i (nproc); Ddof[0] = DomDof[0]; for (i=1;i<nproc;i++) Ddof[i] = Ddof[i-1] + DomDof[i];
  if (myrank == 0) *GDof = 0; else *GDof = Ddof[myrank-1];
  GL = aloc1l (DomDof[myrank]);
  */

  if (PFEM_DEBUG || PFEM_DEBUG_ALL){

    for (i=0;i<DomDof[myrank];i++){
      for (j=0;j<*ndofd;j++){
	if (i != LG[j] - *GDof && j == *ndofd - 1 && GL[i] == 0)
	  PGFEM_printf("There may not be a local mapping for global id %ld[%ld] on [%d]!\n  Checking::LG[0] = %ld\n",i+*GDof,i,myrank,LG[0]);
      }
    }
  }

  for(i=0;i<DomDof[myrank];i++) Ap[i] = ap[GL[i]];

  /******************************/
  /* ASSEMBLE ARRAYS TO BE SENT */
  /******************************/
  
  if (myrank == 0) LI1 = 0; else LI1 = Ddof[myrank-1]; LI2 = Ddof[myrank] - 1;

  for (i=0;i<*ndofd;i++){/* Number of rows per domain to be sent */
    if (LI1 <= LG[i] && LG[i] <= LI2) continue;
    else Nrs++;
  }
  
  ap1 = aloc1l (Nrs);
  
  Nrs = 0;
  for (i=0;i<*ndofd;i++){/* Local row indexes on domain to be sent */
    if (LI1 <= LG[i] && LG[i] <= LI2) continue;
    else {ap1[Nrs] = i; Nrs++;}
  }
  
  /* Determine where to send */
  for (j=0;j<nproc;j++){
    k = 1; if (j == myrank) continue;
    if (j == 0) LI1 = 0; else LI1 = Ddof[j-1]; LI2 = Ddof[j] - 1;
    for (i=0;i<Nrs;i++){
      if (LI1 <= LG[ap1[i]] && LG[ap1[i]] <= LI2) {comm->S[j] = k; k++;}
    }
  }

  AA = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<nproc;i++) {
    if (myrank == i || comm->S[i] == 0) k = 1; else k = comm->S[i];
    AA[i] = (long*) PGFEM_calloc (k,sizeof(long));
  }
  if (AA == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  comm->SLID = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<nproc;i++) {
    if (myrank == i || comm->S[i] == 0) continue;
    else comm->SLID[i] = (long*) PGFEM_calloc (comm->S[i],sizeof(long));
  }
  if (comm->SLID == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  /* Determine what GID/LID to send */
  for (j=0;j<nproc;j++){
    k = 0;
    if (j == myrank) continue;
    if (j == 0) LI1 = 0; else LI1 = Ddof[j-1]; LI2 = Ddof[j] - 1;
    for (i=0;i<Nrs;i++){
      if (LI1 <= LG[ap1[i]] && LG[ap1[i]] <= LI2) {AA[j][k] = comm->SLID[j][k] = ap1[i]; k++;} 
    }
  }

  /* SEND TO WHOM I WILL TALK TO */
  for (i=0;i<nproc;i++){
    if (i == myrank) continue;
    MPI_Send (&comm->S[i],1,MPI_INT,i,myrank,Comm_RN);
  }
  
  /* RECIEVE FROM SENDER */
  for (i=0;i<nproc;i++){
    if (i == myrank) continue;
    MPI_Recv (&j,1,MPI_INT,i,MPI_ANY_TAG,Comm_RN,&stat);
    comm->R[i]=j;
  }
  
  /* Compute how much to send and recieve */
  comm->Ns = 0;
  for (i=0;i<nproc;i++){
    if (i == myrank || comm->S[i] == 0) continue;
    comm->Ns++;
  }
  comm->Nr = 0;
  for (i=0;i<nproc;i++){
    if (i == myrank || comm->R[i] == 0) continue;
    comm->Nr++;
  }

  /* Allocate arrays */
  if (comm->Ns == 0) KK = 1; else KK = comm->Ns; comm->Nss = (long*) PGFEM_calloc (KK,sizeof(long*));
  if (comm->Nr == 0) KK = 1; else KK = comm->Nr; comm->Nrr = (long*) PGFEM_calloc (KK,sizeof(long*));
  
  comm->Ns = 0;
  for (i=0;i<nproc;i++){
    if (i == myrank || comm->S[i] == 0) continue;
    comm->Nss[comm->Ns] = i;
    comm->Ns++;
  }
  
  comm->Nr = 0;
  for (i=0;i<nproc;i++){
    if (i == myrank || comm->R[i] == 0) continue;
    comm->Nrr[comm->Nr] = i;
    comm->Nr++;
  }
  
  /*************************************************************************************/
  /* SEND ALL INFORMATION || Global row number : Number of nonzeros in row : Global ID */
  /*************************************************************************************/  
  
  /* How many numbers I will send */
  for (i=0;i<nproc;i++){
    if (i == myrank) continue;
    comm->AS[i] = 2*comm->S[i];
  }
  
  /* Send allocation information */
  for (i=0;i<comm->Ns;i++){
    KK = comm->Nss[i];
    MPI_Send (&comm->AS[KK],1,MPI_LONG,KK,myrank,Comm_RN);
  }
  
  /* Recieve allocation information */
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    MPI_Recv (&j,1,MPI_LONG,KK,MPI_ANY_TAG,Comm_RN,&stat);
    comm->AR[KK]=j;
  }
  
  /* Allocate status and request fields */
  if (comm->Ns == 0) KK = 1; else KK = comm->Ns; sta_s = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
  if (comm->Nr == 0) KK = 1; else KK = comm->Nr; sta_r = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
  if (comm->Ns == 0) KK = 1; else KK = comm->Ns; req_s = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));
  if (comm->Nr == 0) KK = 1; else KK = comm->Nr; req_r = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));
  
  /* Allocate recieve */
  SEND = (long**) PGFEM_calloc (nproc,sizeof(long*));
  RECI = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<nproc;i++) {
    if (comm->AS[i] == 0) JJ = 1; else JJ = comm->AS[i]; SEND[i] = (long*) PGFEM_calloc (JJ,sizeof(long));
    if (comm->AR[i] == 0) JJ = 1; else JJ = comm->AR[i]; RECI[i] = (long*) PGFEM_calloc (JJ,sizeof(long));
  }
  if (SEND == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}  
  if (RECI == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}  
  
  /*************/
  /* Send data */
  /*************/
  for (i=0;i<comm->Ns;i++){
    KK = comm->Nss[i];
    
    for (j=0;j<comm->S[KK];j++) {SEND[KK][j] = LG[AA[KK][j]]; SEND[KK][comm->S[KK]+j] = ap[AA[KK][j]];}
    MPI_Isend (SEND[KK],comm->AS[KK],MPI_LONG,KK,myrank,Comm_RN,&req_s[i]);
  }
  
  /****************/
  /* Recieve data */
  /****************/
  for (i=0;i<nproc;i++) NRr += comm->R[i];
  GNRr = aloc1l (NRr); ApRr = aloc1l (NRr);
  
  comm->RGID = (long**) PGFEM_calloc (nproc,sizeof(long*));
  comm->RAp = (long**) PGFEM_calloc (nproc,sizeof(long*));

  for (i=0;i<comm->Nr;i++) {
    KK = comm->Nrr[i];
    comm->RGID[KK] = (long*) PGFEM_calloc (comm->R[KK],sizeof(long));
    comm->RAp[KK] = (long*) PGFEM_calloc (comm->R[KK],sizeof(long));
  }
  if (comm->RGID == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}

  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    MPI_Irecv (RECI[KK],comm->AR[KK],MPI_LONG,KK,MPI_ANY_TAG,Comm_RN,&req_r[i]);
  }
  
  /* Wait to complete the communications */
  MPI_Waitall (comm->Ns,req_s,sta_s);
  MPI_Waitall (comm->Nr,req_r,sta_r);
  
  /* Unpack it */
  II = 0;
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    for (j=0;j<comm->R[KK];j++) {GNRr[II+j] = comm->RGID[KK][j] = RECI[KK][j]; ApRr[II+j] = comm->RAp[KK][j] = RECI[KK][comm->R[KK]+j];}
    II += comm->R[KK];
  }
  
  /* Deallocate recieve */
  for (i=0;i<nproc;i++) {free (SEND[i]); free (RECI[i]);} free (SEND); free (RECI); 
  
  /* send */
  comm->SAp = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<comm->Ns;i++) {
    KK = comm->Nss[i];
    comm->SAp[KK] = (long*) PGFEM_calloc (comm->S[KK],sizeof(long));
  }
  if (comm->SAp == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  /* How many numbers I will send */
  for (i=0;i<nproc;i++){comm->AS[i] = 0; comm->AR[i] = 0;}
  for (i=0;i<nproc;i++){
    if (i == myrank) continue;
    k = 0; for (j=0;j<comm->S[i];j++) k += comm->SAp[i][j] = ap[AA[i][j]];
    comm->AS[i] = k;
  }
  
  /* Send allocation information */
  for (i=0;i<comm->Ns;i++){
    KK = comm->Nss[i];
    MPI_Send (&comm->AS[KK],1,MPI_LONG,KK,myrank,Comm_RN);
  }
  
  /* Recieve allocation information */
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    MPI_Recv (&j,1,MPI_LONG,KK,MPI_ANY_TAG,Comm_RN,&stat);
    comm->AR[KK]=j;
  }
  /*************/
  /* Send data */
  /*************/
  comm->SGRId = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<comm->Ns;i++){
    KK = comm-> Nss[i];
    comm->SGRId[KK] = (long*) PGFEM_calloc (comm->AS[KK],sizeof(long));
  }
  if (comm->SGRId == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  /* Allocate recieve */
  RECI = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<nproc;i++) {if (comm->AR[i] == 0) JJ = 1; else JJ = comm->AR[i]; RECI[i] = (long*) PGFEM_calloc (JJ,sizeof(long));}
  if (RECI == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}    
  
  for (i=0;i<comm->Ns;i++){
    KK = comm->Nss[i];
    
    II = 0;
    for (j=0;j<comm->S[KK];j++){
      for (k=0;k<ap[AA[KK][j]];k++){
	comm->SGRId[KK][II] = ID[AA[KK][j]][k];
	II++;
      }
    }
    
    MPI_Isend (comm->SGRId[KK],comm->AS[KK],MPI_LONG,KK,myrank,Comm_RN,&req_s[i]);
    
  }/* end i < comm->Ns */
  
  /****************/
  /* Recieve data */
  /****************/
  if (NRr == 0) KK = 1; else KK = NRr; GIDRr = (long**) PGFEM_calloc (KK,sizeof(long*));
  for (i=0;i<NRr;i++) GIDRr[i]= (long*) PGFEM_calloc (ApRr[i],sizeof(long));
  if (GIDRr == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  comm->RGRId = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<comm->Nr;i++) {
    KK = comm->Nrr[i];
    comm->RGRId[KK] = (long*) PGFEM_calloc (comm->AR[KK],sizeof(long));
  }
  if (comm->RGRId == NULL){PGFEM_printf ("\n Memory is full.\n");  fflush(stdout); PGFEM_Comm_code_abort (Comm_RN,0);}
  
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    MPI_Irecv (RECI[KK],comm->AR[KK],MPI_LONG,KK,MPI_ANY_TAG,Comm_RN,&req_r[i]);
  }
  
  /* Wait to complete the communications */
  MPI_Waitall (comm->Ns,req_s,sta_s);
  MPI_Waitall (comm->Nr,req_r,sta_r);
  
  /* Unpack ir */
  JJ = 0;
  for (i=0;i<comm->Nr;i++){
    KK = comm->Nrr[i];
    II = 0;
    for (j=0;j<comm->R[KK];j++){
      for (k=0;k<ApRr[JJ];k++){
	GIDRr[JJ][k] = comm->RGRId[KK][II] = RECI[KK][II];
	II++;
      }
      JJ++;
    }
  }/* end i < nproc */


  /****************************/
  /* END SEND ALL INFORMATION */
  /****************************/
  
  if (PFEM_DEBUG_ALL || PFEM_PRINT){
  
    PGFEM_fprintf (out,"Process [%d] || Number of local unknowns = %ld :: ndofn = %ld\n\n",myrank,*ndofd,ndofn);
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Number of GLOBAL non-zeros per row\n");
    for (i=0;i<DomDof[myrank];i++){
      PGFEM_fprintf (out,"[%ld-%ld-%ld]: Ap=%d || ",i,GL[i],LG[GL[i]],Ap[i]);
      for (j=0;j<Ap[i];j++) PGFEM_fprintf (out,"%ld ",ID[GL[i]][j]);
      PGFEM_fprintf (out,"\n");
    }
  }

  /* Add recieved rows to Ap and ID */
  for (i=0;i<NRr;i++){
    II = ApRr[i] + Ap[GNRr[i]-*GDof];
      
    send = aloc1l (II);
    JJ = 0; for (j=0;j<Ap[GNRr[i]-*GDof];j++) {send[j] = ID[GL[GNRr[i]-*GDof]][j]; JJ++;}
    for (j=0;j<ApRr[i];j++) send[JJ+j] = GIDRr[i][j];
    
    /* Sort Row IDs */

    qsort (send,II, sizeof(long), compare_long);
    
    /* GLOBAL number of non-zeros in row */
    JJ = 1; for (j=0;j<II-1;j++) if (send[j] < send[j+1]) JJ++; 
    
    Ap[GNRr[i]-*GDof] = JJ;

    free (ID[GL[GNRr[i]-*GDof]]); ID[GL[GNRr[i]-*GDof]]= (long*) PGFEM_calloc (JJ,sizeof(long));
    
    JJ = 0; for (j=0;j<II-1;j++) if (send[j] < send[j+1]) {ID[GL[GNRr[i]-*GDof]][JJ] = send[j]; JJ++;}
    
    ID[GL[GNRr[i]-*GDof]][JJ] = send[II-1];
    
    dealoc1l (send);
  }/* end i <NRr */
  
  if (PFEM_PRINT || PFEM_DEBUG_ALL){
    PGFEM_fprintf (out,"NODES\n");
    for (i=0;i<nn;i++){
      PGFEM_fprintf (out,"[%ld] || LID : %ld %ld %ld %ld || GID : %ld %ld %ld %ld\n",i,node[i].id[0],node[i].id[1],node[i].id[2],node[i].id[3],node[i].Gid[0],node[i].Gid[1],node[i].Gid[2],node[i].Gid[3]);
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Number of non-zeros per row\n");
    for (i=0;i<*ndofd;i++){
      PGFEM_fprintf (out,"[%ld-%ld]: ap=%ld | ID : ",i,LG[i],ap[i]);
      for (j=0;j<ap[i];j++) PGFEM_fprintf (out,"%ld ",ID[i][j]);
      PGFEM_fprintf (out,"\n");
    }
  
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Number of Global Dofs on domain : %ld\n",DomDof[myrank]);
    PGFEM_fprintf (out,"Number of GLOBAL non-zeros per row\n");
    for (i=0;i<DomDof[myrank];i++){
      PGFEM_fprintf (out,"[%ld-%ld-%ld]: Ap=%d || ",i,GL[i],LG[GL[i]],Ap[i]);
      for (j=0;j<Ap[i];j++) PGFEM_fprintf (out,"%ld ",ID[GL[i]][j]);
      PGFEM_fprintf (out,"\n");
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Number of rows to be send from domain [%d] :: %ld\n",myrank,Nrs);
    for (i=0;i<Nrs;i++){
      PGFEM_fprintf (out,"[%ld-%ld]: ap=%ld || ",ap1[i],LG[ap1[i]],ap[ap1[i]]);
      for (j=0;j<ap[ap1[i]];j++) PGFEM_fprintf (out,"%ld ",ID[ap1[i]][j]);
      PGFEM_fprintf (out,"\n");
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"[%d] - send to :\n",myrank);
    for (i=0;i<nproc;i++){
      PGFEM_fprintf (out,"%ld || %ld ::  ",i,comm->S[i]);
      for (j=0;j<comm->S[i];j++) PGFEM_fprintf (out,"%ld-%ld  ",AA[i][j],LG[AA[i][j]]);
      PGFEM_fprintf (out,"\n");
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"[%d] - send to :\n",myrank);
    for (i=0;i<nproc;i++) PGFEM_fprintf (out,"%ld[%ld]  ",comm->S[i],comm->AS[i]);
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"[%d] - recieve from :\n",myrank);
    for (i=0;i<nproc;i++) PGFEM_fprintf (out,"%ld[%ld]  ",comm->R[i],comm->AR[i]);
    PGFEM_fprintf (out,"\n");
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Recieved rows - %ld\n",NRr);
    for (i=0;i<NRr;i++){
      PGFEM_fprintf (out,"%ld : Ap=%ld  ||  ",GNRr[i],ApRr[i]);
      for (j=0;j<ApRr[i];j++) PGFEM_fprintf (out,"%ld ",GIDRr[i][j]);
      PGFEM_fprintf (out,"\n");
    }
    PGFEM_fprintf (out,"\n");
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"FINAL Number of GLOBAL non-zeros per row\n");
    for (i=0;i<DomDof[myrank];i++){
      PGFEM_fprintf (out,"[%ld-%ld-%ld]: Ap=%d || ",i,GL[i],LG[GL[i]],Ap[i]);
      for (j=0;j<Ap[i];j++) PGFEM_fprintf (out,"%ld ",ID[GL[i]][j]);
      PGFEM_fprintf (out,"\n");
    }
    
    fclose (out);
  }
 
  /* Prepare Ap and Ai vectors on domains */
  II = 0; for (i=0;i<DomDof[myrank];i++) II += Ap[i];
  Ai = aloc1i (II);
  
  II = Ap[0]; Ap[0] = 0; k = 0;
  for (i=1;i<DomDof[myrank]+1;i++){
    for (j=0;j<II;j++) {Ai[k] = ID[GL[i-1]][j]; k++;}
    JJ = Ap[i];
    Ap[i] = II + Ap[i-1];
    II = JJ;
  }

  
  /* Deallocate recieve */
  for (i=0;i<nproc;i++) free (RECI[i]); free (RECI);

  for (i=0;i<*ndofd;i++) free (ID[i]); free (ID); for (i=0;i<nproc;i++) free (AA[i]); free (AA); for (i=0;i<NRr;i++) free (GIDRr[i]); free (GIDRr);
  dealoc1l (cnL); dealoc1l (cnG); dealoc1l (cncL); dealoc1l (cncG); dealoc1l (nodc); dealoc1l (nod); dealoc1l (ap); dealoc1l (LG); dealoc1i (Ddof); 
  dealoc1l (ap1); dealoc1l (GL); dealoc1l (GNRr); dealoc1l (ApRr); free (sta_s); free (sta_r); free (req_s); free (req_r);

  if (PFEM_DEBUG || PFEM_DEBUG_ALL){
    if(myrank == 0) PGFEM_printf("\nExiting %s: %s\n",__FILE__,__func__);
  }

  return (Ai);
}
