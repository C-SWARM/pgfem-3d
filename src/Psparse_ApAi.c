/* HEADER */
#include "Psparse_ApAi.h"
#include "string.h"
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

/**
 * Determine what processors talk to each other and allocate
 * appropriate receive space in comm.
 *
 * \param[in/out] comm Allocates Nss and Nrr, sets Ns and Nr
 *
 * \return non-zero on error.
 *
 * Side effects: non-blocking communication between all processes.
 */
static int determine_comm_pattern(COMMUN comm,
				  const MPI_Comm mpi_comm);

/**
 * Communicate the number of rows/columns.
 *
 * \param[in/out] comm Allocates internal space and modifies values
 * \param[out] NRr Returns number of rows to on the domain
 * \param[out] GNRr Allocates space and returns ids of rows that are received
 * \param[out] ApRr Allocates space and returns column ids of communicated entires.
 *
 * \return non-zero on error
 * Side effects: Non-blocking point-to-point communication based on comm
 */
static int communicate_number_row_col(COMMUN comm,
				      long *NRr,
				      long **GNRr,
				      long **ApRr,
				      const long *LG,
				      const long *ap,
				      long **AA,
				      const MPI_Comm mpi_comm);

/**
 * Communicate all the row/col numbers.
 *
 * \param[in,out] comm Allocates internal space and modifies values
 * \param[out] GIDRr Allocates space and initializes values
 * \param[in] AA,ID are not modified
 *
 * \return non-zero on error
 * Side-effects: Non-blocking pt2pt communication based on comm
 */
static int communicate_row_info(COMMUN comm,
				long ***GIDRr,
				const long NRr,
				const long *ApRr,
				const long *ap,
				long **AA,
				long **ID,				
				const MPI_Comm mpi_comm);

int* Psparse_ApAi (int nproc,
		   int myrank,
		   long ne,
		   long n_be,
		   long nn,
		   long ndofn,
		   long ndofd,
		   ELEMENT *elem,
		   BOUNDING_ELEMENT *b_elems,
		   NODE *node,
		   int *Ap,
		   long nce,
		   COEL *coel,
		   long *DomDof,
		   int *GDof,
		   COMMUN comm,
		   MPI_Comm Comm_Orig,
		   const int cohesive)
{
  char jmeno[200];
  FILE *out=NULL;
  long i,j,k,II,JJ,*cnL=NULL,*cnG=NULL,nne,ndofe,*nod=NULL,*ap=NULL;
  long **AA=NULL,*ap1=NULL,**ID=NULL,*LG=NULL,Nrs=0,*GL=NULL;
  long LI1,LI2,*cncL=NULL,*cncG=NULL,ndofc,*nodc=NULL;
  long *send=NULL,NRr=0,*GNRr=NULL,*ApRr=NULL,**GIDRr=NULL;
  int *Ddof = NULL;
  int *Ai = NULL;
  MPI_Status *sta_s=NULL, *sta_r=NULL;
  MPI_Request *req_s=NULL, *req_r=NULL;

  if (PFEM_DEBUG_ALL || PFEM_PRINT){
    sprintf (jmeno,"%s%d.report","ApAi_",myrank);
    if ((out = fopen(jmeno,"w")) == NULL ){
      PGFEM_printf("Output file is not possible to"
		   " open on processor [%d]\n",myrank);
      PGFEM_printf("Check the output file and run program again\n");
      return (0);
    } 
  }
  
  /* Alocation has to be changed for elements with more than 10 nodes
     and 30 degreess of freedom */
  cnL = aloc1l (30);
  cnG = aloc1l (30);
  nod = aloc1l (10);
  cncL = aloc1l (30);
  cncG = aloc1l (30);
  nodc = aloc1l (10);
  ap = aloc1l (ndofd);
  ap1 = aloc1l (ndofd);
  LG = aloc1l (ndofd); 

  comm->S = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->R = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->AS = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->AR = (long*) PGFEM_calloc (nproc,sizeof(long));
  comm->LG = (long*) PGFEM_calloc (ndofd,sizeof(long));
  comm->GL = (long*) PGFEM_calloc (DomDof[myrank],sizeof(long));

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
      if (II < 0)
	continue;
      for (k=0;k<ndofe;k++){/* column */
	JJ = cnL[k]-1;
	if (JJ < 0)  
	  continue;
	ap[II]++;
      }
    }
  }/* end i < ne */
  
  /* COHESIVE ELEMENTS */
  if (cohesive == 1){
    /* ndofc = ndofn
       I have two groups of elements and pressure dof is
       not accounted for (in STABILIZED) in the element loop, when
       cohesive element lives in the other domain.
    */
    ndofc = ndofn; 
    for (i=0;i<nce;i++){
      
      nne = coel[i].toe/2;
      ndofe = coel[i].toe*ndofc;
      for (j=0;j<coel[i].toe;j++)
	nodc[j] = coel[i].nod[j];
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nodc,node,cncL);
      
      for (j=0;j<ndofe;j++){/* row */
	II = cncL[j]-1;
	if (II < 0)
	  continue;
	for (k=0;k<ndofe;k++){/* column */
	  JJ = cncL[k]-1;
	  if (JJ < 0)
	    continue;
	  ap[II]++;
	}
      }
    }/* end i < nce */
  }/* end coh == 1 */

  AA = (long**) PGFEM_calloc (ndofd,sizeof(long*));
  for (i=0;i<ndofd;i++)
    AA[i]= (long*) PGFEM_calloc (ap[i],sizeof(long));
  null_quit((void*) AA,0);
  
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
      for (j=0;j<coel[i].toe;j++)
	nodc[j] = coel[i].nod[j];
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nodc,node,cncL);
      get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nodc,node,cncG);

      for (j=0;j<ndofe;j++){/* row */
	II = cncL[j]-1;
	if (II < 0)  continue;
	LG[II] = cncG[j]-1;
	for (k=0;k<ndofe;k++){/* column */
	  JJ = cncG[k]-1;
	  if (JJ < 0)  continue;
	  AA[II][ap1[II]] = JJ;
	  ap1[II]++;
	}
      }
    }/* end i < nce */
  }/* end coh == 1 */
  
  for (k=0;k<ndofd;k++){/* Sort list of IDs for Aij */
    ap[k] = 0; /* null ap */
    qsort (AA[k], ap1[k], sizeof(long*), compare_long);
  }
  
  ap[0] = 0;
  for (k=0;k<ndofd;k++){ /* Number of non-zeros in rows */
    comm->LG[k] = LG[k];
    for (j=0;j<ap1[k]-1;j++)
      if (AA[k][j] < AA[k][j+1])
	ap[k]++;
    ap[k]++;
  }
  
  ID = (long**) PGFEM_calloc (ndofd,sizeof(long*));
  for (i=0;i<ndofd;i++) {
    ID[i]= (long*) PGFEM_calloc (ap[i],sizeof(long));
    ap[i] = 0; /* null */
  }

  null_quit((void*) ID,0);

  for (k=0;k<ndofd;k++){ /* Global row indexes */
    for (j=0;j<ap1[k]-1;j++) {
      if (AA[k][j] < AA[k][j+1]){
	ID[k][ap[k]] = AA[k][j];
	ap[k]++;
      }
    }
    ID[k][ap[k]] = AA[k][ap1[k]-1];
    ap[k]++;
  }
  
  /* Free unnecessary constants */
  for (i=0;i<ndofd;i++)
    free (AA[i]);
  free (AA);
  dealoc1l (ap1); ap1 = NULL;
  
  Ddof = aloc1i (nproc);
  Ddof[0] = DomDof[0];

  for (i=1;i<nproc;i++)
    Ddof[i] = Ddof[i-1] + DomDof[i];

  if (myrank == 0)
    *GDof = 0;
  else
    *GDof = Ddof[myrank-1];
  GL = aloc1l (DomDof[myrank]);

  /* Sort Global Dof */
  for (i=0;i<DomDof[myrank];i++){
    for (j=0;j<ndofd;j++){
      if (i == LG[j] - *GDof) {
	Ap[i] = ap[j];
	GL[i] = comm->GL[i] = j;
	break;
      }
      if (i != LG[j] - *GDof && j == ndofd - 1 && GL[i] == 0){
	PGFEM_printf("There is no local mapping for "
	       "global id %ld on [%d]!\n",i+*GDof,myrank);
      }
    }
  }

  /******************************/
  /* ASSEMBLE ARRAYS TO BE SENT */
  /******************************/
  
  if (myrank == 0)
    LI1 = 0;
  else
    LI1 = Ddof[myrank-1];
  LI2 = Ddof[myrank] - 1;

  for (i=0;i<ndofd;i++){/* Number of rows per domain to be sent */
    if (LI1 <= LG[i] && LG[i] <= LI2)
      continue;
    else Nrs++;
  }
  
  if(Nrs > 0){
    ap1 = aloc1l (Nrs);

    Nrs = 0;
    for (i=0;i<ndofd;i++){/* Local row indexes on domain to be sent */
      if (LI1 <= LG[i] && LG[i] <= LI2)
	continue;
      else {ap1[Nrs] = i; Nrs++;}
    }
  }
  
  /* Determine where to send */
  for (j=0;j<nproc;j++){
    k = 1;
    if (j == myrank)
      continue;
    if (j == 0)
      LI1 = 0;
    else
      LI1 = Ddof[j-1];
    LI2 = Ddof[j] - 1;
    for (i=0;i<Nrs;i++){
      if (LI1 <= LG[ap1[i]] && LG[ap1[i]] <= LI2) {
	comm->S[j] = k;
	k++;
      }
    }
  }

  AA = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<nproc;i++) {
    if (myrank == i || comm->S[i] == 0) 
      k = 1;
    else
      k = comm->S[i];
    AA[i] = (long*) PGFEM_calloc (k,sizeof(long));
  }

  null_quit((void*) AA,0);  

  comm->SLID = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (i=0;i<nproc;i++) {
    if (myrank == i || comm->S[i] == 0)
      continue;
    else comm->SLID[i] = (long*) PGFEM_calloc (comm->S[i],sizeof(long));
  }

  null_quit((void*) comm->SLID,0);  

  /* Determine what GID/LID to send */
  for (j=0;j<nproc;j++){
    k = 0;
    if (j == myrank)
      continue;
    if (j == 0)
      LI1 = 0;
    else
      LI1 = Ddof[j-1];
    LI2 = Ddof[j] - 1;
    for (i=0;i<Nrs;i++){
      if (LI1 <= LG[ap1[i]] && LG[ap1[i]] <= LI2) {
	AA[j][k] = comm->SLID[j][k] = ap1[i];
	k++;
      } 
    }
  }

  /* Communicate who I am communicating with
   *============================================= */
  determine_comm_pattern(comm,Comm_Orig);
  
  /* Communicate how many rows/columns I am sending
   *============================================= */
  communicate_number_row_col(comm,&NRr,&GNRr,&ApRr,
			     LG,ap,AA,Comm_Orig);

  /* Communicate the row/column information
   *============================================= */
  communicate_row_info(comm,&GIDRr,NRr,ApRr,
		       ap,AA,ID,Comm_Orig);

  /****************************/
  /* END SEND ALL INFORMATION */
  /****************************/
  
  /* PRINT */

  if (PFEM_DEBUG_ALL || PFEM_PRINT) {
    PGFEM_fprintf (out,"Process [%d] || Number of local unknowns"
	     " = %ld :: ndofn = %ld\n\n",myrank,ndofd,ndofn);
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
    JJ = 0;
    for (j=0;j<Ap[GNRr[i]-*GDof];j++) {
      send[j] = ID[GL[GNRr[i]-*GDof]][j];
      JJ++;
    }

    for (j=0;j<ApRr[i];j++) 
      send[JJ+j] = GIDRr[i][j];
    
    /* Sort Row IDs */

    qsort (send,II, sizeof(long), compare_long);
    
    /* GLOBAL number of non-zeros in row */
    JJ = 1; for (j=0;j<II-1;j++) if (send[j] < send[j+1]) JJ++; 
    
    Ap[GNRr[i]-*GDof] = JJ;

    free (ID[GL[GNRr[i]-*GDof]]);
    ID[GL[GNRr[i]-*GDof]]= (long*) PGFEM_calloc (JJ,sizeof(long));
    
    JJ = 0; 
    for (j=0;j<II-1;j++) {
      if (send[j] < send[j+1]) {
	ID[GL[GNRr[i]-*GDof]][JJ] = send[j];
	JJ++;
      }
    }
    
    ID[GL[GNRr[i]-*GDof]][JJ] = send[II-1];
    
    dealoc1l (send);
  }/* end i <NRr */
  
  /* PRINT */  
  if (PFEM_DEBUG_ALL || PFEM_PRINT) {
    PGFEM_fprintf (out,"NODES\n");
    for (i=0;i<nn;i++){
      PGFEM_fprintf (out,"[%ld] || LID : %ld %ld %ld %ld ||"
		     " GID : %ld %ld %ld %ld\n", i,node[i].id[0],
		     node[i].id[1],node[i].id[2],node[i].id[3],
		     node[i].Gid[0],
		     node[i].Gid[1],node[i].Gid[2],node[i].Gid[3]);
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Number of non-zeros per row\n");
    for (i=0;i<ndofd;i++){
      PGFEM_fprintf (out,"[%ld-%ld]: ap=%ld | ID : ",i,LG[i],ap[i]);
      for (j=0;j<ap[i];j++) PGFEM_fprintf (out,"%ld ",ID[i][j]);
      PGFEM_fprintf (out,"\n");
    }
  
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Number of Global Dofs on domain : %ld\n",
		   DomDof[myrank]);
    PGFEM_fprintf (out,"Number of GLOBAL non-zeros per row\n");
    for (i=0;i<DomDof[myrank];i++){
      PGFEM_fprintf (out,"[%ld-%ld-%ld]: Ap=%d || "
		     ,i,GL[i],LG[GL[i]],Ap[i]);
      for (j=0;j<Ap[i];j++) PGFEM_fprintf (out,"%ld ",ID[GL[i]][j]);
      PGFEM_fprintf (out,"\n");
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Number of rows to be send "
	     "from domain [%d] :: %ld\n",myrank,Nrs);
    for (i=0;i<Nrs;i++){
      PGFEM_fprintf (out,"[%ld-%ld]: ap=%ld || ",
		     ap1[i],LG[ap1[i]],ap[ap1[i]]);
      for (j=0;j<ap[ap1[i]];j++) PGFEM_fprintf (out,"%ld ",
						ID[ap1[i]][j]);
      PGFEM_fprintf (out,"\n");
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"[%d] - send to :\n",myrank);
    for (i=0;i<nproc;i++){
      PGFEM_fprintf (out,"%ld || %ld ::  ",i,comm->S[i]);
      for (j=0;j<comm->S[i];j++)
	PGFEM_fprintf (out,"%ld-%ld  ",AA[i][j],LG[AA[i][j]]);
      PGFEM_fprintf (out,"\n");
    }
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"[%d] - send to :\n",myrank);
    for (i=0;i<nproc;i++)
      PGFEM_fprintf (out,"%ld[%ld]  ",comm->S[i],comm->AS[i]);
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"[%d] - recieve from :\n",myrank);
    for (i=0;i<nproc;i++)
      PGFEM_fprintf (out,"%ld[%ld]  ",comm->R[i],comm->AR[i]);
    PGFEM_fprintf (out,"\n");
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"Recieved rows - %ld\n",NRr);
    for (i=0;i<NRr;i++){
      PGFEM_fprintf (out,"%ld : Ap=%ld  ||  ",GNRr[i],ApRr[i]);
      for (j=0;j<ApRr[i];j++)
	PGFEM_fprintf (out,"%ld ",GIDRr[i][j]);
      PGFEM_fprintf (out,"\n");
    }
    PGFEM_fprintf (out,"\n");
    
    PGFEM_fprintf (out,"\n");
    PGFEM_fprintf (out,"FINAL Number of GLOBAL non-zeros per row\n");
    for (i=0;i<DomDof[myrank];i++){
      PGFEM_fprintf (out,"[%ld-%ld-%ld]: Ap=%d || ",
		     i,GL[i],LG[GL[i]],Ap[i]);
      for (j=0;j<Ap[i];j++)
	PGFEM_fprintf (out,"%ld ",ID[GL[i]][j]);
      PGFEM_fprintf (out,"\n");
    }
    
    fclose (out);
  }/* END DEBUGGING */


  /* Prepare Ap and Ai vectors on domains */
  II = 0;
  for (i=0;i<DomDof[myrank];i++)
    II += Ap[i];
  Ai = aloc1i (II);
  
  II = Ap[0]; Ap[0] = 0; k = 0;
  for (i=1;i<DomDof[myrank]+1;i++){
    for (j=0;j<II;j++) {
      Ai[k] = ID[GL[i-1]][j];
      k++;
    }
    JJ = Ap[i];
    Ap[i] = II + Ap[i-1];
    II = JJ;
  }

  
  /* Deallocate recieve */
  dealoc2l(AA,nproc);
  dealoc2l(ID,ndofd);
  dealoc2l(GIDRr,NRr);

  dealoc1l (cnL);
  dealoc1l (cnG);
  dealoc1l (cncL);
  dealoc1l (cncG);
  dealoc1l (nodc);
  dealoc1l (nod);
  dealoc1l (ap);
  dealoc1l (LG);
  dealoc1i (Ddof); 
  dealoc1l (ap1);
  dealoc1l (GL);
  dealoc1l (GNRr);
  dealoc1l (ApRr);
  free (sta_s);
  free (sta_r);
  free (req_s);
  free (req_r);

  return (Ai);
}

static int determine_comm_pattern(COMMUN comm,
				  const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(mpi_comm,&myrank);
  MPI_Comm_size(mpi_comm,&nproc);

  MPI_Status t_sta_r;
  MPI_Status *t_sta_s = NULL;
  MPI_Request *t_req_s = NULL;
  MPI_Request *t_req_r = NULL;
  if(nproc > 1){
    t_sta_s = PGFEM_calloc(nproc-1,sizeof(MPI_Status));
    t_req_s = PGFEM_calloc(nproc-1,sizeof(MPI_Status));
    t_req_r = PGFEM_calloc(nproc-1,sizeof(MPI_Status));
  }

  /* Post recieves from all otehr processes */
  int t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i == myrank)
      continue;
    err += MPI_Irecv(&comm->R[i],1,MPI_LONG,i,MPI_ANY_TAG,
		    mpi_comm,&t_req_r[t_count]);
    t_count++;
  }

  /* Send size to all other processors */
  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank){
      err += MPI_Isend(&comm->S[i],1,MPI_LONG,i,myrank,
		       mpi_comm,&t_req_s[t_count]);
      t_count++;
      if(comm->S[i] > 0) comm->Ns++;
    }
  }

  /* Allocate send space and determine the reduced list of procs to
     send to. */
  {
    long KK = 0;
    if (comm->Ns == 0) KK = 1; else  KK = comm->Ns;
    comm->Nss = (long*) PGFEM_calloc (KK,sizeof(long));
  }

  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->S[i] > 0) comm->Nss[t_count++] = i;
  }

  /* Process received messages as they arrive */
  t_count = 0;
  comm->Nr = 0;
  while (t_count < nproc-1){
    int idx = 0;
    err +=  MPI_Waitany(nproc-1,t_req_r,&idx,&t_sta_r);
    int source = t_sta_r.MPI_SOURCE;
    if (source != myrank && comm->R[source] > 0){
      comm->Nr++;
    }
    t_count++;
  }
  free(t_req_r);


  /* Allocate receive space and determine the reduced list of
     processors to receive from */
  {
    long KK = 0;
    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
    comm->Nrr = (long*) PGFEM_calloc (KK,sizeof(long));
  }

  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->R[i] > 0) comm->Nrr[t_count++] = i;
  }

  /* Wait for send communications to finish */
  err += MPI_Waitall(nproc-1,t_req_s,t_sta_s);

  /* deallocate */
  free(t_req_s);
  free(t_sta_s);

  return err;
}

static int communicate_number_row_col(COMMUN comm,
				      long *NRr,
				      long **GNRr,
				      long **ApRr,
				      const long *LG,
				      const long *ap,
				      long **AA,
				      const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);
  err += MPI_Comm_size(mpi_comm,&nproc);

  /* How many numbers I will send */
  for (int i = 0; i < comm->Ns; i++){
    comm->AS[comm->Nss[i]] = 2*comm->S[comm->Nss[i]];
  }

  /* How many numbers I will receive */
  for(int i = 0; i < comm->Nr; i++){
    comm->AR[comm->Nrr[i]] = 2*comm->R[comm->Nrr[i]];
  }

 
  MPI_Status *sta_s = NULL;
  MPI_Status *sta_r = NULL;
  MPI_Request *req_s = NULL;
  MPI_Request *req_r = NULL;
  long **SEND = NULL;
  long **RECI = NULL;
 
  /* Allocate status and request fields */
  {
    long KK = 0;
    if (comm->Ns == 0) KK = 1; else KK = comm->Ns;
    sta_s = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
    req_s = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));

    if (comm->Nr == 0) KK = 1; else KK = comm->Nr; 
    sta_r = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status));
    req_r = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request));
  }

  if(comm->Ns > 0) SEND = PGFEM_calloc (comm->Ns,sizeof(long*));
  if(comm->Nr > 0) RECI = PGFEM_calloc (comm->Nr,sizeof(long*));
  
  /* =======================================================*/

  /* Post receive */
  for (int i = 0; i < comm->Nr; i++){
    int r_idx = comm->Nrr[i];
    int n_rec = comm->AR[r_idx];

    RECI[i] = PGFEM_calloc(n_rec,sizeof(long));

    err += MPI_Irecv (RECI[i],n_rec,MPI_LONG,r_idx,
		      MPI_ANY_TAG,mpi_comm,&req_r[i]);
  }

  /* Post sends */
  for (int i = 0; i < comm->Ns; i++){
    int s_idx = comm->Nss[i];
    int n_send = comm->AS[s_idx];

    SEND[i] = PGFEM_calloc(n_send,sizeof(long));

    /* populate send buffer */
    for (int j = 0; j < comm->S[s_idx]; j++){
      SEND[i][j] = LG[AA[s_idx][j]]; /* global? row id */
      SEND[i][comm->S[s_idx] + j] = ap[AA[s_idx][j]]; /* # columns */
    }

    err += MPI_Isend (SEND[i],n_send,MPI_LONG,s_idx,
		      myrank,mpi_comm,&req_s[i]);
  }
  
  /* Compute the number of rows to receive */
  {
    long n_rows_recv = 0;
    for (int i = 0; i < comm->Nr;i++)
      n_rows_recv += comm->R[comm->Nrr[i]];

    *NRr = n_rows_recv;
    if(n_rows_recv > 0){
      *GNRr = aloc1l (n_rows_recv);
      *ApRr = aloc1l (n_rows_recv);
    }
  }

  /* THIS IS BAD! RGID/RAp should only be comm->Nr long, however
     propagating this fix is not trivial */
  comm->RGID = (long**) PGFEM_calloc (nproc,sizeof(long*));
  comm->RAp = (long**) PGFEM_calloc (nproc,sizeof(long*));

  for (int i = 0; i < comm->Nr; i++) {
    int KK = comm->Nrr[i];
    comm->RGID[KK] = (long*) PGFEM_calloc (comm->R[KK],sizeof(long));
    comm->RAp[KK] = (long*) PGFEM_calloc (comm->R[KK],sizeof(long));
  }

  
  /* CHANGE TO WAIT ANY */
  /* Wait to complete the communications */
  err += MPI_Waitall (comm->Ns,req_s,sta_s);
  err += MPI_Waitall (comm->Nr,req_r,sta_r);
  

  /* Unpack it */
  int II = 0;
  for (int i = 0; i < comm->Nr; i++){
    int KK = comm->Nrr[i];
    for (int j = 0; j < comm->R[KK]; j++) {
      (*GNRr)[II+j] = comm->RGID[KK][j] = RECI[i][j];
      (*ApRr)[II+j] = comm->RAp[KK][j] = RECI[i][comm->R[KK]+j];
    }
    II += comm->R[KK];
  }
  
  /* Deallocate */
  dealoc2l(SEND,comm->Ns);
  dealoc2l(RECI,comm->Nr);
  free(sta_s);
  free(sta_r);
  free(req_s);
  free(req_r);

  return err;
}

static int communicate_row_info(COMMUN comm,
				long ***GIDRr,
				const long NRr,
				const long *ApRr,
				const long *ap,
				long **AA,
				long **ID,				
				const MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  err += MPI_Comm_rank(mpi_comm,&myrank);
  err += MPI_Comm_size(mpi_comm,&nproc);

  MPI_Status *sta_s = NULL;
  MPI_Status *sta_r = NULL;
  MPI_Request *req_s = NULL;
  MPI_Request *req_r = NULL;

  /* Allocate status and request fields */ 
  {
    int KK = 0;
    if (comm->Ns == 0) KK = 1; else KK = comm->Ns;
    sta_s = PGFEM_calloc (KK,sizeof(MPI_Status));
    req_s = PGFEM_calloc (KK,sizeof(MPI_Request));

    if (comm->Nr == 0) KK = 1; else KK = comm->Nr; 
    sta_r = PGFEM_calloc (KK,sizeof(MPI_Status));
    req_r = PGFEM_calloc (KK,sizeof(MPI_Request));
  }

  /* clear the number of communication */
  memset(comm->AS,0,nproc*sizeof(long));
  memset(comm->AR,0,nproc*sizeof(long));

  /* post receive  for allocation information*/
  for (int i = 0; i < comm->Nr; i++){
    int r_idx = comm->Nrr[i];
    err += MPI_Irecv (&comm->AR[r_idx],1,MPI_LONG,
		      r_idx,MPI_ANY_TAG,mpi_comm,&req_r[i]);
  }

  /* Allocate and compute quantity information to send */
  comm->SAp = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (int i = 0; i < comm->Ns; i++) {
    int idx = comm->Nss[i];
    comm->SAp[idx] = PGFEM_calloc (comm->S[idx],sizeof(long));
  }
  
  for (int i = 0; i < nproc;i++){
    if (i != myrank){
      int ncols = 0;
      for (int j = 0; j < comm->S[i]; j++){
	ncols += comm->SAp[i][j] = ap[AA[i][j]];
      }
      comm->AS[i] = ncols;
    }
  }
  
  /* Send allocation information */
  for (int i = 0; i < comm->Ns; i++){
    int s_idx = comm->Nss[i];
    err += MPI_Isend (&comm->AS[s_idx],1,MPI_LONG,
		      s_idx,myrank,mpi_comm,&req_s[i]);
  }
 
  comm->SGRId = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (int i = 0; i < comm->Ns; i++){
    int KK = comm-> Nss[i];
    comm->SGRId[KK] = (long*) PGFEM_calloc (comm->AS[KK],sizeof(long));
  }

  /* wait for receives to complete */
  err += MPI_Waitall(comm->Nr,req_r,sta_r);

  /* Allocate recieve buffer */
  long **RECI = PGFEM_calloc (nproc,sizeof(long*));
  for (int i = 0; i < nproc; i++) {
    int JJ = 0;
    if (comm->AR[i] == 0) JJ = 1; else JJ = comm->AR[i];
    RECI[i] = PGFEM_calloc (JJ,sizeof(long));
  }

  free(sta_r);
  free(req_r);

  /* Wait for sends to complete */
  err += MPI_Waitall(comm->Ns,req_s,sta_s); 

  /* reallocate status and request fields */
  free(sta_s);
  free(req_s);
  {
    int KK = 0;
    if (comm->Ns == 0) KK = 1; else KK = comm->Ns;
    sta_s = PGFEM_calloc (KK,sizeof(MPI_Status));
    req_s = PGFEM_calloc (KK,sizeof(MPI_Request));

    if (comm->Nr == 0) KK = 1; else KK = comm->Nr; 
    sta_r = PGFEM_calloc (KK,sizeof(MPI_Status));
    req_r = PGFEM_calloc (KK,sizeof(MPI_Request));
  }

  /* post receive */
  for (int i = 0; i < comm->Nr; i++){
    int KK = comm->Nrr[i];
    MPI_Irecv (RECI[KK],comm->AR[KK],MPI_LONG,KK,
	       MPI_ANY_TAG,mpi_comm,&req_r[i]);
  }

  for (int i = 0; i < comm->Ns; i++){
    int KK = comm->Nss[i];
    
    int II = 0;
    for (int j = 0; j < comm->S[KK]; j++){
      for (int k = 0; k < ap[AA[KK][j]]; k++){
	comm->SGRId[KK][II] = ID[AA[KK][j]][k];
	II++;
      }
    }
    
    MPI_Isend (comm->SGRId[KK],comm->AS[KK],MPI_LONG,KK,
	       myrank,mpi_comm,&req_s[i]);
    
  }/* end i < comm->Ns */
  
  /****************/
  /* Recieve data */
  /****************/
  {
    int KK = 0;
    if (NRr == 0) KK = 1; else KK = NRr;
    *GIDRr = (long**) PGFEM_calloc (KK,sizeof(long*));
  }

  for (int i = 0; i < NRr; i++){
    (*GIDRr)[i]= (long*) PGFEM_calloc (ApRr[i],sizeof(long));
  }

  comm->RGRId = (long**) PGFEM_calloc (nproc,sizeof(long*));
  for (int i = 0; i < comm->Nr; i++) {
   int KK = comm->Nrr[i];
    comm->RGRId[KK] = (long*) PGFEM_calloc (comm->AR[KK],sizeof(long));
  }

  /* WAIT ANY and unpack */
  /* Wait to complete the communications */
  MPI_Waitall (comm->Ns,req_s,sta_s);
  MPI_Waitall (comm->Nr,req_r,sta_r);
  
  /* Unpack ir */
  int JJ = 0;
  for (int i = 0; i < comm->Nr;i ++){
    int KK = comm->Nrr[i];
    int II = 0;
    for (int j = 0; j < comm->R[KK]; j++){
      for (int k = 0; k < ApRr[JJ]; k++){
	(*GIDRr)[JJ][k] = comm->RGRId[KK][II] = RECI[KK][II];
	II++;
      }
      JJ++;
    }
  }/* end i < nproc */

  free(sta_s);
  free(sta_r);
  free(req_s);
  free(req_r);

  dealoc2l(RECI,nproc);

  return err;
}
