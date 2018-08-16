#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFEM_io.h"
#include "Psparse_ApAi.h"
#include "allocation.h"
#include "get_dof_ids_on_elem.h"
#include "get_ndof_on_elem.h"
#include "incl.h"
#include "matice.h"
#include "string.h"
#include "utils.h"

#ifdef HAVE_PHOTON
#include "communication/photon/PhotonNetwork.hpp"
#endif

using namespace pgfem3d;
using namespace pgfem3d::net;

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#ifndef PFEM_DEBUG_ALL
#define PFEM_DEBUG_ALL 0
#endif

#ifndef PFEM_PRINT
#define PFEM_PRINT 0
#endif

/*
 * dof = degrees of freedom
 * # procs and # domains is the same
 */


/**
 * Determine what processors talk to each other and allocate
 * appropriate receive space in comm.  Uses ISIR network pattern.
 *
 * \param[in/out] comm Allocates Nss and Nrr, sets Ns and Nr
 *
 * \return non-zero on error.
 *
 * Side effects: non-blocking communication between all processes.
 */
static int determine_comm_pattern_ISIR(CommunicationStructure *com,
				       const int *preSend,
				       const int *preRecv,
				       const int nsend,
				       const int nrecv);

static int determine_comm_pattern_ISIR(CommunicationStructure *com);

/**
 * Determine what processors talk to each other and allocate
 * appropriate receive space in comm.  Uses PWC network pattern.
 *
 * \param[in/out] comm Allocates Nss and Nrr, sets Ns and Nr
 *
 * \return non-zero on error.
 *
 * Side effects: non-blocking communication between all processes.
 */
static int determine_comm_pattern_PWC(CommunicationStructure *com,
				       const int *preSend,
				       const int *preRecv,
				       const int nsend,
				       const int nrecv);

static int determine_comm_pattern_PWC(CommunicationStructure *com);


/**
 * Communicate the number of rows/columns. Uses ISIR network pattern.
 *
 * \param[in/out] comm Allocates internal space and modifies values
 * \param[out] NRr Returns number of rows to on the domain
 * \param[out] GNRr Allocates space and returns ids of rows that are received
 * \param[out] ApRr Allocates space and returns column ids of communicated entires.
 *
 * \return non-zero on error
 * Side effects: Non-blocking point-to-point communication based on comm
 */
static int communicate_number_row_col_ISIR(CommunicationStructure *com,
					   long *NRr,
					   long **GNRr,
					   long **ApRr,
					   const long *LG,
					   const long *ap,
					   long **AA);

/**
 * Communicate all the row/col numbers. Uses ISIR network pattern.
 *
 * \param[in,out] comm Allocates internal space and modifies values
 * \param[out] GIDRr Allocates space and initializes values
 * \param[in] AA,ID are not modified
 *
 * \return non-zero on error
 * Side-effects: Non-blocking pt2pt communication based on comm
 */
static int communicate_row_info_ISIR(CommunicationStructure *com,
				     long ***GIDRr,
				     const long NRr,
				     const long *ApRr,
				     const long *ap,
				     long **AA,
				     long **ID);

/**
 * Communicate the number of rows/columns. Uses PWC network pattern.
 *
 * \param[in/out] comm Allocates internal space and modifies values
 * \param[out] NRr Returns number of rows to on the domain
 * \param[out] GNRr Allocates space and returns ids of rows that are received
 * \param[out] ApRr Allocates space and returns column ids of communicated entires.
 *
 * \return non-zero on error
 * Side effects: Non-blocking point-to-point communication based on comm
 */
static int communicate_number_row_col_PWC(CommunicationStructure *com,
					   long *NRr,
					   long **GNRr,
					   long **ApRr,
					   const long *LG,
					   const long *ap,
					   long **AA);

/**
 * Communicate all the row/col numbers. Uses PWC network pattern.
 *
 * \param[in,out] comm Allocates internal space and modifies values
 * \param[out] GIDRr Allocates space and initializes values
 * \param[in] AA,ID are not modified
 *
 * \return non-zero on error
 * Side-effects: Non-blocking pt2pt communication based on comm
 */
static int communicate_row_info_PWC(CommunicationStructure *com,
				    long ***GIDRr,
				    const long NRr,
				    const long *ApRr,
				    const long *ap,
				    long **AA,
				    long **ID);

int* Psparse_ApAi (long ne,
                   long n_be,
                   long nn,
                   long ndofn,
                   long ndofd,
                   Element *elem,
                   BoundingElement *b_elems,
                   Node *node,
                   long nce,
                   COEL *coel,
		   CommunicationStructure *com,
                   const int cohesive,
                   const int mp_id)
{
  char jmeno[200];
  FILE *out=NULL;
  long i,j,k,II,JJ,*cnL=NULL,*cnG=NULL,nne,ndofe,*nod=NULL,*ap=NULL;
  long **AA=NULL,*ap1=NULL,**ID=NULL,*LG=NULL,Nrs=0,*GL=NULL;
  long LI1,LI2,*cncL=NULL,*cncG=NULL,ndofc{},*nodc=NULL;
  long *send=NULL,NRr=0,*GNRr=NULL,*ApRr=NULL,**GIDRr=NULL;
  int *Ddof = NULL;
  int *Ai = NULL;

  /* pull necessary info from Communication handle */
  int myrank = com->rank;
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  long *DomDof = com->DomDof;
  int *GDof = &(com->GDof);
  int *Ap = com->Ap;
  
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
  cnL = aloc1l (30);                                                          //dof_ids
  cnG = aloc1l (30);
  nod = aloc1l (10);                                                          //node_ids_on_elem
  cncL = aloc1l (30);
  cncG = aloc1l (30);
  nodc = aloc1l (10);
  ap = aloc1l (ndofd);                                                        //used in multiple places
  ap1 = aloc1l (ndofd);                                                       //used in multiple places
  LG = aloc1l (ndofd);                                                        //local to global

  int nsend = 0;
  int nrecv = 0;
  const int *preSend = NULL;
  const int *preRecv = NULL;

  if (com->hints != NULL){                                                         //check if comm hints were provided
    nsend = com->hints->get_nrecv();
    nrecv = com->hints->get_nsend();
    preSend = com->hints->get_recv_list();                                    //returns hints->recv
    preRecv = com->hints->get_send_list();
  }

  comm->S = PGFEM_calloc (long, nproc);       //amount of rows to send to each processor
  comm->R = PGFEM_calloc (long, nproc);       //amount of rows to receive from each processsor
  comm->AS = PGFEM_calloc (long, nproc);      //amount to send
  comm->AR = PGFEM_calloc (long, nproc);      //amount to recieve
  comm->LG = PGFEM_calloc (long, ndofd);                       //local to global (ndofd)
  comm->GL = PGFEM_calloc (long, DomDof[myrank]);              //global to local

  for (i=0;i<ne;i++){/* Number of contributions from elements to Aij */       //loop over elements

    /* Number of element nodes */
    nne = elem[i].toe;                                                        //find out how many nodes are associated with this element
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);                                               //returns the nodes of this element (array of nodes stored in nod)
    /* Element Dof */
    ndofe = get_total_ndof_on_elem(nne,nod,node,b_elems,&elem[i],ndofn);            //get total degrees of freedom on element, put it into ndofe
    /* Id numbers */
    get_all_dof_ids_on_elem(0,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnL,mp_id); //get dof ids for element nodes, elements, and boundary elements


    /*
     * This section seems to count the number of non-prescribed neighbors for each element.
     * The sum continues across elements. Prescribed elements dont get a number (0).
     */
    for (j=0;j<ndofe;j++){/* row */
      II = cnL[j]-1;
      if (II < 0)                                                             //if this node is part of another domain (prescribed),
        continue;                                                             //ignore it
      for (k=0;k<ndofe;k++){/* column */                                      //else, for this degree of freedom, loop through again
        JJ = cnL[k]-1;
        if (JJ < 0)                                                           //and count the number of non-prescribed degrees of freedom which are in this element
          continue;
        ap[II]++;                                                             //add them to the domain-wide total for this node
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

      ndofe = coel[i].toe*ndofc;
      for (j=0;j<coel[i].toe;j++)
        nodc[j] = coel[i].nod[j];
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nodc,node,cncL,mp_id);
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


  //CREATE AA (Aij) matrix

  AA = PGFEM_calloc (long*, ndofd);
  for (i=0;i<ndofd;i++)
    AA[i]= PGFEM_calloc (long, ap[i]);                         //AA is ndofd*ap large
  null_quit((void*) AA,0);

  for (i=0;i<ne;i++){/* List of ID number for Aij */                          // loop over elements in one domain
    /* Number of element nodes */
    nne = elem[i].toe;                                                        //How many nodes each element is associated with
    /* Nodes on element */
    elemnodes (i,nne,nod,elem);                                               //fills nod, an array which contains the ids of the nodes here
    /* Element Dof */
    ndofe = get_total_ndof_on_elem(nne,nod,node,b_elems,&elem[i],ndofn);
    /* Id numbers */
    get_all_dof_ids_on_elem(0,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnL,mp_id); //get local ids for degrees of freedom
    get_all_dof_ids_on_elem(1,nne,ndofe,ndofn,nod,node,b_elems,&elem[i],cnG,mp_id); //get global ids for degrees of freedom


    /*
     * This loop starts the main Aij matrix. It creates a matrix which is
     * (degrees of freedom for element n)*(degrees of freedom for element n)*(number of elements?) large.
     * Also where ap1 (addresses for the columns) is created.
     * This matrix is filled densely thanks to ap1. Once a row has a value in it, ap1, makes AA
     * go to the next row.
     */
    for (j=0;j<ndofe;j++){/* row */
      II = cnL[j]-1;                                                          //II is the address  (dof ID) for this row (-1). 2 purposes:
      if (II < 0)  continue;                                                  //1.by using the ID, we only have 1 entry per column for this element
      LG[II] = cnG[j]-1;                                                      //2.if it is negative (prescribed), skip it.
      for (k=0;k<ndofe;k++){/* column */                                      //else for each element-friend that also belongs here
        JJ = cnG[k]-1;                                                        //set JJ to the global dof id - 1
        if (JJ < 0)  continue;
        AA[II][ap1[II]] = JJ;                                                 //fills the AA matrix with id's
        ap1[II]++;                                                            //count how many things have been put in this row
      }
    }
  }/* end i < ne */


  /* COHESIVE ELEMENTS */
  if (cohesive == 1){
    for (i=0;i<nce;i++){

      ndofe = coel[i].toe*ndofc;
      for (j=0;j<coel[i].toe;j++)
        nodc[j] = coel[i].nod[j];
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nodc,node,cncL,mp_id);
      get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nodc,node,cncG,mp_id);

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
    ap[k] = 0; /* null ap */                                                  //also reset ap (optimization combined for-loops)
    qsort (AA[k], ap1[k], sizeof(long*), compare_long);                       // sort each row(column?) of (Aij?) with respect to the global dof_ids
  }

  ap[0] = 0;
  for (k=0;k<ndofd;k++){ /* Number of non-zeros in rows */                    //loop over max degrees of freedom
    comm->LG[k] = LG[k];                                                      //move degree of freedom id's to comm
    for (j=0;j<ap1[k]-1;j++)                                                  //loop over number of dofs
      if (AA[k][j] < AA[k][j+1])                                              //if the number less than the following number,
        ap[k]++;                                                              //(as they should be since we just sorted) then add 2 to that ap
    ap[k]++;                                                                  //else (if they are equal) then just add 1
  }

  {
    int n_dup = number_of_duplicates(comm->LG,ndofd,sizeof(long),compare_long);
    if(n_dup){
      PGFEM_printerr("[%d]:ERROR comm->LG contains %d duplicate values!\n",
                     myrank,n_dup);
    }
    /* if(n_dup) PGFEM_Abort(); */
  }

  ID = PGFEM_calloc (long*, ndofd);                           //create an ID matrix, approximately the size of AA
  for (i=0;i<ndofd;i++) {                                                     //loop over domain degrees of freedom
    ID[i] = PGFEM_calloc (long, ap[i]);                         //allocate an ID array for each dof (huge)
    ap[i] = 0; /* null */                                                     //ID is around ndofd*ndofd
  }

  null_quit((void*) ID,0);

  for (k=0;k<ndofd;k++){ /* Global row indexes */                             //loop over all dofs in this domain
    for (j=0;j<ap1[k]-1;j++) {
      if (AA[k][j] < AA[k][j+1]){
        ID[k][ap[k]] = AA[k][j];                                                    //fill ID using AA
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

  Ddof = aloc1i (nproc);                                                      //a total (global) count of the degrees of freedom
  Ddof[0] = DomDof[0];

  for (i=1;i<nproc;i++)                                                       //Ddof is size nproc
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
        Ap[i] = ap[j];                                                              //manually sort ap into Ap
        GL[i] = comm->GL[i] = j;                                                    //also global to local
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

  for (i=0;i<ndofd;i++){/* Number of rows per domain to be sent (Nrs)*/     //loop over domain dof
    if (LI1 <= LG[i] && LG[i] <= LI2)                                       //if this dof belongs here,
      continue;                                                             //ignore it
    else Nrs++;                                                             //else add it to the counter
  }

  if(Nrs > 0){
    ap1 = aloc1l (Nrs);                                                     //allocate an array the size of all the dofs which
                                                                            //were not part of this domain
    Nrs = 0;
    for (i=0;i<ndofd;i++){/* Local row indexes on domain to be sent */      //same as before, but this time
      if (LI1 <= LG[i] && LG[i] <= LI2)                                     //write down which nodes are not
        continue;                                                                 //in this domain
      else {ap1[Nrs] = i; Nrs++;}
    }
  }


  //This is it

  /* Determine where to send */
  for (j=0;j<nproc;j++){                                                      //loop over processes
    k = 1;
    if (j == myrank)                                                          //if this is the current process, theres nothing to be sent
      continue;
    if (j == 0)                                                               //(first) process is special
      LI1 = 0;
    else
      LI1 = Ddof[j-1];                                                        //using the prefix sum, we can categorize which
    LI2 = Ddof[j] - 1;                                                        //nproc a particular dof is in.
    for (i=0;i<Nrs;i++){                                                      //loop over all rows to be sent
      if (LI1 <= LG[ap1[i]] && LG[ap1[i]] <= LI2) {                           //if this dof belongs to proc j,
        comm->S[j] = k;                                                       //write down how many things will be sent there
        k++;
      }
    }
  }

  AA = PGFEM_calloc (long*, nproc);
  for (i=0;i<nproc;i++) {                                                     //loop over procs
    if (myrank == i || comm->S[i] == 0)                                       //if this is my rank, or if there I dont share with this proc,
      k = 1;                                                                  //just allocate 1 space
    else
      k = comm->S[i];                                                         // allocate enough memory to send the appropriate amount of info
    AA[i] = PGFEM_calloc (long, k);
  }

  null_quit((void*) AA,0);

  comm->SLID = PGFEM_calloc (long*, nproc);
  for (i=0;i<nproc;i++) {                                                     //loop over nproc
    if (myrank == i || comm->S[i] == 0)
      continue;
    else comm->SLID[i] = PGFEM_calloc (long, comm->S[i]);      //allocate more memory for
  }                                                                           //"local ID of communicated rows"

  null_quit((void*) comm->SLID,0);

  /* Determine what GID/LID to send */
  for (j=0;j<nproc;j++){                                                      //loop over nproc
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
        AA[j][k] = comm->SLID[j][k] = ap1[i];                                       //same as before, but this time fill in
        k++;                                                                        //the local ID of communicated nodes
      }
    }
  }

  /* Communicate who I am communicating with
   *============================================= */
  switch(com->net->type()) {
  case NET_ISIR:
    if (com->hints == NULL)                                          //checks if comm hints weren't provided
      determine_comm_pattern_ISIR(com);
    else
      determine_comm_pattern_ISIR(com, preSend, preRecv, nsend, nrecv);
    
    /* Communicate how many rows/columns I am sending
     *============================================= */
    communicate_number_row_col_ISIR(com,&NRr,&GNRr,&ApRr,
				    LG,ap,AA);
    
    /* Communicate the row/column information
     *============================================= */
    communicate_row_info_ISIR(com,&GIDRr,NRr,ApRr,
			      ap,AA,ID);
    break;
  case NET_PWC:
    if (com->hints == NULL)                                          //checks if comm hints weren't provided
      determine_comm_pattern_PWC(com);
    else
      determine_comm_pattern_PWC(com, preSend, preRecv, nsend, nrecv);
    
    /* Communicate how many rows/columns I am sending
     *============================================= */
    communicate_number_row_col_PWC(com,&NRr,&GNRr,&ApRr,
				    LG,ap,AA);
    
    /* Communicate the row/column information
     *============================================= */
    communicate_row_info_PWC(com,&GIDRr,NRr,ApRr,
			      ap,AA,ID);
    break;
  default:
    PGFEM_printerr("[%d]ERROR: Unknown network type", com->rank);
    PGFEM_Abort();
  }
    
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
    ID[GL[GNRr[i]-*GDof]]= PGFEM_calloc (long, JJ);

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
                     " GID : %ld %ld %ld %ld\n", i,node[i].id_map[mp_id].id[0],
                     node[i].id_map[mp_id].id[1],node[i].id_map[mp_id].id[2],node[i].id_map[mp_id].id[3],
                     node[i].id_map[mp_id].Gid[0],
                     node[i].id_map[mp_id].Gid[1],node[i].id_map[mp_id].Gid[2],node[i].id_map[mp_id].Gid[3]);
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

  return (Ai);
}


//If comm hints were provided
static int determine_comm_pattern_ISIR(CommunicationStructure *com,
				       const int *preSend,
				       const int *preRecv,
				       const int nsend,
				       const int nrecv)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(com->net);
  
  int err = 0;
  int myrank = com->rank;
  int nproc = com->nproc;
  //mpi_status contains 3 things: the rank of the sender, tag of the message, and the length of the message
  Status t_sta_r;
  Status *t_sta_s = NULL;
  Request *t_req_s = NULL;
  Request *t_req_r = NULL;

  SparseComm *comm = com->spc;
  
  if(nproc > 1){
    net->allocStatusArray(nproc-1, &t_sta_s);
    net->allocRequestArray(nproc-1, &t_req_s);
    net->allocRequestArray(nproc-1, &t_req_r);
  }

  int recvFrom;
  int t_count;
  for (int i = 0; i < nsend; i++){                                               //prepares mailboxes to receive from all nodes
    recvFrom = preSend[i];
    try {
      net->irecv(&comm->R[recvFrom], 1, NET_DT_LONG, recvFrom, NET_ANY_TAG,  //put received info in comm->R
		      com->comm, t_req_r+i);                              //save info of proc from which things came
    } catch(...) {
      err++;
    }
  }
  //can be changed to smaller comm
  /* Send size to all other processors */

  comm->Ns = nrecv;
  for (int i = 0; i < nrecv; i++){                                              //send size to all

    int sendTo = preRecv[i];
    try {
      net->isend(&comm->S[sendTo], 1, NET_DT_LONG, sendTo, myrank,
		      com->comm, &t_req_s[i]);                                         //t_req_s is required for each non-blocking call
    } catch(...) {
      err++;
    }

    if (comm->S[sendTo] == 0) comm->Ns--;                                            //calculate number of procs that I sent to
  }

  /* Allocate send space and determine the reduced list of procs to
     send to. */
  {
    long KK = 0;
    if (comm->Ns == 0) KK = 1; else  KK = comm->Ns;                             //if theres neighbors, set KK equal to the number of procs I sent to
    comm->Nss = PGFEM_calloc (long, KK);
  }

  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->S[i] > 0) comm->Nss[t_count++] = i;                //not sure why comm->S > 0 check is there since line 407
  }                                                                             //guarantees that

  /* Process received messages as they arrive */
  t_count = 0;
  comm->Nr = nsend;
  while (t_count < nsend) //wait until ive heard back from everyone
  {
    int idx = 0;
    try {
      net->waitany(nsend, t_req_r, &idx, &t_sta_r);                         //listen for messages
    } catch(...) {
      err++;
    }
      
    int source = t_sta_r.NET_SOURCE;

    if(comm->R[source]==0)
      comm->Nr--;                                                               //write down how many non-empty letters I received
    
    t_count++;
  }
  delete [] t_req_r;


  /* Allocate receive space and determine the reduced list of
     processors to receive from */
  {
    long KK = 0;
    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
    comm->Nrr = PGFEM_calloc (long, KK);                         //allocate memory for number of responses
  }

  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->R[i] > 0) comm->Nrr[t_count++] = i;                //write down number of responses for each domain
  }                                                                             //that I talked to

  /* Wait for send communications to finish */
  try {
    net->waitall(nrecv, t_req_s, t_sta_s);                                  //wait until all expected messages arrive
  } catch(...) {
    err++;
  }

  /* deallocate */
  delete [] t_req_s;
  delete [] t_sta_s;

  return err;
}


//If comm hints were not provided
static int determine_comm_pattern_ISIR(CommunicationStructure *com)
{
  int err = 0;
  int myrank = com->rank;
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(com->net);
  
  Status t_sta_r;
  Status *t_sta_s = NULL;
  Request *t_req_s = NULL;
  Request *t_req_r = NULL;
  if(nproc > 1){
    net->allocStatusArray(nproc-1, &t_sta_s);
    net->allocRequestArray(nproc-1, &t_req_s);
    net->allocRequestArray(nproc-1, &t_req_r);
  }

  /* Post recieves from all otehr processes */
  int t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i == myrank)
      continue;
    try {
      net->irecv(&comm->R[i], 1, NET_DT_LONG, i, NET_ANY_TAG,
		 com->comm, &t_req_r[t_count]);
    } catch(...) {
      err++;
    }
    t_count++;
  }

  /* Send size to all other processors */
  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank){
      try {
	net->isend(&comm->S[i], 1, NET_DT_LONG, i, myrank,
		   com->comm, &t_req_s[t_count]);
      } catch(...) {
	err++;
      }
      t_count++;
      if(comm->S[i] > 0) comm->Ns++;
    }
  }
  
  /* Allocate send space and determine the reduced list of procs to
     send to. */
  {
    long KK = 0;
    if (comm->Ns == 0) KK = 1; else  KK = comm->Ns;
    comm->Nss = PGFEM_calloc (long, KK);
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
    try {
      net->waitany(nproc-1, t_req_r, &idx, &t_sta_r);
    } catch(...) {
      err++;
    }
    int source = t_sta_r.NET_SOURCE;
    if (source != myrank && comm->R[source] > 0){
      comm->Nr++;
    }
    t_count++;
  }
  delete [] t_req_r;

  /* Allocate receive space and determine the reduced list of
     processors to receive from */
  {
    long KK = 0;
    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
    comm->Nrr = PGFEM_calloc (long, KK);
  }

  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->R[i] > 0) comm->Nrr[t_count++] = i;
  }

  /* Wait for send communications to finish */
  try {
    net->waitall(nproc-1, t_req_s, t_sta_s);
  } catch(...) {
    err++;
  }
  
  /* deallocate */
  delete [] t_req_s;
  delete [] t_sta_s;
  
  return err;
}

static int determine_comm_pattern_PWC(CommunicationStructure *com,
				       const int *preSend,
				       const int *preRecv,
				       const int nsend,
				       const int nrecv)
{
  int myrank = com->rank;
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  PWCNetwork *net = static_cast<PWCNetwork*>(com->net);
  CID lid = 0xcafebabe;

  comm->Ns = nrecv;
  
  // send the local S value to each pre-determined rank
  for (int p = 0; p < nrecv; p++) {
    int sendTo = preRecv[p];
    CID rid = (CID)comm->S[sendTo];
    net->pwc(sendTo, 0, 0, 0, lid, rid);
    if (comm->S[sendTo] == 0) comm->Ns--;
  }
  
  /* Allocate send space and determine the reduced list of procs to
     send to. */
  {
    long KK = 0;
    if (comm->Ns == 0) KK = 1; else  KK = comm->Ns;
    comm->Nss = PGFEM_calloc (long, KK);
  }

  int t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->S[i] > 0)
      comm->Nss[t_count++] = i;
  }

  // wait for local PWC completions from above
  net->wait_n_id(nrecv, lid);
  
  // Process received message as they arrive
  t_count = 0;
  comm->Nr = nsend;
  while (t_count < nsend) {
    int flag;
    CID val;
    Status stat;
    net->probe(&flag, &val, &stat, 0);
    if (flag) {
      int p = stat.NET_SOURCE;
      // the long values exchanged are encoded in the PWC RIDs
      comm->R[p] = (long)val;
      if (comm->R[p] == 0){
	comm->Nr--;
      }
      t_count++;
    }
  }

  /* Allocate receive space and determine the reduced list of
     processors to receive from */
  {
    long KK = 0;
    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
    comm->Nrr = PGFEM_calloc (long, KK);
  }

  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->R[i] > 0)
      comm->Nrr[t_count++] = i;
  }
  
  return 0;
}

static int determine_comm_pattern_PWC(CommunicationStructure *com)
{
  int myrank = com->rank;
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  PWCNetwork *net = static_cast<PWCNetwork*>(com->net);
  CID lid = 0xcafebabe;
  
  // send the local S[p] value to each associated rank
  for (int p = 0; p < nproc; p++) {
    if (p == myrank) continue;
    CID rid = (CID)comm->S[p];
    net->pwc(p, 0, 0, 0, lid, rid);
    if(comm->S[p] > 0) comm->Ns++;
  }
  
  /* Allocate send space and determine the reduced list of procs to
     send to. */
  {
    long KK = 0;
    if (comm->Ns == 0) KK = 1; else  KK = comm->Ns;
    comm->Nss = PGFEM_calloc (long, KK);
  }

  int t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->S[i] > 0)
      comm->Nss[t_count++] = i;
  }

  // wait for local PWC completions from above
  net->wait_n_id(nproc-1, lid);
  
  // Process received message as they arrive
  t_count = 0;
  while (t_count < nproc-1) {
    int flag;
    CID val;
    Status stat;
    net->probe(&flag, &val, &stat, 0);
    if (flag) {
      int p = stat.NET_SOURCE;
      // the long values exchanged are encoded in the PWC RIDs
      comm->R[p] = (long)val;
      if (p != myrank && comm->R[p] > 0){
	comm->Nr++;
      }
      t_count++;
    }
  }

  /* Allocate receive space and determine the reduced list of
     processors to receive from */
  {
    long KK = 0;
    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
    comm->Nrr = PGFEM_calloc (long, KK);
  }

  t_count = 0;
  for (int i = 0; i < nproc; i++){
    if (i != myrank && comm->R[i] > 0)
      comm->Nrr[t_count++] = i;
  }
  
  return 0;
}

static int communicate_number_row_col_ISIR(CommunicationStructure *com,
					   long *NRr,
					   long **GNRr,
					   long **ApRr,
					   const long *LG,
					   const long *ap,
					   long **AA)
{
  int err = 0;
  int myrank = com->rank;
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(com->net);
  
  /* How many numbers I will send */
  for (int i = 0; i < comm->Ns; i++){                                           //loop over number of procs to send to
    comm->AS[comm->Nss[i]] = 2*comm->S[comm->Nss[i]];                           //number to send was in comm->S
  }

  /* How many numbers I will receive */
  for(int i = 0; i < comm->Nr; i++){                                            //loop over number to receive
    comm->AR[comm->Nrr[i]] = 2*comm->R[comm->Nrr[i]];                           //number to receive was in comm->R
  }

  Status *sta_s = NULL;                                                         //storing send message info
  Status *sta_r = NULL;                                                         //storing receive message info
  Request *req_s = NULL;                                                        //required for
  Request *req_r = NULL;                                                        //nonblocking comms
  long **SEND = NULL;
  long **RECI = NULL;

  /* Allocate status and request fields */
  {
    long KK = 0;
    if (comm->Ns == 0) KK = 1; else KK = comm->Ns;                              //if number to send is 0 then allocate 1 space
    net->allocStatusArray(KK, &sta_s);
    net->allocRequestArray(KK, &req_s);
    
    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;                              //if number to receive is 0 then allocate 1 space
    net->allocStatusArray(KK, &sta_r);
    net->allocRequestArray(KK, &req_r);
  }

  if(comm->Ns > 0) SEND = PGFEM_calloc (long*, comm->Ns);                       //similar with SEND
  if(comm->Nr > 0) RECI = PGFEM_calloc (long*, comm->Nr);                       //and RECI

  /* =======================================================*/

  /* Post receive */
  for (int i = 0; i < comm->Nr; i++){                                           //loop over number to receive
    int r_idx = comm->Nrr[i];                                                   //idx is # being received at i
    int n_rec = comm->AR[r_idx];                                                //rec is the amount to receive for idx things

    RECI[i] = PGFEM_calloc(long, n_rec);                                        //create a mailbox big enough for AR

    try {
      net->irecv(RECI[i], n_rec, NET_DT_LONG, r_idx,                       //post mailboxes
		      NET_ANY_TAG, com->comm, &req_r[i]);
    } catch(...) {
      err++;
    }
  }

  /* Post sends */
  for (int i = 0; i < comm->Ns; i++){
    int s_idx = comm->Nss[i];                                                   //idx is # being sent to i
    int n_send = comm->AS[s_idx];                                               //send is the amount to send for idx things

    SEND[i] = PGFEM_calloc(long, n_send);                                       //create envelopes

    /* populate send buffer */
    for (int j = 0; j < comm->S[s_idx]; j++){
      SEND[i][j] = LG[AA[s_idx][j]]; /* global? row id */                       //write letters
      SEND[i][comm->S[s_idx] + j] = ap[AA[s_idx][j]]; /* # columns */           //
    }

    try {
      net->isend(SEND[i], n_send, NET_DT_LONG, s_idx,                      //send letters
                      myrank, com->comm, &req_s[i]);
    } catch(...) {
      err++;
    }
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
  comm->RGID = PGFEM_calloc (long*, nproc);
  comm->RAp = PGFEM_calloc (long*, nproc);

  for (int i = 0; i < comm->Nr; i++) {
    int KK = comm->Nrr[i];
    comm->RGID[KK] = PGFEM_calloc (long, comm->R[KK]);
    comm->RAp[KK] = PGFEM_calloc (long, comm->R[KK]);
  }


  /* CHANGE TO WAIT ANY */
  /* Wait to complete the communications */
  try {
    net->waitall(comm->Ns, req_s, sta_s);
    net->waitall(comm->Nr, req_r, sta_r);
  } catch(...) {
    err++;
  }
  
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
  delete [] sta_s;
  delete [] sta_r;
  delete [] req_s;
  delete [] req_r;

  return err;
}

static int communicate_number_row_col_PWC(CommunicationStructure *com,
					   long *NRr,
					   long **GNRr,
					   long **ApRr,
					   const long *LG,
					   const long *ap,
					   long **AA)
{
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  pwc::PhotonNetwork *net = static_cast<pwc::PhotonNetwork*>(com->net);
  CID lid = 0xcafebabe;

  /* How many numbers I will send */
  for (int i = 0; i < comm->Ns; i++){                                           //loop over number of procs to send to
    comm->AS[comm->Nss[i]] = 2*comm->S[comm->Nss[i]];                           //number to send was in comm->S
  }

  /* How many numbers I will receive */
  for(int i = 0; i < comm->Nr; i++){                                            //loop over number to receive
    comm->AR[comm->Nrr[i]] = 2*comm->R[comm->Nrr[i]];                           //number to receive was in comm->R
  }

  long **SEND = NULL;
  long **RECI = NULL;

  if(comm->Ns > 0) SEND = PGFEM_calloc (long*, comm->Ns);                       //similar with SEND
  RECI = PGFEM_calloc (long*, nproc);                       //and RECI

  /* =======================================================*/

  /* Allocate recieve buffer */
  Buffer *rbuffers = PGFEM_calloc (Buffer, nproc);
  for (int i = 0; i < nproc; i++) {
    int n_rec = 0;
    if (comm->AR[i] == 0) n_rec = 1; else n_rec = comm->AR[i];

    RECI[i] = PGFEM_calloc_pin (long, n_rec,
				net, &rbuffers[i].key);
    rbuffers[i].addr = reinterpret_cast<uintptr_t> (RECI[i]);
    rbuffers[i].size = sizeof(long)*n_rec;
  }

  /* Exchange receive buffers */
  for (int i = 0; i < nproc; i++) {
    net->gather(&rbuffers[i], sizeof(Buffer), NET_DT_BYTE,
        net->wbuf, sizeof(Buffer), NET_DT_BYTE, i, com->comm);
  }

  /* Post sends */
  Buffer *sbuffers = PGFEM_calloc (Buffer, nproc);
  for (int i = 0; i < comm->Ns; i++){
    int s_idx = comm->Nss[i];                                                   //idx is # being sent to i
    int n_send = comm->AS[s_idx];                                               //send is the amount to send for idx things

    SEND[i] = PGFEM_calloc_pin (long, n_send,
				net, &sbuffers[s_idx].key);
    sbuffers[s_idx].addr = reinterpret_cast<uintptr_t> (SEND[i]);
    sbuffers[s_idx].size = sizeof(long)*n_send;

    /* populate send buffer */
    for (int j = 0; j < comm->S[s_idx]; j++){
      SEND[i][j] = LG[AA[s_idx][j]]; /* global? row id */                       //write letters
      SEND[i][comm->S[s_idx] + j] = ap[AA[s_idx][j]]; /* # columns */           //
    }

    CID rid = (CID)n_send;
    net->pwc(s_idx, sbuffers[s_idx].size, &sbuffers[s_idx], &net->wbuf[s_idx], lid, rid);
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
  comm->RGID = PGFEM_calloc (long*, nproc);
  comm->RAp = PGFEM_calloc (long*, nproc);

  for (int i = 0; i < comm->Nr; i++) {
    int KK = comm->Nrr[i];
    comm->RGID[KK] = PGFEM_calloc (long, comm->R[KK]);
    comm->RAp[KK] = PGFEM_calloc (long, comm->R[KK]);
  }

  /* Wait to complete the communications */
  net->wait_n_id(comm->Ns, lid);
  
  /* Unpack it */
  int t_count = 0;
  while (t_count < comm->Nr) {
    int flag;
    CID val;
    Status stat;
    net->probe(&flag, &val, &stat, 0);
    if (flag) {
      t_count++;
    }
  }

  int II = 0;
  for (int i = 0; i < comm->Nr; i++){
    int KK = comm->Nrr[i];
    for (int j = 0; j < comm->R[KK]; j++) {
      (*GNRr)[II+j] = comm->RGID[KK][j] = RECI[KK][j];
      (*ApRr)[II+j] = comm->RAp[KK][j] = RECI[KK][comm->R[KK]+j];
    }
    II += comm->R[KK];
  }
  
  net->barrier(com->comm);

  /* Deallocate */
  for (int i = 0; i < nproc; i++) {
    net->unpin(reinterpret_cast<void *> (rbuffers[i].addr), rbuffers[i].size);
  }

  for (int i = 0; i < comm->Ns; i++){
    int s_idx = comm->Nss[i];
    net->unpin(reinterpret_cast<void *> (sbuffers[s_idx].addr), sbuffers[s_idx].size);
  }

  dealoc1l(rbuffers);
  dealoc1l(sbuffers);
  dealoc2l(SEND,comm->Ns);
  dealoc2l(RECI,nproc);

  return 0;
}

static int communicate_row_info_ISIR(CommunicationStructure *com,
				    long ***GIDRr,
				    const long NRr,
				    const long *ApRr,
				    const long *ap,
				    long **AA,
				    long **ID)
{
  int err = 0;
  int myrank = com->rank;
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(com->net);
  
  Status *sta_s = NULL;
  Status *sta_r = NULL;
  Request *req_s = NULL;
  Request *req_r = NULL;

  /* Allocate status and request fields */
  {
    int KK = 0;
    if (comm->Ns == 0) KK = 1; else KK = comm->Ns;
    net->allocStatusArray(KK, &sta_s);
    net->allocRequestArray(KK, &req_s);

    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
    net->allocStatusArray(KK, &sta_r);
    net->allocRequestArray(KK, &req_r);
  }

  /* clear the number of communication */
  memset(comm->AS,0,nproc*sizeof(long));
  memset(comm->AR,0,nproc*sizeof(long));

  /* post receive for allocation information*/
  for (int i = 0; i < comm->Nr; i++){
    int r_idx = comm->Nrr[i];
    try {
      net->irecv(&comm->AR[r_idx], 1 ,NET_DT_LONG,
		 r_idx, NET_ANY_TAG, com->comm, &req_r[i]);
    } catch (...) {
      err++;
    }
  }

  /* Allocate and compute quantity information to send */
  comm->SAp = PGFEM_calloc (long*, nproc);
  for (int i = 0; i < comm->Ns; i++) {
    int idx = comm->Nss[i];
    comm->SAp[idx] = PGFEM_calloc (long, comm->S[idx]);
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
    try {
      net->isend(&comm->AS[s_idx], 1 ,NET_DT_LONG,
		 s_idx, myrank, com->comm, &req_s[i]);
    } catch(...) {
      err++;
    }
  }

  comm->SGRId = PGFEM_calloc (long*, nproc);
  for (int i = 0; i < comm->Ns; i++){
    int KK = comm-> Nss[i];
    comm->SGRId[KK] = PGFEM_calloc (long, comm->AS[KK]);
  }

  /* wait for receives to complete */
  try {
    net->waitall(comm->Nr, req_r, sta_r);
  } catch(...) {
    err++;
  }
  
  /* Allocate recieve buffer */
  long **RECI = PGFEM_calloc (long*, nproc);
  for (int i = 0; i < nproc; i++) {
    int JJ = 0;
    if (comm->AR[i] == 0) JJ = 1; else JJ = comm->AR[i];
    RECI[i] = PGFEM_calloc (long, JJ);
  }

  delete [] sta_r;
  delete [] req_r;

  /* Wait for sends to complete */
  try {
    net->waitall(comm->Ns, req_s, sta_s);
  } catch(...) {
    err++;
  }

  /* reallocate status and request fields */
  delete [] sta_s;
  delete [] req_s;
  {
    int KK = 0;
    if (comm->Ns == 0) KK = 1; else KK = comm->Ns;
    net->allocStatusArray(KK, &sta_s);
    net->allocRequestArray(KK, &req_s);

    if (comm->Nr == 0) KK = 1; else KK = comm->Nr;
    net->allocStatusArray(KK, &sta_r);
    net->allocRequestArray(KK, &req_r);
  }

  /* post receive */
  for (int i = 0; i < comm->Nr; i++){
    int KK = comm->Nrr[i];
    net->irecv(RECI[KK], comm->AR[KK], NET_DT_LONG, KK,
	       NET_ANY_TAG, com->comm, &req_r[i]);
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

    net->isend(comm->SGRId[KK], comm->AS[KK], NET_DT_LONG, KK,
	       myrank, com->comm, &req_s[i]);
    
  }/* end i < comm->Ns */

  /****************/
  /* Recieve data */
  /****************/
  {
    int KK = 0;
    if (NRr == 0) KK = 1; else KK = NRr;
    *GIDRr = PGFEM_calloc (long*, KK);
  }

  for (int i = 0; i < NRr; i++){
    (*GIDRr)[i]= PGFEM_calloc (long, ApRr[i]);
  }

  comm->RGRId = PGFEM_calloc (long*, nproc);
  for (int i = 0; i < comm->Nr; i++) {
    int KK = comm->Nrr[i];
    comm->RGRId[KK] = PGFEM_calloc (long, comm->AR[KK]);
  }

  /* WAIT ANY and unpack */
  /* Wait to complete the communications */
  net->waitall(comm->Ns, req_s, sta_s);
  net->waitall(comm->Nr,req_r,sta_r);

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


  delete [] sta_s;
  delete [] sta_r;
  delete [] req_s;
  delete [] req_r;

  dealoc2l(RECI,nproc);

  return 0;
}

static int communicate_row_info_PWC(CommunicationStructure *com,
				     long ***GIDRr,
				     const long NRr,
				     const long *ApRr,
				     const long *ap,
				     long **AA,
				     long **ID)
{
  int myrank = com->rank;
  int nproc = com->nproc;
  SparseComm *comm = com->spc;
  pwc::PhotonNetwork *net = static_cast<pwc::PhotonNetwork*>(com->net);
  CID lid = 0xcafebabe;

  /* clear the number of communication */
  memset(comm->AS,0,nproc*sizeof(long));
  memset(comm->AR,0,nproc*sizeof(long));

  /* Allocate and compute quantity information to send */
  comm->SAp = PGFEM_calloc (long*, nproc);
  for (int i = 0; i < comm->Ns; i++) {
    int idx = comm->Nss[i];
    comm->SAp[idx] = PGFEM_calloc (long, comm->S[idx]);
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
    CID rid = (CID)comm->AS[s_idx];
    net->pwc(s_idx, 0, 0, 0, lid, rid);
  }

  comm->SGRId = PGFEM_calloc (long*, nproc);
  Buffer *sbuffers = PGFEM_calloc (Buffer, nproc);
  for (int i = 0; i < comm->Ns; i++){
    int KK = comm-> Nss[i];
    comm->SGRId[KK] = PGFEM_calloc_pin (long, comm->AS[KK],
				        net, &sbuffers[KK].key);
    sbuffers[KK].addr = reinterpret_cast<uintptr_t> (comm->SGRId[KK]);
    sbuffers[KK].size = sizeof(long)*comm->AS[KK];
    //net->pin(reinterpret_cast<void *> (sbuffers[KK].addr), sbuffers[KK].size, &sbuffers[KK].key);
  }

  /* wait for local PWC completions from above */
  net->wait_n_id(comm->Ns, lid);
 
  int t_count = 0;
  while (t_count < comm->Nr) {
    int flag;
    CID val;
    Status stat;
    net->probe(&flag, &val, &stat, 0);
    if (flag) {
      int p = stat.NET_SOURCE;
      // the long values exchanged are encoded in the PWC RIDs
      comm->AR[p] = (long)val;
      t_count++;
    }
  }
 
  /* Allocate recieve buffer */
  long **RECI = PGFEM_calloc (long*, nproc);
  Buffer *rbuffers = PGFEM_calloc (Buffer, nproc);
  for (int i = 0; i < nproc; i++) {
    int JJ = 0;
    if (comm->AR[i] == 0) JJ = 1; else JJ = comm->AR[i];
    RECI[i] = PGFEM_calloc_pin (long, JJ,
				net, &rbuffers[i].key);
    //RECI[i] = PGFEM_calloc (long, JJ);
    rbuffers[i].addr = reinterpret_cast<uintptr_t> (RECI[i]);
    rbuffers[i].size = sizeof(long)*JJ;
    //net->pin(reinterpret_cast<void *> (rbuffers[i].addr), rbuffers[i].size, &rbuffers[i].key);
  }

  /* Exchange receive buffers */
  for (int i = 0; i < nproc; i++) {
    net->gather(&rbuffers[i], sizeof(Buffer), NET_DT_BYTE,
	net->wbuf, sizeof(Buffer), NET_DT_BYTE, i, com->comm);
  }

  /* Send data with pwc */
  for (int i = 0; i < comm->Ns; i++){
    int KK = comm->Nss[i];
    
    int II = 0;
    for (int j = 0; j < comm->S[KK]; j++){
      for (int k = 0; k < ap[AA[KK][j]]; k++){
        comm->SGRId[KK][II] = ID[AA[KK][j]][k];
        II++;
      }
    }
    CID rid = (CID)comm->AS[KK];
    net->pwc(KK, sbuffers[KK].size, &sbuffers[KK], &net->wbuf[KK], lid, rid);
  }/* end i < comm->Ns */

  /****************/
  /* Recieve data */
  /****************/
  {
    int KK = 0;
    if (NRr == 0) KK = 1; else KK = NRr;
    *GIDRr = PGFEM_calloc (long*, KK);
  }

  for (int i = 0; i < NRr; i++){
    (*GIDRr)[i]= PGFEM_calloc (long, ApRr[i]);
  }

  comm->RGRId = PGFEM_calloc (long*, nproc);
  for (int i = 0; i < comm->Nr; i++) {
    int KK = comm->Nrr[i];
    comm->RGRId[KK] = PGFEM_calloc (long, comm->AR[KK]);
  }

  /* WAIT ANY and unpack */
  /* Wait to complete the communications */
  net->wait_n_id(comm->Ns, lid);

  t_count = 0;
  while (t_count < comm->Nr) {
    int flag;
    CID val;
    Status stat;
    net->probe(&flag, &val, &stat, 0);
    if (flag) {
      t_count++;
    }
  }

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
  }

  //net->barrier(com->comm);
  for (int i = 0; i < nproc; i++) {
    net->unpin(reinterpret_cast<void *> (rbuffers[i].addr), rbuffers[i].size);
  }

  for (int i = 0; i < comm->Ns; i++){
    int KK = comm->Nss[i];
    net->unpin(reinterpret_cast<void *> (sbuffers[KK].addr), sbuffers[KK].size);
  }

  dealoc1l(rbuffers);
  dealoc1l(sbuffers);
  dealoc2l(RECI,nproc);

  return 0;
}
