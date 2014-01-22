/* HEADER */

#include "bounding_element.h"
#include <stdlib.h>

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

typedef long idx_type;

int read_bounding_elements(FILE *ifile,
			   const int ndof_be,
			   int *nbe,
			   BOUNDING_ELEMENT **b_elems,
			   MPI_Comm mpi_comm)
{
  int err = 0;
  int myrank =0;
  MPI_Comm_rank(mpi_comm,&myrank);

  /* read number of bounding elements */
  fscanf(ifile,"%d",nbe);

  /* read the number of nodes per element and build the elements */

  /*=== NOTE: Currently only support triangles, but coding for mixed
    meshes ===*/

  int *nn_per_elem;
  if(*nbe == 0){ /* no elements, explicitly set NULL */
    nn_per_elem = NULL;
    (*b_elems) = NULL;
    return err;
  } else {
    nn_per_elem = PGFEM_calloc (*nbe,sizeof(idx_type));
    for(int i=0; i< (*nbe); i++){
      fscanf(ifile,"%d",&nn_per_elem[i]);
    }
    err += construct_bounding_elements(*nbe,nn_per_elem,b_elems);
  }

  /* read the element connectivity and bounding element id from the
     file. Note that all of the material property information etc is
     inherited from the volumetric element it bounds and is obtained
     through the bounding element id*/
  for(int i=0; i<(*nbe); i++){

    /* read element id and associated vol_element id */
    fscanf(ifile,"%d %d",&(*b_elems)[i].id,
	   &(*b_elems)[i].vol_elem_id);

    /* read connectivity */
    switch ((*b_elems)[i].nnodes){
    case 3: /* triangles */
      fscanf(ifile,"%ld %ld %ld",&(*b_elems)[i].nodes[0],
	     &(*b_elems)[i].nodes[1],&(*b_elems)[i].nodes[2]);
      break;

    default: /* all others */
      for(int j=0; j<(*b_elems)[i].nnodes; j++){
	fscanf(ifile,"%ld",&(*b_elems)[i].nodes[j]);
      }
      break;
    } /* end switch on nnodes */

    /* read normal */
    /* fscanf(ifile,"%lf %lf %lf",&(*b_elems)[i].normal[0], */
    /* 	   &(*b_elems)[i].normal[1],&(*b_elems)[i].normal[2]); */

    /* read periodic */
    fscanf(ifile,"%d",&(*b_elems)[i].periodic);
    if((*b_elems)[i].periodic){
      fscanf(ifile,"%d %d %d %d",
	     &(*b_elems)[i].master_dom,&(*b_elems)[i].master_be_id,
	     &(*b_elems)[i].slave_dom,&(*b_elems)[i].slave_be_id);
      if((*b_elems)[i].master_dom == myrank
	 && (*b_elems)[i].master_be_id == i){
	(*b_elems)[i].master = 1;
      } else {
	(*b_elems)[i].master = 0;
      }
      (*b_elems)[i].n_dofs = ndof_be;
      (*b_elems)[i].other_val = 0;
      (*b_elems)[i].G_dof_ids = PGFEM_calloc((*b_elems)[i].n_dofs,sizeof(long));
      (*b_elems)[i].L_dof_ids = PGFEM_calloc((*b_elems)[i].n_dofs,sizeof(long));
    } else {
      (*b_elems)[i].master = 0;
      (*b_elems)[i].n_dofs = 0;
      (*b_elems)[i].other_val = 0;
      (*b_elems)[i].master_dom = (*b_elems)[i].master_be_id = -1;
      (*b_elems)[i].slave_dom = (*b_elems)[i].slave_be_id = -1;
      (*b_elems)[i].G_dof_ids = NULL;
      (*b_elems)[i].L_dof_ids = NULL;
    }

  }/* end for each element */

  free(nn_per_elem);

  return err;
}

int read_bounding_elements_fname(char *filename,
				 const int ndof_be,
				 int *nbe,
				 BOUNDING_ELEMENT **b_elems,
				 MPI_Comm mpi_comm)
{
  int err = 0;

  FILE *ifile;
  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);

  /* open file to read */
  if((ifile = fopen(filename,"r")) == NULL){
    /* PGFEM_printf("[%d] Could not open file %s (%s)\n",myrank,filename,__func__); */
    /* PGFEM_Abort(); */

    /* if file does not exist, assume there are no bounding
       elements. Do not abort, but return control to calling function
       with err = 1 */
    *nbe = 0;
    (*b_elems) = NULL;
    err = 1;
    return err;
  }

  err += read_bounding_elements(ifile,ndof_be,nbe,b_elems,mpi_comm);

  /* close the file */
  fclose(ifile);

  return err;
}/* read_bounding_elements */

int write_bounding_elements(FILE *ofile,
			    int nbe,
			    BOUNDING_ELEMENT *b_elems)
{
  int err = 0;
  PGFEM_fprintf(ofile,"%d\n\n",nbe);
  for(int i=0; i<nbe; i++){
    PGFEM_fprintf(ofile,"%d ",b_elems[i].nnodes);
  }
  PGFEM_fprintf(ofile,"\n\n");

  for(int i=0; i<nbe; i++){
    /* print id numbers */
    PGFEM_fprintf(ofile,"%5d %5d\t",
	    b_elems[i].id,
	    b_elems[i].vol_elem_id);

    /* print connectivity */
    for(int j=0; j<b_elems[i].nnodes-1; j++){
      PGFEM_fprintf(ofile,"%5ld ",b_elems[i].nodes[j]);
    }
    PGFEM_fprintf(ofile,"%5ld\t",b_elems[i].nodes[b_elems[i].nnodes-1]);

    /* print normal */
    /* PGFEM_fprintf(ofile,"%19.12e %19.12e %19.12e\n",b_elems[i].normal[0], */
    /* 	    b_elems[i].normal[1],b_elems[i].normal[2]); */

    /* print periodicity */
    PGFEM_fprintf(ofile,"%d %d %d\t",b_elems[i].periodic,
	    b_elems[i].master_dom,b_elems[i].master_be_id);

    if(b_elems[i].periodic){
      for(int j=0; j<b_elems[i].n_dofs; j++){
	PGFEM_fprintf(ofile,"%5ld::%-5ld ",b_elems[i].L_dof_ids[j],
		b_elems[i].G_dof_ids[j]);
      }
    }
    PGFEM_fprintf(ofile,"\n");
  } /* for each element */

  return err;
}/* write_bounding_elements */


int write_bounding_elements_fname(char *filename,
				  int nbe,
				  BOUNDING_ELEMENT *b_elems)
{
  int err = 0;
  FILE *ofile = NULL;

  /* open file to read */
  if((ofile = fopen(filename,"w")) == NULL){
    PGFEM_printerr("Could not open file %s (%s)\n",filename,__func__);
    err ++;
    return err;
  }

  err += write_bounding_elements(ofile,nbe,b_elems);

  PGFEM_fclose(ofile);
  return err;
}/* write_bounding_elements_fname */

int construct_bounding_elements(int nbe,
				int *nn_per_elem,
				BOUNDING_ELEMENT **b_elems)
{
  int err = 0;

  /* allocate the overall structure */
  (*b_elems) = PGFEM_calloc (nbe,sizeof(BOUNDING_ELEMENT));

  for(int i=0; i<nbe; i++){
    (*b_elems)[i].nnodes = nn_per_elem[i];
    /* (*b_elems)[i].normal = PGFEM_calloc(3,sizeof(double)); */
    (*b_elems)[i].nodes = PGFEM_calloc ((*b_elems)[i].nnodes,sizeof(idx_type));
    (*b_elems)[i].loc_nodes = PGFEM_calloc ((*b_elems)[i].nnodes,sizeof(idx_type));
  }

  return err;
}/* construct_bounding_elements */

int destroy_bounding_elements(int nbe,
			      BOUNDING_ELEMENT *b_elems)
{
  int err = 0;

  if(b_elems != NULL){
    for(int i=0; i<nbe; i++){
      if(b_elems[i].nodes != NULL){   
	free(b_elems[i].nodes);
      }
      if(b_elems[i].loc_nodes != NULL){   
	free(b_elems[i].loc_nodes);
      }
      /* if(b_elems[i].normal != NULL){ */
      /* 	free(b_elems[i].normal); */
      /* } */
      free(b_elems[i].G_dof_ids);
      free(b_elems[i].L_dof_ids);
    }
    free(b_elems);
  }

  return err;
}/* destroy_bounding_elements */


/** A short driver function to test the functionality of the
    "class".
=================================================================
Compile with the following command:

mpicc -DBOUNDING_ELEMENT_DRIVER -I../include bounding_element.c -std=c99 -Wall -g -o be_driver

Run the driver in a directory containing a valid boundary element
input file named test.in by issuing the command:

mpirun -np 1 ./be_driver

=================================================================
*/

#ifdef BOUNDING_ELEMENT_DRIVER
int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);

  int nbe;
  BOUNDING_ELEMENT *b_elems;

  /* read/allocate the bounding elements */
  read_bounding_elements_fname("test.in",&nbe,&b_elems);

  /* write the elements to stdout */
  write_bounding_elements(stdout,nbe,b_elems);

  /* deallocate the elements */
  destroy_bounding_elements(nbe,b_elems);

  MPI_Finalize();
}
#endif
