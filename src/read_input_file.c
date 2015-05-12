/* HEADER */
#include "read_input_file.h"

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef IN_H
#include "in.h"
#endif

static const int ndim = 3;

int read_input_file(const PGFem3D_opt *opts,
		    MPI_Comm comm,
		    long *nn,
		    long *Gnn,
		    long *ndofn,
		    long *ne,
		    long *lin_maxit,
		    double *lin_err,
		    double *lim_zero,
		    long *nmat,
		    long *n_concentrations,
		    long *n_orient,
		    NODE **node,
		    ELEMENT **elem,
		    MATERIAL **material,
		    MATGEOM *matgeom,
		    SUPP *sup,
		    long *nln,
		    ZATNODE **znod,
		    long *nel_s,
		    ZATELEM **zelem_s,
		    long *nel_v,
		    ZATELEM **zelem_v)
{
  int err = 0;
  int myrank = 0;
  MPI_Comm_rank(comm,&myrank);

  /* compute filename and open file */
  char *filename = PGFEM_calloc(500,sizeof(char));
  sprintf(filename,"%s/%s%d.in",opts->ipath,opts->ifname,myrank);
  FILE *in = PGFEM_fopen(filename,"r");

  /* read header lines */
  fscanf (in,"%ld %ld %ld",nn,ndofn,ne);
  fscanf (in,"%ld %lf %lf",lin_maxit,lin_err,lim_zero);
  fscanf (in,"%ld %ld %ld",nmat,n_concentrations,n_orient);

  /* Set ndofn according to analysis type */
  switch(opts->analysis_type){
  case STABILIZED: case MINI: case MINI_3F: *ndofn = 4; break;
  default: *ndofn = ndim; break;
  }

  (*node) = build_node(*nn,*ndofn);
  (*elem) = build_elem(in,*ne,opts->analysis_type);
  (*material) = PGFEM_calloc(*nmat,sizeof(MATERIAL));
  (*matgeom) = build_matgeom(*n_concentrations,*n_orient);

  *Gnn = read_nodes(in,*nn,*node,opts->legacy,comm);
  /* NOTE: Supports assume only ndim supported dofs per node! */
  *sup = read_supports(in,*nn,ndim,*node);

  read_elem(in,*ne,*elem,*sup,opts->legacy);

  for(int i=0, e=*nmat; i<e; i++){
    if ( read_material(in,*nmat,*material,opts->legacy) ){
      PGFEM_Abort();
    }
  }

  read_matgeom(in,*n_concentrations,*n_orient,*matgeom);

  /* NOTE: Node/Element loading assumes forces only in ndim
     directions */
  /* node */
  fscanf(in,"%ld",nln);
  *znod = build_zatnode (ndim,*nln);
  read_nodal_load (in,*nln,ndim,*znod);
  /* surface */
  fscanf (in,"%ld",nel_s);
  *zelem_s = build_zatelem (ndim,*nel_s);
  read_elem_surface_load (in,*nel_s,ndim,*elem,*zelem_s);
  /* volume */
  fscanf (in,"%ld",nel_v);
  *zelem_v = build_zatelem (ndim,*nel_v);

  /* check the ferror bit */
  if(ferror(in)) err++;

  /* free local memory and close file */
  free(filename);
  fclose(in);
  return err;
}
