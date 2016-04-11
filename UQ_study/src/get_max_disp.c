#include "allocation.h"
#include "homogen.h"
#include "utils.h"
#include "constitutive_model.h"

#include "read_input_file.h"
#include "post_processing.h"
#include "enumerations.h"
#include "PGFem3D_to_VTK.hpp"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define DEBUG_PRINT 0

/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/

double sphere_on_four_pt(double *p)
{
  int id1[] = {1,2,3};
  int id2[] = {0,1,0};
    
  double x[4],y[4],z[4];
  Matrix(double) pt;
  pt.m_row = 4; pt.m_col = 3; pt.m_pdata = p;
  
  for(int a=0; a<4; a++)
  {
    x[a] = Mat_v(pt,a+1,1);
    y[a] = Mat_v(pt,a+1,2);
    z[a] = Mat_v(pt,a+1,3);
  }
  
  Matrix(double) A,AI,B,cp,dr;
  Matrix_construct_redim(double,A, 3,3);
  Matrix_construct_redim(double,B, 3,1);
  Matrix_construct_redim(double,AI,3,3);
  Matrix_construct_redim(double,cp,3,1);
  Matrix_construct_redim(double,dr,4,3);
  
  for(int a=0; a<3; a++)
  {
    Mat_v(A,a+1,1) = x[id1[a]] - x[id2[a]];
    Mat_v(A,a+1,2) = y[id1[a]] - y[id2[a]];
    Mat_v(A,a+1,3) = z[id1[a]] - z[id2[a]];
    Mat_v(B,a+1,1) = x[id1[a]]*x[id1[a]] - x[id2[a]]*x[id2[a]]
                   + y[id1[a]]*y[id1[a]] - y[id2[a]]*y[id2[a]]
                   + z[id1[a]]*z[id1[a]] - z[id2[a]]*z[id2[a]];
  }

  Matrix_inv(A,AI);
  Matrix_AxB(cp,0.5,0.0,AI,0,B,0);

  for(int a=0; a<4; a++)
  {
    Mat_v(dr,a+1,1) = (x[a] - Mat_v(cp,1,1))*(x[a] - Mat_v(cp,1,1));
    Mat_v(dr,a+1,2) = (y[a] - Mat_v(cp,2,1))*(y[a] - Mat_v(cp,2,1));
    Mat_v(dr,a+1,3) = (z[a] - Mat_v(cp,3,1))*(z[a] - Mat_v(cp,3,1));
  }    
  
  double r = 0;
  for(int a=0; a<4; a++)
  {
    double r_n = sqrt(Mat_v(dr,a+1,1) + Mat_v(dr,a+1,2) +Mat_v(dr,a+1,3));
    r += r_n;
  }      

  r = r/4.0;

  Matrix_cleanup(A);
  Matrix_cleanup(B);
  Matrix_cleanup(AI);
  Matrix_cleanup(cp);
  Matrix_cleanup(dr);
  return r;
}

int find_grid_coord_max_min(double *Xmin, double *Xmax, int nn, NODE *node)
{
  int err = 0;
  int nsd = 3;      
  double x_max[nsd], x_min[nsd], x[nsd];
  
  for(int a=0; a<nsd; a++)
  {
     x_max[a] = -1.0e-15;
     x_min[a] =  1.0e-15;       
  }
    
  for(int a=0; a<nn; a++)
  {
    x[0] = node[a].x1_fd;
    x[1] = node[a].x2_fd;
    x[2] = node[a].x3_fd;
    
    for(int n=0; n<nsd; n++)
    {
      x_max[n] = MAX(x_max[n], x[n]);
      x_min[n] = MIN(x_min[n], x[n]);      
    }
  }
  
  MPI_Allreduce(x_max,Xmax,nsd,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(x_min,Xmin,nsd,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);  
      
  return err;  
}

int compute_min_element_size(double *min_r, int nn, int ne, NODE *node, ELEMENT *elem)
{
  int err = 0;
  int nsd = 3;      
  
  double x[3];

  double L_min_r = 1.0e+15;
  *min_r = 1.0e+15;
  
  for(int e=0; e<ne; e++)
  {
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe, e, elem, node, 1,1);
    double r = sphere_on_four_pt(fe.node_coord.m_pdata);
    L_min_r = MIN(r,L_min_r);        
    FEMLIB_destruct(&fe);    
  }
  
  MPI_Allreduce(&L_min_r,min_r,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);  
          
  return err;
}


int find_grid_point_interest(int *is_pt_in_here, int *nid, double *pt, int nn, NODE *node, int myrank, double tol)
{
  int err = 0;
  int nsd = 3;      
  
  double r_min = tol;
  int node_id = -1;
  for(int a=0; a<nn; a++)
  {
    double x = node[a].x1_fd - pt[0];
    double y = node[a].x2_fd - pt[1];
    double z = node[a].x3_fd - pt[2];
    double r = sqrt(x*x + y*y + z*z);
    if(r<r_min)
    {
      node_id = a;
      r_min = r;
    }
  }
  
  double R_min = 1.0e+15;
  double r_min_loc = r_min;
  
  MPI_Allreduce(&r_min,&R_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

  if(fabs(R_min-r_min_loc)<(tol*0.1))
  {  
    *is_pt_in_here = 1;
    *nid = node_id;
  }
  else
  {
    *is_pt_in_here = 0;
    *nid = -1;    
  }  
  return err;
}

int main(int argc,char *argv[])
{
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int myrank = 0;
  int nproc = 0;

  /*=== END INITIALIZATION === */
  MPI_Init (&argc,&argv);
  MPI_Comm_rank (mpi_comm,&myrank);
  MPI_Comm_size (mpi_comm,&nproc);
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen = 0;  
  MPI_Get_processor_name (processor_name,&namelen);
  PGFEM_initialize_io(NULL,NULL);  
    
  PGFem3D_opt options;
  if (argc <= 2){
    if(myrank == 0){
      print_usage(stdout);
    }
    exit(0);
  }
  set_default_options(&options);
  re_parse_command_line(myrank,2,argc,argv,&options);

  long nn = 0;
  long Gnn = 0;
  long ndofn = 0;
  long ne = 0;
  long ni = 0;
  double err = 0.0;
  double limit = 0.0;  
  long nmat = 0;  
  long nc = 0;
  long np = 0;
  NODE *node = NULL;  
  ELEMENT *elem = NULL; 
  MATERIAL *mater = NULL;   
  MATGEOM matgeom = NULL;
  SUPP sup = NULL;
  long nln = 0;  
  ZATNODE *znod = NULL;   
  long nle_s = 0;
  ZATELEM *zele_s = NULL;
  long nle_v = 0;
  ZATELEM *zele_v = NULL;
  Model_parameters *param_list = NULL;   
  
  int in_err = 0;
  in_err = read_input_file(&options,mpi_comm,&nn,&Gnn,&ndofn,
         &ne,&ni,&err,&limit,&nmat,&nc,&np,&node,
         &elem,&mater,&matgeom,&sup,&nln,&znod,
         &nle_s,&zele_s,&nle_v,&zele_v);
  if(in_err){
    PGFEM_printerr("[%d]ERROR: incorrectly formatted input file!\n",
      myrank);
    PGFEM_Abort();
  }   
      
  HOMMAT *hommat = NULL;
  Mat_3D_orthotropic (nmat,mater,options.analysis_type);
  
  long ***a = NULL;
  a = aloc3l (nmat,nmat,nc);
  long nhommat = list(a,ne,nmat,nc,elem); 
  
  /*  alocation of the material matrices  */
  hommat = build_hommat(nhommat);

  hom_matrices(a,ne,nmat,nc,elem,mater,matgeom,
    hommat,matgeom->SH,options.analysis_type);

  dealoc3l(a,nmat,nmat);
  free(mater);  
  
  EPS *eps = NULL;    
  eps = build_eps_il(ne,elem,options.analysis_type);
    
  if (options.analysis_type == CM) {
    /* parameter list and initialize const. model at int points.
     * NOTE: should catch/handle returned error flag...
     */
    char *cm_filename = NULL;
    alloc_sprintf(&cm_filename,"%s/model_params.in",options.ipath);
    FILE *cm_in = PGFEM_fopen(cm_filename, "r");
    read_model_parameters_list(&param_list, nhommat, hommat, cm_in);
    free(cm_filename);
    fclose(cm_in);
    init_all_constitutive_model(eps,ne,elem,nhommat,param_list);
  } 
  
  /////////////////////////////////////////////////////////////////////////////////////
  // read solver file  
  FILE *in_st;
  char filename[1024];

  if(options.override_solver_file)
  {
    if(myrank == 0)
	    printf("Overriding the default solver file with:\n%s\n", options.solver_file);

    in_st = fopen(options.solver_file,"r");
  }
  else 
  {
    sprintf (filename,"%s/%s%d.in.st",options.ipath,options.ifname,myrank);
    in_st = fopen(filename,"r");  
  }
  
  double nor_min;
  long iter_max, npres, FNR, nt, ostepno;

  
  fscanf (in_st,"%lf %ld %ld %ld",&nor_min,&iter_max,&npres,&FNR);
  fscanf (in_st,"%ld",&nt);
    
  Matrix(double) t;
  Matrix_construct_redim(double, t, nt+1, 1);
  
  for(int a=0; a<nt+1; a++)
    fscanf(in_st,"%lf", (t.m_pdata) + a);

  fscanf(in_st, "%ld", &ostepno);
  
  Matrix(long) ostep;
  Matrix_construct_redim(long, ostep, ostepno, 1);

  for(int a=0; a<ostepno; a++)
    fscanf(in_st,"%ld", ostep.m_pdata + a);
  
  fclose(in_st);  

  /////////////////////////////////////////////////////////////////////////////////////
  // read inputs  
  double xmax[3], xmin[3], min_r;
  int is_pt_in_here_up   = 0;
  int is_pt_in_here_down = 0;
    
  int nid_up   = -1;
  int nid_down = -1;  
  
  FILE *fp = fopen("displacement.out", "w");
  FILE *fpVec = fopen("displacement.out.vec", "w");

  for(int A=0; A<ostepno; A++)
  {
    sprintf(filename,"%s/VTK/STEP_%.5ld/%s_%d_%ld.vtu",options.opath,ostep.m_pdata[A],options.ofname,myrank,ostep.m_pdata[A]);   
  
    double *u = aloc1(nn*ndofn);
    double *P, *V;

    int nVol = 1;

    if(options.analysis_type==TF)
    {
      if(elem[0].toe==10 && ndofn==3)
      {  
        npres = 1;
        nVol = 1;
      }
    
      P = (double *) malloc(sizeof(double)*ne*npres);
      V = (double *) malloc(sizeof(double)*ne*nVol);
      read_VTK_file4TF(filename, u, P, V);        
      for (int e=0;e<ne;e++)
      {
        if(npres==1)
        {
      	  eps[e].d_T   = (double *) PGFEM_calloc(3,sizeof(double));
      	  eps[e].d_T[0] = P[e];
        }
    
     	  eps[e].T   = (double *) PGFEM_calloc(nVol*3,sizeof(double));	
  		  eps[e].T[0] = V[e];
      }
                              
      free(P);
      free(V);          
    }
    else
      read_VTK_file(filename, u);
  
    // find interest point id (node id) for the future use:
    // two points are found (up and down) from a pt to compute QOI (displacement)
    // using two point, displacement will be interpolated. 
    if(A==0)
    {  
      find_grid_coord_max_min(xmin,xmax,nn,node);
      compute_min_element_size(&min_r,nn,ne,node,elem);
      double pt[3];
      pt[0] = 0.0;
      pt[1] = xmax[1];
      pt[2] = xmax[2]/2.0-min_r;
      
      find_grid_point_interest(&is_pt_in_here_down,&nid_down,pt,nn,node,myrank,min_r*2.0);
      
      pt[0] = 0.0;
      pt[1] = xmax[1];
      pt[2] = xmax[2]/2.0+min_r;
      
      find_grid_point_interest(&is_pt_in_here_up,&nid_up,pt,nn,node,myrank,min_r*2.0);
            
      if(myrank==0)
      {  
        printf("xmax = %e %e %e\nxmin = %e %e %e\n", xmax[0],xmax[1],xmax[2],xmin[0],xmin[1],xmin[2]);
        printf("smallest element size = %e\n", min_r); 
      }
    }

    // Each process will have or not have point where QOI is computed. 
    // Communicates to share point information (up and down points)
    double L_disp_down[5] = {0};
    double G_disp_down[5] = {0};

    if(is_pt_in_here_down)
    {
      L_disp_down[0] = node[nid_down].x3_fd;
      L_disp_down[1] = u[nid_down*3+0];
      L_disp_down[2] = u[nid_down*3+1];
      L_disp_down[3] = u[nid_down*3+2];
      L_disp_down[4] = 1.0;
    }                   
    MPI_Allreduce(L_disp_down,G_disp_down,5,MPI_DOUBLE,MPI_SUM,mpi_comm);

    double L_disp_up[5] = {0};
    double G_disp_up[5] = {0};

    if(is_pt_in_here_up)
    {
      L_disp_up[0] = node[nid_up].x3_fd;
      L_disp_up[1] = u[nid_up*3+0];
      L_disp_up[2] = u[nid_up*3+1];
      L_disp_up[3] = u[nid_up*3+2];
      L_disp_up[4] = 1.0;      
    }                   
    MPI_Allreduce(L_disp_up,G_disp_up,5,MPI_DOUBLE,MPI_SUM,mpi_comm);
    
    // Multiple process possibily share one point by domain decomposition
    // In this case, values are averaged and recover to original values
    for(int a=0; a<4; a++)
    {
      G_disp_down[a] /= G_disp_down[4];
        G_disp_up[a] /=   G_disp_up[4];      
    }
    
    // linear interpolation of displacement
    double z = xmax[2]/2.0;
    double dz = G_disp_up[0] - G_disp_down[0];
    double disp[3];
    for(int a=1; a<=3; a++)
    {
      disp[a-1] = (G_disp_up[a] - G_disp_down[a])/dz*(z - G_disp_down[0]) + G_disp_down[a];
    }
        
    if(myrank==0){
      fprintf(fp, "%e %e %e %e\n",t.m_pdata[ostep.m_pdata[A]], disp[0],disp[1],disp[2]); 
      fprintf(fpVec, "%e\n", disp[0]);
      fprintf(fpVec, "%e\n", disp[1]);
      fprintf(fpVec, "%e\n", disp[2]);   
    }
    free(u); 
  }
  
  fclose(fp);  
  fclose(fpVec);
  Matrix_cleanup(t);
  Matrix_cleanup(ostep);
 
      
  destroy_zatnode(znod,nln);
  destroy_zatelem(zele_s,nle_s);
  destroy_zatelem(zele_v,nle_v);
  destroy_matgeom(matgeom,np);
  destroy_hommat(hommat,nhommat);
  destroy_model_parameters_list(nhommat,param_list);  
  destroy_eps_il(eps,elem,ne,options.analysis_type);
  destroy_supp(sup);
  destroy_elem(elem,ne);
  destroy_node(nn,node);

  /*=== FINALIZE AND EXIT ===*/
  PGFEM_finalize_io();
  MPI_Finalize(); 
  return(0);
}
