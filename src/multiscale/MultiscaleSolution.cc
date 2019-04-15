#include "msnet/Multiscale.hpp"
#include "pgfem3d/MultiscaleSolution.hpp"
#include "pgfem3d/MultiscaleCommon.hpp"
#include "allocation.h"
#include "constitutive_model.h"
#include "elem3d.h"
#include "incl.h"
#include "initialize_damage.h"
#include "set_fini_def.h"
#include "utils.h"
#include <search.h>
#include <cassert>

using namespace pgfem3d;

namespace pgfem3d {
  
  /**
   * \brief Private data type for storing common solution buffers.
   *
   * This portion is allocated and held by COMMON_MULTISCALE and
   * MULTISCALE_SOLUTION simply holds pointers to the various buffers.
   */
  typedef struct MULTISCALE_SOLUTION_BUFFERS {
    /* local vectors */
    double *f;
    double *d_r;
    double *rr;
    double *D_R;
    double *R;
    double *RR;
    double *f_u;
    double *f_defl;
    double *RRn;
    double *U;
    double *DK;
    double *dR;
    
    /* global vectors */
    double *BS_f;
    double *BS_x;
    double *BS_RR;
    double *BS_f_u;
    double *BS_d_r;
    double *BS_D_R;
    double *BS_rr;
    double *BS_R;
    double *BS_U;
    double *BS_DK;
    double *BS_dR;
  } MULTISCALE_SOLUTION_BUFFERS;
  
  /*==== API FUNCTIONS ====*/
  static inline int
  sort_first(const void *a, const void *b)
  {
    return *((const int*)a) - *((const int*) b);
  }
  
  static inline int
  sort_second(const void *a, const void *b)
  {
    return *((const int*)a+1) - *((const int*)b+1);
  }
  
  void
  sol_idx_map_build(sol_idx_map *map, const size_t size)
  {
    /* map stores id : idx pairs */
    map->size = size;
    map->map = new int[2*size];
    for(size_t i=0; i<2*size; i += 2){
      map->map[i] = -1;
      map->map[i+1] = i/2;
    }
  }
  
  void sol_idx_map_destroy(sol_idx_map *map)
  {
    delete [] map->map;
    map->map = NULL;
    map->size = 0;
  }
  
  void sol_idx_map_sort_id(sol_idx_map *map)
  {
    qsort(map->map,map->size,sizeof(*(map->map)),sort_first);
  }
  
  void sol_idx_map_sort_idx(sol_idx_map *map)
  {
    qsort(map->map,map->size,sizeof(*(map->map)),sort_second);
  }
  
  int sol_idx_map_id_get_idx(const sol_idx_map *map,
			     const int id)
  {
    int val[2] = {0,0}; val[0] = id;
    size_t len = map->size;
    int *ptr = ((int *) lfind(val,map->map,&len,2*sizeof(*(map->map)),sort_first) + 1);
    return (ptr == NULL)? -1 : *ptr;
  }
  
  int sol_idx_map_idx_get_id(const sol_idx_map *map,
			     const int idx)
  {
    int val[2] = {0,0}; val[1] = idx;
    size_t len = map->size;
    int *ptr = ((int *) lfind(val,map->map,&len,2*sizeof(*(map->map)),sort_second));
    return (ptr == NULL)? -1 : *ptr;
  }
  
  int sol_idx_map_get_idx_reset_id(sol_idx_map *map,
				   const int cur_id,
				   const int new_id)
  {
    int val[2] = {0,0}; val[0] = cur_id;
    size_t len = map->size;
    
    /* get pointer to matching pair */
    int *ptr = static_cast<int*>(lfind(val, map->map, &len, sizeof(val),
				       sort_first));
    assert(ptr != NULL);
    ptr[0] = new_id;
    return ptr[1];
  }
  
  void sol_idx_map_idx_set_id(sol_idx_map *map,
			      const int idx,
			      const int id)
  {
    int val[2] = {0,0}; val[1] = idx;
    size_t len = map->size;
    int *ptr = static_cast<int*>(lfind(val, map->map, &len, sizeof(val),
				       sort_second));
    assert(ptr != NULL);
    *ptr = id;
  }
  
  void initialize_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol)
  {
    /* stress/strain/state information */
    sol->sig_e = NULL;
    sol->sig_n = NULL;
    sol->eps = NULL;
    sol->statv_list = NULL;
    sol->crpl = NULL;
    sol->npres = 0;
    
    /* solution information */
    /* local vectors */
    sol->r = NULL;
    sol->f = NULL;
    sol->d_r = NULL;
    sol->rr = NULL;
    sol->D_R = NULL;
    sol->R = NULL;
    sol->f_defl = NULL;
    sol->RR = NULL;
    sol->f_u = NULL;
    sol->RRn = NULL;
    sol->U = NULL;
    sol->DK = NULL;
    sol->dR = NULL;
    
    /* global vectors */
    sol->BS_f = NULL;
    sol->BS_f_u = NULL;
    sol->BS_x = NULL;
    sol->BS_RR = NULL;
    sol->BS_d_r = NULL;
    sol->BS_D_R = NULL;
    sol->BS_rr = NULL;
    sol->BS_R = NULL;
    sol->BS_U = NULL;
    sol->BS_DK = NULL;
    sol->BS_dR = NULL;
    
    /* convergence info */
    sol->dt = 0.0;
    sol->times = NULL;
    sol->tim = 0;
    sol->p_tim = 0;
    sol->NORM = 0.0;
    
    /* failure flag */
    sol->failed = 0;
  }
  
  void build_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
				 const MultiscaleCommon *common,
				 const int analysis)
  {
    int myrank = common->rank;
    const long local_len = common->ndofd;

    sol->sig_e = build_sig_il(common->ne,analysis,common->elem);
    
    if (analysis == CM)
      {
	int n_state_variables = 0;
	for(int eid=0; eid<common->ne; eid++)
	  {
	    int nne = common->elem[eid].toe;
	    long nint = 0;
	    int_point(nne,&nint);
	    n_state_variables += nint;
	  }
	if (n_state_variables > 0)
	  sol->statv_list = (State_variables *) malloc(sizeof(State_variables)*n_state_variables);
      }
    
    sol->eps = build_eps_il(common->ne,common->elem,analysis,&(sol->statv_list));
    initialize_damage(common->ne,common->elem,common->hommat,
		      sol->eps,analysis);
    
    if (analysis == CM) {
      init_all_constitutive_model(sol->eps, common->ne,
				  common->elem, common->nhommat,
				  common->hommat, myrank);
    }
    
    /* initialize state variable buffer at macro time (n) */
    /* length of the solution vector */
    sol->packed_state_var_len = local_len*sizeof(double);
    
    /* length of the EPS list */
    sol->packed_state_var_len += sizeof_eps_list(sol->eps,
						 common->ne,
						 common->elem,
						 analysis);
    
    /* length of the state variables stored in COEL */
    sol->packed_state_var_len += coel_list_get_state_length_bytes(common->nce,
								  common->coel);
    
    /* length of NORM */
    sol->packed_state_var_len += sizeof(sol->NORM);
    
    /* length of dt */
    sol->packed_state_var_len += sizeof(sol->dt);
    
    /* length of failed */
    sol->packed_state_var_len += sizeof(sol->failed);
    
    /* allocate the packed state buffer */
    sol->packed_state_var_n = PGFEM_calloc(char, sol->packed_state_var_len);
    
    /* initialize buffer */
    {
      /* pack eps after end of sol vector buffer */
      size_t pos = local_len*sizeof(double);
      pack_eps_list(sol->eps,common->ne,common->elem,analysis,
		    sol->packed_state_var_n,&pos);
    }
    /* need to figure out elem/coel_state_info indexing */
    
    switch(analysis){
    case STABILIZED: case MINI: case MINI_3F: sol->npres = 4; break;
    case DISP: case CM: sol->npres = 0; break;
    default: sol->npres = 1; break;
    }
    build_pressure_nodes(common->ne,sol->npres,common->elem,
			 sol->sig_e,sol->eps,analysis);
    set_fini_def (common->ne,sol->npres,common->elem,
		  sol->eps,sol->sig_e,analysis);
    
    /*=== crystal plasticity is not currently supported ===*/
    /* if (analysis == FS_CRPL) { */
    /*   sol->crpl = PGFEM_calloc (common->nhommat,sizeof(CRPL)); */
    /*   read_cryst_plast (in1,common->nhommat,sol->crpl,plc); */
    /*   build_crystal_plast (common->ne,common->elem,sol->sig_e,sol->eps, */
    /*             sol->crpl,analysis,plc); */
    /*   set_fini_def_pl(common->ne,sol->npres,common->elem, */
    /*            sol->eps,sol->sig_e,sol->crpl,analysis,plc); */
    /* } */
    
    /* local solution vectors */
    sol->r = PGFEM_calloc(double, local_len);
    
    /* Get pointers to the shared solution workspace */
    {
      MULTISCALE_SOLUTION_BUFFERS *buff =
	(MULTISCALE_SOLUTION_BUFFERS *) common->solution_buffer;
      
      sol->f  = buff->f     ;
      sol->d_r    = buff->d_r   ;
      sol->rr = buff->rr    ;
      sol->D_R    = buff->D_R   ;
      sol->R  = buff->R     ;
      sol->f_defl = buff->f_defl;
      sol->RR = buff->RR    ;
      sol->f_u    = buff->f_u   ;
      sol->RRn    = buff->RRn   ;
      sol->U  = buff->U     ;
      sol->DK = buff->DK    ;
      sol->dR = buff->dR    ;
      
      sol->BS_f   = buff->BS_f  ;
      sol->BS_f_u = buff->BS_f_u;
      sol->BS_x   = buff->BS_x  ;
      sol->BS_RR  = buff->BS_RR ;
      sol->BS_d_r = buff->BS_d_r;
      sol->BS_D_R = buff->BS_D_R;
      sol->BS_rr  = buff->BS_rr ;
      sol->BS_R   = buff->BS_R  ;
      sol->BS_U   = buff->BS_U  ;
      sol->BS_DK  = buff->BS_DK ;
      sol->BS_dR  = buff->BS_dR ;
    }
    
    sol->dt = 0.0;
    sol->times = PGFEM_calloc(double, 3);
    sol->tim = 0;
    sol->NORM = 0.0;
    
    sol->failed = 0;
  }

  void destroy_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
				   const MultiscaleCommon *common,
				   const int analysis)
  {
    free(sol->r);
    free(sol->packed_state_var_n);
    free(sol->times);
    
    destroy_sig_il(sol->sig_e,common->elem,common->ne,analysis);
    /* destroy sig_n */
    
    destroy_eps_il(sol->eps,common->elem,common->ne,analysis);
    if(sol->statv_list != NULL)
      free(sol->statv_list);
    
    //destroy_crpl
  }
  
  void build_MULTISCALE_SOLUTION_BUFFERS(void **buffer,
					 const int local_len,
					 const int global_len)
  {
    // allocate space for the solution
    *buffer = PGFEM_calloc(MULTISCALE_SOLUTION_BUFFERS, 1);

    MULTISCALE_SOLUTION_BUFFERS *buff = (MULTISCALE_SOLUTION_BUFFERS*) *buffer;
    
    /* local solution vectors */
    buff->f = PGFEM_calloc(double, local_len);
    buff->d_r = PGFEM_calloc(double, local_len);
    buff->rr = PGFEM_calloc(double, local_len);
    buff->D_R = PGFEM_calloc(double, local_len);
    buff->R = PGFEM_calloc(double, local_len);
    buff->f_defl = PGFEM_calloc(double, local_len);
    buff->RR = PGFEM_calloc(double, local_len);
    buff->f_u = PGFEM_calloc(double, local_len);
    buff->RRn = PGFEM_calloc(double, local_len);
    buff->U = PGFEM_calloc(double, local_len);
    buff->DK = PGFEM_calloc(double, local_len);
    buff->dR = PGFEM_calloc(double, local_len);
    
    /* global solution vectors */
    buff->BS_f = PGFEM_calloc(double, global_len);
    buff->BS_f_u = PGFEM_calloc(double, global_len);
    buff->BS_x = PGFEM_calloc(double, global_len);
    buff->BS_RR = PGFEM_calloc(double, global_len);
    buff->BS_d_r = PGFEM_calloc(double, global_len);
    buff->BS_D_R = PGFEM_calloc(double, global_len);
    buff->BS_rr = PGFEM_calloc(double, global_len);
    buff->BS_R = PGFEM_calloc(double, global_len);
    buff->BS_U = PGFEM_calloc(double, global_len);
    buff->BS_DK = PGFEM_calloc(double, global_len);
    buff->BS_dR = PGFEM_calloc(double, global_len);
  }
  
  void destroy_MULTISCALE_SOLUTION_BUFFERS(void *buffer)
  {
    MULTISCALE_SOLUTION_BUFFERS *buff = (MULTISCALE_SOLUTION_BUFFERS*) buffer;
    free(buff->f);
    free(buff->d_r);
    free(buff->rr);
    free(buff->D_R);
    free(buff->R);
    free(buff->f_defl);
    free(buff->RR);
    free(buff->f_u);
    free(buff->RRn);
    free(buff->U);
    free(buff->DK);
    free(buff->dR);
    
    free(buff->BS_f);
    free(buff->BS_f_u);
    free(buff->BS_x);
    free(buff->BS_RR);
    free(buff->BS_d_r);
    free(buff->BS_D_R);
    free(buff->BS_rr);
    free(buff->BS_R);
    free(buff->BS_U);
    free(buff->BS_DK);
    free(buff->BS_dR);
  }

  int reset_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
				const MultiscaleCommon *micro)
  {
    int err = 0;
    int myrank = micro->rank;
    const int loc_ndof = micro->ndofd;
    const int g_ndof = micro->DomDof[myrank];
    size_t pos = 0;
    
    /* reset displacement (solution) vector */
    unpack_data(sol->packed_state_var_n,
		sol->r,&pos,loc_ndof,sizeof(*(sol->r)));
    
    /* null all of the other local vectors */
    nulld(sol->f,loc_ndof);
    nulld(sol->d_r,loc_ndof);
    nulld(sol->rr,loc_ndof);
    nulld(sol->D_R,loc_ndof);
    nulld(sol->R,loc_ndof);
    nulld(sol->RR,loc_ndof);
    nulld(sol->f_u,loc_ndof);
    nulld(sol->f_defl,loc_ndof);
    nulld(sol->RRn,loc_ndof);
    nulld(sol->U,loc_ndof);
    nulld(sol->DK,loc_ndof);
    nulld(sol->dR,loc_ndof);
    
    /* null all of the "global" vectors */
    nulld(sol->BS_f,g_ndof);
    nulld(sol->BS_x,g_ndof);
    nulld(sol->BS_RR,g_ndof);
    nulld(sol->BS_f_u,g_ndof);
    nulld(sol->BS_d_r,g_ndof);
    nulld(sol->BS_D_R,g_ndof);
    nulld(sol->BS_rr,g_ndof);
    nulld(sol->BS_R,g_ndof);
    nulld(sol->BS_U,g_ndof);
    nulld(sol->BS_DK,g_ndof);
    nulld(sol->BS_dR,g_ndof);
    
    /* reset state variables */
    unpack_eps_list(sol->eps,
		    micro->ne,
		    micro->elem,
		    micro->opts->analysis_type,
		    sol->packed_state_var_n,
		    &pos);
    
    /* reset coel state variables */
    coel_list_unpack_state(micro->nce,
			   micro->coel,
			   micro->co_props,
			   sol->packed_state_var_n,
			   &pos);
    
    /* reset NORM */
    unpack_data(sol->packed_state_var_n,&sol->NORM,
		&pos,1,sizeof(sol->NORM));
    
    /* reset dt */
    unpack_data(sol->packed_state_var_n,&sol->dt,
		&pos,1,sizeof(sol->dt));
    
    /* reset failed flag */
    unpack_data(sol->packed_state_var_n,&sol->failed,
		&pos,1,sizeof(sol->failed));
    
    assert(pos == sol->packed_state_var_len);
    if(pos != sol->packed_state_var_len) err++;
    return err;
  }/* reset_MULTISCALE_SOLUTION */
  
  int update_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
				 const MultiscaleCommon *micro)
  {
    int err = 0;
    const int loc_ndof = micro->ndofd;
    size_t pos = 0;
    
    /* copy r -> rn  */
    pack_data(sol->r,sol->packed_state_var_n,&pos,
	      loc_ndof,sizeof(*(sol->r)));
    
    /* leave other solution vectors alone */
    
    /* update state variables */
    pack_eps_list(sol->eps,
		  micro->ne,
		  micro->elem,
		  micro->opts->analysis_type,
		  sol->packed_state_var_n,
		  &pos);
    
    /* update cohesive state variables */
    coel_list_pack_state(micro->nce,
			 micro->coel,
			 sol->packed_state_var_n,
			 &pos);
    
    /* pack NORM */
    pack_data(&sol->NORM,sol->packed_state_var_n,
	      &pos,1,sizeof(sol->NORM));
    
    /* pack dt */
    pack_data(&sol->dt,sol->packed_state_var_n,
	      &pos,1,sizeof(sol->dt));
    
    /* pack failed flag */
    pack_data(&sol->failed,sol->packed_state_var_n,
	      &pos,1,sizeof(sol->failed));
    
    assert(pos == sol->packed_state_var_len);
    if(pos != sol->packed_state_var_len) err++;
    
    return err;
    
  }/* update_MULTISCALE_SOLUTION */
  
  int dump_MULTISCALE_SOLUTION_state(const MULTISCALE_SOLUTION *sol,
				     FILE *out)
  {
    int err = 0;
    size_t n_write = fwrite(sol->packed_state_var_n,sizeof(char),
			    sol->packed_state_var_len,out);
    assert(n_write == sol->packed_state_var_len);
    if(n_write != sol->packed_state_var_len) err++;
    return err;
  }
  
  int read_MULTISCALE_SOLUTION_state(MULTISCALE_SOLUTION *sol,
				     FILE *in)
  {
    int err = 0;
    size_t n_read = fread(sol->packed_state_var_n,sizeof(char),
			  sol->packed_state_var_len,in);
    assert(n_read == sol->packed_state_var_len);
    if(n_read != sol->packed_state_var_len) err++;
    return err;
  }
  
} // end namespace pgfem3d
