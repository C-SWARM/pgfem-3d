#ifndef PGFEM3D_MULTISCALE_SOLUTION_H
#define PGFEM3D_MULTISCALE_SOLUTION_H

#include "eps.h"
#include "crpl.h"
#include "sig.h"

namespace pgfem3d {
  
  /** This structure contains the information for the
      history-dependent solution. */
  struct MULTISCALE_SOLUTION{
    /* stress/strain/state information */
    SIG *sig_e;
    SIG *sig_n;
    EPS *eps;
    State_variables *statv_list;    /// object to store element variables
    CRPL *crpl;
    long npres;
    
    /* solution information */
    double *r; /**< current solution r at macro time n+1 */
    
    /* State variables at time n */
    size_t packed_state_var_len;
    char *packed_state_var_n;
    
    /** The following are pointers to a shared buffer elsewhere! The
	buffers are used purely as a workspace. */
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
    
    /* convergence info */
    double dt;
    double *times;
    int tim;
    int p_tim;
    double NORM;
    
    /**
     * Flag denoting cell has failed. Do not compute
     * solutions/response when set to non-zero value
     */
    int failed;
  };
  
  struct SolutionIndexMap {
    size_t size = 0;
    int *map = nullptr;
  };

  using sol_idx_map = SolutionIndexMap;
  
  void sol_idx_map_build(sol_idx_map *map,
			 const size_t size);
  void sol_idx_map_destroy(sol_idx_map *map);
  void sol_idx_map_sort_id(sol_idx_map *map);
  void sol_idx_map_sort_idx(sol_idx_map *map);
  
  /**
   * Get solution idx from job id. Returns -1 if job id not found in
   * the map.
   */
  int sol_idx_map_id_get_idx(const sol_idx_map *map,
			     const int id);
  
  /**
   * Get the idx from the job id and then assign a new id. Aborts if
   * the current id is not valid.
   */
  int sol_idx_map_get_idx_reset_id(sol_idx_map *map,
				   const int cur_id,
				   const int new_id);
  
  /**
   * Get the job id from the idx. Return -1 if not found.
   */
  int sol_idx_map_idx_get_id(const sol_idx_map *map,
			     const int idx);

  /**
   * Set the id at the given idx. Aborts if the idx is not found.
   */
  void sol_idx_map_idx_set_id(sol_idx_map *map,
			      const int idx,
			      const int id);

  void build_MULTISCALE_SOLUTION_BUFFERS(void **buffer,
					 const int local_len,
					 const int global_len);
  
  void destroy_MULTISCALE_SOLUTION_BUFFERS(void *buffer);
  
  /**
   * Dump the solution state vector to a binary file. Returns non-zero
   * if there is a problem writing the file.
   */
  int dump_MULTISCALE_SOLUTION_state(const MULTISCALE_SOLUTION *sol,
				     FILE *out);
  
  /**
   * Read a dumped binary state file. Returns non-zero if there is a
   * problem reading the file.
   */
  int read_MULTISCALE_SOLUTION_state(MULTISCALE_SOLUTION *sol,
				     FILE *in);
  
} // end namespace pgfem3d

#endif
