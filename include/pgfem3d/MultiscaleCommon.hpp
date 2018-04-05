#ifndef PGFEM3D_MULTISCALE_COMMON_H
#define PGFEM3D_MULTISCALE_COMMON_H

/// @brief This file defines the multiscale common classes
#include "pgfem3d/MultiscaleCommunication.hpp"
#include "pgfem3d/MultiscaleSolution.hpp"
#include "PGFem3D_options.h"
#include "cohesive_element.h"
#include "element.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "supp.h"

namespace pgfem3d {

/** This is the structure of microscale information that is
    identical for all microstructures. */
class MultiscaleCommon : public CommunicationStructure {
public:
  MultiscaleCommon();
  ~MultiscaleCommon();
  
  void initialize(const int argc,
		  char **argv,
		  const CommunicationStructure *com,
		  const int mp_id);

  /** build n solutions to compute on the scale */
  void build_solutions(const int n_solutions);
  
  /* options /solver information */
  pgfem3d::solvers::SparseSystem *SOLVER;
  double lin_err;
  double lim_zero;
  int maxit_nl;

  long ndofd;     /* dof on dom */
  
  /* mesh info */
  long nn;        /**< no. nodes */
  long ne;        /**< no. elements */
  long nce;       /**< no. cohesive elements */
  long ndofn;     /**< no. dof on each node */
  long npres;     /**< no. pressure nodes on elem */
  double VVolume; /**< original volume */
  Node *node;
  Element *elem;  /* NOTE: state/solution information is copied from
		     solution structure */
  COEL *coel;     /* NOTE: state/solution information is copied from
		     solution structure */
  long n_orient;
  MATGEOM matgeom; /**< !pointer */
  long nhommat;
  HOMMAT *hommat;
  Ensight *ensight;
  SUPP supports;   /**< !pointer */
  int n_co_props;
  cohesive_props *co_props;

  /* mixed tangent matrices. */
  void *K_01;
  void *K_10;
  
  PGFem3D_opt *opts;
  MULTISCALE_SOLUTION *sol;

  /* buffer for solution stuff */
  void *solution_buffer;
  sol_idx_map idx_map;
  
private:
  void build_common(const CommunicationStructure *com,
		    const int mp_id);
  void destroy_common();
};

class Microscale : public MultiscaleCommon {};
class Macroscale : public MultiscaleCommon {};

// External solution methods
void initialize_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol);
void build_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
			       const pgfem3d::MultiscaleCommon *common,
			       const int analysis);
void destroy_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
				 const pgfem3d::MultiscaleCommon *common,
				 const int analysis);
  
/** resets a single MULTISCALE_SOLUTION to time (n) */
int reset_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
			      const pgfem3d::MultiscaleCommon *msc);

/** updates a single MULTISCALE_SOLUTION to time (n+1) */
int update_MULTISCALE_SOLUTION(MULTISCALE_SOLUTION *sol,
			       const pgfem3d::MultiscaleCommon *msc);

} // end namespace pgfem3d

#endif
