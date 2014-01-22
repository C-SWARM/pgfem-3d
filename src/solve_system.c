/*HEADER*/

#include "solve_system.h"
#include "PGFEM_mpi.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef PGFEM_HYPRE_PRECOND_H
#include "PGFEM_HYPRE_preconditioners.h"
#endif

#ifndef DEBUG_SOLVE
#define DEBUG_SOLVE 0
#endif

/* create aliases for HYPRE functions so we can set and run with fewer
   switches/code */

typedef int (*My_HYPRE_ptr_solve_func)(HYPRE_Solver solver,
				       HYPRE_ParCSRMatrix A,
				       HYPRE_ParVector b,
				       HYPRE_ParVector x);

typedef int (*My_HYPRE_ptr_set_pre_func)(HYPRE_Solver solver,
					 My_HYPRE_ptr_solve_func pre_solve,
					 My_HYPRE_ptr_solve_func pre_setup,
					 HYPRE_Solver precond);

typedef int (*My_HYPRE_ptr_check_pre_func)(HYPRE_Solver solver,
					   HYPRE_Solver *precond);

typedef int (*My_HYPRE_ptr_iter_func)(HYPRE_Solver solver,
				      int *num_it);

typedef int (*My_HYPRE_ptr_norm_func)(HYPRE_Solver solver,
				      double *norm);

/* This function handles the setup and solve for the various solver
   libraries (currently only HYPRE and BlockSolve95). The important
   solver information is passed out of the function through the
   SOLVER_INFO structure defined in the header file */

double solve_system(const PGFem3D_opt *opts,
		    double *loc_rhs,
		    double *loc_sol,
		    const int tim,
		    const int iter,
		    const long *DomDof,
		    SOLVER_INFO *info,
		    PGFEM_HYPRE_solve_info *PGFEM_hypre,
		    BSprocinfo *BSinfo,
		    BSspmat *k,
		    BSpar_mat **pk,
		    BSpar_mat **f_pk,
		    MPI_Comm mpi_comm)
{
  double func_time = -MPI_Wtime();
  int myrank=0;
  MPI_Comm_rank(mpi_comm,&myrank);

  switch(opts->solverpackage){
  case BLOCKSOLVE:
    {
      /* Set symmetry */
      BSset_mat_symmetric (k,FALSE);
	
      /* and storage */
      BSset_mat_icc_storage (k,FALSE);
	
      /* permute the matrix */
      static int BS_perm = 0;
      if (tim == 0 && BS_perm == 0) {
	*pk = BSmain_perm (BSinfo,k);
	CHKERRN(0);
	BS_perm = 1;
      }else{
	BSmain_reperm (BSinfo,k,*pk);
	CHKERRN(0);
      }
	
      /* diagonally scale the matrix */
      if(BSinfo->scaling) {
	BSscale_diag (*pk,(*pk)->diag,BSinfo);
	CHKERRN(0);
      }
	
      /* set up the communication structure for triangular matrix
	 solution */
      BScomm *kcomm = BSsetup_forward (*pk,BSinfo);
      CHKERRN(0);
	
      /* get a copy of the sparse matrix */
      *f_pk = BScopy_par_mat (*pk);
      CHKERRN(0);
	
      /* set up a communication structure for factorization */
      BScomm *k_comm = BSsetup_factor (*f_pk,BSinfo);
      CHKERRN(0);
	
      /* shifted_diag is the initial diagonal */
      double BS_diag = 1.0;
	
      /* factor the matrix until successful */
      while (BSfactor (*f_pk,k_comm,BSinfo) != 0) {
	CHKERRN(0);
	/* recopy the nonzeroes */
	BScopy_nz (*pk,*f_pk);
	CHKERRN(0);
	/* increment the diagonal shift */
	BS_diag += 0.1;
	BSset_diag (*f_pk,BS_diag,BSinfo);
	CHKERRN(0);
      }
      CHKERRN(0);
	
      /* Solve system of equations */
      info->n_iter = BSpar_solve (*pk,*f_pk,kcomm,loc_rhs,loc_sol,
				  &(info->res_norm),BSinfo);
      CHKERRN(0);

      /* tell it to print coloring, reordering and linear system
	 solution options */
      BSctx_set_pr (BSinfo,FALSE);
      CHKERRN(0);
    }
    break;
  case HYPRE:
    {
      /* Assemble the rhs and solution vector */
      nulld(loc_sol,DomDof[myrank]);
      HYPRE_IJVectorSetValues(PGFEM_hypre->hypre_rhs,DomDof[myrank],
			      &PGFEM_hypre->grows[0],loc_rhs);
      HYPRE_IJVectorSetValues(PGFEM_hypre->hypre_sol,DomDof[myrank],
			      &PGFEM_hypre->grows[0],loc_sol);
      HYPRE_IJVectorAssemble(PGFEM_hypre->hypre_rhs);
      HYPRE_IJVectorAssemble(PGFEM_hypre->hypre_sol);

      /* set the solver and preconditioner function pointers*/
      My_HYPRE_ptr_solve_func precond_solve = NULL;
      My_HYPRE_ptr_solve_func precond_setup = NULL;
      My_HYPRE_ptr_solve_func solver_solve = NULL;
      My_HYPRE_ptr_solve_func solver_setup = NULL;

      My_HYPRE_ptr_set_pre_func set_precond = NULL;
      My_HYPRE_ptr_check_pre_func get_precond = NULL;
      My_HYPRE_ptr_iter_func get_num_iter = NULL;
      My_HYPRE_ptr_norm_func get_res_norm = NULL;

      switch(PGFEM_hypre->precond_type){
      case PARA_SAILS:
	if(iter <= 1){
	  HYPRE_ParaSailsSetReuse (PGFEM_hypre->hypre_pc,0);
	} else {
	  HYPRE_ParaSailsSetReuse (PGFEM_hypre->hypre_pc,1);
	}
	precond_solve = HYPRE_ParCSRParaSailsSolve;
	precond_setup = HYPRE_ParCSRParaSailsSetup;
	break;

      case PILUT:
	precond_solve = HYPRE_ParCSRPilutSolve;
	precond_setup = HYPRE_ParCSRPilutSetup;
	break;

      case EUCLID:
	precond_solve = HYPRE_EuclidSolve;
	precond_setup = HYPRE_EuclidSetup;
	break;

      case BOOMER:
	precond_solve = HYPRE_BoomerAMGSolve;
	precond_setup = HYPRE_BoomerAMGSetup;
	break;

      case DIAG_SCALE:
	precond_solve = PGFEM_HYPRE_ScaleDiagSolve;
	precond_setup = PGFEM_HYPRE_ScaleDiagSetup;
	break;

      case JACOBI:
	precond_solve = PGFEM_HYPRE_JacobiSolve;
	precond_setup = PGFEM_HYPRE_JacobiSetup;
	break;

      case NONE:
	break;

      default:
	info->err = UNREC_PRECOND;
	return 0.0;
      }
	
      switch(PGFEM_hypre->solver_type){
      case HYPRE_GMRES:
	set_precond = HYPRE_ParCSRGMRESSetPrecond;
	get_precond = HYPRE_ParCSRGMRESGetPrecond;
	solver_setup = HYPRE_ParCSRGMRESSetup;
	solver_solve = HYPRE_ParCSRGMRESSolve;
	get_num_iter = HYPRE_ParCSRGMRESGetNumIterations;
	get_res_norm = HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm;
	break;

      case HYPRE_BCG_STAB:
	set_precond = HYPRE_ParCSRBiCGSTABSetPrecond;
	get_precond = HYPRE_ParCSRBiCGSTABGetPrecond;
	solver_setup = HYPRE_ParCSRBiCGSTABSetup;
	solver_solve = HYPRE_ParCSRBiCGSTABSolve;
	get_num_iter = HYPRE_ParCSRBiCGSTABGetNumIterations;
	get_res_norm = HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm;
	break;

      case HYPRE_AMG:
	/* no preconditioner */
	solver_setup = HYPRE_BoomerAMGSetup;
	solver_solve = HYPRE_BoomerAMGSolve;
	get_num_iter = HYPRE_BoomerAMGGetNumIterations;
	get_res_norm = HYPRE_BoomerAMGGetFinalRelativeResidualNorm;
	get_num_iter = HYPRE_ParCSRFlexGMRESGetNumIterations;
	get_res_norm = HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm;
	break;

      case HYPRE_FLEX:
	set_precond = HYPRE_ParCSRFlexGMRESSetPrecond;
	get_precond = HYPRE_ParCSRFlexGMRESGetPrecond;
	solver_setup = HYPRE_ParCSRFlexGMRESSetup;
	solver_solve = HYPRE_ParCSRFlexGMRESSolve;
	get_num_iter = HYPRE_ParCSRFlexGMRESGetNumIterations;
	get_res_norm = HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm;
	break;

      case HYPRE_HYBRID:
	set_precond = HYPRE_ParCSRHybridSetPrecond;
	solver_setup = HYPRE_ParCSRHybridSetup;
	solver_solve = HYPRE_ParCSRHybridSolve;
	get_num_iter = HYPRE_ParCSRHybridGetNumIterations;
	get_res_norm = HYPRE_ParCSRHybridGetFinalRelativeResidualNorm;
	break;

      default:
	info->err = UNREC_SOLVER_TYPE;
	return 0.0;
      }

      /* do actual setup and solve */
      if(PGFEM_hypre->precond_type != NONE 
	 && PGFEM_hypre->solver_type != HYPRE_AMG){
	set_precond(PGFEM_hypre->hypre_solver,precond_solve,precond_setup,
		    PGFEM_hypre->hypre_pc);
	if(opts->solver != HYPRE_HYBRID) {/* does not have get function */
	  get_precond(PGFEM_hypre->hypre_solver,
		      &PGFEM_hypre->hypre_pc_gotten);
	  if(PGFEM_hypre->hypre_pc_gotten != PGFEM_hypre->hypre_pc) {
	    info->err = BAD_PRECOND;
	    return 0.0;
	  }
	}
      }

      solver_setup(PGFEM_hypre->hypre_solver,PGFEM_hypre->hypre_pk,
		   PGFEM_hypre->hypre_prhs,PGFEM_hypre->hypre_psol);
      info->err = solver_solve(PGFEM_hypre->hypre_solver,
			       PGFEM_hypre->hypre_pk,
			       PGFEM_hypre->hypre_prhs,
			       PGFEM_hypre->hypre_psol);
      get_num_iter(PGFEM_hypre->hypre_solver,&(info->n_iter));
      get_res_norm(PGFEM_hypre->hypre_solver,&(info->res_norm));
      HYPRE_IJVectorGetValues(PGFEM_hypre->hypre_sol,DomDof[myrank],
			      &PGFEM_hypre->grows[0],loc_sol);

      /* reset hypre error flag so we can restart cleanly based on OUR
	 error tolerances etc. */
      hypre__global_error = 0;

    }
    break;
  default:
    info->err = UNREC_SOLVER;
    return 0.0;
  }
  /* update timer and return */
  func_time += MPI_Wtime();
  return func_time;
}

int solve_system_check_error(FILE *out,
			     const SOLVER_INFO info)
{
  switch(info.err){
  case SOLVE_SUCCESS:
    return 0;
  case SOLVE_ERROR:
    PGFEM_fprintf(out,"solve system returned error code SOLVE_ERROR.\n");
    return 1;
  case UNREC_SOLVER:
    PGFEM_fprintf(out,"solve system returned error code UNREC_SOLVER.\n");
    return 1;
  case UNREC_SOLVER_TYPE:
    PGFEM_fprintf(out,"solve system returned error code UNREC_SOLVER_TYPE.\n");
    return 1;
  case UNREC_PRECOND:
    PGFEM_fprintf(out,"solve system returned error code UNREC_PRECOND.\n");
    return 1;
  case BAD_PRECOND:
    PGFEM_fprintf(out,"solve system returned error code BAD_PRECOND.\n");
    return 1;
  default:
    PGFEM_fprintf(out,"solve system returned unrecognized error code (%d).\n",
	    info.err);
    return 1;
  }
}
