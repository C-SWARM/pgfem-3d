/* HEADER */
#pragma once
#ifndef PGFEM3D_SUBDIVISION_H
#define PGFEM3D_SUBDIVISION_H

#include "PGFem3D_data_structure.h"
#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "sig.h"
#include "supp.h"

class SubdivisionScheme{
  public:
  int step_size;              /// step size
  int step_id;                /// stepping id with subdivision
  int decellerate;            /// decelleration parameter
  int accellerate;            /// accelleration parameter
  double dt_0;                /// saved time step size to accellerate or decellerate
  int reset_variables;        /// if 1, reset variables
  int is_subdivided;          /// 1 if subdivided
  int need_to_update_loading; /// if 1, new time step is updated. Need to update loads too,
                              /// in order to increase loads based on the new time step size
  double loading_factor;      /// loading increments factor
  bool no_subdivision_limits; /// if ture: no limit of subdivisions
                              /// default is false. There is a limit
                              
  const int max_subdivision_allowed = 10000;   /// maximum number of subdivision allowed
  const int min_dt_allowed          = 1.0e-15; /// minimun time step size allowed in the subdivision scheme 
  
  SubdivisionScheme(){
    no_subdivision_limits = false;    
    step_id = decellerate = accellerate = 0;
    step_size = 1;
    dt_0 = 0.0;
  }
  SubdivisionScheme(const bool subdivision_limit){
    no_subdivision_limits = subdivision_limit;
    step_id = decellerate = accellerate = 0;
    step_size = 1;
    dt_0 = 0.0;    
  }
  /// subdevide time step size
  void do_subdivision(int was_NR_ok,
                      double *dt,
                      double *times,
                      long tim,
                      int iter,
                      int max_iter,
                      double alpha,
                      const pgfem3d::CommunicationStructure *com);
};

double subdiv_arc (long INFO,
                   double *dt,
                   double dt0,
                   long *STEP,
                   long *DIV,
                   long tim,
                   double *times,
                   long *ST,
                   long ne,
                   long ndofd,
                   long npres,
                   Element *elem,
                   CRPL *crpl,
                   EPS *eps,
                   SIG *sig,
                   SUPP sup,
                   double *sup_defl,
                   double *rr,
                   double *d_r,
                   double *D_R,
                   double *f_defl,
                   double *f,
                   long *GAMA,
                   double *DT,
                   long *OME,
                   double stab,
                   double dAL0,
                   double dAL,
                   double dALMAX,
                   double nor_min,
                   double dlm0,
                   long *ITT,
                   long iter,
                   long iter_max,
                   long TYPE,
                   const pgfem3d::CommunicationStructure *com,
                   const int analysis);

#endif /* #define PGFEM3D_SUBDIVISION_H */
