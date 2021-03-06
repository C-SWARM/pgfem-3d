#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "pgfem3d/Communication.hpp"
#include "subdivision.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "PGFEM_io.h"
#include "constitutive_model.h"
#include "elem3d.h"
#include "enumerations.h"
#include "matice.h"
#include "res_fini_def.h"
#include "stabilized.h"
#include <cmath>

using namespace pgfem3d;
using namespace multiscale::net;
                       
/// determine time step size based on the previous time step convergence history
/// and physics based evolution rate.
///
/// \param[in] was_NR_ok if 1, previous iteration was successful
/// \param[in, out] dt time step size, new value will be updated
/// \param[in, out] times time at t(tim-1), t(tim), and t(tim+1), times[tim] will be updated
/// \param[in] tim time step id
/// \param[in] iter number of iteration taken in the iterative solver
/// \param[in] maximum number of iteration defined in the iterative solver
/// \param[in] alpha physics based evolution parameters
/// \return non-zero on internal error
void 
SubdivisionScheme::do_subdivision(int was_NR_ok,
                                  double *dt,
                                  double *times,
                                  long tim,
                                  int iter,
                                  int max_iter,
                                  double alpha,
                                  const CommunicationStructure *com)
{
  int var_step_size   = this->step_size;
  int var_step_id     = this->step_id;
  int var_decellerate = this->decellerate;
  int var_accellerate = this->accellerate;
  double var_dt_0     = this->dt_0;

  int is_subdivided_when_step_size_eq_1 = 0;

  this->reset_variables        = 0;
  this->is_subdivided          = 0;
  this->need_to_update_loading = 0;

  long i;

  double nor,gama{};

  int myrank = com->rank;

  if(was_NR_ok == 0) // Step converged
  {
    var_decellerate = 0; // reset decelleration flag
    if(var_step_id == 0) //  1st subdivision
    {
      if(tim == 0)
        var_step_size = 1;
      else
      {
        // not 1st step, compute subdivision based on previous step
        // compute step_size based on previous convergence properties
        // (times[tim-1] is modified by subdivision procedure
        var_step_size +=((times[tim+1] - times[tim])/(times[tim] - times[tim-1]) - 1);

        if(var_step_size <= 0)
          var_step_size = 1;

        *dt = (times[tim+1] - times[tim])/(var_step_size);
      }
    }
    else
    {
      // var_step_id != 0
      //===== ACCELERATE =====
      // compute the acceleration scaling parameter. NOTE: smaller
      // gama --> larger acceleration
      if(alpha <= 0.5) // physics allows large acceleration
      {
        if(iter > max_iter)       // max_iter reached, accel. slowly
        {
          var_dt_0 = *dt;
          gama = 0.91;            // slow acceleration
          var_accellerate = 0;        // do not over accelerate if catch next step
        }
        else if(var_accellerate == 0) // first acceleration
        {
          var_dt_0 = *dt;
          gama = 0.75;
          var_accellerate = 1;
        }
        else
        {
          gama = 1./exp(var_accellerate*1./2.);
          var_accellerate += 1;
        }
      }
      else if(alpha <= 0.8)       // physics nearing optimum
      {
        var_dt_0 = *dt;
        gama = 0.8;
        var_accellerate = 0;
      }
      else if(alpha > 0.8)        // physics very near optimum, asymptotically accelerate
      {
        var_dt_0 = *dt;
        gama = alpha;             // note alpha > 1 automatically decellerates!
        (var_accellerate) = 0;
      }

      *dt = (times[tim+1] - times[tim])/(var_step_size)*var_step_id;
      times[tim] += *dt;

      this->need_to_update_loading = 1;
      this->loading_factor = 1.0/var_step_size*var_step_id;

      // compute new number of steps
      *dt = var_dt_0/gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      var_step_size = round1 (nor);
      if(var_step_size <= 0)
        var_step_size = 1;
      *dt = (times[tim+1] - times[tim])/(var_step_size);
      if(var_step_size == 1)
        is_subdivided_when_step_size_eq_1 = 1;
      var_step_id = 0;
    }
  }
  else
  {
    if(var_step_id == 0) // 1st subdivision
    {
      if(var_decellerate == 0) // 1st load stepping
      {
        // decellerate based on physics. If would
        // decellerate more based on physics, then
        // do it, otherwise use exponential
        if(alpha > 4./3.)
        {
          var_dt_0 = *dt;
          gama = 1.0/alpha;
          var_decellerate = 0;
        }
        else
        {
          var_dt_0 = *dt;
          gama = 0.75;
          var_decellerate = 1;
        }
      }
      else
      {
        gama = 1./exp(var_decellerate*1./2.);
        var_decellerate += 1;
      }

      // compute number of steps
      *dt = var_dt_0*gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      i = round1 (nor);
      if(i <= var_step_size)
        var_step_size += 1;
      else
        var_step_size = i;
      if(var_step_size <= 0)
        var_step_size = 1;
      *dt = (times[tim+1] - times[tim])/(var_step_size);
      var_accellerate = 0;
    }
    else
    {
      // load is subdivided and successfully completed at least
      // one step. Update time, compute remaining load and subdivide again
      var_dt_0 = *dt;
      gama = 0.75;
      var_decellerate = 1;
      *dt = (times[tim+1] - times[tim])/(var_step_size)*var_step_id;
      times[tim] += *dt;

      this->need_to_update_loading = 1;
      this->loading_factor = 1.0/var_step_size*var_step_id;

      // compute number of steps
      *dt = var_dt_0*gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      i = round1 (nor);
      if(i <= (var_step_size-var_step_id))
        var_step_size += 1;
      if(var_step_size <= 0)
        var_step_size = 1;
      *dt = (times[tim+1] - times[tim])/(var_step_size);
      var_step_id = 0;
      var_accellerate = 0;
    }

    if(myrank == 0)
    {
      PGFEM_printf("\n[%ld] Sub. steps = %ld :: gama = %e ||"
                   " Time %e | dt = %e\n", tim,var_step_size,gama,times[tim+1],*dt);
    }

    this->reset_variables = 1;

    if(gama < 1.0 / this->max_subdivision_allowed || 
       var_step_size > this->max_subdivision_allowed || 
       *dt  < this->min_dt_allowed)
    {
      if(this->no_subdivision_limits){
        if(myrank == 0)
          PGFEM_printf ("Subdivision reaches maximum limits, but keep proceeding.\n");
      } else{
        if(myrank == 0)
          PGFEM_printf ("Error in Subdivision routine\n");
        PGFEM_Comm_code_abort (com,i);
      }
    }
  }

  if(var_step_size > 1 || is_subdivided_when_step_size_eq_1)
    this->is_subdivided = 1;

  // update new time steps
  this->step_size   = var_step_size;
  this->step_id     = var_step_id;
  this->decellerate = var_decellerate;
  this->accellerate = var_accellerate;
  this->dt_0        = var_dt_0;

}

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
           const CommunicationStructure *com,
                   const int analysis)
{
  long i;
  double nor,gama,DLM0,I1,I2;

  static int last_it  = 1;

  int myrank = com->rank;

  DLM0 = 0.0;

  if (tim == 0 && iter == 0 && INFO == 0)
    return ((*dt/dt0)*dAL0);

  if (TYPE == 0){/* Time subdivision scheme */
    if (INFO == 0) {
      if (*DIV == 0){
        nor = (times[tim+1] - times[tim])/(times[tim] - times[tim-1]);
        *STEP =  round1 (nor) - 1;
        if (*STEP <= 0)
          *STEP = 1;
        *dt = (times[tim+1] - times[tim])/(*STEP);

        /*
          DDT = (times[tim+1] - times[tim])/(times[tim] - times[tim-1]);
          I1 = *ITT-1;
          I2 = (iter_max+2)/2;
          if (I1 <= 0)
          I1 = 1;
          *dt = sqrt (I2/I1)*DDT;
          nor = (times[tim+1] - times[tim])/(*dt);
          *STEP = round1 (nor);
          if (*STEP <= 0)
          *STEP = 1;
          *dt = (times[tim+1] - times[tim])/(*STEP);
          PGFEM_printf ("NS =  %ld || I1 = %ld || I2 = %ld ::"
          " Time %f | dt = %10.10f\n",*STEP,*ITT-1,
          (iter_max+2)/2,times[tim+1],*dt);
        */
      }/* end DIV == 0 */
      else{
        if (*OME == 0) {
          *DT = *dt;
          gama = 0.75;
          *OME = 1;
        }
        else {
          gama = 1./exp(*OME*1./2.);
          *OME += 1;
        }
        *dt = (times[tim+1] - times[tim])/(*STEP)**DIV;
        times[tim] += *dt;

        *dt = *DT/gama;
        nor = (times[tim+1] - times[tim])/(*dt);
        *STEP = round1 (nor);
        if (*STEP <= 0)
          *STEP = 1;

        /* Restart  */
        if (fabs(dlm0) < nor_min) {
          if (myrank == 0) PGFEM_printf ("Trying increment length restart\n");
          *STEP = 1;
        }
        *dt = (times[tim+1] - times[tim])/(*STEP);
        if (*STEP == 1)
          *ST = 1;
        *DIV = 0;

        /* DDT = *dt; */
        /* *dt = (times[tim+1] - times[tim])/(*STEP)**DIV; */
        /* times[tim] += *dt; */
        /* I1 = *ITT-1; */
        /* I2 = (iter_max+2)/2; */
        /* if (I1 <= 0) */
        /*   I1 = 1; */
        /* *dt = sqrt (I2/I1)*DDT; */

        /* nor = (times[tim+1] - times[tim])/(*dt); */
        /* *STEP = round1 (nor); */
        /* if (*STEP <= 0) */
        /*   *STEP = 1; */

        /* /\* Restart *\/ */
        /* if (fabs(dlm0) < nor_min) { */
        /*   PGFEM_printf ("Trying increment length restart\n"); */
        /*   *STEP = 1; */
        /*   } */

        /* *dt = (times[tim+1] - times[tim])/(*STEP); */
        /* if (*STEP == 1) */
        /*   *ST = 1; */
        /* *DIV = 0; */
        /* PGFEM_printf ("STEP = %ld :: NS =  %ld || I1 = %ld || I2 = %ld" */
        /*         " :: Time %f | dt = %10.10f\n",*DIV,*STEP,*ITT-1, */
        /*    (iter_max+2)/2,times[tim+1],*dt); */

      }
    }/* end INFO = 0 */
    else{
      if (*DIV == 0){
        if (*GAMA == 0) {
          *DT = *dt;
          gama = 0.75;
          *GAMA = 2;
        } else {
          gama = 1./exp(*GAMA*1./2.);
          *GAMA += 2;
        }

        *dt = *DT*gama;
        nor = (times[tim+1] - times[tim])/(*dt);
        i = round1 (nor);
        if (i <= *STEP)
          *STEP += 1;
        else
          *STEP = i;
        if (*STEP <= 0)
          *STEP = 1;

        /* Limit long overflow */
        if (*STEP >= 200000000) {
          *STEP = 1;
          *DIV = 0;
          *GAMA = 0;
          gama = 0.0;
        }

        *dt = (times[tim+1] - times[tim])/(*STEP);
        *OME = 0;
      }/* end DIV = 0 */
      else{
        *DT = *dt;
        gama = 0.75;
        *GAMA = 2;

        *dt = (times[tim+1] - times[tim])/(*STEP)**DIV;
        times[tim] += *dt;
        *dt = *DT*gama;
        nor = (times[tim+1] - times[tim])/(*dt);
        i = round1 (nor);
        if (i <= (*STEP-*DIV))
          *STEP += 1;
        if (*STEP <= 0)
          *STEP = 1;

        /* Limit long overflow */
        if (*STEP >= 200000000) {
          *STEP = 1;
          *DIV = 0;
          *GAMA = 0;
          gama = 0.0;
        }

        *dt = (times[tim+1] - times[tim])/(*STEP);
        *DIV = 0;
        *OME = 0;
      }/* DIV != 0 */

      if (myrank == 0)
        PGFEM_printf ("\n[%ld] Sub. steps = %ld :: gama = %8.8e ||"
                      " Time %f | dt = %10.10f\n",tim,*STEP,
                      gama,times[tim+1],*dt);

      /* Reset variables */
      switch(analysis){
       case FS_CRPL:
       case FINITE_STRAIN:
        res_fini_def (ne,npres,elem,eps,sig,crpl,analysis);
        break;
       case STABILIZED:
        res_stab_def (ne,npres,elem,eps,sig,stab);
        break;
       case MINI:
        MINI_reset(elem,ne,npres,sig);
        break;
       case MINI_3F:
        MINI_3f_reset(elem,ne,npres,4,sig,eps);
        break;
       default: break;
      }

      /* reset damage variables */
      for(i=0; i<ne; i++){
        long n_ip;
        int_point(elem[i].toe,&n_ip);
        for(int j=0; j<n_ip; j++){
          reset_damage(&eps[i].dam[j]);
        }
        if(analysis == STABILIZED){/* get pressure terms too */
          int_point(10,&n_ip);
          for(int j=0; j<n_ip; j++){
            reset_damage(&eps[i].st[j].dam);
          }
        }
      }

      for (i=0;i<ndofd;i++) {
        rr[i] = d_r[i] = D_R[i] = f_defl[i] = f[i] = 0.0;
      }

      if (gama < nor_min && (*dt/dt0)*dAL0 < nor_min*nor_min) {
        if (myrank == 0)
          PGFEM_printf ("Error in Subdivision rutine\n");
        PGFEM_Comm_code_abort(com, 0);
      }
    }/* end INFO = 1 */

    DLM0 = (*dt/dt0)*dAL0;
  }/* end TYPE =  0 */
  /************************************************************/
  /************************************************************/
  if (TYPE == 1){/* Arc-length subdivision procedure */

    I1 = (double) *ITT;
    I2 = (double) iter_max; //round1((iter_max+2)/2.);

    if (INFO == 0) {
      if (tim == 0){
        DLM0 = dAL0;
      } else {
        if (I1 <= 0) I1 = 1;
        if ( I1 <= last_it && last_it < iter_max){
          I1 = (I1<last_it)?I1:(last_it - 1);
        }
        DLM0 = sqrt (I2/I1)*dAL;
        if (DLM0 > dALMAX) DLM0 = dALMAX;
      }
    }/* end INFO = 0 */
    else{
      if (*GAMA == 0) {
        gama = 0.75;
        *GAMA = 2;
      } else {
        gama = 1./exp(*GAMA*1./2.);
        *GAMA += 2;
      }

      DLM0 = dAL*gama;

      if (myrank == 0)
        PGFEM_printf ("\n[%ld] gama = %8.8e || Time %f | dt = %10.10f\n",
                      tim,gama,times[tim+1],*dt);

      /* Reset variables */
      switch(analysis){
       case FS_CRPL:
       case FINITE_STRAIN:
        res_fini_def (ne,npres,elem,eps,sig,crpl,analysis);
        break;
       case STABILIZED:
        res_stab_def (ne,npres,elem,eps,sig,stab);
        break;
       case MINI:
        MINI_reset(elem,ne,npres,sig);
        break;
       case MINI_3F:
        MINI_3f_reset(elem,ne,npres,4,sig,eps);
        break;
       default: break;
      }

      /* reset damage variables */
      for(i=0; i<ne; i++){
        long n_ip;
        int_point(elem[i].toe,&n_ip);
        for(int j=0; j<n_ip; j++){
          reset_damage(&eps[i].dam[j]);
        }
        if(analysis == STABILIZED){/* get pressure terms too */
          int_point(10,&n_ip);
          for(int j=0; j<n_ip; j++){
            reset_damage(&eps[i].st[j].dam);
          }
        }
      }


      for (i=0;i<ndofd;i++) {
        rr[i] = d_r[i] = D_R[i] = f_defl[i] = f[i] = 0.0;
      }

      if (DLM0 < nor_min*nor_min) {
        if (myrank == 0)
          PGFEM_printf ("Error in Subdivision rutine\n");
        PGFEM_Comm_code_abort(com, 0);
      }
    }/* end INFO = 1 */

    last_it = *ITT; // Save last number of iterations in static variale

    *STEP = 1; *DIV = 0;

    if (myrank == 0)
      PGFEM_printf ("I1 = %f || I2 = %f :: dAL = %12.12f\n",I1,I2,dAL);
  }/* end TYPE = 1 */

  return (DLM0);
}
