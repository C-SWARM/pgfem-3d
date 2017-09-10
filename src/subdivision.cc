#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

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

static constexpr int MAX_STEP = 10000;
static constexpr double MIN_D_TIME = 1.0e-15;
static constexpr int periodic = 0;

/// determine time step size based on the previous time step convergence history
/// and physics based evolution rate.
///
/// \param[in] was_NR_ok if 1, previous iteration was successful
/// \param[in, out] sp container of subdivision parameters
/// \param[in, out] dt time step size, new value will be updated
/// \param[in, out] times time at t(tim-1), t(tim), and t(tim+1), times[tim] will be updated
/// \param[in] tim time step id
/// \param[in] iter number of iteration taken in the iterative solver
/// \param[in] maximum number of iteration defined in the iterative solver
/// \param[in] alpha physics based evolution parameters
/// \return non-zero on internal error
int subdivision_scheme(int was_NR_ok,
                       SUBDIVISION_PARAM *sp,
                       double *dt,
                       double *times,
                       long tim,
                       int iter,
                       int max_iter,
                       double alpha,
                       MPI_Comm mpi_comm)
{
  int err = 0;

  int step_size   = sp->step_size;
  int step_id     = sp->step_id;
  int decellerate = sp->decellerate;
  int accellerate = sp->accellerate;
  double dt_0     = sp->dt_0;

  int is_subdivided_when_step_size_eq_1 = 0;

  sp->reset_variables        = 0;
  sp->is_subdivided          = 0;
  sp->need_to_update_loading = 0;

  long i;

  double nor,gama{};

  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);

  if(was_NR_ok == 0) // Step converged
  {
    decellerate = 0; // reset decelleration flag
    if(step_id == 0) //  1st subdivision
    {
      if(tim == 0)
        step_size = 1;
      else
      {
        // not 1st step, compute subdivision based on previous step
        // compute step_size based on previous convergence properties
        // (times[tim-1] is modified by subdivision procedure
        step_size +=((times[tim+1] - times[tim])/(times[tim] - times[tim-1]) - 1);

        if(step_size <= 0)
          step_size = 1;

        *dt = (times[tim+1] - times[tim])/(step_size);
      }
    }
    else
    {
      // step_id != 0
      //===== ACCELERATE =====
      // compute the acceleration scaling parameter. NOTE: smaller
      // gama --> larger acceleration
      if(alpha <= 0.5) // physics allows large acceleration
      {
        if(iter > max_iter)       // max_iter reached, accel. slowly
        {
          dt_0 = *dt;
          gama = 0.91;            // slow acceleration
          accellerate = 0;        // do not over accelerate if catch next step
        }
        else if(accellerate == 0) // first acceleration
        {
          dt_0 = *dt;
          gama = 0.75;
          accellerate = 1;
        }
        else
        {
          gama = 1./exp(accellerate*1./2.);
          accellerate += 1;
        }
      }
      else if(alpha <= 0.8)       // physics nearing optimum
      {
        dt_0 = *dt;
        gama = 0.8;
        accellerate = 0;
      }
      else if(alpha > 0.8)        // physics very near optimum, asymptotically accelerate
      {
        dt_0 = *dt;
        gama = alpha;             // note alpha > 1 automatically decellerates!
        (accellerate) = 0;
      }

      *dt = (times[tim+1] - times[tim])/(step_size)*step_id;
      times[tim] += *dt;

      sp->need_to_update_loading = 1;
      sp->loading_factor = 1.0/step_size*step_id;

      // compute new number of steps
      *dt = dt_0/gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      step_size = round1 (nor);
      if(step_size <= 0)
        step_size = 1;
      *dt = (times[tim+1] - times[tim])/(step_size);
      if(step_size == 1)
        is_subdivided_when_step_size_eq_1 = 1;
      step_id = 0;
    }
  }
  else
  {
    accellerate = 0; // reset acceleration flag
    if(step_id == 0) // 1st subdivision
    {
      if(decellerate == 0) // 1st load stepping
      {
        // decellerate based on physics. If would
        // decellerate more based on physics, then
        // do it, otherwise use exponential
        if(alpha > 4./3.)
        {
          dt_0 = *dt;
          gama = 1.0/alpha;
          decellerate = 0;
        }
        else
        {
          dt_0 = *dt;
          gama = 0.75;
          decellerate = 1;
        }
      }
      else
      {
        gama = 1./exp(decellerate*1./2.);
        decellerate += 1;
      }

      // compute number of steps
      *dt = dt_0*gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      i = round1 (nor);
      if(i <= step_size)
        step_size += 1;
      else
        step_size = i;
      if(step_size <= 0)
        step_size = 1;
      *dt = (times[tim+1] - times[tim])/(step_size);
      accellerate = 0;
    }
    else
    {
      // load is subdivided and successfully completed at least
      // one step. Update time, compute remaining load and subdivide again
      dt_0 = *dt;
      gama = 0.75;
      decellerate = 1;
      *dt = (times[tim+1] - times[tim])/(step_size)*step_id;
      times[tim] += *dt;

      sp->need_to_update_loading = 1;
      sp->loading_factor = 1.0/step_size*step_id;

      // compute number of steps
      *dt = dt_0*gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      i = round1 (nor);
      if(i <= (step_size-step_id))
        step_size += 1;
      if(step_size <= 0)
        step_size = 1;
      *dt = (times[tim+1] - times[tim])/(step_size);
      step_id = 0;
      accellerate = 0;
    }

    if(myrank == 0)
    {
      PGFEM_printf("\n[%ld] Sub. steps = %ld :: gama = %e ||"
                   " Time %e | dt = %e\n", tim,step_size,gama,times[tim+1],*dt);
    }

    sp->reset_variables = 1;

    if(gama < 1.0 / MAX_STEP || step_size > MAX_STEP || *dt  < MIN_D_TIME)
    {
      if(myrank == 0)
        PGFEM_printf ("Error in Subdivision routine\n");
      PGFEM_Comm_code_abort (mpi_comm,i);
    }
  }

  if(step_size > 1 || is_subdivided_when_step_size_eq_1)
    sp->is_subdivided = 1;

  // update new time steps
  sp->step_size   = step_size;
  sp->step_id     = step_id;
  sp->decellerate = decellerate;
  sp->accellerate = accellerate;
  sp->dt_0        = dt_0;

  return err;
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
                   MPI_Comm mpi_comm,
                   const int analysis)
{
  long i;
  double nor,gama,DLM0,I1,I2;

  static int last_it  = 1;

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

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
        PGFEM_Comm_code_abort(mpi_comm,0); }
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
        PGFEM_Comm_code_abort(mpi_comm,0);
      }
    }/* end INFO = 1 */

    last_it = *ITT; // Save last number of iterations in static variale

    *STEP = 1; *DIV = 0;

    if (myrank == 0)
      PGFEM_printf ("I1 = %f || I2 = %f :: dAL = %12.12f\n",I1,I2,dAL);
  }/* end TYPE = 1 */

  return (DLM0);
}
