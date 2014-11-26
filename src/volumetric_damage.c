/* HEADER */
#include "volumetric_damage.h"
#include "utils.h"
#include "PGFEM_io.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef VD_DEBUG
#define VD_DEBUG 0
#endif

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

static double DAMAGE_THRESH = 0.9999;

/*========================================================*/
/*               Static Function Declarations             */
/*========================================================*/

static void damage_initialize_values(damage *dam);

static int set_damage_function(damage *dam,
			       vd_fun function);

static int set_damage_evolution(damage *dam,
				vd_fun evolution);

static int set_damage_parameters(damage *dam,
				 const double *params,
				 const int len);

static int set_damage_evolution_rate(damage *dam,
				     vd_fun evolution_rate);

/*========================================================*/
/*                 Main interface                         */
/*========================================================*/
void init_damage(damage *dam,
		 vd_fun function,
		 vd_fun evolution,
		 vd_fun evolution_rate,
		 const double *params,
		 const int len)
{
  int err = 0;
  damage_initialize_values(dam);
  if(set_damage_function(dam,function) != 0){
    PGFEM_printerr("set_damage_function returned error in %s.\n",__func__);
    err = 1;
  }
  if(set_damage_evolution(dam,evolution) != 0){
    PGFEM_printerr("set_damage_evolution returned error in %s.\n",__func__);
    err = 1;
  }
  if(set_damage_parameters(dam,params,len) != 0){
    PGFEM_printerr("set_damage_parameters returned error in %s.\n",__func__);
    err = 1;
  }
  if(err != 0){
    abort();
  }
}

void init_damagef(damage *dam,
		  const int eq_flag,
		  const double *params,
		  const int len)
{
  int perr = 0;

  damage_initialize_values(dam);
  dam->eq_flag = eq_flag;

  if(set_damage_parameters(dam,params,len) != 0){
    PGFEM_printerr("set_damage_parameters returned error in %s.\n",__func__);
    perr = 1;
  }
  
  /* damage functions */
  reset_damage_functions(dam,eq_flag);

  if(perr){
    abort();
  }
}

void reset_damage_functions(damage *dam,
			    const int eq_flag)
{
  int ferr = 0;
  int eerr = 0;
  int ererr = 0;

  switch(eq_flag){
  default:
    ferr = set_damage_function(dam,weibull_function);
    eerr = set_damage_evolution(dam,weibull_evolution);
    ererr = set_damage_evolution_rate(dam,weibull_evolution_rate);
    break;
  }

  if(ferr != 0){
    PGFEM_printerr("set_damage_function returned error in %s.\n",__func__);
    ferr = 1;
  }
  if(eerr != 0){
    PGFEM_printerr("set_damage_evolution returned error in %s.\n",__func__);
    eerr = 1;
  }
  if(ererr != 0){
    PGFEM_printerr("set_damage_evolution_rate returned error in %s.\n",__func__);
    ererr = 1;
  }
  if(ferr+eerr+ererr != 0){
    abort();
  }
}

void copy_damage(damage *restrict dest,
		 const damage *restrict src)
{
  if(dest == src) return;
  memcpy(dest,src,sizeof(*src));
  /* may not be necessary, but do it anyhow */
  reset_damage_functions(dest,dest->eq_flag);
}

size_t sizeof_damage(const damage *dam)
{
  size_t s = 0;
  s += sizeof(dam->params);
  s += 9*sizeof(double) + 4*sizeof(int);
  return s;
}

void pack_damage(const damage *src,
		 char *buffer,
		 size_t *pos)
{
  pack_data(&(src->params),buffer,pos,1,sizeof(src->params));
  pack_data(&(src->wn),buffer,pos,1,sizeof(src->wn));
  pack_data(&(src->w),buffer,pos,1,sizeof(src->w));
  pack_data(&(src->Hn),buffer,pos,1,sizeof(src->Hn));
  pack_data(&(src->H),buffer,pos,1,sizeof(src->H));
  pack_data(&(src->Hpn),buffer,pos,1,sizeof(src->Hpn));
  pack_data(&(src->Hp),buffer,pos,1,sizeof(src->Hp));
  pack_data(&(src->Xn),buffer,pos,1,sizeof(src->Xn));
  pack_data(&(src->X),buffer,pos,1,sizeof(src->X));
  pack_data(&(src->dmu),buffer,pos,1,sizeof(src->dmu));
  pack_data(&(src->damaged_n),buffer,pos,1,sizeof(src->damaged_n));
  pack_data(&(src->damaged),buffer,pos,1,sizeof(src->damaged));
  pack_data(&(src->broken),buffer,pos,1,sizeof(src->broken));
  pack_data(&(src->eq_flag),buffer,pos,1,sizeof(src->eq_flag));
}

void unpack_damage(damage *dest,
		   const char *buffer,
		   size_t *pos)
{
  unpack_data(buffer,&(dest->params),pos,1,sizeof(dest->params));
  unpack_data(buffer,&(dest->wn),pos,1,sizeof(dest->wn));
  unpack_data(buffer,&(dest->w),pos,1,sizeof(dest->w));
  unpack_data(buffer,&(dest->Hn),pos,1,sizeof(dest->Hn));
  unpack_data(buffer,&(dest->H),pos,1,sizeof(dest->H));
  unpack_data(buffer,&(dest->Hpn),pos,1,sizeof(dest->Hpn));
  unpack_data(buffer,&(dest->Hp),pos,1,sizeof(dest->Hp));
  unpack_data(buffer,&(dest->Xn),pos,1,sizeof(dest->Xn));
  unpack_data(buffer,&(dest->X),pos,1,sizeof(dest->X));
  unpack_data(buffer,&(dest->dmu),pos,1,sizeof(dest->dmu));
  unpack_data(buffer,&(dest->damaged_n),pos,1,sizeof(dest->damaged_n));
  unpack_data(buffer,&(dest->damaged),pos,1,sizeof(dest->damaged));
  unpack_data(buffer,&(dest->broken),pos,1,sizeof(dest->broken));
  unpack_data(buffer,&(dest->eq_flag),pos,1,sizeof(dest->eq_flag));
  /* may not be necessary, but do it anyhow */
  reset_damage_functions(dest,dest->eq_flag);
}

void reset_damage(damage *dam){
  dam->X = dam->Xn;
  dam->w = dam->wn;
  dam->H = dam->Hn;
  dam->Hp = dam->Hn;
  dam->damaged = dam->damaged_n;
}

void update_damage(damage *dam){
  dam->Xn = dam->X;
  dam->wn = dam->w;
  dam->Hn = dam->H;
  dam->Hpn = dam->Hp;
  dam->damaged_n = dam->damaged;
}

double damage_int_alg(damage *dam,
		      const double Ybar,
		      const double dt)
{
  dam->dmu = dt*dam->params.mu;
  double G = dam->function(Ybar,&(dam->params));
  double g = G - dam->Xn;

  if (g > 0.0){ /* Damage propagation */
    /* flag integration point */
    dam->damaged = 1;

    /* update damage parameter */
    dam->w = dam->wn + dam->dmu/(1+dam->dmu)*g;

    /* update softening parameter */
    dam->X = MAX(dam->Xn,(dam->Xn + dam->dmu*G)/(1+dam->dmu));

    /* update evolution parameter */
    dam->H = dam->evolution(Ybar,&(dam->params));

    /* update evolution rate parameter */
    dam->Hp = dam->evolution_rate(Ybar,&(dam->params));

  } else { /* no damage propagation */
    dam->damaged = 0;
    dam->w = dam->wn;
    dam->X = dam->Xn;
    dam->H = 0.0;
    dam->Hp = 0.0;
  }

  return g;
}

/*========================================================*/
/*                     Damage Equations  (G)              */
/*========================================================*/
double weibull_function(const double Ybar,
			const damage_params *params)
{
  if(Ybar <= params->Yin) return 0.0;
  return (DAMAGE_THRESH - DAMAGE_THRESH
	  *exp(-pow((Ybar-params->Yin)/(params->p1*params->Yin),
		    params->p2)
	       )
	  );

}

/*========================================================*/
/*               Damage Evolution Equations  (H)          */
/*========================================================*/
double weibull_evolution(const double Ybar,
			 const damage_params *params)
{
  if(Ybar <= params->Yin) return 0.0;
  return (DAMAGE_THRESH*params->p2/(params->p1*params->Yin)
	  *exp(-pow((Ybar - params->Yin)/(params->p1*params->Yin),params->p2))
	  *pow((Ybar - params->Yin)/(params->p1*params->Yin),params->p2-1.0)
	  );
}

/*========================================================*/
/*         Damage Rate of Evolution Equations  (Hp)       */
/*========================================================*/	      
double weibull_evolution_rate(const double Ybar,
			      const damage_params *params)
{
  if(Ybar <= params->Yin) return 0.0;
  return (
	  -DAMAGE_THRESH/pow(Ybar-params->Yin,2.)
	  *exp(-pow((Ybar-params->Yin)/(params->p1*params->Yin),params->p2))
	  *params->p2*(1.
		       + params->p2*(-1 + pow((Ybar-params->Yin)
					      /(params->p1*params->Yin)
					      ,params->p2)))
	  *pow((Ybar-params->Yin)/(params->p1*params->Yin),params->p2)
	  );
}

/*========================================================*/
/*               Static Function Definitions              */
/*========================================================*/
static void damage_params_initialize_vals(damage_params *params)
{
  params->Yin = params->p1 = params->p2 = params->mu = 0.0;
}

static void damage_initialize_values(damage *dam)
{
  /* intialize damage parameters */
  damage_params_initialize_vals(&(dam->params));

  /* initialize internal variables */
  dam->wn = dam->w = 0.0;
  dam->Hn = dam->H = 0.0;
  dam->Hpn = dam->Hp = 0.0;
  dam->Xn = dam->X = 0.0;
  dam->dmu = 0.0;
  dam->damaged_n = dam->damaged = 0;
  dam->broken = 0;
  dam->eq_flag = 0;
  dam->function = NULL;
  dam->evolution = NULL;
  dam->evolution_rate = NULL;
}

static int set_damage_function(damage *dam,
			       vd_fun function)
{
  dam->function = function;
  if(dam->function != NULL) return 0;
  else return 1;
}

static int set_damage_evolution(damage *dam,
				vd_fun evolution)
{
  dam->evolution = evolution;
  if(dam->evolution != NULL) return 0;
  else return 1;
}

static int set_damage_evolution_rate(damage *dam,
				     vd_fun evolution_rate)
{
  dam->evolution_rate = evolution_rate;
  if(dam->evolution_rate != NULL) return 0;
  else return 1;
}

static int set_damage_parameters(damage *dam,
				 const double *params,
				 const int len)
{
  /* return error if not enough parameters given */
  if(len < sizeof(damage_params)/sizeof(double)){
    return 1;
  }

  /* parameters taken in order listed */
  dam->params.Yin = params[0];
  dam->params.p1  = params[1];
  dam->params.p2  = params[2];
  dam->params.mu  = params[3];

  return 0;
}
