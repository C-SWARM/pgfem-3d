/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  Aaron Howell, [1], <ahowell3@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFEM_io.h"
#include "PGFem3D_data_structure.h"
#include "allocation.h"
#include "constitutive_model.h"
#include "constitutive_model_3f.h"
#include "cm_placeholder_functions.h"
#include "cm_iso_viscous_damage.h"
#include "cm_j2_plasticity.h"
#include "cm_uqcm.h"
#include "cm_poro_viscoplasticity.h"
#include "cm_iso_viscous_damage_split.h"
#include "cm_MMS.h"
#include "dynamics.h"
#include "hommat.h"
#include "supp.h"
#include "elem3d.h"
#include "femlib.h"
#include "get_dof_ids_on_elem.h"
#include "hommat.h"
#include "hyperelasticity.h"                 // <= constitutive model elasticity
#include "index_macros.h"
#include "material_properties.h"    // <= constitutive model material properties
#include "plasticity_model_none.h"
#include "plasticity_model.h"
#include "plasticity_model_BPA.h"
#include "supp.h"
#include "utils.h"
#include <ttl/ttl.h>
#include <algorithm>
#include "condense.h"
#include "new_potentials.h"

using namespace pgfem3d;

static constexpr int       DIM_3 = 3;
static constexpr int     DIM_3x3 = 9;
static constexpr int   DIM_3x3x3 = 27;
static constexpr int DIM_3x3x3x3 = 81;
static double I[DIM_3x3] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};

//ttl declarations
namespace {
  template<int R, int D = 3, class S = double>
  using Tensor = ttl::Tensor<R, D, S>;

  template<int R, int D = 3, class S = double *>
  using TensorA = ttl::Tensor<R, D, S>;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'o'> o;
  static constexpr ttl::Index<'p'> p;

  template<class T1, class T2>
  inline int inv(T1 &A, T2 &AI) {
    int err = inv3x3(A.data, AI.data);
    return err;
  }
  template<class T1, class T2>
  inline int symm(T1 &A, T2 &Asym)
  {
    Asym(i,j) = 0.5*(A(i,j)+A(j,i));
    return 0;
  }
}

/// print constitutive model interface infomation
void print_constitutive_model_info(FILE *out)
{
  const CMVariableNames variable_names[] = {{HYPER_ELASTICITY,         "Hyper Elasticity         "},        
                                            {CRYSTAL_PLASTICITY,       "Crystal Plasticity       "},      
                                            {BPA_PLASTICITY,           "BPA Plasticity           "},          
                                            {ISO_VISCOUS_DAMAGE,       "Iso Viscous Damage       "},      
                                            {J2_PLASTICITY_DAMAGE,     "J2 Plasticity with Damage"},    
                                            {POROVISCO_PLASTICITY,     "Poro Visco Plasticity    "},    
                                            {ISO_VISCOUS_SPLIT_DAMAGE, "Splited Iso Viscous Dmage"},
                                            {MANUFACTURED_SOLUTIONS,   "Method of Manufactured Solutions (use only for MMS study)"}
                                          };

  for(int ia=0; ia<NUM_MODELS; ++ia)
    PGFEM_fprintf(out, "\t\t%s %d\n", variable_names[ia].name, variable_names[ia].id);
}


// define identity
static TensorA<2> delta_ij(I);

/// compute nodal temperature in an element
///
/// if time step is subdivided, temperature needs to be subdivided
/// dT   = T(t(n+1)) - T(t(n))
/// Tn_e = T(t(n)) + dT*factor_n
/// Tn_e = T(t(n)) + dT*factor_n
///
/// T(t=n)---------Tn_e----Tnp1_e-----T(t=np1)
/// |<--factor_n-->|       |
/// |<-------factor_np1--->|
///
/// \param[in] fe container of finite element resources for an element
/// \param[in] grid an object containing all mesh info
/// \param[in] fv_h field variable object for thermal
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[out] Tnp1_e computed nodal temperature for current element at t(n+1)
/// \param[out] Tn_e computed nodal temperature for current element at t(n+1)
/// \param[out] Tnm1_e computed nodal temperature for current element at t(n+1)
/// \param[in] factor_np1 factor of computing temperature for t(n+1)
/// \param[in] factor_n factor of computing temperature for t(n)
/// \return non-zero on internal error
int get_nodal_temperatures(FEMLIB *fe,
                           Grid *grid,
                           FieldVariables *fv_h,
                           LoadingSteps *load,
                           int mp_id,
                           double *Tnp1_e,
                           double *Tn_e,
                           double *Tnm1_e,
                           double factor_np1,
                           double factor_n)
{
  int err = 0;

  int nne = fe->nne;
  long *cn = (long *) aloc1l(nne);
  long *nod = (fe->node_id).m_pdata;
  int ndofn = fv_h->ndofn;
  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cn,mp_id);

  double T0    = fv_h->u0;
  double *T    = fv_h->u_np1;
  double *dT   = fv_h->d_u;
  double *Tn   = fv_h->temporal->u_n;
  double *Tnm1 = fv_h->temporal->u_nm1;

  SUPP sup = load->sups[mp_id];

  for(int ia=0; ia<nne; ia++)
  {
    Tnm1_e[ia] = Tnm1[nod[ia]];
    const int id = cn[ia];
    const int aid = abs(id) - 1;

    double T_n = Tn[nod[ia]];
    double T_np1 = 0.0;

    if (id == 0)
      T_np1 = T0;
    else if(id > 0)
      T_np1 = T[aid] + dT[aid];
    else
    {
      T_np1 = T0 + sup->defl[aid] + sup->defl_d[aid];
    }

      Tn_e[ia] = T_n + (T_np1 - T_n)*factor_n;
    Tnp1_e[ia] = T_n + (T_np1 - T_n)*factor_np1;
  }

  free(cn);
  return err;
}

/// compute temperature dependent variables at the integration point
///
/// \param[in] fe container of finite element resources for an element
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] T0 initial temperature
/// \param[in] Tnp1 nodal temperature array at t(n+1)
/// \param[in] Tn nodal temperature array at t(n)
/// \param[in] Tnm1 nodal temperature array at t(n-1)
/// \param[out] hFnp1 computed deformation gradient at t(n+1)
/// \param[out] hFn computed deformation gradient at t(n)
/// \param[out] hFnm1 computed deformation gradient at t(n-1)
/// \return non-zero on internal error
int compute_temperature_at_ip(FEMLIB *fe,
                              const Grid *grid,
                              const MaterialProperty *mat,
                              const double T0,
                              double *Tnp1,
                              double *Tn,
                              double *Tnm1,
                              double *hFnp1,
                              double *hFn,
                              double *hFnm1)
{
  int err = 0;
  double T     = 0.0;
  double dTn   = 0.0;
  double dTnp1 = 0.0;
  double dTnm1 = 0.0;
  double *N = (fe->N).m_pdata;

  for(int ia=0; ia<fe->nne; ia++)
  {
    T     += N[ia]*Tnp1[ia];
    dTnp1 += N[ia]*(Tnp1[ia] - T0);
    dTn   += N[ia]*(Tn[ia]   - T0);
    dTnm1 += N[ia]*(Tnm1[ia] - T0);
  }

  const int eid = fe->curt_elem_id;
  const int hmat_id = (grid->element[eid]).mat[2];
  const int mat_id = (mat->hommat[hmat_id]).mat_id;

  const double ax = mat->mater[mat_id].ax;
  const double ay = mat->mater[mat_id].ay;
  const double az = mat->mater[mat_id].az;

  hFnp1[1] = hFnp1[2] = hFnp1[3] = hFnp1[5] = hFnp1[6] = hFnp1[7] = 0.0;
    hFn[1] =   hFn[2] =   hFn[3] =   hFn[5] =   hFn[6] =   hFn[7] = 0.0;
  hFnm1[1] = hFnm1[2] = hFnm1[3] = hFnm1[5] = hFnm1[6] = hFnm1[7] = 0.0;

  hFnp1[0] = 1.0 + ax*dTnp1;
  hFnp1[4] = 1.0 + ay*dTnp1;
  hFnp1[8] = 1.0 + az*dTnp1;

  hFn[0] = 1.0 + ax*dTn;
  hFn[4] = 1.0 + ay*dTn;
  hFn[8] = 1.0 + az*dTn;

  hFnm1[0] = 1.0 + ax*dTnm1;
  hFnm1[4] = 1.0 + ay*dTnm1;
  hFnm1[8] = 1.0 + az*dTnm1;

  return err;
}

/* add the macroscopic deformation gradient to the _TOTAL_ deformation
   gradient */
static void cm_add_macro_F(const SUPP sup,
                           double * __restrict F)
{
  const double * __restrict F0 = sup->F0;
  for (int i = 0; i < DIM_3x3; i++) F[i] += F0[i];
}

int construct_model_context(void **ctx,
                            const int type,
                            double *F,
                            const double dt,
                            const double alpha,
                            double *eFnpa,
                            int npa);

int construct_model_context_with_thermal(void **ctx,
                                         const int type,
                                         double *F,
                                         const double dt,
                                         const double alpha,
                                         double *eFnpa,
                                         double *hFn,
                                         double *hFnp1,
                                         int npa);                                                           

// Nodal Temerature for transient
class NodalTemerature
{
  public:
    double T0;
    Matrix<double> np1, n, nm1;

    NodalTemerature(){};
    NodalTemerature(int nne){ initialization(nne);};
    void initialization(int nne)
    {
      T0 = 0.0;
      np1.initialization(nne, 1, 0.0);
        n.initialization(nne, 1, 0.0);
      nm1.initialization(nne, 1, 0.0);
    };
    void get_temperature(FEMLIB *fe,
                         Grid *grid,
                         FieldVariables *fv,
                         LoadingSteps *load,
                         const Multiphysics& mp,
                         int mp_id,
                         int is_it_couple_w_thermal)
    {
      FieldVariables *fv_h = fv->fvs[is_it_couple_w_thermal];
      T0 = fv_h->u0;
      // get temperature for this element for each nodal point

      int mp_cp_id = mp.coupled_ids[mp_id][is_it_couple_w_thermal+1];
      get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                             np1.m_pdata,n.m_pdata,nm1.m_pdata,
                             fv->subdivision_factor_np1,fv->subdivision_factor_n);
    };
};

// deformation gradients for steady state
class DeformationGradient_ss
{
  public:
    double Jn, Jnp1;
    Tensor<2> n, np1;

    DeformationGradient_ss(){ initialization(); };

    void initialization(void)
    {
      Jn = Jnp1 = 1.0;
      n   = delta_ij(i,j);
      np1 = delta_ij(i,j);
    };
    void update_thermal_part(FEMLIB *fe,
                             Grid *grid,
                             MaterialProperty *mat,
                             FieldVariables *fv,
                             NodalTemerature &T)

    {
      // compute thermal part at integration point
        Tensor<2> hFnm1;
        compute_temperature_at_ip(fe,grid,mat,T.T0,
                                  T.np1.m_pdata,T.n.m_pdata,T.nm1.m_pdata,
                                  np1.data,n.data,hFnm1.data);
    };
};

// deformation gradients for transient
class DeformationGradient_tr : public DeformationGradient_ss
{
  public:
    double Jnm1;
    Tensor<2> nm1;
    Tensor<2> npa;
    Tensor<2> nma;

    DeformationGradient_tr(){ initialization(); };

    void initialization(void)
    {
      Jn = Jnp1 = Jnm1 = 1.0;
      nm1 = delta_ij(i,j);
      n   = delta_ij(i,j);
      np1 = delta_ij(i,j);
      npa = delta_ij(i,j);
      nma = delta_ij(i,j);
    };
    void set_npa_nma(double alpha)
    {
      /// F_nma = (1-alpha)*F(n-1) + alpha*F(n)
      /// F_npa = (1-alpha)*F(n)   + alpha*F(n+1)
      nma = (1.0-alpha)*nm1(i,j) + alpha*n(i,j);
      npa = (1.0-alpha)*n(i,j) + alpha*np1(i,j);
    };
    void update_thermal_part(FEMLIB *fe,
                             Grid *grid,
                             MaterialProperty *mat,
                             FieldVariables *fv,
                             NodalTemerature &T)

    {
      // compute thermal part at integration point
      compute_temperature_at_ip(fe,grid,mat,T.T0,
                                T.np1.m_pdata,T.n.m_pdata,T.nm1.m_pdata,
                                np1.data,n.data,nm1.data);
    };
};

template <class T> class DeformationGradient
{
  public:
    Tensor<2> Fr, M;
    T F, eF, pF;
    T *hF;
    DeformationGradient()
    {
      hF = NULL;
      initialization();
    };
    DeformationGradient(T *hF_in){ hF = hF_in;};
    void initialization(void)
    {
      Fr = delta_ij(i,j);
      M  = delta_ij(i,j);
      F.initialization();
      eF.initialization();
      pF.initialization();
      if(hF != NULL)
        hF->initialization();
    };
    void update_total_deformation_gradient(FEMLIB *fe,
                                           FieldVariables *fv,
                                           SUPP sup,
                                           Matrix<double> &u,
                                           const int total_Lagrangian)
    {
      int eid = fe->curt_elem_id;
      int ip = fe->curt_itg_id - 1;
      Constitutive_model *m = (fv->eps[eid]).model + ip;      
      fe->update_deformation_gradient(fv->ndofn,u.m_pdata,Fr.data,(m->param)->pFI);

      if(total_Lagrangian)
      {
        if(sup->multi_scale)
          cm_add_macro_F(sup,Fr.data);

        // Total Lagrangian formulation Fn = 1, Fnp1 = Fr
        F.np1 = Fr(i,j);
      }
      else
      {
        if(sup->multi_scale)
        {
          PGFEM_printerr("Multi-scale formulation does not support UL!\n");
          PGFEM_Abort();
        }
        m->param->get_pF(m, F.n.data, 1);
        F.np1 = Fr(i,k)*F.n(k,j); // Fn+1 = Fr*Fn
      }
    };
    int compute_deformation_gradient_for_ss(FEMLIB *fe,
                                            FieldVariables *fv,
                                            const int total_Lagrangian,
                                            const int is_it_couple_w_thermal,
                                            const int is_it_couple_w_chemical)
    {
      int err = 0;
      int eid = fe->curt_elem_id;
      int ip = fe->curt_itg_id - 1;

      Constitutive_model *m = (fv->eps[eid]).model + ip;
      err += m->param->get_pF(m,pF.np1.data,2);

      if(total_Lagrangian)
      {
        if(is_it_couple_w_thermal>=0)
        {
          Tensor<2> pFnp1_I, hFnp1_I;
          err += inv(pF.np1, pFnp1_I);
          err += inv(hF->np1, hFnp1_I);
          M = hFnp1_I(i,k)*pFnp1_I(k,j);
        }
        else
          err += inv(pF.np1,M);

        eF.np1 = F.np1(i,k)*M(k,j);
      }
      else
      {
        Tensor<2> pFnp1_I, pFn_I;
        err += m->param->get_pF(m, pF.n.data, 1);
        err += m->param->get_F(m,  F.n.data,  1);

        err += inv(pF.n, pFn_I);
        err += inv(pF.np1, pFnp1_I);

        if(is_it_couple_w_thermal>=0)
        {
          Tensor<2> hFn_I, hFnp1_I;
          err += inv(hF->n,   hFn_I);
          err += inv(hF->np1, hFnp1_I);

          M = pF.n(i,k)*hF->n(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);
          eF.n = F.n(i,k)*hFn_I(k,l)*pFn_I(l,j);
        }
        else
        {
          M = pF.n(i,k)*pFnp1_I(k,j);
          eF.n = F.n(i,k)*pFn_I(k,j);
        }

        eF.np1 = Fr(i,k)*eF.n(k,l)*M(l,j);
      }
      return err;
    };
    int compute_deformation_gradient_for_tr(FEMLIB *fe,
                                            FieldVariables *fv,
                                            const int total_Lagrangian,
                                            const int is_it_couple_w_thermal,
                                            const int is_it_couple_w_chemical)
    {
      int err = 0;
      return err;
    }
};

class Integrate
{
  public:
    FieldVariables *fv;
    double *du;
    int eid;
    bool run_integration_algorithm;
    bool run_for_update;
    Integrate(){
      fv = NULL;
      du = NULL;
      eid = -1;
      run_integration_algorithm = false;
      run_for_update = false;
    };

    int inner_loop(Matrix<double> &dMdu,
                   Matrix<double> &dMdt,
                   CM_ThreeField &cm_tf);

    int out_loop(double *out,
                 int nne,
                 int nsd,
                 int Pno,
                 int Vno);
};

class IntegrateThreeFieldStiffness : public Integrate
{
  public:
  ThreeFieldStiffness K;

  IntegrateThreeFieldStiffness(FEMLIB *fe,
                               int Vno,
                               int Pno)
  {
    fv = NULL;
    du = NULL;
    eid = fe->curt_elem_id;
    run_integration_algorithm = false;
    run_for_update = false;
    K.initialization(fe, Vno, Pno);
  };

  int inner_loop(Matrix<double> &dMdu,
                 Matrix<double> &dMdt,
                 CM_ThreeField &cm_tf)
  {
    int err = 0;
    err += K.compute_stiffness(cm_tf, dMdu, dMdt);
    return err;
  };
  int out_loop(double *out,
               int nne,
               int nsd,
               int Pno,
               int Vno)
  {
    int err = 0;
    K.Kpt.trans(K.Ktp);
    K.Kpu.trans(K.Kup);

    err += condense_K_3F_to_1F(out, nne, nsd, Pno, Vno,
                               K.Kuu.m_pdata, K.Kut.m_pdata, K.Kup.m_pdata,
                               K.Ktu.m_pdata, K.Ktt.m_pdata, K.Ktp.m_pdata,
                               K.Kpu.m_pdata, K.Kpt.m_pdata, NULL);

    // check diagonal for zeros/nans
    for (int a = 0; a < nne; a++) {
      for (int b = 0; b < nsd; b++) {
        if ( !isnormal(out[idx_K(a,b,a,b,nne,nsd)]) ) err++;
      }
    }

    return err;
  };

};

class IntegrateThreeFieldResidual : public Integrate
{
  public:
  ThreeFieldStiffness K;
  ThreeFieldResidual R;

  IntegrateThreeFieldResidual(FEMLIB *fe,
                              int Vno,
                              int Pno)
  {
    fv = NULL;
    du = NULL;
    eid = fe->curt_elem_id;
    run_integration_algorithm = true;
    run_for_update = false;
    K.initialization(fe, Vno, Pno, true);
    R.initialization(fe, Vno, Pno);
  };

  int inner_loop(Matrix<double> &dMdu,
                 Matrix<double> &dMdt,
                 CM_ThreeField &cm_tf)
  {
    int err = 0;
    err += R.compute_residual(cm_tf);
    err += K.compute_stiffness(cm_tf, dMdu, dMdt);
    return err;
  };
  int out_loop(double *out,
               int nne,
               int nsd,
               int Pno,
               int Vno)
  {
    int err = 0;
    K.Kpt.trans(K.Ktp);

    err += condense_F_3F_to_1F(out, nne, nsd, Pno, Vno,
                               R.Ru.m_pdata,  R.Rt.m_pdata,  R.Rp.m_pdata,
                               K.Kut.m_pdata, K.Kup.m_pdata, K.Ktp.m_pdata, K.Ktt.m_pdata, K.Kpt.m_pdata);
    return err;
  };

};

class IntegrateThreeFieldUpdate : public Integrate
{
  public:

  ThreeFieldStiffness K;
  ThreeFieldResidual R;

  IntegrateThreeFieldUpdate(FEMLIB *fe,
                            int Vno,
                            int Pno)
  {
    fv = NULL;
    du = NULL;
    eid = fe->curt_elem_id;
    run_integration_algorithm = false;
    run_for_update = true;
    K.initialization(fe, Vno, Pno, true);
    R.initialization(fe, Vno, Pno, true);
  };

  int inner_loop(Matrix<double> &dMdu,
                 Matrix<double> &dMdt,
                 CM_ThreeField &cm_tf)
  {
    int err = 0;
    err += R.compute_residual(cm_tf);
    err += K.compute_stiffness(cm_tf, dMdu, dMdt);
    return err;
  };
  int out_loop(double *out,
               int nne,
               int nsd,
               int Pno,
               int Vno)
  {
    int err = 0;
    Matrix<double> Kpu(Pno,nne*nsd);
    Kpu.trans(K.Kup);
    K.Kpt.trans(K.Ktp);

    Matrix<double> d_theta(Vno, 1), dP(Pno, 1);
    err += compute_d_theta_dP(d_theta.m_pdata, dP.m_pdata, du,
                              nne, nsd, Pno, Vno,
                              R.Ru.m_pdata, R.Rt.m_pdata, R.Rp.m_pdata,
                              Kpu.m_pdata, K.Ktu.m_pdata, K.Ktp.m_pdata, K.Ktt.m_pdata, K.Kpt.m_pdata);

    for(int ia=1; ia<=Pno; ia++)
      fv->tf.ddP(eid+1,ia) = dP(ia);


    for(int ia=1; ia<=Vno; ia++)
      fv->tf.ddV(eid+1,ia) = d_theta(ia);
    return err;
  };

};

template <class CM> class ConstitutiveModelIntregrate
{
  public:
    /// perform element integration in quasi steady state
    ///
    /// Updated Lagrangian and total Lagrangian based. When thermal
    /// is couled, temperature is assumed constant.
    ///
    /// \param[in] fe finite element helper object
    /// \param[out] lk computed element stiffness matrix
    /// \param[in] r_e nodal variabls(displacements) on the current element
    /// \param[in] grid a mesh object
    /// \param[in] mat a material object
    /// \param[in] fv object for field variables
    /// \param[in] sol object for solution scheme
    /// \param[in] load object for loading
    /// \param[in] opts structure PGFem3D option
    /// \param[in] mp mutiphysics object
    /// \param[in] mp_id mutiphysics id
    /// \param[in] dt time step size
    /// \return non-zero on internal error
    int integrate_ss(FEMLIB *fe,
                     double *out,
                     double *r_e,
                     Grid *grid,
                     MaterialProperty *mat,
                     FieldVariables *fv,
                     int sol_run_integration_algorithm,
                     LoadingSteps *load,
                     const PGFem3D_opt *opts,
                     const Multiphysics& mp,
                     int mp_id,
                     double dt)
    {
      int err = 0;
      double alpha = -1.0; // if alpha < 0, no inertia

      int total_Lagrangian = 1;
      if(opts->cm==UPDATED_LAGRANGIAN)
        total_Lagrangian = 0;

      int is_it_couple_w_thermal  = -1;
      int is_it_couple_w_chemical = -1;

      for(int ia=0; ia<fv->n_coupled; ia++)
      {
        if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
          is_it_couple_w_thermal = ia;
        if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
          is_it_couple_w_chemical = ia;
      }

      int eid = fe->curt_elem_id;
      int nsd = fe->nsd;
      int nne = fe->nne;
      int ndofn = fv->ndofn;
      int Pno   = fv->npres;
      int Vno   = fv->nVol;
      SUPP sup = load->sups[mp_id];

      Matrix<double> u(nne*nsd, 1), P(Pno, 1);
      Matrix<double> dMdu(DIM_3x3*nne*nsd,1);
      Matrix<double> dMdt(DIM_3x3*Vno,1);

      for(int a=0;a<nne;a++)
      {
        for(int b=0; b<nsd;b++)
          u.m_pdata[a*nsd+b] = r_e[a*ndofn+b];

        if(Pno==nne)
          P.m_pdata[a] = r_e[a*ndofn+nsd];
      }
      if(Pno==1)
      {
        P(1) = fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1);
      }

      DeformationGradient_ss *hF = NULL;

      Tensor<2> eSd = {};
      Tensor<4>  Ld = {};

      NodalTemerature *T = NULL;

      if(is_it_couple_w_thermal >= 0)
      {
        hF = new DeformationGradient_ss;
        T = new NodalTemerature;
        T->initialization(fe->nne);
        T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
      }

      if(is_it_couple_w_chemical >=0)
      {}

      DeformationGradient<DeformationGradient_ss> Fs(hF);

      CM cm_method(fe, Vno, Pno);
      Matrix<double> du;
      if(cm_method.run_for_update)
      {
        cm_method.fv = fv;

        int ndofe = nne*ndofn;
        Matrix<long> cn(ndofe,1);
        long *nod = fe->node_id.m_pdata;
        Matrix<double> dr_e(ndofe,1);
        get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,grid->node,cn.m_pdata,mp_id);
        def_elem(cn.m_pdata,ndofe,fv->dd_u,grid->element,grid->node,dr_e.m_pdata,sup,2);

        du.initialization(nne*nsd, 1, 0.0);
        for(int a=0;a<nne;a++)
        {
          for(int b=0; b<nsd;b++)
            du.m_pdata[a*nsd+b] = dr_e.m_pdata[a*ndofn+b];
        }
        cm_method.du = du.m_pdata;
      }

      Matrix<double> Nt(Vno,1), Np(Pno,1);

      for(int ip = 1; ip<=fe->nint; ip++)
      {
        fe->elem_basis_V(ip);
        fe->update_shape_tensor();

        fe->elem_shape_function(ip,Pno, Np.m_pdata);
        fe->elem_shape_function(ip,Vno, Nt.m_pdata);

        double theta_r = 0.0;
        double theta_n = 0.0;
        double Pnp1    = 0.0;

        for(int ia=1; ia<=Pno; ia++)
          Pnp1 += Np(ia)*P(ia);

        for(int ia=1; ia<=Vno; ia++)
        {
          theta_r += (fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia))*Nt(ia);
          theta_n += fv->tf.V_n(  eid+1, ia)*Nt(ia);
        }
        if(total_Lagrangian)
          theta_n = 1.0;

        Fs.initialization();
        if(is_it_couple_w_thermal>=0)
          Fs.hF->update_thermal_part(fe, grid, mat, fv, *T);
        Fs.update_total_deformation_gradient(fe, fv, sup, u, total_Lagrangian);

        Constitutive_model *m = &(fv->eps[eid].model[ip-1]);
        if(cm_method.run_integration_algorithm && sol_run_integration_algorithm){
          double tJ = det(Fs.F.np1);
          if(tJ<0)
            ++err;
          else{
            double tf_factor = pow(theta_r*theta_n/tJ, 1.0/3.0);
            err += m->run_integration_algorithm(Fs.F.np1.data,hF->n.data,hF->np1.data,
                                                dt,alpha,fe->x_ip.m_pdata, 0.0, is_it_couple_w_thermal, tf_factor);
          }
        }
        if(err>0)
          return err;

        err += Fs.compute_deformation_gradient_for_ss(fe, fv, total_Lagrangian,
                                                      is_it_couple_w_thermal,
                                                      is_it_couple_w_chemical);

        // get a shortened pointer for simplified CM function calls
        const Model_parameters *mp = m->param;

        void *ctx = NULL;
        if(is_it_couple_w_thermal>=0)
          err += construct_model_context_with_thermal(&ctx, m->param->type, Fs.F.np1.data,dt,alpha, NULL,
                                                      hF->n.data,hF->np1.data,-1);
        else
          err += construct_model_context(&ctx, m->param->type, Fs.F.np1.data,dt,alpha, NULL,-1);

        err += mp->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu.m_pdata);
        err += mp->compute_dMdt(m, ctx, fe->ST, Vno, dMdt.m_pdata);
        err += mp->destroy_ctx(&ctx);

        // <-- update plasticity part
        err += mp->update_elasticity_dev(m, Fs.eF.np1.data, Ld.data, eSd.data, -1, 0, dt, 1);

        double hJnp1 = 1.0;
        double pJnp1 = det(Fs.pF.np1);
        double eJn   = 1.0;
        double Jn    = 1.0;
        double MJ    = det(Fs.M);

        if(is_it_couple_w_thermal>=0)
          hJnp1 = det(hF->np1);

        if(!total_Lagrangian)
        {
          eJn = det(Fs.eF.n);
           Jn = det(Fs.F.n);
        }

        double tJn = Jn;
        double theta_e = theta_r*eJn*MJ;

        double dU  = mp->compute_dudj(  m, theta_e, -1, 0);
        double ddU = mp->compute_d2udj2(m, theta_e, -1, 0);
        // --> update elasticity part

        if(err!=0)
          break;

        Jn = Jn/pJnp1/hJnp1;

        CM_ThreeField cm_tf;
        cm_tf.set_femlib(fe,Vno,Pno,Nt.m_pdata,Np.m_pdata);

        cm_tf.set_tenosrs(Fs.Fr.data, Fs.eF.n.data, Fs.M.data, Fs.pF.np1.data, eSd.data, Ld.data);
        cm_tf.set_scalars(theta_r, theta_n, tJn, Jn, Pnp1, dU, ddU);

        err += cm_method.inner_loop(dMdu, dMdt, cm_tf);
      }

      if(is_it_couple_w_thermal >= 0){
        delete T;
        delete hF;
      }

      err += cm_method.out_loop(out, nne, nsd, Pno, Vno);

      return err;
    };
};

/// this is a wrapper function for the switch that was copy/pasted
/// everywhere. It is no big deal to keep adding to this private
/// function's argument list. Just put everything any model might need
/// and the switch will handle what is actually used.
int construct_model_context(void **ctx,
                            const int type,
                            double *F,
                            const double dt,
                            const double alpha,
                            double *eFnpa,
                            int npa)
{
  int err = 0;
  switch(type) {
  case TESTING:
  case HYPER_ELASTICITY:
    err += plasticity_model_none_ctx_build(ctx, F, eFnpa, NULL, NULL, 0);
    break;
  case CRYSTAL_PLASTICITY:
    err += plasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, NULL, NULL, 0);
    break;
  case BPA_PLASTICITY:
    err += plasticity_model_BPA_ctx_build(ctx, F, dt);
    break;
  case ISO_VISCOUS_DAMAGE:
    err += iso_viscous_damage_model_ctx_build(ctx, F, dt,alpha, eFnpa, NULL, NULL, 0,npa);
    break;
  case J2_PLASTICITY_DAMAGE:
    err += j2d_plasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, NULL, NULL, 0,npa);
    break;
  case POROVISCO_PLASTICITY:
    err += poro_viscoplasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, NULL, NULL, 0,npa);
    break;
  case ISO_VISCOUS_SPLIT_DAMAGE:
    err += iso_viscous_damage_model_split_ctx_build(ctx, F, dt,alpha, eFnpa, NULL, NULL, 0,npa);
    break;
  case MANUFACTURED_SOLUTIONS:
  {
    double x[3] = {};
    err += cm_mms_ctx_build(ctx, F, alpha, eFnpa, npa, x, 0);
    break;
  }
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n", type);
    err++;
    break;
  }
  assert (err == 0);
  return err;
}

/// constructor of constitutive model context
///
/// this is a wrapper function for the switch that was copy/pasted
/// everywhere for the coupled problem. Keep add function's argument list to this private.
/// In adding more constitutive model, just put everything any
/// and the switch will handle what is actually used.
///
/// \param[out] ctx  constructed context based on the model type
/// \param[in] type  constitutive model type
/// \param[in] F     total deformation gradient
/// \param[in] dt    time step size
/// \param[in] alpha mid point alpha
/// \param[in] eFnpa mid point elastic part of deformation gradient
/// \param[in] hFn   thermal part of deformation gradient at t(n)
/// \param[in] hFnp1 thermal part of deformation gradient at t(n+1)
/// \param[in] npa   mid point rule: if npa==0: v_npa = (1-alpha)*v(n-1) + alpha*v(n)
///                                  if npa==1: v_npa = (1-alpha)*v(n)   + alpha*v(n+1)
/// \return non-zeoro on internal error
int construct_model_context_with_thermal(void **ctx,
                                         const int type,
                                         double *F,
                                         const double dt,
                                         const double alpha,
                                         double *eFnpa,
                                         double *hFn,
                                         double *hFnp1,
                                         int npa)
{
  int err = 0;
  switch(type) {
  case TESTING:
  case HYPER_ELASTICITY:
    err += plasticity_model_none_ctx_build(ctx, F, eFnpa, hFn, hFnp1, 1);
    break;
  case CRYSTAL_PLASTICITY:
    err += plasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, hFn, hFnp1, 1);
    break;
  case ISO_VISCOUS_DAMAGE:
    err += iso_viscous_damage_model_ctx_build(ctx, F, dt,alpha, eFnpa, hFn, hFnp1, 1,npa);
    break;
  case J2_PLASTICITY_DAMAGE:
    err += j2d_plasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, hFn, hFnp1, 1,npa);
    break;
  case POROVISCO_PLASTICITY:
    err += poro_viscoplasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, hFn, hFnp1, 1,npa);
    break;
  case ISO_VISCOUS_SPLIT_DAMAGE:
    err += iso_viscous_damage_model_split_ctx_build(ctx, F, dt,alpha, eFnpa, hFn, hFnp1, 1,npa);
    break;
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n", type);
    err++;
    break;
  }
  assert (err == 0);
  return err;
}

/// initialize model parameter
///
/// \param[in, out] a Model_parameters, p[model_id] will be initialized
/// \return non-zeoro on internal error
int
Constitutive_model::initialization(const Model_parameters *p)
{
  int err = 0;
  if (p == NULL)
    err++;
  else
  {
    this->param = p;
    Model_var_info info;
    this->param->get_var_info(info);

    err += this->vars_list[0][this->model_id].initialization(info.n_Fs,
                                                             info.n_vars,
                                                             info.n_flags);
    this->param->set_init_vals(this);
  }
  return err;
}

/// Construct a Model_parameters. Model type will be determined at runtime.
///
/// \param[out] **p        a Model_parameters, p[model_id] will be constructed
/// \param[in]  model_id   index of p
/// \param[in]  model_type constitutive model type
/// \return non-zero on error.
int construct_Model_parameters(Model_parameters **p, int model_id, int model_type)
{
  int err = 0;

  if(p[model_id] != NULL)
  {
    delete p[model_id];
    p[model_id] = NULL;
  }

  switch(model_type)
  {
    case TESTING:
    case HYPER_ELASTICITY:
      p[model_id] = new HE_PARAM;
      break;
    case CRYSTAL_PLASTICITY:
      p[model_id] = new CP_PARAM;
      break;
    case BPA_PLASTICITY:
      p[model_id] = new BPA_PARAM;
      break;
    case ISO_VISCOUS_DAMAGE:
      p[model_id] = new CM_IVD_PARAM;
      break;
    case J2_PLASTICITY_DAMAGE:
      p[model_id] = new CM_J2P_PARAM;
      break;
    case POROVISCO_PLASTICITY:
      p[model_id] = new CM_PVP_PARAM;
      break;
    case ISO_VISCOUS_SPLIT_DAMAGE:
      p[model_id] = new CM_IVDS_PARAM;
      break;
  case MANUFACTURED_SOLUTIONS:
      p[model_id] = new CM_MMS_PARAM;
    break;      
    default:
      PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",model_type);
      err++;
      return err;
  }

  if(p[model_id] != NULL)
    p[model_id]->set_nulls();

  return err;
}
/// Initialize the Model_parameters object. The object may be used
/// after calling this function. Calling this function on an already
/// initialized object is undefined.
///
/// \param[in] p_hmat material property object
/// \param[in] type   constitutive model type
/// \return non-zero on error.
int
Model_parameters::initialization(const HOMMAT *p_hmat, const size_t type)
{
  int err = 0;

  switch(type)
  {
    case TESTING:
    case HYPER_ELASTICITY:
    case CRYSTAL_PLASTICITY:
    case BPA_PLASTICITY:
    case ISO_VISCOUS_DAMAGE:
    case J2_PLASTICITY_DAMAGE:
    case POROVISCO_PLASTICITY:
    case ISO_VISCOUS_SPLIT_DAMAGE:
    case MANUFACTURED_SOLUTIONS:
      break; // no action
    default:
      PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",type);
      err++;
      return err;
  }

  this->p_hmat = p_hmat;
  this->type = type;

  auto* cm_mat = PGFEM_malloc<MATERIAL_CONSTITUTIVE_MODEL>();
  auto* mat_e = PGFEM_malloc<MATERIAL_ELASTICITY>();
  auto* elast = PGFEM_malloc<ELASTICITY>();

  set_properties_using_E_and_nu(mat_e,p_hmat->E,p_hmat->nu);
  mat_e->m01 = p_hmat->m01;
  mat_e->m10 = p_hmat->m10;
  mat_e->G   = p_hmat->G;
  mat_e->kappa = p_hmat->E/(3.0*(1.0-2.0*p_hmat->nu));
  mat_e->devPotFlag = p_hmat->devPotFlag;
  mat_e->volPotFlag = p_hmat->volPotFlag;

  set_properties_constitutive_model(cm_mat,mat_e,NULL);
  construct_elasticity(elast, mat_e, 1);

  this->cm_mat   = cm_mat;
  this->cm_elast = elast;
  err += this->model_dependent_initialization();

  return err;
}

/// destory members in a Model_parameter
/// \return non-zero on error.
int
Model_parameters::finalization()
{
  /* drop pointer to material (material free'd elsewhere) */
  this->p_hmat = NULL;
  delete (this->cm_mat)->mat_e;
  delete this->cm_mat;
  destruct_elasticity(this->cm_elast);
  delete this->cm_elast;

  /* reset counters/flags */
  this->type = -1;

  /* free model constants */
  this->n_param = -1;
  delete this->model_param;
  this->model_param = NULL;
  this->n_param_index = -1;
  if(this->model_param_index !=NULL)
    delete this->model_param_index;
  this->model_param_index = NULL;
  return 0;
}

/// User defined function that returns the size of the data to be
/// packed/unpacked.
/// Does not modify the CM object or any of the data it holds.
///
/// \param[in] m, CM object with internal data set from the buffer
/// \return size in bytes of the pack/unpack data
int
Constitutive_model::get_size()
{
  return this->vars_list[0][this->model_id].state_variables_get_packed_size();
}

/// User defined function to pack the CM data into a buffer (see pack_data).
/// Does not modify the CM object or any of the data it holds.
///
/// \param[in,out] buffer, a buffer to insert data to
///
/// \param[in,out] pos,    insert position in the buffer. Upon exit - next
///                        insertion position.
/// \return non-zero on error.
int
Constitutive_model::pack(char *buffer, size_t *pos)
{
  return this->vars_list[0][this->model_id].state_variables_pack(buffer, pos);
}

/// User defined function to unpack CM data from a buffer (see also
/// usr_pack, unpack_data).
///
/// \param[in]     buffer, the buffer to read data from
/// \param[in,out] pos,    the position in buffer to begin reading from.
///                        Upon exit - position for next read.
/// \return        non-zero on error.
int
Constitutive_model::unpack(const char *buffer, size_t *pos)
{
  return this->vars_list[0][this->model_id].state_variables_unpack(buffer, pos);
}

/// This function is running integration algorithm such that it modifies
/// the internal state to contain the updated values upon exit, i.e.,
/// subsequent calls to get_F functions will return the correct
/// values without re-integrating the constitutive model.
/// \param[in]
/// \param[in] *Fnp1  deformation gradient at t(n+1)
/// \param[in] *hFn   thermal part of deformation gradient at t(n)
/// \param[in] *hFnp1 thermal part of deformation gradient at t(n+1)
/// \param[in] dt                     time step size
/// \param[in] alpha                  mid point rule alpha
/// \param[in] is_it_couple_w_thermal checking coupling with thermal
///                                   if > 0: apply thermal expansitions. default = 0
/// \param[in] tf_factor              (theta/J)^(1/3) used for computing true Fnp1. default = 1.0
/// \return non-zero on internal error that should be handled by the calling function.
int
Constitutive_model::run_integration_algorithm(double *tFnp1_in,
                                              double *hFn,
                                              double *hFnp1,
                                              double dt,
                                              double alpha,
                                              const double *x,
                                              const double t,                                              
                                              int is_it_couple_w_thermal,
                                              double tf_factor)
{
  int err = 0;
  void *ctx = NULL;
  TensorA<2> tFnp1(tFnp1_in);
  Tensor<2> Fnp1 = tf_factor*tFnp1(i,j);
  if(is_it_couple_w_thermal>=0)
    err += construct_model_context_with_thermal(&ctx, param->type, Fnp1.data,dt,alpha, NULL,
                                                 hFn,hFnp1, -1);
  else{
    if(param->type == MANUFACTURED_SOLUTIONS)
      err += cm_mms_ctx_build(&ctx, Fnp1.data, alpha, NULL, -1, x, t);
    else
      err += construct_model_context(&ctx, param->type, Fnp1.data,dt,alpha, NULL,-1);
  }
  err += param->integration_algorithm(this,ctx); // perform integration algorithm
  err += param->destroy_ctx(&ctx);
  return err;
  
}

/// compute PK2 and elasticity tensor
///
/// \param[in]  m                 constitutive model object
/// \param[in]  eF                elastic part of deformation gradient
/// \param[out] L                 4th order elasticity tensor
/// \param[out] S                 computed PK2 tensor
/// \param[in]  compute_stiffness if 0, no compute L (elasticity tensor)
/// \return     non-zero on error.
int constitutive_model_default_update_elasticity(const Constitutive_model *m,
                                                 const double *eF,
                                                 double *L,
                                                 double *S,
                                                 const int compute_stiffness)
{
  int err = 0;

  ELASTICITY *elast = (m->param)->cm_elast; // get elasticity handle
  double *tempS = elast->S; // temporal pointer to update *L, and *S using elast
  double *tempL = elast->L;
  elast->S = S;

  double F[DIM_3x3];
  memcpy(F,eF,sizeof(double)*DIM_3x3);
  if(compute_stiffness)
    elast->L = L;
  else
    elast->L = NULL;

  elast->update_elasticity(elast,F,compute_stiffness);

  elast->S = tempS;
  elast->L = tempL;
  return err;
}

static int compare_mat_id(const void *a, const void *b)
{
  return ((const HOMMAT *)a)->mat_id - ((const HOMMAT *)b)->mat_id;
}

int read_model_parameters_list(const int n_mat,
                               HOMMAT *hmat_list,
                               FILE *in)
{
  /* see issue #28 */
  /* File format:
     ------
     num_entries
     material_id model_type
     { # begin model info
     ... # specified by model
     } # end model info
     ...
     ------
     caveats:
     - Comments can only go at the end of a line after all of the data
       on that line has been specified
     - No support for sub-braces in model sections (yet)
     - Undefined behavior for duplicate entires (will try to
       re-initialize the object)
   */
  int err = 0;

  if (n_mat <= 0) return 1;
  int *is_set = (int *) calloc(n_mat, sizeof(int));

  int num_entries = -1;
  HOMMAT *key = (HOMMAT *) calloc(1, sizeof(*key));
  err += scan_for_valid_line(in);
  CHECK_SCANF(in, "%d", &num_entries);

  int i = 0;
  for (i = 0; i < num_entries; i++) {
    int model_type = -1;
    int model_type_uq = 0;

    HOMMAT *p_hmat = NULL;
    err += scan_for_valid_line(in);
    if (feof(in)) break;

    CHECK_SCANF(in, "%d %d", &(key->mat_id), &model_type);
    err += scan_for_valid_line(in);

    int brace = fgetc(in);
    assert(brace == '{' && "Expect opening brace as next valid entry");

    /*
     * NOTE: The material ID in the input files is not necessarily the
     * index of the HOMMAT material. The hmat_list is the reduced set
     * of material properties that are actually used on the
     * domain. Therefore, we need to search for the matching
     * mat_id. Futheremore, no warning is issued if no match is found,
     * but we perform checks to ensire that all materials are
     * specified.
     */

    if(key->mat_id>=0)
    {
      /* search for matching pointer in hmat_list (assume unique) */
      p_hmat = static_cast<HOMMAT*>(bsearch(key, hmat_list, n_mat,
                                            sizeof(*hmat_list), compare_mat_id));

      /* check for match */
      if (p_hmat != NULL)
      {
        int idx = p_hmat - hmat_list;
        if(is_set[idx])
        {
          if(-1 == model_type)
          {
            if(hmat_list[idx].param->pF==NULL)
            {
              hmat_list[idx].param->pF  = (double *) malloc(sizeof(double)*DIM_3x3);
              hmat_list[idx].param->pFI = (double *) malloc(sizeof(double)*DIM_3x3);
            }
            double *pF  = hmat_list[idx].param->pF;
            double *pFI = hmat_list[idx].param->pFI; 
            err += scan_for_valid_line(in);
            int match = fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", pF+0, pF+1, pF+2,
                                                                          pF+3, pF+4, pF+5,
                                                                          pF+6, pF+7, pF+8);
            TensorA<2> pF0(pF), pF0I(pFI);
            inv(pF0, pF0I);

            if(match != DIM_3x3)
            {
              ++err;
              assert(match == DIM_3x3 && "Did not read expected number of parameters");
            }

            err += scan_for_valid_line(in);
            // not expecting EOF, check and return error if encountered
            if (feof(in)) err ++;
              assert(!feof(in) && "EOF reached prematurely");
          }
        }
        else{

          is_set[idx] = 1;
          
          // construct and initialize this object
          if(model_type==CM_UQCM)
          {
            model_type_uq = 1;
            err += scan_for_valid_line(in);
            CHECK_SCANF(in, "%d", &model_type);
          }

          err += construct_Model_parameters(&(hmat_list[idx].param), 0, model_type);
          hmat_list[idx].param->uqcm = model_type_uq;

          err += hmat_list[idx].param->initialization(hmat_list + idx, model_type);
          err += hmat_list[idx].param->read_param(in);
          hmat_list[idx].param->mat_id = key->mat_id;
        }
      }
    }

    /* scan to closing brace and continue on to the next entry */
    while(fgetc(in) != '}' && !feof(in)){}
    if (feof(in)) break;
  }

  if (feof(in) && i != num_entries) {
    err++;
    assert(0 && "Prematurely reached EOF");
  }

  int sum = 0;
  for (int i = 0; i < n_mat; i++){
    sum += is_set[i];
  }
  if (sum != n_mat) err++;
  assert(sum == n_mat && "require that all model params are set");

  free(key);
  free(is_set);
  return err;
}

int init_all_constitutive_model(EPS *eps,
                                const int ne,
                                const Element *elem,
                                const int n_mat,
                                const HOMMAT *hmat_list,
				int myrank)
{
  int err = 0;
  if (ne <= 0) return 1;

  for(int i = 0; i < ne; i++) {
    /* aliases */
    EPS *p_eps = eps + i;
    const Element *p_el = elem + i;
    const Model_parameters *p_param = hmat_list[p_el->mat[2]].param;

    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++)
      p_eps->model[j].initialization(p_param);
  }

  plasticity_model_set_orientations(eps, ne, elem, n_mat, hmat_list, myrank); // nothing will happen if there is no use of the crystal plasticity model
  return err;
}

int constitutive_model_reset_state(EPS *eps,
                                   const int ne,
                                   const Element *elem)
{
  int err = 0;
  if (ne <= 0) return 1;

  for (int i = 0; i < ne; i++) {
    const Element *p_el = elem + i;
    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++) {
      Constitutive_model *m = &(eps[i].model[j]);
      m->param->reset_state_vars(m);
    }
  }
  return err;
}

/// save state variables
///
/// save state variables(t(n-1) and t(n)) to temporal variables(t(n-1) and t(n))
/// in order to use when solution step is failed and requires go to initial step
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_save_state_vars_to_temporal(FieldVariables *fv,
                                                   Grid *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;
  const Element *elem = grid->element;
  for(int eid=0; eid<grid->ne; eid++)
  {
    long nint = 1;
    int_point(elem[eid].toe,&nint);
    for (int ip = 0; ip < nint; ip++)
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->save_state_vars_to_temporal(m, var+m->model_id);
    }
  }
  return err;
}

/// update state variables(t(n+1))
///
/// temporary save state variables(t(n+1))
/// when coupled physics calls dependent physics,
/// original soultions at t(n-1) and t(n) and updated solutions t(n+1) are needed.
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_update_np1_state_vars_to_temporal(FieldVariables *fv,
                                                         Grid *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;

  const Element *elem = grid->element;

  for(int eid=0; eid<grid->ne; eid++)
  {
    long nint = 1;
    int_point(elem[eid].toe,&nint);
    for (int ip = 0; ip < nint; ip++)
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->update_np1_state_vars_to_temporal(m, var+m->model_id);
    }
  }
  return err;
}

/// reset state variables using priori stored values
///
/// reset state variables(t(n-1) and t(n)) using temporal variables(t(n-1) and t(n))
/// when solution step is failed and requires go to initial step
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_reset_state_using_temporal(FieldVariables *fv,
                                                  Grid *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;

  const Element *elem = grid->element;

  for(int eid=0; eid<grid->ne; eid++)
  {
    long nint = 1;
    int_point(elem[eid].toe,&nint);
    for (int ip = 0; ip < nint; ip++)
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->reset_state_vars_using_temporal(m, var+m->model_id);
    }
  }
  return err;
}

int constitutive_model_update_time_steps(const Element *elem,
                                         Node *node,
                                         EPS *eps,
                                         const int ne,
                                         const int nn,
                                         const int ndofn,
                                         const double* r,
                                         const double dt,
                                         const int total_Lagrangian,
                                         const int mp_id)
{
  int nsd = 3;
  int err = 0;
  if (ne <= 0) return 1;

  for (int i = 0; i < ne; i++)
  {
    const Element *p_el = elem + i;
    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++)
    {
      Constitutive_model *m = &(eps[i].model[j]);
      m->param->update_state_vars(m);
    }
  }

  /*********************/
  /* Coordinate update */
  /*********************/
  if(total_Lagrangian) {
    for(int n = 0;n<nn; n++) {
      for(int a=0;a<nsd;a++) {
        int II = node[n].id_map[mp_id].id[a];
        if (II != 0){
          if (a == 0)      node[n].x1 = node[n].x1_fd + r[n*ndofn + a];
          else if (a == 1) node[n].x2 = node[n].x2_fd + r[n*ndofn + a];
          else if (a == 2) node[n].x3 = node[n].x3_fd + r[n*ndofn + a];
        }
      }
    }/* end n < nn */
  } else {
    for(int n = 0;n<nn; n++) {
      for(int a=0;a<nsd;a++) {
        int II = node[n].id_map[mp_id].id[a];
        if (II != 0){
          if (a == 0)      node[n].x1 += r[n*ndofn + a];
          else if (a == 1) node[n].x2 += r[n*ndofn + a];
          else if (a == 2) node[n].x3 += r[n*ndofn + a];
        }
      }
    }/* end n < nn */
  }

  return err;
}

int constitutive_model_test(const HOMMAT *hmat, double *L_in, int Print_results)
{
  int err = 0;
  return err;
}

/// Common part of computing stiffness Matrix
///
/// \param[out] lk computed element stiffness matrix
/// \param[in] fe finite element helper object
/// \param[in] Fr 2nd order tensor Fr
/// \param[in] eFnMT 2nd order tensor (eFn*M)'
/// \param[in] eFn 2nd order tensor eFn
/// \param[in] M 2nd order tensor M
/// \param[in] FrTFr 2nd order tensor Fr'*Fr
/// \param[in] eFnM 2nd order tensor eFn*M
/// \param[in] S 2nd order tensor S
/// \param[in] L 4th order elasticity tensor L
/// \param[in] dMdu_all list of 2nd order dMdu tensors
/// \param[in] Jn det(Fn)
/// \return non-zero on internal error
int compute_stiffness_matrix(double *lk,
                             const FEMLIB *fe,
                             double *Fr_in,
                             double *eFnMT_in,
                             double *eFn_in,
                             double *M_in,
                             double *FrTFr_in,
                             double *eFnM_in,
                             double *S_in,
                             double *L_in,
                             double *dMdu_all,
                             const double Jn)
{
  int err = 0;
  const int nne = fe->nne;
  const int nsd = fe->nsd;

  TensorA<2> Fr(Fr_in);
  TensorA<2> eFnMT(eFnMT_in);
  TensorA<2> eFn(eFn_in);
  TensorA<2> M(M_in);
  TensorA<2> FrTFr(FrTFr_in);
  TensorA<2> eFnM(eFnM_in);
  TensorA<2> S(S_in);
  TensorA<4> L(L_in);
  Tensor<2> MI;
  
  inv(M,MI);

  Tensor<2> MTeFnT_sAA_eFn, MTeFnT_sAA_eFnM,
            sBB,sCC, dCdu;

  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const int id_ab = idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd);
      TensorA<2> ST_ab((fe->ST)+id_ab);

      Tensor<2> AA =  Fr(k,i)*ST_ab(k,j);
      Tensor<2> sAA = 0.5*(AA(i,j) + AA(j,i));

      MTeFnT_sAA_eFn(i,j) = eFnMT(i,k) * sAA(k,l) * eFn(l,j);
      MTeFnT_sAA_eFnM(i,j) = MTeFnT_sAA_eFn(i,k) * M(k,j);

      for(int w=0; w<nne; w++)
      {
        for(int g=0; g<nsd; g++)
        {
          const int id_wg = idx_4_gen(w,g,0,0,nne,nsd,nsd,nsd);
          TensorA<2> ST_wg((fe->ST)+id_wg);

          TensorA<2> dMdu(dMdu_all + id_wg);

          Tensor<2> BB =  Fr(k,i)*ST_wg(k,j);
          Tensor<2> sBB = 0.5 * (BB(i,j) + BB(j,i));

          Tensor<2> CC =  ST_ab(k,i) * ST_wg(k,j);
          sCC = 0.5 * (CC(i,j) + CC(j,i));

          // compute dCdu
          Tensor<2> MTeFnT_FrTFreFndMdu = eFnMT(i,k)*FrTFr(k,l)*eFn(l,o)*dMdu(o,j);
          dCdu(i,j) = 0.5 * (MTeFnT_FrTFreFndMdu(i,j)
                                    + MTeFnT_FrTFreFndMdu(j,i)) + eFnMT(i,k)*sBB(k,l)*eFnM(l,j);


          // compute MTeFnT_sAA_eFnM:L:dCdu
          Tensor<2> L_dCdu = L(i,j,k,l) * dCdu(k,l);
          double MTeFnT_sAA_eFnM_L_dCdu = MTeFnT_sAA_eFnM(i,j) * L_dCdu(i,j);

          // compute MTeFnT_sCC_eFnM
          Tensor<2> MTeFnT_sCC_eFnM = eFnMT(i,k) * sCC(k,l) * eFnM(l,j);

          // compute MTeFnT_sCC_eFnM:S
          double MTeFnT_sCC_eFnM_S = MTeFnT_sCC_eFnM(i,j) * S(i,j);

          // compute MTeFnT_sAA_eFndMdu
          Tensor<2> MTeFnT_sAA_eFndMdu = MTeFnT_sAA_eFn(i,k) * dMdu(k,j);
          // compute MTeFnT_sAA_eFndMdu:S
          double sMTeFnT_sAA_eFndMdu_S = 0.5*(MTeFnT_sAA_eFndMdu(i,j)
                                                       + MTeFnT_sAA_eFndMdu(j,i))* S(i,j);
                                                       
          double MTeFnT_sAA_eFnMS = MTeFnT_sAA_eFnM(i,j)*S(i,j);
          double DpJ = -MI(j,i)*dMdu(i,j);
          const int lk_idx = idx_K(a,b,w,g,nne,nsd);

          lk[lk_idx] += 1.0/Jn*fe->detJxW*(MTeFnT_sAA_eFnM_L_dCdu 
                                           + 2.0*sMTeFnT_sAA_eFndMdu_S
                                           + MTeFnT_sCC_eFnM_S
                                           + DpJ*MTeFnT_sAA_eFnMS);
        }
      }
    }
  }

  return err;
}

/// Common part of computing residual vector
///
/// \param[out] f computed residual vector
/// \param[in] fe finite element helper object
/// \param[in] Fr 2nd order tensor Fr
/// \param[in] eFnMT 2nd order tensor (eFn*M)'
/// \param[in] eFnM 2nd order tensor eFn*M
/// \param[in] S 2nd order tensor S
/// \param[in] Jn det(Fn)
/// \return non-zero on internal error
int compute_residual_vector(double *f,
                            const FEMLIB *fe,
                            double *Fr_in,
                            double *eFnMT_in,
                            double *eFnM_in,
                            double *S_in,
                            const double Jn)
{
  int err = 0;
  const int nne = fe->nne;
  const int nsd = fe->nsd;

  TensorA<2> Fr(Fr_in);
  TensorA<2> eFnMT(eFnMT_in);
  TensorA<2> eFnM(eFnM_in);
  TensorA<2> S(S_in);

  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const int id_ab = idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd);
      TensorA<2> ST_ab((fe->ST)+id_ab);

      Tensor<2> AA = Fr(k,i)*ST_ab(k,j);
      Tensor<2> sAA = 0.5*(AA(i,j)+AA(j,i));
      double MTeFnT_sAA_eFnM_S = eFnMT(i,k)*sAA(k,l)*eFnM(l,j)*S(i,j);

      int fe_id = a*nsd + b;
      f[fe_id] += 1.0/Jn*fe->detJxW*MTeFnT_sAA_eFnM_S;
    }
  }

  return err;
}

/// compute ouput variables e.g. effective stress and strain
///
/// Visit each element and compute output variables according to the element model type.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] load object for loading
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] alpha mid point rule alpha
/// \return non-zero on internal error
int constitutive_model_update_output_variables(Grid *grid,
                                               MaterialProperty *mat,
                                               FieldVariables *FV,
                                               LoadingSteps *load,
                                               PGFem3D_opt *opts,
                                               const Multiphysics& mp,
                                               int mp_id,
                                               const double dt,
                                               double alpha)
{
  int err = 0;

  int total_Lagrangian = 1;
  if(opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  FieldVariables *fv = FV+mp_id;
  SIG *sig = fv->sig;
  EPS *eps = fv->eps;
  Node *node = grid->node;
  Element *elem = grid->element;

  int idx[6];
  idx[0] = idx_2(0,0); //XX
  idx[1] = idx_2(1,1); //YY
  idx[2] = idx_2(2,2); //ZZ
  idx[3] = idx_2(1,2); //YZ
  idx[4] = idx_2(0,2); //XZ
  idx[5] = idx_2(0,1); //XY

  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  /* deformation gradient */
  Tensor<2> F, FI, eF, pF, pFI, eSd, eS, S, S_bar, eC, eCI, sigma, E, b, bI;

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;


  for (int eid = 0; eid < grid->ne; eid++)
  { 
    FEMLIB fe;   
    Constitutive_model *m = (fv->eps[eid]).model;
    fe.initialization(eid,elem,node,0,total_Lagrangian,m->param->pFI);
        
    int nne   = fe.nne;
    int ndofn = fv->ndofn;

    memset(sig[eid].el.o,0,6*sizeof(double));
    memset(eps[eid].el.o,0,6*sizeof(double));

    NodalTemerature *T = NULL;
    long *nod = fe.node_id.m_pdata;

    if(is_it_couple_w_thermal >= 0)
    {
      T = new NodalTemerature;
      T->initialization(fe.nne);
      T->get_temperature(&fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
    }

    if(is_it_couple_w_chemical >= 0){}

    double V = 0.0;

    for(int ip=1; ip<=fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      V += fe.detJxW;

      Constitutive_model *m = &(eps[eid].model[ip-1]);
      const Model_parameters *func = m->param;
      err += func->get_F(m,  F.data,1);

      double theta = 0.0;
      double J = det(F);
      if(opts->analysis_type==CM3F)
      {
        int Vno   = fv->nVol;
        Matrix<double> Nt(Vno, 1, 0.0);
        fe.elem_shape_function(ip,Vno, Nt.m_pdata);
        for(int ia=0; ia<Vno; ia++)
        {
          double Via = 0;
          if(Vno==nne){
            long id = nod[ia]*ndofn + ia;
            Via = fv->u_n[id];
          }
          if(Vno==1)
            Via = fv->tf.V_n(eid+1, 1);

          theta += Nt(ia+1)*Via;
        }
        double factor = pow(theta/J, 1.0/3.0);
        for(int ia=0; ia<DIM_3x3; ia++)
          F.data[ia] *= factor;
      }
      else
        theta = J;
        
      inv(F, FI);        

      err += func->get_pF(m,pF.data,1);
      inv(pF,pFI);

      double hJ  = 1.0;
      double pJ  = det(pF);

      if(is_it_couple_w_thermal>=0)
      {
        // compute temperature at the integration point
        double hFnm1[9],hFn[9];
        Tensor<2> hFnp1;
        err += compute_temperature_at_ip(&fe,grid,mat,T->T0,
                                         T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                         hFnp1.data,hFn,hFnm1);
        hJ = det(hFnp1);
        Tensor<2> hFI;
        inv(hFnp1, hFI);
        eF(i,j) = F(i,k)*hFI(k,l)*pFI(l,j);
      }
      else
        eF(i,j) = F(i,k)*pFI(k, j);

      double theta_e = theta/pJ/hJ;
      double Up = 0.0;

      S_bar = {};
      if((m->param)->uqcm)
      {
        double *x_ip = (fe.x_ip).m_pdata;

        ELASTICITY *elast = (m->param)->cm_elast;
        mat_e_in = elast->mat;
        err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
        elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation

        err += func->update_elasticity_dev(m, eF.data, NULL, S_bar.data, -1, 0, dt, 0);
        Up  = func->compute_dudj(m, theta_e, -1, 0);

        elast->mat = mat_e_in;
      }
      else
      {
        err += func->update_elasticity_dev(m, eF.data, NULL, S_bar.data, -1, 0, dt, 0);
        Up  = func->compute_dudj(m, theta_e, -1, 0);
      }
      // <-- update elasticity part
      eC = eF(k,i)*eF(k,j);
      inv(eC, eCI);
      eS = S_bar(i,j) + Up*theta_e*eCI(i,j);

      // Compute Cauchy stress (theta)^-1 F S F'
      sigma(i,j) = 1.0/theta_e*eF(i,k)*eS(k,l)*eF(j,l);
      
      // compute PKII
      S(i,j) = theta*FI(i,k)*sigma(k,l)*FI(j,l);

      // Elastic Green Lagrange strain
      E = 0.5*(F(k,i)*F(k,j) - delta_ij(i,j));

      // Compute the logarithmic strain e = 1/2(I - inv(FF'))
      b = F(i,k)*F(j,k);
      inv(b, bI);
      Tensor<2> e = 0.5*(delta_ij(i,j) - bI(i,j));

      for(int ia=0; ia<6; ia++)
      {
        int id = idx[ia];
        sig[eid].il[ip-1].o[ia] = S.data[id];
        eps[eid].il[ip-1].o[ia] = E.data[id];

        sig[eid].el.o[ia] += fe.detJxW*sigma.data[id];

        // engineering strain
        if(ia<3)
          eps[eid].el.o[ia] += fe.detJxW*e.data[id];
        else
          eps[eid].el.o[ia] += 2.0*fe.detJxW*e.data[id];
      }

      /* store total deformation */
      memcpy(eps[eid].il[ip-1].F, F.data, DIM_3x3 * sizeof(double));

      /* store the hardening parameter */
      err += func->get_hardening(m, &eps[eid].dam[ip-1].wn,2);

      /* compute/store the plastic strain variable */
      err += func->get_plast_strain_var(m, &eps[eid].dam[ip-1].Xn);
    }
    for(int ia=0; ia<6; ia++)
    {
      sig[eid].el.o[ia] = sig[eid].el.o[ia]/V;
      eps[eid].el.o[ia] = eps[eid].el.o[ia]/V;
    }
    if(is_it_couple_w_thermal >=0)
      delete T;
  }

  return err;
}

/// compute element residual for mid point rule
///
/// compute residual(n+alpha)
///
/// \param[out] f computed element residual
/// \param[in] m constitutive model object
/// \param[in] ii element id
/// \param[in] ndofn number of degree freedom on node
/// \param[in] pFnp1 2nd order tenosr pF(n+1)
/// \param[in] pFn   2nd order tenosr pF(n)
/// \param[in] Fnp1  2nd order tenosr F(n+1)
/// \param[in] Fn    2nd order tenosr F(n)
/// \param[in] hFnp1 2nd order tenosr hF(n+1)
/// \param[in] hFn   2nd order tenosr hF(n)
/// \param[in] is_it_couple_w_chemical flag for coupling with thermal
/// \param[in] alpha mid point rule alpha
/// \param[in] dt_alpha_1_minus_alpha -dt*(1-alpha) for mid point btw t(n) t(n+1)
///                                   -dt*alpha     for mid point btw t(n-1) t(n)
/// \param[in] fe finite element helper object
/// \param[in] npa   mid point rule: if npa==0: v_npa = (1-alpha)*v(n-1) + alpha*v(n)
///                                  if npa==1: v_npa = (1-alpha)*v(n)   + alpha*v(n+1)
/// \return non-zero on internal error
int residuals_el_constitutive_model_n_plus_alpha(double *f,
                                                 const Constitutive_model *m,
                                                 const int ii,
                                                 const int ndofn,
                                                 double *pFnp1_in,
                                                 double *pFn_in,
                                                 double *Fnp1_in,
                                                 double *Fn_in,
                                                 double *hFnp1_in,
                                                 double *hFn_in,
                                                 const int is_it_couple_w_thermal,
                                                 const double alpha,
                                                 const double dt_alpha_1_minus_alpha,
                                                 FEMLIB *fe,
                                                 const int npa)
{
  // Total Lagrangian based
  int err = 0;
  const int nsd = fe->nsd;

  Tensor<2> M = {},eFnpa = {},pFnpa,pFnpa_I,hFnpa,Fnpa,S = {},MT;
  TensorA<2> pFnp1(pFnp1_in),pFn(pFn_in),
                        Fnp1(Fnp1_in),Fn(Fn_in),
                        hFnp1(hFnp1_in),hFn(hFn_in);


  int compute_stiffness = 0;

  mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, nsd*nsd);
  hFnpa = hFnp1(i,j);
//  mid_point_rule(hFnpa.data, hFn.data, hFnp1.data, alpha, nsd*nsd);
  mid_point_rule( Fnpa.data,  Fn.data,  Fnp1.data, alpha, nsd*nsd);

  if(is_it_couple_w_thermal>=0)
  {
    Tensor<2> hFnpa_I;

    err += inv(hFnpa, hFnpa_I);
    err += inv(pFnpa, pFnpa_I);
    M = hFnpa_I(i,k)*pFnpa_I(k,j);
  }
  else
    err += inv(pFnpa, M);

  eFnpa = Fnpa(i,k)*M(k,j);

  MT = M(i,j).to(j,i);

  {
    // check that deformation is invertible -> J > 0
    int terr = 0;
    getJacobian(Fnpa.data, ii, &terr);
    err += terr;
  }

  void *ctx;

  err += construct_model_context(&ctx, m->param->type, Fnp1.data,0.0,alpha,eFnpa.data,npa);
  err += (m->param)->update_elasticity(m,ctx,NULL,S.data,compute_stiffness);
  err += m->param->destroy_ctx(&ctx);

  if(err==0)
  {
    double Jn = 1.0;
    err += compute_residual_vector(f,fe,Fnpa.data,MT.data,M.data,S.data,Jn/dt_alpha_1_minus_alpha);
  }

  return err;
}

int cm_get_subdivision_parameter(double *subdiv_param,
                                 const int ne,
                                 const Element *elem,
                                 const EPS *eps,
                                 const double dt)
{
  int err = 0;
  *subdiv_param = 0.0;
  double cur_val = 0.0;
  double max_val = 0.0;
  for (int i = 0; i < ne; i++) {
    long n_ip = 0;
    int_point(elem[i].toe, &n_ip);
    for (int ip = 0; ip < n_ip; ip++) {
      err += eps[i].model[ip].param->get_subdiv_param(&(eps[i].model[ip]), &cur_val, dt);
      max_val = std::max(max_val, cur_val);
    }
  }

  *subdiv_param = max_val;
  return err;
}

/// compute element stiffness matrix in transient for single field formulation
///
/// Compute element stiffness matrix based on the mid point rule. When thermal
/// is couled, temperature is assumed constant. Currenlty total Lagrangian is
/// active (Updated lagrangian is implemented, but accelleration term is not
/// fully implemented for updated lagrangian)
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_el_constitutive_model_w_inertia_1f(FEMLIB *fe,
                                                 double *lk,
                                                 double *r_e,
                                                 Grid *grid,
                                                 MaterialProperty *mat,
                                                 FieldVariables *fv,
                                                 Solver *sol,
                                                 LoadingSteps *load,
                                                 CRPL *crpl,
                                                 const PGFem3D_opt *opts,
                                                 const Multiphysics& mp,
                                                 int mp_id,
                                                 double dt)
{
  int err = 0;
  double alpha = sol->alpha;
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

// when updated Lagrangian is used should be un-commented
//  if(opts->cm != 0)
//    total_Lagrangian = 1;

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;

  double *u = (double *) malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = (double *) malloc(sizeof(*dMdu_all)*DIM_3x3*nne*nsd);
  memset(dMdu_all,0,DIM_3x3*nne*nsd*sizeof(double));

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2> Fn,Fr,Fnp1,pFn,pFnp1,S = {},eFnpa = {},pFnpa,pFnpa_I,
            eFn,M = {},eFnM,eFnMT,FrTFr = {},hFnp1,hFn,hFnpa,hFnpa_I;

  Tensor<4,3,double> L;

  NodalTemerature *T = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    T = new NodalTemerature;
    T->initialization(fe->nne);
    T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
  }

  if(is_it_couple_w_chemical >= 0)
  {}

  int compute_stiffness = 1;

  for(int ip = 1; ip<=fe->nint; ip++)
  {

    double Jn = 1.0; // det(F2[Fn], Jn);
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    
    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);    
    fe->update_deformation_gradient(ndofn,u,Fnp1.data,(m->param)->pFI);

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    hFnpa_I = delta_ij(i,j);
    if(is_it_couple_w_thermal >= 0)
    {
      double hFnm1[9];
      err += compute_temperature_at_ip(fe,grid,mat,T->T0,
                                       T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1);

      hFnpa = hFnp1(i,j);
      hJ = det(hFnp1);
      //mid_point_rule(hFnpa.data, hFn.data, hFnp1.data, alpha, DIM_3x3);
      inv(hFnpa,hFnpa_I);
    }

    err += func->get_pF(m,pFnp1.data,2);
    pJ = det(pFnp1);

    err += func->get_pF(m,pFn.data,1);
    err += func->get_F(m,  Fn.data,1);


    mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, DIM_3x3);
    mid_point_rule(   Fr.data,  Fn.data,  Fnp1.data, alpha, DIM_3x3);

    err += inv(pFnpa, pFnpa_I);

    M = hFnpa_I(i,k)*pFnpa_I(k,j);
    eFnpa = Fr(i,k)*M(k,j);
    FrTFr = Fr(k,i)*Fr(k,j);

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dt,alpha, eFnpa.data,
                                                  hFn.data,hFnp1.data,1);
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dt,alpha, eFnpa.data,1);

    err += m->param->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu_all);

    // --> update elasticity part
    memset(L.data, 0, sizeof(double)*DIM_3x3x3x3);
    memset(S.data, 0, sizeof(double)*DIM_3x3);

    err += (m->param)->update_elasticity(m,ctx,L.data,S.data,compute_stiffness);
    // <-- update elasticity part
    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;

    // --> start computing tagent
    // total Lagrangian
    eFn = delta_ij(i,j);
    eFnM = M(i,j);
    eFnMT(j,i) = M(i,j).to(j,i);

    Jn = Jn/pJ/hJ;
    err += compute_stiffness_matrix(lk,fe,
                                    Fr.data,eFnMT.data,eFn.data,M.data,FrTFr.data,eFnM.data,S.data,
                                    L.data,dMdu_all,Jn);
  }

  if(is_it_couple_w_thermal >=0)
    delete T;

  free(u);
  free(dMdu_all);

  return err;
}

/// compute element stiffness matrix in transient for three-field mixed method formulation
///
/// Compute element stiffness matrix based on the mid point rule. When thermal
/// is couled, temperature is assumed constant. Currenlty total Lagrangian is
/// active (Updated lagrangian is implemented, but accelleration term is not
/// fully implemented for updated lagrangian)
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] re_np1 nodal variables at t(n+1) in the current element
/// \param[in] re_npa nodal variables at (1-alpha)r(n) + alpha*r(n+1) in the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_el_constitutive_model_w_inertia_3f(FEMLIB *fe,
                                                 double *lk,
                                                 double *re_np1,
                                                 double *re_npa,
                                                 Grid *grid,
                                                 MaterialProperty *mat,
                                                 FieldVariables *fv,
                                                 Solver *sol,
                                                 LoadingSteps *load,
                                                 CRPL *crpl,
                                                 const PGFem3D_opt *opts,
                                                 const Multiphysics& mp,
                                                 int mp_id,
                                                 double dt)
{
  int err = 0;
  double alpha = sol->alpha;

  double alpha_1 = 1.0 - alpha;
  double alpha_2 = alpha;
  double dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;

  if(opts->cm==UPDATED_LAGRANGIAN)
  {
    PGFEM_printf("Constitutive model with inertia are not supported for updated Lagrangian \n");
    abort();
  }

  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  SUPP sup = load->sups[mp_id];

  Matrix<double> r_np1(nne*nsd, 1), r_npa(nne*nsd, 1), u_npa(nne*nsd, 1), P(Pno, 1);
  Matrix<double> dMdu(DIM_3x3*nne*nsd,1);
  Matrix<double> dMdt(DIM_3x3*Vno,1);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      r_np1.m_pdata[a*nsd+b] = re_np1[a*ndofn+b];
      r_npa.m_pdata[a*nsd+b] = re_npa[a*ndofn+b];      
    }

    if(Pno==nne)
      P.m_pdata[a] = re_npa[a*ndofn+nsd];
  }
  if(Pno==1)
  {
    P(1) = alpha_1*fv->tf.P_n(eid+1,1) +
           alpha_2*(fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1));
  }

  // define xFnp1
  Tensor<2>  Fnp1, pFnp1, hFnp1;

  // define xFnpa
  Tensor<2>  Fr_npa, eFnpa = {}, pFnpa;

  // define xFn
  Tensor<2> Fn, hFn, pFn;

  // define other F and M
  Tensor<2> Fr, eSd = {}, M_npa;
  Tensor<4> Ld = {};

  NodalTemerature *T = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    T = new NodalTemerature;
    T->initialization(fe->nne);
    T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
  }

  if(is_it_couple_w_chemical >= 0)
  {}

  ThreeFieldStiffness K(fe, Vno, Pno);
  Matrix<double> Nt(Vno,1), Np(Pno,1);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double hJnp1 = 1.0;
    double pJnp1 = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    
    Tensor<2> Fr_temp, Fr_npa_temp;
    
    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);        
    fe->update_deformation_gradient(ndofn,r_np1.m_pdata,Fr.data,(m->param)->pFI);
    fe->update_deformation_gradient(ndofn,r_npa.m_pdata,Fr_npa.data,(m->param)->pFI);

    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);


    double theta_npa = 0.0;
    double Pnpa      = 0.0;

    for(int ia=1; ia<=Pno; ia++)
      Pnpa += Np(ia)*P(ia);

    for(int ia=1; ia<=Vno; ia++)
    {
      theta_npa += (alpha_1*fv->tf.V_n(eid+1, ia) +
                    alpha_2*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia)))*Nt(ia);
    }
    

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    if(is_it_couple_w_thermal >= 0)
    {
      double hFnm1[DIM_3x3];
      err += compute_temperature_at_ip(fe,grid,mat,T->T0,
                                       T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1);
      hJnp1 = det(hFnp1);
    }

    // compute deformation gradients
    err += func->get_pF(m, pFnp1.data, 2);
    err += func->get_pF(m, pFn.data,   1);

    if(sup->multi_scale)
      cm_add_macro_F(sup,Fr.data);

    Fnp1 = Fr(i,j);
    mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, DIM_3x3);

    if(is_it_couple_w_thermal>=0)
    {
      Tensor<2> pFnpa_I, hFnp1_I;
      err += inv(pFnpa, pFnpa_I);
      err += inv(hFnp1, hFnp1_I);
      M_npa = hFnp1_I(i,k)*pFnpa_I(k,j);
    }
    else
      err += inv(pFnpa, M_npa);

    eFnpa = Fr_npa(i,k)*M_npa(k,j);

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, func->type, Fnp1.data,dt,alpha,eFnpa.data,
                                                  hFn.data,hFnp1.data,1);
    else
      err += construct_model_context(&ctx, func->type, Fnp1.data,dt,alpha, eFnpa.data,1);

    err += func->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu.m_pdata);
    err += func->compute_dMdt(m, ctx, fe->ST, Vno, dMdt.m_pdata);
    err += func->destroy_ctx(&ctx);


    // <-- update plasticity part
    err += func->update_elasticity_dev(m, eFnpa.data, Ld.data, eSd.data, 1, alpha, dt, 1);

    double MJ_npa = ttl::det(M_npa);
    double theta_e_npa = theta_npa*MJ_npa;

    double dU  = func->compute_dudj(  m, theta_e_npa, 1, alpha);
    double ddU = func->compute_d2udj2(m, theta_e_npa, 1, alpha);
    // --> update elasticity part

    if(err!=0)
      break;

    CM_ThreeField cmtf;
    cmtf.set_femlib(fe,Vno,Pno,Nt.m_pdata,Np.m_pdata);

    double Jn = 1.0/hJnp1/pJnp1;
    cmtf.set_tenosrs(Fr_npa.data, delta_ij.data, M_npa.data, pFnp1.data, eSd.data, Ld.data);
    cmtf.set_scalars(theta_npa, 1.0, 1.0, Jn, Pnpa, dU, ddU, dt_alpha_1_minus_alpha);

    err += K.compute_stiffness(cmtf, dMdu, dMdt);
  }

  if(is_it_couple_w_thermal >=0)
    delete T;

  K.Kpt.trans(K.Ktp);
  K.Kpu.trans(K.Kup);

  err += condense_K_3F_to_1F(lk, nne, nsd, Pno, Vno,
                             K.Kuu.m_pdata, K.Kut.m_pdata, K.Kup.m_pdata,
                             K.Ktu.m_pdata, K.Ktt.m_pdata, K.Ktp.m_pdata,
                             K.Kpu.m_pdata, K.Kpt.m_pdata, NULL);

  // check diagonal for zeros/nans
  for (int a = 0; a < nne; a++) {
    for (int b = 0; b < nsd; b++) {
      if ( !isnormal(lk[idx_K(a,b,a,b,nne,nsd)]) ) err++;
    }
  }

  return err;
}

/// compute element stiffness matrix in transient
///
/// Compute element stiffness matrix based on the mid point rule. When thermal
/// is couled, temperature is assumed constant. Currenlty total Lagrangian is
/// active (Updated lagrangian is implemented, but accelleration term is not
/// fully implemented for updated lagrangian)
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] re_np1 nodal variables at t(n+1) in the current element
/// \param[in] re_npa nodal variables at (1-alpha)r(n) + alpha*r(n+1) in the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_el_constitutive_model_w_inertia(FEMLIB *fe,
                                              double *lk,
                                              double *re_np1,
                                              double *re_npa,
                                              Grid *grid,
                                              MaterialProperty *mat,
                                              FieldVariables *fv,
                                              Solver *sol,
                                              LoadingSteps *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              const Multiphysics& mp,
                                              int mp_id,
                                              double dt)
{
  int err = 0;
  if(opts->analysis_type==CM)
    err += stiffness_el_constitutive_model_w_inertia_1f(fe, lk, re_np1,
                                                        grid, mat, fv, sol, load, crpl,
                                                        opts, mp, mp_id, dt);
  if(opts->analysis_type==CM3F)
    err += stiffness_el_constitutive_model_w_inertia_3f(fe, lk,  re_np1, re_npa,
                                                        grid, mat, fv, sol, load, crpl,
                                                        opts, mp, mp_id, dt);
  return err;
}

/// compute element stiffness matrix in quasi steady state
///
/// Updated Lagrangian and total Lagrangian based. When thermal
/// is couled, temperature is assumed constant.
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_el_constitutive_model_1f(FEMLIB *fe,
                                       double *lk,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       const Multiphysics& mp,
                                       int mp_id,
                                       double dt)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia

  int total_Lagrangian = 1;
  if(opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];

  double *u = (double *) malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = (double *) malloc(sizeof(*dMdu_all)*DIM_3x3*nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2> Fr, Fnp1 = {}, pFnp1, S,
            eFn = {}, M = {}, eFnM = {},eFnMT, FrTFr = {},
            hFn, hFnp1, hFnp1_I;

  Tensor<4> L;

  NodalTemerature *T = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    T = new NodalTemerature;
    T->initialization(fe->nne);
    T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
  }

  if(is_it_couple_w_chemical >=0)
  {}

  int compute_stiffness = 1;

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn = 1.0;
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);    
    fe->update_deformation_gradient(ndofn,u,Fr.data,(m->param)->pFI);

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    // --> update deformations due to coupled physics
    if(is_it_couple_w_thermal >= 0)
    {
      // compute temperature at the integration point
      Tensor<2> hFnm1;
      err += compute_temperature_at_ip(fe,grid,mat,T->T0,
                                       T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1.data);
      hJ = det(hFnp1);
      inv(hFnp1,hFnp1_I);
    }

    // --> update plasticity part
    if(total_Lagrangian)
    {
      if(sup->multi_scale)
        cm_add_macro_F(sup,Fr.data);

      // Total Lagrangian formulation Fn = 1, Fnp1 = Fr
      eFn = delta_ij(i,j);
      Fnp1 = Fr(i,j);
      hFn = delta_ij(i,j);
    }
    else
    {
      if(sup->multi_scale)
      {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }

      Tensor<2> Fn;
      err += m->param->get_F(m, Fn.data,1);
      Fnp1 = Fr(i,k)*Fn(k,j); // Fn+1 = Fr*Fn
      Jn = det(Fn);
    }
    
    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,
                                                  hFn.data,hFnp1.data,-1);
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,-1);

    err += func->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu_all);
    err += func->get_pF(m,pFnp1.data,2);

    pJ = det(pFnp1);


    if(total_Lagrangian) // Total Lagrangian formulation, all xFn = 1
    {
      if(is_it_couple_w_thermal >= 0)
      {
        Tensor<2> pFnp1_I;
        inv(pFnp1, pFnp1_I);
        M = hFnp1_I(i,k)*pFnp1_I(k,j);
      }
      else
        inv(pFnp1, M);

      eFnM = M(i,j);
      eFnMT(j,i) = M(i,j).to(j,i);
    }
    else
    {
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*hFn(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);

        Tensor<2> hFn_I;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        inv(hFn, hFn_I);
        err += m->param->get_eF_of_hF(m,eFn.data,hFn_I.data,stepno);
      }
      else
      {
        Tensor<2> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*pFnp1_I(k,j);
        err += m->param->get_eF(m,eFn.data,1);
      }
      eFnM = eFn(i,k)*M(k,j);
      eFnMT(j,i) = eFnM(i,j).to(j,i);
    }
    // <-- update plasticity part

    FrTFr = Fr(k,i)*Fr(k,j);

    // --> update elasticity part
    memset(L.data, 0, sizeof(double)*DIM_3x3x3x3);
    memset(S.data, 0, sizeof(double)*DIM_3x3);

    if((m->param)->uqcm)
    {
      double *x_ip = (fe->x_ip).m_pdata;
      ELASTICITY *elast = (m->param)->cm_elast;
      mat_e_in = elast->mat;
      err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
      elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation
      err += (m->param)->update_elasticity(m,ctx,L.data,S.data,compute_stiffness);
      elast->mat = mat_e_in;
    }
    else
      err += (m->param)->update_elasticity(m,ctx,L.data,S.data,compute_stiffness);
    // <-- update elasticity part

    err += func->destroy_ctx(&ctx);
    if(err!=0)
      break;      

    // start computing tagent
    Jn = Jn/hJ/pJ;
    err += compute_stiffness_matrix(lk,fe,
                                    Fr.data,eFnMT.data,eFn.data,M.data,FrTFr.data,eFnM.data,S.data,
                                    L.data,dMdu_all,Jn);
  }
  free(u);

  if(is_it_couple_w_thermal >=0)
    delete T;

  /* check diagonal for zeros/nans */
  for (int a = 0; a < nne; a++) {
    for (int b = 0; b < nsd; b++) {
      if ( !isnormal(lk[idx_K(a,b,a,b,nne,nsd)]) ) err++;
    }
  }

  free(dMdu_all);

  return err;
}

/// compute element stiffness matrix in quasi steady state
///
/// Updated Lagrangian and total Lagrangian based. When thermal
/// is couled, temperature is assumed constant.
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error

int stiffness_el_constitutive_model_3f(FEMLIB *fe,
                                       double *lk,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       const Multiphysics& mp,
                                       int mp_id,
                                       double dt)
{
  ConstitutiveModelIntregrate<IntegrateThreeFieldStiffness> cm3f_stiffness;
  return cm3f_stiffness.integrate_ss(fe,lk,r_e,grid,mat,fv,sol->run_integration_algorithm,load,opts,mp,mp_id,dt);
}

/// compute element stiffness matrix in quasi steady state
///
/// Updated Lagrangian and total Lagrangian based. When thermal
/// is couled, temperature is assumed constant.
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_el_constitutive_model(FEMLIB *fe,
                                    double *lk,
                                    double *r_e,
                                    Grid *grid,
                                    MaterialProperty *mat,
                                    FieldVariables *fv,
                                    Solver *sol,
                                    LoadingSteps *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    const Multiphysics& mp,
                                    int mp_id,
                                    double dt)
{
  int err = 0;
  if(opts->analysis_type==CM)
    err += stiffness_el_constitutive_model_1f(fe,lk,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

  if(opts->analysis_type==CM3F)
    err += stiffness_el_constitutive_model_3f(fe,lk,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

  return err;
}

/// compute element residual vector in transient for single field formulation
///
/// redual = residual(n+alpha) + residual(n-1+alpha)
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm.
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1)
///                                                dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] t  time
/// \return non-zero on internal error
int residuals_el_constitutive_model_w_inertia_1f(FEMLIB *fe,
                                                 double *f,
                                                 double *r_e,
                                                 Grid *grid,
                                                 MaterialProperty *mat,
                                                 FieldVariables *fv,
                                                 Solver *sol,
                                                 LoadingSteps *load,
                                                 CRPL *crpl,
                                                 const PGFem3D_opt *opts,
                                                 const Multiphysics& mp,
                                                 const double *dts,
                                                 int mp_id,
                                                 double dt,
                                                 const double t)
{
  int err = 0;
  double alpha = sol->alpha;
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

// when updated Lagrangian is used should be un-commented
//  if(opts->cm != 0)
//    total_Lagrangian = 1;

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;

  double *u       = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_npa   = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_nm1pa = (double *) malloc(sizeof(double)*nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2> Fnp1,Fn,Fnm1,pFnp1,pFn,pFnm1,
            hFnp1,hFn,hFnm1;

  NodalTemerature *T = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    T = new NodalTemerature;
    T->initialization(fe->nne);
    T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
  }

  if(is_it_couple_w_chemical >=0)
  {}

  memset(f_npa, 0, sizeof(double)*nne*ndofn);
  memset(f_nm1pa, 0, sizeof(double)*nne*ndofn);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    
    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);    
    fe->update_deformation_gradient(ndofn,u,Fnp1.data,(m->param)->pFI);

    if(is_it_couple_w_thermal >= 0)
    {
      err += compute_temperature_at_ip(fe,grid,mat,T->T0,
                                       T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1.data);
      hJ = det(hFnp1);
    }

    // perform integration algorithm
    if(sol->run_integration_algorithm)
      err += m->run_integration_algorithm(Fnp1.data,hFn.data,hFnp1.data,dts[DT_NP1],alpha,fe->x_ip.m_pdata,t,is_it_couple_w_thermal);

    if(err!=0)
      break;

    err += m->param->get_pF(m, pFnp1.data,2);
    err += m->param->get_pF(m,   pFn.data,1);
    err += m->param->get_pF(m, pFnm1.data,0);
    err += m->param->get_F(m,     Fn.data,1);
    err += m->param->get_F(m,   Fnm1.data,0);

    pJ = det(pFnp1);

    double dt_1_minus_alpha = -dts[DT_NP1]*(1.0-alpha)*pJ*hJ;
    err += residuals_el_constitutive_model_n_plus_alpha(f_npa,m,eid,ndofn,
                                                        pFnp1.data,pFn.data,Fnp1.data,Fn.data,
                                                        hFnp1.data,hFn.data,
                                                        is_it_couple_w_thermal,
                                                        alpha, dt_1_minus_alpha,fe,1);

    double dt_alpha = -dts[DT_N]*alpha*pJ*hJ;

    err += residuals_el_constitutive_model_n_plus_alpha(f_nm1pa,m,eid,ndofn,
                                                        pFn.data,pFnm1.data,Fn.data,Fnm1.data,
                                                        hFnp1.data,hFn.data,
                                                        is_it_couple_w_thermal,
                                                        alpha, dt_alpha,fe,0);
  }

  if(is_it_couple_w_thermal >=0)
    delete T;

  if(err==0)
  {
    for(int a=0; a<nne*nsd; a++)
      f[a] += f_npa[a] + f_nm1pa[a];
  }

  // free memory for thermal part

  free(u);
  free(f_npa);
  free(f_nm1pa);

  return err;
}


/// compute element residual vector in transient for three-field mixed method formulation
///
/// redual = residual(n+alpha) + residual(n-1+alpha)
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm.
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] re_np1 nodal variabls(displacements) on the current element
/// \param[in] re_npa nodal variables at (1-alpha)r(n)   + alpha*r(n+1) in the current element
/// \param[in] re_nma nodal variables at (1-alpha)r(n-1) + alpha*r(n)   in the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1)
///                                                dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] t  time
/// \return non-zero on internal error
int residuals_el_constitutive_model_w_inertia_3f(FEMLIB *fe,
                                                 double *f,
                                                 double *re_np1,
                                                 double *re_npa,
                                                 double *re_nma,
                                                 Grid *grid,
                                                 MaterialProperty *mat,
                                                 FieldVariables *fv,
                                                 Solver *sol,
                                                 LoadingSteps *load,
                                                 CRPL *crpl,
                                                 const PGFem3D_opt *opts,
                                                 const Multiphysics& mp,
                                                 const double *dts,
                                                 int mp_id,
                                                 double dt,
                                                 const double t)
{
  int err = 0;
  double alpha = sol->alpha;

  double alpha_1 = 1.0 - alpha;
  double alpha_2 = alpha;
  double dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;

  if(opts->cm==UPDATED_LAGRANGIAN)
  {
    PGFEM_printf("Constitutive model with inertia are not supported for updated Lagrangian \n");
    abort();
  }

  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  SUPP sup = load->sups[mp_id];

  Matrix<double> r_np1(nne*nsd, 1), r_npa(nne*nsd, 1), r_nma(nne*nsd, 1), P_npa(Pno, 1), P_nma(Pno, 1);
  Matrix<double> dMdu(DIM_3x3*nne*nsd,1);
  Matrix<double> dMdt(DIM_3x3*Vno,1);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      r_np1.m_pdata[a*nsd+b] = re_np1[a*ndofn+b];
      r_npa.m_pdata[a*nsd+b] = re_npa[a*ndofn+b];
      r_nma.m_pdata[a*nsd+b] = re_nma[a*ndofn+b];
    }

    if(Pno==nne)
    {
      P_npa.m_pdata[a] = re_npa[a*ndofn+nsd];
      P_nma.m_pdata[a] = re_nma[a*ndofn+nsd];
    }
  }
  if(Pno==1)
  {
    P_npa(1) = alpha_1*fv->tf.P_n(eid+1,1) +
               alpha_2*(fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1));
    P_nma(1) = alpha_1*fv->tf.P_nm1(eid+1,1) + alpha_2*fv->tf.P_n(eid+1,1);
  }

  // define xFnp1
  Tensor<2>  Fnp1, pFnp1, hFnp1;

  // define xFnpa
  Tensor<2>  Fr_npa, eFnpa, pFnpa;

  // define xFnma
  Tensor<2> Fr_nma, eFnma, pFnma;

  // define xFn
  Tensor<2> hFn, pFn;

  // define xFnm1
  Tensor<2> pFnm1;

  // define other F and M
  Tensor<2> Fr, S_npa = {}, S_nma = {}, M_npa, M_nma;
  Tensor<4> Ld_npa = {}, Ld_nma;

  NodalTemerature *T = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    T = new NodalTemerature;
    T->initialization(fe->nne);
    T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
  }

  if(is_it_couple_w_chemical >= 0)
  {}

  ThreeFieldStiffness K(fe, Vno, Pno, true);
  ThreeFieldResidual R_npa(fe, Vno, Pno), R_nma(fe, Vno, Pno);

  Matrix<double> Nt(Vno,1), Np(Pno,1);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double hJnp1 = 1.0;
    double pJnp1 = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    
    Tensor<2> Fr_temp, Fr_npa_temp, Fr_nma_temp;

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);    
    fe->update_deformation_gradient(ndofn,r_np1.m_pdata,    Fr.data,(m->param)->pFI);
    fe->update_deformation_gradient(ndofn,r_npa.m_pdata,Fr_npa.data,(m->param)->pFI);
    fe->update_deformation_gradient(ndofn,r_nma.m_pdata,Fr_nma.data,(m->param)->pFI);    
    
    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);

    double theta     = 0.0;
    double theta_npa = 0.0;
    double theta_nma = 0.0;
    double Pnpa      = 0.0;
    double Pnma      = 0.0;

    for(int ia=1; ia<=Pno; ia++)
    {
      Pnpa += Np(ia)*P_npa(ia);
      Pnma += Np(ia)*P_nma(ia);
    }
    for(int ia=1; ia<=Vno; ia++)
    {
      theta     += (fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia))*Nt(ia);
      theta_npa += (alpha_1*fv->tf.V_n(eid+1, ia) +
                    alpha_2*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia)))*Nt(ia);

      theta_nma += (alpha_1*fv->tf.V_nm1(eid+1, ia) + alpha_2*fv->tf.V_n(eid+1, ia))*Nt(ia);
    }

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    if(is_it_couple_w_thermal >= 0)
    {
      double hFnm1[DIM_3x3];
      err += compute_temperature_at_ip(fe,grid,mat,T->T0,
                                       T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1);
      hJnp1 = det(hFnp1);
    }

    if(sup->multi_scale)
      cm_add_macro_F(sup,Fr.data);

    Fnp1 = Fr(i,j);

    // perform integration algorithm
    if(sol->run_integration_algorithm)
    {
      double tJ = det(Fnp1);
      if(tJ<0)
        ++err;
      else{
        double tf_factor = pow(theta/tJ, 1.0/3.0);
        err += m->run_integration_algorithm(Fnp1.data,hFn.data,hFnp1.data,dts[DT_NP1],alpha,fe->x_ip.m_pdata,t,is_it_couple_w_thermal, tf_factor);
      }
    }

    // compute deformation gradients
    err += func->get_pF(m, pFnp1.data, 2);
    err += func->get_pF(m, pFn.data,   1);
    err += func->get_pF(m, pFnm1.data, 0);

    mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, DIM_3x3);
    mid_point_rule(pFnma.data, pFnm1.data, pFn.data, alpha, DIM_3x3);

    if(is_it_couple_w_thermal>=0)
    {
      Tensor<2> pFnpa_I, pFnma_I, hFnp1_I, hFn_I;
      err += inv(pFnpa, pFnpa_I);
      err += inv(pFnma, pFnma_I);
      err += inv(hFnp1, hFnp1_I);
      err += inv(hFn,   hFn_I);
      M_npa = hFnp1_I(i,k)*pFnpa_I(k,j);
      M_nma = hFn_I(i,k)*pFnma_I(k,j);
    }
    else
    {
      err += inv(pFnpa, M_npa);
      err += inv(pFnma, M_nma);
    }

    eFnpa = Fr_npa(i,k)*M_npa(k,j);
    eFnma = Fr_nma(i,k)*M_nma(k,j);

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, func->type, Fnp1.data,dt,alpha,eFnpa.data,
                                                  hFn.data,hFnp1.data,1);
    else
      err += construct_model_context(&ctx, func->type, Fnp1.data,dt,alpha, eFnpa.data,1);

    err += func->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu.m_pdata);
    err += func->compute_dMdt(m, ctx, fe->ST, Vno, dMdt.m_pdata);
    err += func->destroy_ctx(&ctx);

    // <-- update plasticity part
    err += func->update_elasticity_dev(m, eFnpa.data, Ld_npa.data, S_npa.data, 1, alpha, dts[DT_NP1], 1);
    err += func->update_elasticity_dev(m, eFnma.data, Ld_nma.data, S_nma.data, 0, alpha, dts[DT_N],   1);

    double MJ_npa = ttl::det(M_npa);
    double MJ_nma = ttl::det(M_nma);
    double theta_e_npa = theta_npa*MJ_npa;
    double theta_e_nma = theta_nma*MJ_nma;

    double dU_npa  = func->compute_dudj(  m, theta_e_npa, 1, alpha);
    double ddU_npa = func->compute_d2udj2(m, theta_e_npa, 1, alpha);

    double dU_nma  = func->compute_dudj(  m, theta_e_nma, 0, alpha);
    double ddU_nma = func->compute_d2udj2(m, theta_e_nma, 0, alpha);
    // --> update elasticity part

    if(err!=0)
      break;

    CM_ThreeField cmtf_npa, cmtf_nma;
    cmtf_npa.set_femlib(fe,Vno,Pno,Nt.m_pdata,Np.m_pdata);
    cmtf_nma.set_femlib(fe,Vno,Pno,Nt.m_pdata,Np.m_pdata);

    double Jn = 1.0/hJnp1/pJnp1;
    cmtf_npa.set_tenosrs(Fr_npa.data, delta_ij.data, M_npa.data, pFnpa.data, S_npa.data, Ld_npa.data);
    cmtf_npa.set_scalars(theta_npa, 1.0, 1.0, Jn, Pnpa, dU_npa, ddU_npa, dt_alpha_1_minus_alpha);

    cmtf_nma.set_tenosrs(Fr_nma.data, delta_ij.data, M_nma.data, pFnma.data, S_nma.data, Ld_nma.data);
    cmtf_nma.set_scalars(theta_nma, 1.0, 1.0, Jn, Pnma, dU_nma, ddU_nma);

    err += K.compute_stiffness(cmtf_npa, dMdu, dMdt);
    err += R_npa.compute_residual(cmtf_npa);
    err += R_nma.compute_residual(cmtf_nma);
  }

  if(is_it_couple_w_thermal >=0)
    delete T;

  if(err==0)
  {
    ThreeFieldResidual R(fe, Vno, Pno);

    double dt_1_minus_alpha = -dts[DT_NP1]*(1.0-alpha);
    double dt_alpha = -dts[DT_N]*alpha;
    for(int ia=0; ia<nne*nsd; ia++)
      R.Ru.m_pdata[ia] = dt_1_minus_alpha*R_npa.Ru.m_pdata[ia] +  dt_alpha*R_nma.Ru.m_pdata[ia];

    for(int ia=0; ia<Pno; ia++)
      R.Rp.m_pdata[ia] = dt_1_minus_alpha*R_npa.Rp.m_pdata[ia] +  dt_alpha*R_nma.Rp.m_pdata[ia];

    for(int ia=0; ia<Vno; ia++)
      R.Rt.m_pdata[ia] = dt_1_minus_alpha*R_npa.Rt.m_pdata[ia] +  dt_alpha*R_nma.Rt.m_pdata[ia];

    K.Kpt.trans(K.Ktp);

    err += condense_F_3F_to_1F(f, nne, nsd, Pno, Vno,
                               R.Ru.m_pdata, R.Rt.m_pdata, R.Rp.m_pdata,
                               K.Kut.m_pdata, K.Kup.m_pdata, K.Ktp.m_pdata, K.Ktt.m_pdata, K.Kpt.m_pdata);
  }

  return err;
}

/// compute element residual vector in transient
///
/// redual = residual(n+alpha) + residual(n-1+alpha)
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm.
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] re_np1 nodal variables at t(n+1) in the current element
/// \param[in] re_npa nodal variables at (1-alpha)r(n)   + alpha*r(n+1) in the current element
/// \param[in] re_nma nodal variables at (1-alpha)r(n-1) + alpha*r(n)   in the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1)
///                                                dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] t  time
/// \return non-zero on internal error
int residuals_el_constitutive_model_w_inertia(FEMLIB *fe,
                                              double *f,
                                              double *re_np1,
                                              double *re_npa,
                                              double *re_nma,
                                              Grid *grid,
                                              MaterialProperty *mat,
                                              FieldVariables *fv,
                                              Solver *sol,
                                              LoadingSteps *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              const Multiphysics& mp,
                                              const double *dts,
                                              int mp_id,
                                              const double dt,
                                              const double t)
{
  int err = 0;

  if(opts->analysis_type==CM)
    err += residuals_el_constitutive_model_w_inertia_1f(fe, f, re_np1,
                                                        grid, mat, fv, sol, load, crpl,
                                                        opts, mp, dts, mp_id, dt, t);
  if(opts->analysis_type==CM3F)
    err += residuals_el_constitutive_model_w_inertia_3f(fe, f, re_np1, re_npa, re_nma,
                                                        grid, mat, fv, sol, load, crpl,
                                                        opts, mp, dts, mp_id, dt, t);
  return err;
}
/// compute element residual vector in quasi steady state
///
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm.
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int residuals_el_constitutive_model_1f(FEMLIB *fe,
                                       double *f,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       const Multiphysics& mp,
                                       int mp_id,
                                       double dt)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia
  int total_Lagrangian = 1;
  if(opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;
    
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];

  double *u = PGFEM_malloc<double>(nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2> Fr, Fnp1 = {},pFnp1,
            S = {}, M = {}, eFnM = {},eFnMT,
            hFn,hFnp1,hFnp1_I;

  NodalTemerature *T = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    T = new NodalTemerature;
    T->initialization(fe->nne);
    T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
  }

  if(is_it_couple_w_chemical >=0)
  {}

  int compute_stiffness = 0;

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn = 1.0; // if upated Lagrangian, Jn = det(Fn), later updated
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);    
    fe->update_deformation_gradient(ndofn,u,Fr.data,(m->param)->pFI);

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    // --> update deformations due to coupled physics
    if(is_it_couple_w_thermal >= 0)
    {
      // compute temperature at the integration point
      double hFnm1[9];
      err += compute_temperature_at_ip(fe,grid,mat,T->T0,
                                       T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1);
      hJ = det(hFnp1);
      inv(hFnp1, hFnp1_I);
      if(total_Lagrangian)
        hFn = delta_ij(i,j);
    }

    // --> update plasticity part
    if(total_Lagrangian)
    {
      if (sup->multi_scale) {
        cm_add_macro_F(sup,Fr.data);
      }

      // TOTAL LAGRANGIAN FORMULATION Fn = 1, Fnp1 = Fr
      Fnp1 = Fr(i,j);
    }
    else
    {
      if (sup->multi_scale) {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }

      Tensor<2> Fn;
      err += m->param->get_F(m, Fn.data,1);
      Fnp1 = Fr(i,k)*Fn(k,j);  // compute Fnp1 = Fr*Fn
      Jn = det(Fn);
    }

    {
      /* check that deformation is invertible -> J > 0 */
      int terr = 0;
      getJacobian(Fnp1.data, eid, &terr);
      err += terr;
    }

    // perform integration algorithm
    if(sol->run_integration_algorithm)
      err += m->run_integration_algorithm(Fnp1.data,hFn.data,hFnp1.data,dt,alpha,fe->x_ip.m_pdata,0.0,is_it_couple_w_thermal);

    if(err!=0)
      break;

    err += func->get_pF(m,pFnp1.data,2);
    pJ = det(pFnp1);

    if(total_Lagrangian)
    {
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2> pFnp1_I;
        inv(pFnp1, pFnp1_I);
        M = hFnp1_I(i,k)*pFnp1_I(k,j);
      }
      else
        inv(pFnp1, M);

      eFnM = M(i,j);
      eFnMT(j,i) = M(i,j).to(j,i);

    }
    else
    {
      Tensor<2> eFn;
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*hFn(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);

        Tensor<2> hFn_I;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        inv(hFn, hFn_I);
        err += m->param->get_eF_of_hF(m,eFn.data,hFn_I.data,stepno);

      }
      else
      {
        Tensor<2> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*pFnp1_I(k,j);
        err += m->param->get_eF(m,eFn.data,1);
      }
      eFnM = eFn(i,k)*M(k,j);
      eFnMT(j,i) = eFnM(i,j).to(j,i);
    }

    // <-- update plasticity part

    // --> update elasticity part

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,
                                                  hFn.data,hFnp1.data,-1);
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,-1);
    if((m->param)->uqcm)
    {
      double *x_ip = (fe->x_ip).m_pdata;
      ELASTICITY *elast = (m->param)->cm_elast;
      mat_e_in = elast->mat;
      err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
      elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation
      err += (m->param)->update_elasticity(m,ctx,NULL,S.data,compute_stiffness);
      elast->mat = mat_e_in;
    }
    else
      err += (m->param)->update_elasticity(m,ctx,NULL,S.data,compute_stiffness);
    // <-- update elasticity part
    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;

    Jn = Jn/pJ/hJ;
    err += compute_residual_vector(f,fe,Fr.data,eFnMT.data,eFnM.data,S.data,Jn);    
  }

  if(is_it_couple_w_thermal >=0)
    delete T;

  free(u);

  return err;
}

/// compute element residual vector in quasi steady state
///
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm.
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int residuals_el_constitutive_model_3f(FEMLIB *fe,
                                       double *f,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       const Multiphysics& mp,
                                       int mp_id,
                                       double dt)
{
  ConstitutiveModelIntregrate<IntegrateThreeFieldResidual> cm3f_residual;
  return cm3f_residual.integrate_ss(fe,f,r_e,grid,mat,fv,sol->run_integration_algorithm,load,opts,mp,mp_id,dt);
}

/// compute element residual vector in quasi steady state
///
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm.
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int residuals_el_constitutive_model(FEMLIB *fe,
                                    double *f,
                                    double *r_e,
                                    Grid *grid,
                                    MaterialProperty *mat,
                                    FieldVariables *fv,
                                    Solver *sol,
                                    LoadingSteps *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    const Multiphysics& mp,
                                    int mp_id,
                                    double dt)
{
  int err = 0;

  if(opts->analysis_type==CM)
    err += residuals_el_constitutive_model_1f(fe,f,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

  if(opts->analysis_type==CM3F)
    err += residuals_el_constitutive_model_3f(fe,f,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

  return err;
}

int constitutive_model_update_NR_w_inertia_3f(FEMLIB *fe,
                                              double *re_np1,
                                              double *re_npa,
                                              double *re_nma,
                                              double *du,
                                              Grid *grid,
                                              MaterialProperty *mat,
                                              FieldVariables *fv,
                                              LoadingSteps *load,
                                              const PGFem3D_opt *opts,
                                              const Multiphysics& mp,
                                              const double *dts,
                                              int mp_id,
                                              double alpha)
{
  int err = 0;

  double alpha_1 = 1.0 - alpha;
  double alpha_2 = alpha;
  double dt_alpha_1_minus_alpha = dts[DT_NP1]*alpha_1*alpha_2;

  if(opts->cm==UPDATED_LAGRANGIAN)
  {
    PGFEM_printf("Constitutive model with inertia are not supported for updated Lagrangian \n");
    abort();
  }

  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  SUPP sup = load->sups[mp_id];

  Matrix<double> r_np1(nne*nsd, 1), r_npa(nne*nsd, 1), r_nma(nne*nsd, 1), P_npa(Pno, 1), P_nma(Pno, 1);
  Matrix<double> dMdu(DIM_3x3*nne*nsd,1);
  Matrix<double> dMdt(DIM_3x3*Vno,1);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      r_np1.m_pdata[a*nsd+b] = re_np1[a*ndofn+b];
      r_npa.m_pdata[a*nsd+b] = re_npa[a*ndofn+b];
      r_nma.m_pdata[a*nsd+b] = re_nma[a*ndofn+b];            
    }

    if(Pno==nne)
    {
      P_npa.m_pdata[a] = re_npa[a*ndofn+nsd];
      P_nma.m_pdata[a] = re_nma[a*ndofn+nsd];
    }
  }
  if(Pno==1)
  {
    P_npa(1) = alpha_1*fv->tf.P_n(eid+1,1) +
               alpha_2*(fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1));
    P_nma(1) = alpha_1*fv->tf.P_nm1(eid+1,1) + alpha_2*fv->tf.P_n(eid+1,1);
  }

  // define xFnp1
  Tensor<2>  Fnp1, pFnp1, hFnp1;

  // define xFnpa
  Tensor<2>  Fr_npa, eFnpa, pFnpa;

  // define xFnma
  Tensor<2> Fr_nma, eFnma, pFnma;

  // define xFn
  Tensor<2> hFn, pFn;

  // define xFnm1
  Tensor<2> pFnm1;

  // define other F and M
  Tensor<2> Fr, S_npa = {}, S_nma = {}, M_npa, M_nma;
  Tensor<4> Ld_npa = {}, Ld_nma;

  NodalTemerature *T = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    T = new NodalTemerature;
    T->initialization(fe->nne);
    T->get_temperature(fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
  }

  if(is_it_couple_w_chemical >= 0)
  {}

  ThreeFieldStiffness K(fe, Vno, Pno, true);
  ThreeFieldResidual R_npa(fe, Vno, Pno, true), R_nma(fe, Vno, Pno, true);

  Matrix<double> Nt(Vno,1), Np(Pno,1);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double hJnp1 = 1.0;
    double pJnp1 = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    
    Tensor<2> Fr_temp, Fr_npa_temp, Fr_nma_temp;
    
    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);
    fe->update_deformation_gradient(ndofn,r_np1.m_pdata,Fr.data,(m->param)->pFI);
    fe->update_deformation_gradient(ndofn,r_npa.m_pdata,Fr_npa.data,(m->param)->pFI);
    fe->update_deformation_gradient(ndofn,r_nma.m_pdata,Fr_nma.data,(m->param)->pFI);    

    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);


    double theta_npa = 0.0;
    double theta_nma = 0.0;
    double Pnpa      = 0.0;
    double Pnma      = 0.0;

    for(int ia=1; ia<=Pno; ia++)
    {
      Pnpa += Np(ia)*P_npa(ia);
      Pnma += Np(ia)*P_nma(ia);
    }
    for(int ia=1; ia<=Vno; ia++)
    {
      theta_npa += (alpha_1*fv->tf.V_n(eid+1, ia) +
                    alpha_2*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia)))*Nt(ia);

      theta_nma += (alpha_1*fv->tf.V_nm1(eid+1, ia) + alpha_2*fv->tf.V_n(eid+1, ia))*Nt(ia);
    }

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    if(is_it_couple_w_thermal >= 0)
    {
      double hFnm1[DIM_3x3];
      err += compute_temperature_at_ip(fe,grid,mat,T->T0,
                                       T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1);
      hJnp1 = det(hFnp1);
    }

    // compute deformation gradients
    err += func->get_pF(m, pFnp1.data, 2);
    err += func->get_pF(m, pFn.data,   1);
    err += func->get_pF(m, pFnm1.data, 0);

    if(sup->multi_scale)
      cm_add_macro_F(sup,Fr.data);

    Fnp1 = Fr(i,j);
    mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, DIM_3x3);
    mid_point_rule(pFnma.data, pFnm1.data, pFn.data, alpha, DIM_3x3);

    if(is_it_couple_w_thermal>=0)
    {
      Tensor<2> pFnpa_I, pFnma_I, hFnp1_I, hFn_I;
      err += inv(pFnpa, pFnpa_I);
      err += inv(pFnma, pFnma_I);
      err += inv(hFnp1, hFnp1_I);
      err += inv(hFn,   hFn_I);
      M_npa = hFnp1_I(i,k)*pFnpa_I(k,j);
      M_nma = hFn_I(i,k)*pFnma_I(k,j);
    }
    else
    {
      err += inv(pFnpa, M_npa);
      err += inv(pFnma, M_nma);
    }

    eFnpa = Fr_npa(i,k)*M_npa(k,j);
    eFnma = Fr_nma(i,k)*M_nma(k,j);

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, func->type, Fnp1.data,dts[DT_NP1],alpha,eFnpa.data,
                                                  hFn.data,hFnp1.data,1);
    else
      err += construct_model_context(&ctx, func->type, Fnp1.data,dts[DT_NP1],alpha, eFnpa.data,1);

    err += func->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu.m_pdata);
    err += func->compute_dMdt(m, ctx, fe->ST, Vno, dMdt.m_pdata);
    err += func->destroy_ctx(&ctx);

    // <-- update plasticity part
    err += func->update_elasticity_dev(m, eFnpa.data, Ld_npa.data, S_npa.data, 1, alpha, dts[DT_NP1], 1);
    err += func->update_elasticity_dev(m, eFnma.data, Ld_nma.data, S_nma.data, 0, alpha, dts[DT_N],   1);

    double MJ_npa = ttl::det(M_npa);
    double MJ_nma = ttl::det(M_nma);
    double theta_e_npa = theta_npa*MJ_npa;
    double theta_e_nma = theta_nma*MJ_nma;

    double dU_npa  = func->compute_dudj(  m, theta_e_npa, 1, alpha);
    double ddU_npa = func->compute_d2udj2(m, theta_e_npa, 1, alpha);

    double dU_nma  = func->compute_dudj(  m, theta_e_nma, 0, alpha);
    double ddU_nma = func->compute_d2udj2(m, theta_e_nma, 0, alpha);
    // --> update elasticity part

    if(err!=0)
      break;

    CM_ThreeField cmtf_npa, cmtf_nma;
    cmtf_npa.set_femlib(fe,Vno,Pno,Nt.m_pdata,Np.m_pdata);
    cmtf_nma.set_femlib(fe,Vno,Pno,Nt.m_pdata,Np.m_pdata);

    double Jn = 1.0/hJnp1/pJnp1;
    cmtf_npa.set_tenosrs(Fr_npa.data, delta_ij.data, M_npa.data, pFnp1.data, S_npa.data, Ld_npa.data);
    cmtf_npa.set_scalars(theta_npa, 1.0, 1.0, Jn, Pnpa, dU_npa, ddU_npa, dt_alpha_1_minus_alpha);


    cmtf_nma.set_tenosrs(Fr_nma.data, delta_ij.data, M_nma.data, pFnp1.data, S_nma.data, Ld_nma.data);
    cmtf_nma.set_scalars(theta_nma, 1.0, 1.0, Jn, Pnma, dU_nma, ddU_nma);

    err += K.compute_stiffness(cmtf_npa, dMdu, dMdt);
    err += R_npa.compute_residual(cmtf_npa);
    err += R_nma.compute_residual(cmtf_nma);
  }

  if(err==0)
  {
    ThreeFieldResidual R(fe, Vno, Pno, true);

    double dt_1_minus_alpha = -dts[DT_NP1]*(1.0-alpha);
    double dt_alpha = -dts[DT_N]*alpha;

    for(int ia=0; ia<Pno; ia++)
      R.Rp.m_pdata[ia] = dt_1_minus_alpha*R_npa.Rp.m_pdata[ia] +  dt_alpha*R_nma.Rp.m_pdata[ia];

    for(int ia=0; ia<Vno; ia++)
      R.Rt.m_pdata[ia] = dt_1_minus_alpha*R_npa.Rt.m_pdata[ia] +  dt_alpha*R_nma.Rt.m_pdata[ia];

    Matrix<double> Kpu(Pno,nne*nsd);
    Kpu.trans(K.Kup);
    K.Kpt.trans(K.Ktp);

    Matrix<double> d_theta(Vno, 1), dP(Pno, 1);
    err += compute_d_theta_dP(d_theta.m_pdata, dP.m_pdata, du,
                              nne, nsd, Pno, Vno,
                              R.Ru.m_pdata, R.Rt.m_pdata, R.Rp.m_pdata,
                              Kpu.m_pdata, K.Ktu.m_pdata, K.Ktp.m_pdata, K.Ktt.m_pdata, K.Kpt.m_pdata);

    for(int ia=1; ia<=Pno; ia++)
      fv->tf.ddP(eid+1,ia) = dP(ia);

    for(int ia=1; ia<=Vno; ia++)
      fv->tf.ddV(eid+1,ia) = d_theta(ia);
  }

  return err;
}

/// compute ouput variables e.g. effective stress and strain
///
/// Visit each element and compute output variables according to the element model type.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] alpha mid point rule alpha
/// \return non-zero on internal error
int constitutive_model_update_NR(Grid *grid,
                                 MaterialProperty *mat,
                                 FieldVariables *fv,
                                 LoadingSteps *load,
                                 const PGFem3D_opt *opts,
                                 const Multiphysics& mp,
                                 int mp_id,
                                 const double *dts,
                                 double alpha)
{
  int err = 0;

  if(opts->analysis_type==CM) // nothing to do for CM model
    return err;

  int total_Lagrangian = 1;
  if(opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  Node *node = grid->node;
  Element *elem = grid->element;

  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];

  for (int eid = 0; eid < grid->ne; eid++)
  {
    FEMLIB fe(eid,elem,node,0,total_Lagrangian);
    int nsd   = fe.nsd;
    int nne   = fe.nne;
    int ndofe = nne*ndofn;

    Matrix<long> cn(ndofe,1);
    long *nod = fe.node_id.m_pdata;

    Matrix<double> r_e(ndofe, 1), dr_e(ndofe,1), du(nne*nsd,1), u(nne*nsd, 1);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,grid->node,cn.m_pdata,mp_id);

    // get the deformation on the element
    def_elem_total(cn.m_pdata,ndofe,fv->u_np1,fv->d_u,elem,node,sup,r_e.m_pdata);
    def_elem(cn.m_pdata,ndofe,fv->dd_u,elem,node,dr_e.m_pdata,sup,2);


    for(int ia=0;ia<nne;ia++)
    {
      for(int ib=0; ib<nsd;ib++)
        du.m_pdata[ia*nsd+ib] = dr_e.m_pdata[ia*ndofn+ib];
    }

    const int mat_id = grid->element[eid].mat[2];
    double rho = mat->hommat[mat_id].density;
    // make a decision to include ineria
    bool include_inertia = true;
    if(fabs(rho)<MIN_DENSITY)
      include_inertia = false;

    if(include_inertia)
    {
      Matrix<double> re_npa(nne*nsd, 1, 0.0), re_nma(nne*nsd, 1, 0.0);
      Matrix<double> u_n(nne*nsd, 1, 0.0), u_nm1(nne*nsd, 1, 0.0);
      for(int ia=0;ia<nne;ia++)
      {
        for(int ib=0; ib<nsd; ib++)
        {
            u_n.m_pdata[ia*nsd + ib] =   fv->u_n[nod[ia]*ndofn + ib];
          u_nm1.m_pdata[ia*nsd + ib] = fv->u_nm1[nod[ia]*ndofn + ib];
        }
      }

      mid_point_rule(re_nma.m_pdata, u_nm1.m_pdata, u_n.m_pdata, alpha, nne*nsd);
      mid_point_rule(re_npa.m_pdata,   u_n.m_pdata, r_e.m_pdata, alpha, nne*nsd);

      err += constitutive_model_update_NR_w_inertia_3f(&fe, r_e.m_pdata, re_npa.m_pdata, re_nma.m_pdata, du.m_pdata,
                                                       grid, mat, fv, load, opts, mp, dts, mp_id, alpha);
    }
    else
    {
      ConstitutiveModelIntregrate<IntegrateThreeFieldUpdate> cm3f_update;
      err += cm3f_update.integrate_ss(&fe,NULL,r_e.m_pdata,grid,mat,fv,0,load,opts,mp,mp_id,dts[DT_NP1]);
    }
  }

  return err;
}


int cm_write_tensor_restart(FILE *fp, const double *tensor)
{
  int err = 0;
  fprintf(fp, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          tensor[0], tensor[1], tensor[2],
          tensor[3], tensor[4], tensor[5],
          tensor[6], tensor[7], tensor[8]);
  return err;
}

int cm_read_tensor_restart(FILE *fp, double *tensor)
{
  int err = 0;
  fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &tensor[0], &tensor[1], &tensor[2],
         &tensor[3], &tensor[4], &tensor[5],
         &tensor[6], &tensor[7], &tensor[8]);
  return err;
}

/// compute and set initial conditions for three field mixed method
///
/// \param[in] grid  an object containing all mesh info
/// \param[in] mat   a material object
/// \param[in,out]   fv array of field variable object
/// \param[in] load  object for loading
/// \param[in] mp    mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] analysis_type type of analysis such as cm, cm3f
/// \return non-zero on internal error
void compute_cm_initial_conditions(Grid *grid,
                                   const MaterialProperty *mat,
                                   FieldVariables *fv,                                     
                                   LoadingSteps *load,
                                   const Multiphysics &mp,
                                   const int mp_id,
                                   const int analysis_type){

  const int total_Lagrangian = 1;
  const int intg_order = 0;
     
  EPS *eps = fv->eps;
  Node *node = grid->node;
  Element *elem = grid->element;
  
  int is_it_couple_w_thermal  = -1;

  for(int ia=0; ia<fv->n_coupled; ia++){
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
  }    

  int ndofn = fv->ndofn;      

  for (int eid=0;eid<grid->ne;eid++)
  {
    FEMLIB fe;
    Constitutive_model *m = (fv->eps[eid]).model;
    fe.initialization(eid,elem,node,intg_order,total_Lagrangian,m->param->pFI);
    
    int nne   = fe.nne;
    int nsd   = fe.nsd;
    
    NodalTemerature *T = NULL;
    long *nod = fe.node_id.m_pdata;

    if(is_it_couple_w_thermal >= 0)
    {
      T = new NodalTemerature;
      T->initialization(fe.nne);
      T->get_temperature(&fe, grid, fv, load, mp, mp_id, is_it_couple_w_thermal);
    }    

    Matrix<double> u_nm1(nne*nsd, 1), u_n(nne*nsd, 1);
    
    for (long I=0;I<nne;I++){
      for(long J=0; J<nsd; J++){
        u_n.m_pdata[I*ndofn + J] =   fv->u_n[nod[I]*ndofn + J];
        u_nm1.m_pdata[  I*ndofn + J] = fv->u_nm1[nod[I]*ndofn + J];
      }
    }

    double eJ_n    = 0.0;
    double eJ_nm1  = 0.0;
    double Up_n   = 0.0;
    double Up_nm1 = 0.0;
    double V = 0.0;
    
    for(int ip=1; ip<=fe.nint; ip++){

      fe.elem_basis_V(ip);

      const Model_parameters *func =eps[eid].model[ip-1].param;

      Tensor<2> Fn, Fnm1;
      
      fe.update_shape_tensor();
      fe.update_deformation_gradient(nsd,u_n.m_pdata,    Fn.data);
      fe.update_deformation_gradient(nsd,u_nm1.m_pdata,Fnm1.data);
      
      func->set_F(m, Fn.data, 1);
      func->set_F(m, Fnm1.data, 0);
      
      if(analysis_type == CM)
        continue;
      
      Tensor<2> eFnm1, eFn, pFnm1, pFn, pFnm1I, pFnI;
        
      func->get_pF(m,pFn.data,  1);
      func->get_pF(m,pFnm1.data,0);
      
      inv(pFn,   pFnI);
      inv(pFnm1, pFnm1I);          
      
      if(is_it_couple_w_thermal>=0){
        // compute temperature at the integration point
        double hFnp1[9];
        Tensor<2> hFn, hFnI, hFnm1, hFnm1I;
        compute_temperature_at_ip(&fe,grid,mat,T->T0,
                                         T->np1.m_pdata,T->n.m_pdata,T->nm1.m_pdata,
                                         hFnp1,hFn.data,hFnm1.data);
        inv(hFn,   hFnI);
        inv(hFnm1, hFnm1I);        
        eFn(i,j)   = Fn(i,k)*hFnI(k,l)*pFnI(l,j);
        eFnm1(i,j) = Fnm1(i,k)*hFnm1I(k,l)*pFnm1I(l,j);        
      }
      else{
        eFn(i,j) = Fn(i,k)*pFnI(k, j);
        eFnm1(i,j) = Fnm1(i,k)*pFnm1I(k, j);        
      }
                  
      double eJip_n   = ttl::det(eFn);
      double eJip_nm1 = ttl::det(eFnm1);
        
      V      += fe.detJxW;
      eJ_n   += fe.detJxW*eJip_n;
      eJ_nm1 += fe.detJxW*eJip_nm1;
      Up_n   += fe.detJxW*func->compute_dudj(m, eJip_n,   -1, 0);
      Up_nm1 += fe.detJxW*func->compute_dudj(m, eJip_nm1, -1, 0);
    }
    if(analysis_type == CM3F){
      if(fv->npres == 1){
        fv->tf.P_np1(eid+1, 1) = fv->tf.P_n(eid+1, 1) = Up_n/V;
        fv->tf.P_nm1(eid+1, 1) = Up_nm1/V;
      }
    
      fv->tf.V_np1(eid+1, 1) = fv->tf.V_n(eid+1, 1) = eJ_n/V;
      fv->tf.V_nm1(eid+1, 1) = eJ_nm1/V;
    }
  }
}