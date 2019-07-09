//////////////////////////////////////////////////////////////////////
/// Define the PGFem3D data structs
///
/// Authors:
///   Sangmin Lee, University of Notre Dame <slee43 [at] nd.edu>
//////////////////////////////////////////////////////////////////////
#ifndef PGFEM3D_DATA_STRUCTURE_H
#define PGFEM3D_DATA_STRUCTURE_H

#include "bounding_element.h"
#include "cohesive_element.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "material.h"
#include "matgeom.h"
#include "mesh_load.h"
#include "node.h"
#include "sig.h"
#include "state_variables.h"
#include "pgfem3d/Solver.hpp"
#include "pgfem3d/Communication.hpp"
//#include "pgfem3d/MultiscaleCommon.hpp"

#include <sstream>
#include <vector>

static const constexpr int MAX_BND_ELEMENT_NO = 6; // maximum number of boundary element 
                                                    // for PGFem3D is 6 (HEXAHEDRAL)

/// Carrying Neumann boundary (NB) conditions for Multiphysics problems
class NeumannBoundaryElement{
  public:
    int nbe_no;               //! number of elements for NB
    int feature_no;           //! number of features for NB
    Matrix<int> element_ids;  //!< list of elements for NB
    Matrix<int> features;     //!< geometry features for NB
    Matrix<std::string> load; //!< NBC value
    Matrix<int> load_type;    //!< NBC type,    1: for load = t.N
                              //             nsd : load = [tx, ty, tz] nsd = [1,3]
    Matrix<Matrix<int>> bnd_elements; //!< list of boundary element for element in element_ids
    Matrix<double> load_values_np1; //!< load values compted from load at t(n+1)
    Matrix<double> load_values_n;   //!< load values compted from load at t(n) 
    Matrix<double> load_values_nm1; //!< load values compted from load at t(n-1)
    Matrix<int>    element_check;   //!< check for element loop 
    
    // set all zero
    NeumannBoundaryElement(){
      nbe_no     = 0;
      feature_no = 0;
    }
    
    // set size of feature related objects
    void set_feature_size(const int feature_no_in,
                          const int dof){
      feature_no = feature_no_in;
      features.initialization(feature_no, 2, 0);
      load.initialization(feature_no, dof);
      load_type.initialization(feature_no, 1);
      load_values_np1.initialization(feature_no, dof);
      load_values_n.initialization(  feature_no, dof);
      load_values_nm1.initialization(feature_no, dof);
    }
    
    // get T3D feature
    int bnd_feature(int id){return features(id, 0);}

    // get T3D feature ID
    int bnd_feature_id(int id){return features(id, 1);}

    // onstruct list of elemnets contributing for NBC
    void construct_list_of_boundary_elements(Element *elem,
                                             const int elemno,
                                             const int nsd);
    // print feature values
    void print_feature(void);
    
    // print element values
    void print_bnd_elements(void);
        
    // read math expression for load values
    void read_loads(const std::string &in_dir,
                    const std::string &load_name, 
                    const int feature_id,
                    const int load_id);
    // compute load values at t(n+1), t(n), and t(n-1_
    void compute_load_values(double tnp1,
                             double tn,
                             double tnm1);
};

/// Time stepping struct
/// Has time stepping information
struct TimeStepping {
  long nt;          //!< total number of times
  long tim;         //!< current time step number
  double *times;    //!< list of time
  double dt_n;      //!< dt at n
  double dt_np1;    //!< dt at n+1
  long *print;      //!< step numbers to be printed
  double *tns;      //!< time at n for multiple physics
};

/// Mesh
/// Has all mesh data
struct Grid {
  long Gnn;                  //!< global number of nodes
  long Gne;                  //!< global number of elements
  long Gnbndel;              //!< global number of boundary elements
  long Gn_be;                //!< global number of bounding elements
  long ne;                   //!< number of elements
  long nn;                   //!< number of nodes
  int nsd;                   //!< number of spatial dimension
  int n_be;                  //!< number of bounding elements
  long nce;                  //!< number of COEL (cohesive elements)
  long Gnce;                 //!< global number of cohesive elements
  Node *node;                //!< list of node
  Element *element;          //!< list of element
  BoundingElement *b_elems;  //!< list of bounding element
  COEL *coel;                //!< list of cohesive elements
  Matrix<NeumannBoundaryElement> NBE; //!< list of Neumann boundary condition element
                                      //!< as many as number of physics
};

/// Variables for three-field mixed method
class ThreeFieldVariables
{
  public:
  /// volume
  bool is_for_temporal;
  Matrix<double> V_np1;
  Matrix<double> V_nm1;
  Matrix<double> V_n;
  Matrix<double> dV, ddV;

  /// pressure
  Matrix<double> P_np1;
  Matrix<double> P_nm1;
  Matrix<double> P_n;
  Matrix<double> dP, ddP;

  ThreeFieldVariables(){is_for_temporal=false;};

  /// constructor
  void construct(const int ne,
                 const int Vno,
                 const int Pno,
                 bool is_temporal = false)
  {
    V_nm1.initialization(ne,Vno,1.0); // must be initialized to 1.0
      V_n.initialization(ne,Vno,1.0); // must be initialized to 1.0

    P_nm1.initialization(ne,Pno,0.0);
      P_n.initialization(ne,Pno,0.0);

    if(is_temporal == false)
    {
      V_np1.initialization(ne,Vno,1.0); // must be initialized to 1.0
         dV.initialization(ne,Vno,0.0);
        ddV.initialization(ne,Vno,0.0);

      P_np1.initialization(ne,Pno,0.0);
         dP.initialization(ne,Pno,0.0);
        ddP.initialization(ne,Pno,0.0);
    }
  };

  /// update variable at n-1 from n
  ///                 at n from n+1
  /// and reset increments to zeros
  void update_for_next_time_step(bool updated_Lagrangian = false)
  {
    if(is_for_temporal)
      return;

    for(int ia=0; ia<V_np1.m_row*V_np1.m_col; ia++)
    {
      V_nm1.m_pdata[ia] = V_n.m_pdata[ia];
      if(updated_Lagrangian)
        V_n.m_pdata[ia] *= V_np1.m_pdata[ia];
      else
        V_n.m_pdata[ia] = V_np1.m_pdata[ia];

      dV.m_pdata[ia] = 0.0;
      ddV.m_pdata[ia] = 0.0;
    }

    for(int ia=0; ia<P_np1.m_row*P_np1.m_col; ia++)
    {
      P_nm1.m_pdata[ia] = P_n.m_pdata[ia];
      P_n.m_pdata[ia] = P_np1.m_pdata[ia];
      dP.m_pdata[ia] = 0.0;
      ddP.m_pdata[ia] = 0.0;
    }
  };

  /// update dP and dV (increments) from ddP and ddV(increments of increments)
  /// \param[in] gamma if not converged smaller gamma will be applied, default = 1.0
  void update_increments_from_NR(double gamma = 1.0)
  {
    if(is_for_temporal)
      return;

    for(int ia=0; ia<V_np1.m_row*V_np1.m_col; ia++)
    {
      ddV.m_pdata[ia] *= gamma;
      dV.m_pdata[ia] += ddV.m_pdata[ia];
    }

    for(int ia=0; ia<P_np1.m_row*P_np1.m_col; ia++)
    {
      ddP.m_pdata[ia] *= gamma;
      dP.m_pdata[ia] += ddP.m_pdata[ia];
    }
  };

  /// update P_np1 and V_np1 from NR solution
  void update_np1_from_increments(void)
  {
    if(is_for_temporal)
      return;

    for(int ia=0; ia<V_np1.m_row*V_np1.m_col; ia++)
      V_np1.m_pdata[ia] += dV.m_pdata[ia];

    for(int ia=0; ia<P_np1.m_row*P_np1.m_col; ia++)
      P_np1.m_pdata[ia] += dP.m_pdata[ia];
  };

  /// if not converged, return to 0
  void reset(void)
  {
    if(is_for_temporal)
      return;

    for(int ia=0; ia<V_np1.m_row*V_np1.m_col; ia++)
      ddV.m_pdata[ia] = dV.m_pdata[ia] = 0.0;

    for(int ia=0; ia<P_np1.m_row*P_np1.m_col; ia++)
      ddP.m_pdata[ia] = dP.m_pdata[ia] = 0.0;
  };
};

/// struct for field variables
struct FieldVariablesTemporal {
  double *u_nm1;           //!< displacement at n-1
  double *u_n;             //!< displacement at n

  ThreeFieldVariables tf;  //!< variables for three-field mixed method

  int element_variable_no; //!< number of element variables
  State_variables *var;    //!< object to store element variables
};

/// struct for field variables
struct FieldVariables {
  bool   apply_initial_velocity; //!< if true, initial velocity is applied at the first time steps
                                 //!< Each process can have DIFFERNT value.
  double u0;      //!< reference value of field variables
  long Gndof;     //!< total number of degree freedom
  long ndofn;     //!< number of degree of freedom on a node
  long ndofd;     //!< number of degree of freedom in the domain
  long npres;     //!< number of pressure per element
  long nVol;      //!< number of volume per element
  long n_concentrations; //!< number of concentrations
  double *u_np1;  //!< displacement at n+1
  double *u_n;    //!< displacement at n
  double *u_nm1;  //!< displacement at n-1
  double *d_u;    //!< workspace for local increment of the solution n->n+1
  double *dd_u;   //!< workspace for local _iterative_ increment of the solution
  double *f;      //!< workspace for local residual
  double *R;      //!< [in] vector of Neumann loads (Incramental forces)
  double *f_defl; //!< workspace for the load vector due to derichlet conditions
  double *RR;     //!< [out] total Neumann load (Total forces for subdivided increment)
  double *f_u;    //!< workspace for load due to body force
  double *RRn;    //!< [in] Neumann load to time n (Total force after equiblirium)
  double pores;   //!< [out] opening volume of failed cohesive interfaces
  double *BS_x;   //!< workspace for the locally owned part of the global solution 'rr'
  double *BS_f;   //!< Global part of 'f'
  double *BS_f_u; //!< Global part of 'f_u'
  double *BS_RR;  //!< Global part of 'RR'
  double NORM;    //!< [out] residual of first iteration (tim = 0, iter = 0).
  SIG *sig;       //!< pointer for the stress
  EPS *eps;       //!< pointer for strain
  SIG *sig_n;     //!< smoothed stress
  int n_coupled;  //!< number of coupled physics
  int *coupled_physics_ids;     //!< array of phyiscs ids to be coupled
                                //!< it tells physics e.g.) fv.coupled_physics_ids[ib] == MULTIPHYSICS_MECHANICAL
                                //!<                        fv.coupled_physics_ids[ib] == MULTIPHYSICS_THERMAL
                                //!<                               :                         :
  FieldVariables **fvs; //!< array of FieldVariables pointers for multiphysics coupling
  FieldVariablesTemporal *temporal; //!< temporal space for transient time stepping
  State_variables *statv_list;        //!< list of state variables for constitutive model interface
  double subdivision_factor_n;        //!< use for linearly map subdivided parameters at t(n)
                                      //!<   v = v_n + (v_np1 - v_n)*subdivision_factor_n
  double subdivision_factor_np1;      //!< use for linearly map subdivided parameters at t(n+1)
                                      //!<   v = v_n + (v_np1 - v_n)*subdivision_factor_np1
  ThreeFieldVariables tf;  //!< variables for three-field mixed method
};

/// struct for field variables
struct FieldVariablesThermal {
  long Gndof;    //!< total number of degree freedom
  long ndofn;    //!< number of degree of freedom on a node
  long ndofd;    //!< number of degree of freedom in the domain
  double *T_np1; //!< displacement at n+1
  double *T_n;   //!< displacement at n
  double *T_nm1; //!< displacement at n-1
  double *dT;    //!< workspace for local increment of the temperature solution n->n+1
  double *dd_T;  //!< workspace for local _iterative_ increment of the solution
  double NORM;   //!< [out] residual of first iteration (tim = 0, iter = 0).
};

/// struct for material properties
struct MaterialProperty {
  double           *density; //!< list of material density
  Material           *mater; //!< list of material properites (Mechanical)
  MaterialThermal  *thermal; //!< list of material properites (Thermal)

  HOMMAT *hommat;  //!< list of homogeneous material properites
  MATGEOM matgeom; //!< information related to material geometry (for crystal plasticity)
  long nhommat;    //!< number of homogeneous materials
  long nmat;       //!< number of materials
  long n_orient;   //!< number of orientations
  int n_co_props;  //!< number of cohesive material properites
  cohesive_props *co_props; //!< list of cohesive material properites
};

/// struct for the boundary conditions
struct LoadingSteps {
  SUPP *sups;
  double **sup_defl;  //!< sum of Dirichlet BC increments to step n
  long nln;           //!< number of nodes with loads
  long nle_s;         //!< number of surface element with loads
  long nle_v;         //!< number of volume element with loads
  ZATNODE *znod;      //!< list of nodes with loads
  ZATELEM *zele_s;    //!< list of surface element with loads
  ZATELEM *zele_v;    //!< list of volume element with loads
  long **tim_load;    //!< list of time steps to be saved
  FILE **solver_file; //!< file pointer for reading loads increments
};

/// for setting physics ids
enum MultiphysicsAnalysis {
  MULTIPHYSICS_MECHANICAL,
  MULTIPHYSICS_THERMAL,
  MULTIPHYSICS_CHEMICAL,
  MULTIPHYSICS_NO
};

/// struct for setting multiphysics
struct Multiphysics {
  int physicsno;      //!< number of physics
  char **physicsname; //!< physics names
  int *physics_ids;   //!< physics ids
  int *ndim;          //!< degree of feedom of the physics
  int *write_no;      //!< number of variables to be written as results
  int total_write_no; //!< total number of variables to be written as results
  std::vector<std::vector<int>> write_ids;    //!< index of physical variables to be written
  std::vector<std::vector<int>> coupled_ids;  //!< coupled physics id
};

/// initialize time stepping variable
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int time_stepping_initialization(TimeStepping *ts);

/// destruct time stepping variable
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int destruct_time_stepping(TimeStepping *ts);

/// initialize mesh object
///
/// \param[in, out] grid an object containing all mesh data
/// \return non-zero on internal error
int grid_initialization(Grid *grid);

/// destruct of mesh
///
/// \param[in, out] grid an object containing all mesh data
/// \param[in] mp multiphysics object
/// \return non-zero on internal error
int destruct_grid(Grid *grid,
                  const PGFem3D_opt *opts,
                  const Multiphysics& mp);

/// initialize field variables
///
/// \param[in, out] fv an object containing all field variables
/// \return non-zero on internal error
int field_varialbe_initialization(FieldVariables *fv);

/// construct field variables
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \param[in] com an object for communication
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id physics id
/// \return non-zero on internal error
int construct_field_varialbe(FieldVariables *fv,
                             Grid *grid,
                             pgfem3d::CommunicationStructure *com,
                             const PGFem3D_opt *opts,
                             const Multiphysics& mp,
                             int myrank,
                             int mp_id);

/// destruct field variables
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id physics id
/// \return non-zero on internal error
int destruct_field_varialbe(FieldVariables *fv,
                            Grid *grid,
                            const PGFem3D_opt *opts,
                            const Multiphysics& mp,
                            int mp_id);

/// initialize field variables thermal part
///
/// \param[in, out] fv an object containing all field variables for thermal
/// \return non-zero on internal error
int thermal_field_varialbe_initialization(FieldVariablesThermal *fv);

/// prepare temporal varialbes for staggering Newton Raphson iterations
///
/// Before call this function, physics coupling should be defined in
/// fv->n_coupled and fv->coupled_physics_ids
///
/// \param[in, out] fv an object containing all field variables for thermal
/// \param[in] grid an object containing all mesh data
/// \param[in] is_for_Mechanical if yes, prepare constitutive models
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int prepare_temporal_field_varialbes(FieldVariables *fv,
                                     Grid *grid,
                                     int is_for_Mechanical,
                                     const PGFem3D_opt *opts);

/// destory temporal varialbes for staggering Newton Raphson iterations
///
/// should be called before destroying fv
///
/// \param[in, out] fv an object containing all field variables for thermal
/// \param[in] is_for_Mechanical if yes, prepare constitutive models
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int destory_temporal_field_varialbes(FieldVariables *fv,
                                     int is_for_Mechanical,
                                     const PGFem3D_opt *opts);

/// initialize material properties
///
/// \param[in, out] mat an object containing all material parameters
/// \return non-zero on internal error
int material_initialization(MaterialProperty *mat);

/// destruct material properties
///
/// \param[in, out] mat an object containing all material parameters
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int destruct_material(MaterialProperty *mat, const PGFem3D_opt *opts);

/// initialize iterative solver object
///
/// \param[in, out] sol an object containing data for linear solver
/// \return non-zero on internal error
int solution_scheme_initialization(pgfem3d::Solver *sol);

/// initialize loading steps object
///
/// \param[in, out] load an object containing boundary increments
/// \return non-zero on internal error
int loading_steps_initialization(LoadingSteps *load);

/// construct loading steps object
///
/// \param[in, out] load an object containing boundary increments
/// \param[in] mp multiphysics object
/// \return non-zero on internal error
int construct_loading_steps(LoadingSteps *load, const Multiphysics& mp);

/// destruct loading steps object
///
/// \param[in, out] load an object containing boundary increments
/// \param[in] mp multiphysics object
/// \return non-zero on internal error
int destruct_loading_steps(LoadingSteps *load, const Multiphysics& mp);

/// initialize communication structures
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int communication_structure_initialization(pgfem3d::CommunicationStructure *com);

/// destruct communication structures
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int destruct_communication_structure(pgfem3d::CommunicationStructure *com);

/// initialize multiphysics object
///
/// \param[in, out] mp an object for multiphysics stepping
/// \return non-zero on internal error
int multiphysics_initialization(Multiphysics& mp);

/// construct multiphysics object
///
/// \param[in, out] mp an object for multiphysics stepping
/// \param[in] physicsno number of physics
/// \return non-zero on internal error
int construct_multiphysics(Multiphysics& mp,
                           int physicsno);

/// destruct multiphysics object
///
/// \param[in, out] mp an object for multiphysics stepping
/// \return non-zero on internal error
int destruct_multiphysics(Multiphysics& mp);

/// read and construct multiphysics
///
/// \param[in, out] mp an object for multiphysics stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_multiphysics_settings(Multiphysics& mp,
                               const PGFem3D_opt *opts,
                               int myrank);

#endif // #define PGFEM3D_DATA_STRUCTURE_H
