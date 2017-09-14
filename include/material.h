#ifndef PGFEM3D_MATERIAL_H
#define PGFEM3D_MATERIAL_H

/// Structure of material properties
struct Material {
  double Ex,Ey,Ez,Gyz,Gxz,Gxy,nyz,nxz,nxy,ax,ay,az,sig;
  /*Elastic stiffness*/
  double L[9];                                  //!< Local coordinate system
  double M[9];                                  //!< Local coordinate system

  /// potential function flags
  int devPotFlag;
  int volPotFlag;
};

struct MaterialThermal {
  double k[9];                //!< heat conductivity
  double cp;                  //!< heat capacity
  double FHS_MW;              //!< fraction_of_heat_sorce_due_to_mechanical_work
};

#endif // #ifndef PGFEM3D_MATERIAL_H
