#ifndef PGFEM3D_ENSIGHT_H
#define PGFEM3D_ENSIGHT_H

struct Ensight {
  long ncn = 0,
    *Sm = nullptr,
    *Sp = nullptr,
    *No = nullptr,
    NVp = 0,
    NCp = 0,
    *Vp = nullptr,
    *Cp = nullptr;

  ~Ensight();
};

#endif // #define PGFEM3D_ENSIGHT_H
