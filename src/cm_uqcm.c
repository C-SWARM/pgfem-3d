#include "cm_uqcm.h"
#include "material_properties.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif
int material_properties_elasticity_at_ip(MATERIAL_ELASTICITY *m_in, MATERIAL_ELASTICITY *m_out, 
                                         double x, double y, double z)
{
  int err = 0;
  double values = 1.0 + 0.5*sin(y*2.0*M_PI);
  double E = (m_in->E)*values;
  err += set_properties_using_E_and_nu(m_out, E, m_in->nu);
  
  m_out->m01 = (m_in->m01)*values;
  m_out->m10 = (m_in->m10)*values;
  m_out->volPotFlag = m_in->volPotFlag;
  m_out->devPotFlag = m_in->devPotFlag;
  return err;
}