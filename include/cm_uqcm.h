/* HEADER */
/**
 * This file defines the interface for the UQ study interface through constitutive mode interface
 * Sangmin Lee, University of Notre Dame, <slee43@nd.edu>
 */
#pragma once
#ifndef H__H__UQ_THROUGH_CM_H
#define H__H__UQ_THROUGH_CM_H

typedef struct MATERIAL_ELASTICITY MATERIAL_ELASTICITY;

/**
 * Construct and initialize the model context for calling functions
 * through the constitutive modeling interface.
 *
 * \param[in] m_in - handle to input material properties
 * \param[out] m_out - handle to output material properties
 * \param[in] x,y,z - coordinate at the integration point
 * \return non-zero on internal error.
 */
int material_properties_elasticity_at_ip(MATERIAL_ELASTICITY *m_in, MATERIAL_ELASTICITY *m_out, 
                                         double x, double y, double z);

#endif
