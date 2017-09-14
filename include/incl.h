/* HEADER */
#pragma once
#ifndef PGFEM3D_INCL_H
#define PGFEM3D_INCL_H

#include "data_structure.h"
#include "element.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "allocation.h"

void build_elem_inelas (long ne,
                        Element *elem);

void build_pressure_nodes (long ne,
                           long npres,
                           Element *elem,
                           SIG *sig,
                           EPS *eps,
                           const int analysis);

void build_crystal_plast (long ne,
                          Element *elem,
                          SIG *sig,
                          EPS *eps,
                          CRPL *crpl,
                          const int analysis,
                          const int plc);

void nulld (double *a,
            long n);

void nulld2 (double **a,
             long m,
             long n);

#endif // PGFEM3D_INCL_H
