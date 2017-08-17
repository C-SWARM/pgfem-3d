/* HEADER */
#pragma once
#ifndef INCL_H
#define INCL_H

#include "element.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "allocation.h"

void build_elem_inelas (long ne,
            ELEMENT *elem);

void build_pressure_nodes (long ne,
               long npres,
               ELEMENT *elem,
               SIG *sig,
               EPS *eps,
               const int analysis);

void build_crystal_plast (long ne,
              ELEMENT *elem,
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

#endif
