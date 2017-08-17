/********************
 *  Solvers         *
 *  Jaroslav Kruis  *
 ********************/

#ifndef RESICE_H
#define RESICE_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

void reseni_rovnic(double *a,
           double *x,
           double *y,
           long n,
           long m,
           long as);

int resic_sky (double *a,
           double *x,
           double *y,
           long *adr,
           long n,
           long tc,
           const int analysis);

void lokalizace_scr (double *a,
             double *b,
             long *lcn,
             long *adr,
             long *ci,
             long n);

void strci_scr (long *adr,
        long ne,
        long ndofn,
        ELEMENT *elem,
        NODE *node,
        const int mp_id);

void aci_scr (long *ci,
          long *adr,
          long ne,
          long ndofn,
          ELEMENT *elem,
          NODE *node,
          const int mp_id);

void sort_cr (long *ci,
          long *adrb,
          long *adre,
          long ndof);

long size (long *adrb,
       long *adre,
       long n);

void colindex_cr (long *ci,
          long *adr,
          long *aci,
          long *adre,
          long ndof);

void mv_scr (double *a,
         double *b,
         double *c,
         long *adr,
         long *ci,
         long n);

void minimize_cr (double *a,
          long *b,
          long *adrn,
          long *adro,
          long n,
          double limit);

double energie_scr (double *a,
            double *x,
            double *b,
            long *adr,
            long *ci,
            long n);

void cg_scr (double *a,
         double *x,
         double *y,
         long *adr,
         long *ci,
         long n,
         long ni,
         double err,
         long *ani,
         double *ares,
         double limit,
         long iv);

long nonzero_sky (double *a,
          long *adr,
          double limit,
          long n);

void sky_scr (double *sky,
          long *adrs,
          double *scr,
          long *adrc,
          long *ci,
          double limit,
          long n);

#endif
