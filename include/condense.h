#ifndef _H_CONDENSE_H_
#define _H_CONDENSE_H_

void condense_K_out(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp);
void condense_F_out(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, double *Kut, double *Kup, double *Ktp, double *Ktt,double *Kpt);

#endif
