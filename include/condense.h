#ifndef _H_CONDENSE_H_
#define _H_CONDENSE_H_

void condense_Fupt_to_Fup(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, double *Kut, double *Ktt, double *Kpt);
void condense_Fupt_to_Fu(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, double *Kup, double *Ktp, double *Ktt,double *Kpt);
void condense_K2_to_K1(double *K11, int nne, int nsd, int npres,
                   double *Kuu, double *Kup,
                   double *Kpu, double *Kpp);
void condense_K3_to_K2(double *K11, double *K12, double *K21, double *K22, 
                   int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp);
void condense_Kupt_to_Ku(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp); 
void condense_Kupt_to_Kup(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp);
void condense_K_out(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp);
#endif