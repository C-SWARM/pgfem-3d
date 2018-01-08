#ifndef _H_CONDENSE_H_
#define _H_CONDENSE_H_

void condense_K_out(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp);
void condense_F_out(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, double *Kut, double *Kup, double *Ktp, double *Ktt,double *Kpt);


int condense_K_3F_to_1F(double *Ks, int nne, int nsd, int Pno, int Vno,
                        double *Kuu_in, double *Kut_in, double *Kup_in,
                        double *Ktu_in, double *Ktt_in, double *Ktp_in,
                        double *Kpu_in, double *Kpt_in, double *Kpp_in);

int condense_F_3F_to_1F(double *fe, int nne, int nsd, int Pno, int Vno,
                        double *fu_in, double *ft_in, double *fp_in, 
                        double *Kut_in, double *Kup_in, double *Ktp_in, double *Ktt_in,double *Kpt_in);

#endif
