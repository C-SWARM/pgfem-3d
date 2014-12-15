#ifndef _THREE_FIELD_ELEMENT_H_
#define _THREE_FIELD_ELEMENT_H_

#ifdef __cplusplus
extern "C" {
#endif

void add_3F_Kuu_ip(double *K,
        int nne, double *ST, double *F, double jj, double wt, double Pn, double Tn,
        double dt, double alpha);
        
void add_3F_Kut_ip(double *K,
        int nne, int nVol, double *ST, double *F, double jj, double wt, double *Nt,
        double dt, double alpha);
        
void add_3F_Kup_ip(double *K,
        int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
        double dt, double alpha);
        
void add_3F_Ktu_ip(double *K,
        int nne, int nVol,
        double *ST, double *F, double jj, double wt, double *Nt,
        double dt, double alpha);

void add_3F_Ktt_ip(double *K, int nVol, double jj, double wt, double *Nt, double Upp,
        double dt, double alpha);
        
void add_3F_Ktp_ip(double *K,
        int nVol, int npres, double jj, double wt, double *Nt, double *Np,
        double dt, double alpha);

void add_3F_Kpu_ip(double *K,
        int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
        double dt, double alpha);

void add_3F_Kpt_ip(double *K,
        int nVol, int npres, double jj, double wt, double *Nt, double *Np,
        double dt, double alpha);

void resid_w_inertia_Ru_ip(double *fu,
        int nne, double *ST, double *F, double jj, double wt, double Pn);

void resid_w_inertia_Rt_ip(double *ft, int nVol, double jj, double wt, double *Nt, double Pn, double Up);

void resid_w_inertia_Rp_ip(double *fp, int npres, double *F, double jj, double wt, double *Np, double Tn);

#ifdef __cplusplus
}
#endif

#endif