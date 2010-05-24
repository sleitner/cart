#ifndef __ANALYSIS_XRAY_H__
#define __ANALYSIS_XRAY_H__

/* Chandra parameters for X-ray spectroscopic fitting (Vikhlinin 2006) */
#define alpha_xray      (0.875)
#define beta_xray       (1.0)
#define delta_xray_1    (0.19)
#define delta_xray_2    (0.25)

/* parameters for ASCA Lbol and Tsl (Bartlemann & Steinmetz 1996, Mazzotta 2004)*/
#define Tsl_alpha	(0.75)
#define xray_cj		(2.42e-24)
#define eband_min0	(0.7)
#define eband_max0	(7.0)

extern double eband_min;
extern double eband_max;

void init_xray_tables();
void set_xray_redshift( double aexp );
void xray_calibration( double Tg, double *cT, double *lambda, double *fT );
double xray_calibrated_line_temperature( double avgE );

#endif /* ANALYSIS_XRAY */
