#ifndef __EXT_XRAYS_H__
#define __EXT_XRAYS_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/* Chandra parameters for X-ray spectroscopic fitting (Vikhlinin 2006) */
#define alpha_xray      (0.875)
#define beta_xray       (1.0)
#define delta_xray_1    (0.19)
#define delta_xray_2    (0.25)

void init_xray_tables();
void set_xray_redshift( double aexp );
void xray_calibration( double Tg, double *cT, double *lambda, double *fT );
double xray_calibrated_line_temperature( double avgE );

#endif /* EXT_XRAYS */
