#ifndef __ANALYSIS_XRAY_H__
#define __ANALYSIS_XRAY_H__

void init_xray_tables();
void set_xray_redshift( double aexp );
void xray_calibration( double Tg, double *cT, double *lambda, double *fT );
double xray_calibrated_line_temperature( double avgE );

#endif /* ANALYSIS_XRAY */
