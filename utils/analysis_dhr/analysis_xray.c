#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "analysis_xray.h"
#include "auxiliary.h"
#include "units.h"
#include "constants.h"

#ifdef ANALYSIS

#define xnltmax		500
#define xnlzmax		500

int xnlt, xnlz;
double xtlmin, xtlmax, xdlt; /* temperature */
double xzlmin, xzlmax, xdlz; /* redshift */
double xdlzi, xdlti;

double eband_min;
double eband_max;

double cT_table[xnlzmax][xnltmax];
double cT_fixed[xnltmax];
double lambda_table[xnlzmax][xnltmax];
double lambda_fixed[xnltmax];
double fT_table[xnlzmax][xnltmax];
double fT_fixed[xnltmax];

void init_xray_tables() {
	int i;
	FILE *data;
	int ilz, ilt;
	double z, T;

	cart_debug("loading chandra calibration...");
	
	data = fopen("chandra_calibration.dat", "r");
	if ( data == NULL ) {
		cart_error("Unable to open chandra_calibration.dat");
	}

	fscanf( data, "%lf %lf %lf %u\n", &xzlmin, &xzlmax, &xdlz, &xnlz );
	fscanf( data, "%lf %lf %lf %u\n", &xtlmin, &xtlmax, &xdlt, &xnlt );

	xdlti = 1.0 / xdlt;
	xdlzi = 1.0 / xdlz;

	for ( ilz = 0; ilz < xnlz; ilz++ ) {
		for ( ilt = 0; ilt < xnlt; ilt++ ) {
			fscanf( data, "%lf %lf %lf %lf %lf\n", &z, &T, 
					&cT_table[ilz][ilt], 
					&lambda_table[ilz][ilt],
					&fT_table[ilz][ilt] );
		}
	}

	fclose(data);
}

void set_xray_redshift( double aexp ) {
	double rs, rs1, rs2, ac, bc;
	int irs, irs1, irs2, ilt;

	rs = max( 1.0 / aexp - 1.0, 0.0 );

	eband_min = eband_min0 / aexp;
	eband_max = eband_max0 / aexp;

	/* find redshift bin */
	if ( log10(rs) < xzlmin ) {
		irs1 = irs2 = 0;
	} else {
		irs = (int)((log10(rs) - xzlmin)*xdlzi);
		irs1 = max(irs,0);
		irs1 = min(irs1, xnlz-1);
		irs2 = min(irs+1,xnlz-1);
		irs2 = max(irs2,0);
	}

	if ( irs1 == irs2 ) {
		/* just copy over specific redshift */
		for ( ilt = 0; ilt < xnlt; ilt++ ) {
			cT_fixed[ilt] = cT_table[irs1][ilt];	
			lambda_fixed[ilt] = lambda_table[irs1][ilt];
			fT_fixed[ilt] = fT_table[irs1][ilt];
		}
	} else {
		/* need to interpolate to specific redshift */
		rs1 = xzlmin + xdlz*(double)irs1;
		rs2 = xzlmin + xdlz*(double)irs2;

		for ( ilt = 0; ilt < xnlt; ilt++ ) {
			ac = ( cT_table[irs2][ilt] - cT_table[irs1][ilt] ) / (rs2 - rs1);
			bc = cT_table[irs1][ilt] - ac * rs1;

			cT_fixed[ilt] = ac*log10(rs) + bc;
			cart_assert( cT_fixed[ilt] >= 0.0 );

			ac = ( lambda_table[irs2][ilt] - lambda_table[irs1][ilt] ) / (rs2 - rs1);
			bc = lambda_table[irs1][ilt] - ac * rs1;

			lambda_fixed[ilt] = ac*log10(rs) + bc;
			cart_assert( lambda_fixed[ilt] >= 0.0 );

			ac = ( fT_table[irs2][ilt] - fT_table[irs1][ilt] ) / (rs2 - rs1);
			bc = fT_table[irs1][ilt] - ac * rs1;

			fT_fixed[ilt] = ac*log10(rs) + bc;
			cart_assert( fT_fixed[ilt] >= 0.0 );
		}
	}
}

void xray_calibration( double Tg, double *cT, double *lambda, double *fT ) {
	double Tlog;
	int it1, it2;
	double td;
	double t1,d1;

	/* compute temperature bin */
	Tlog = log10(Tg);

	if ( Tlog < xtlmin || Tlog > xtlmax ) {
		*cT = 0.0;
		*lambda = 0.0;
		*fT = 0.0;
	} else { 
		it1 = (int)((Tlog - xtlmin)*xdlti);
		it2 = it1 + 1;

		it1 = max(it1,0);
		it1 = min(it1,xnlt-1);
		it2 = max(it2,0);
		it2 = min(it2,xnlt-1);

		/* set up interpolation variables */
		td = xtlmin + xdlt * (double)(it1+1);
		t1 = (td - Tlog) * xdlti;
		d1 = 1.0 - t1;

		*cT = t1*cT_fixed[it1] + d1*cT_fixed[it2];
		*lambda = t1*lambda_fixed[it1] + d1*cT_fixed[it2];
		*fT = t1*fT_fixed[it1] + d1*fT_fixed[it2];
	}
}

double xray_calibrated_line_temperature( double avgE ) {
	int it1, it2;
	double logT;

	/* invert fT(T) to get Tline */
	it2 = 1;
	while ( it2 < xnlt && fT_fixed[it2] < avgE ) {
		it2++;
	}
	it1 = it2-1;

	logT = ( avgE - fT_fixed[it1] ) * xdlt / ( fT_fixed[it2] - fT_fixed[it1] ) + (xtlmin + (double)(it1)*xdlt);
	return pow( 10.0, logT );
}

#endif /* ANALYSIS */
