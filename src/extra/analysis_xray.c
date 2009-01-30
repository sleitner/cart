#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "analysis_xray.h"
#include "auxiliary.h"
#include "units.h"
#include "constants.h"

#ifdef ANALYSIS

#define nltmax		500
#define nlzmax		500

int nlt, nlz;
double tlmin, tlmax, dlt; /* temperature */
double zlmin, zlmax, dlz; /* redshift */
double dlzi, dlti;

double cT_table[nlzmax][nltmax];
double cT_fixed[nltmax];
double lambda_table[nlzmax][nltmax];
double lambda_fixed[nltmax];
double fT_table[nlzmax][nltmax];
double fT_fixed[nltmax];

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

	fscanf( data, "%lf %lf %lf %u\n", &zlmin, &zlmax, &dlz, &nlz );
	fscanf( data, "%lf %lf %lf %u\n", &tlmin, &tlmax, &dlt, &nlt );

	cart_debug("redshift binning: %f %f %f %u", zlmin, zlmax, dlz, nlz );
	cart_debug("temperature binning: %f %f %f %u", tlmin, tlmax, dlt, nlt );

	dlti = 1.0 / dlt;
	dlzi = 1.0 / dlz;

	for ( ilz = 0; ilz < nlz; ilz++ ) {
		for ( ilt = 0; ilt < nlt; ilt++ ) {
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
	int irs, irs1, irs2, ilz, ilt;

	rs = max( 1.0 / aexp - 1.0, 0.0 );

	/* find redshift bin */
	irs = (int)((log10(rs) - zlmin)*dlzi);
	irs1 = max(irs,0);
	irs1 = min(irs1, nlz-1);
	irs2 = min(irs+1,nlz-1);
	irs2 = max(irs2,0);

	if ( irs1 == irs2 ) {
		/* just copy over specific redshift */
		for ( ilt = 0; ilt < nlt; ilt++ ) {
			cT_fixed[ilt] = cT_table[irs1][ilt];	
			lambda_fixed[ilt] = lambda_table[irs1][ilt];	
			fT_fixed[ilt] = fT_table[irs1][ilt];
		}
	} else {
		/* need to interpolate to specific redshift */
		rs1 = zlmin + dlz*(double)irs1;
		rs2 = zlmin + dlz*(double)irs2;

		for ( ilt = 0; ilt < nlt; ilt++ ) {
			ac = ( cT_table[irs2][ilt] - cT_table[irs1][ilt] ) / (rs2 - rs1);
			bc = cT_table[irs1][ilt] - ac * rs1;

			cT_fixed[ilt] = ac*rs + bc;

			ac = ( lambda_table[irs2][ilt] - lambda_table[irs1][ilt] ) / (rs2 - rs1);
			bc = lambda_table[irs1][ilt] - ac * rs1;

			lambda_fixed[ilt] = ac*rs + bc;

			ac = ( fT_table[irs2][ilt] - fT_table[irs1][ilt] ) / (rs2 - rs1);
			bc = fT_table[irs1][ilt] - ac * rs1;

			fT_fixed[ilt] = ac*rs + bc;
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

	if ( Tlog < tlmin ) {
		*cT = 0.0;
		*lambda = 0.0;
		*fT = 0.0;
	} else { 
		it1 = (int)((Tlog - tlmin)*dlti);
		it2 = it1 + 1;

		it1 = max(it1,0);
		it1 = min(it1,nlt-1);
		it2 = max(it2,0);
		it2 = min(it2,nlt-1);

		/* set up interpolation variables */
		td = tlmin + dlt * (double)(it1+1);
		t1 = (td - Tlog) * dlti;
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
	while ( it2 < nlt && fT_fixed[it2] < avgE ) {
		it2++;
	}
	it1 = it2-1;

	logT = ( avgE - fT_fixed[it1] ) * dlt / ( fT_fixed[it2] - fT_fixed[it1] ) + (tlmin + (double)(it1)*dlt);
	return pow( 10.0, logT );
}

#endif /* ANALYSIS */
