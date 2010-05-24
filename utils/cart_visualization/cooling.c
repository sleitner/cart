#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "auxiliary.h"
#include "cooling.h"
#include "units.h"
#include "constants.h"
#include "io.h"

#ifdef COOLING

#define smallrate	(1e-30)
#define nltmax		71
#define nldmax		11
#define nlzmax		9
#define nrsmax		21

int nlt, nlz, nld, nrs;
double tlmin, tlmax, dlt; /* temperature */
double dlmin, dlmax, dld; /* density */
double Zlmin, Zlmax, dlZ; /* metallicity */
double rsmin, rsmax, drs; /* scale factor */
double dlti, dldi, dlZi, drsi;

double coolcl[nrsmax][nlzmax][nldmax][nltmax];
double ccl_rs[nlzmax][nldmax][nltmax];
double f_ion[nrsmax][nlzmax][nldmax][nltmax];

void init_cooling() {
	int i;
	FILE *data;
	int irs, ilz, ild, ilt;
	double d[9];
	double cdum, hdum, ct_crit;

#ifdef CLOUDY_COOLING

	data = fopen("clcool.dat", "r");
	if ( data == NULL ) {
		cart_error("Unable to open clcool.dat");
	}

	fscanf( data, "%lf %lf %lf %u\n", &tlmin, &tlmax, &dlt, &nlt );
	fscanf( data, "%lf %lf %lf %u\n", &dlmin, &dlmax, &dld, &nld );
	fscanf( data, "%lf %lf %lf %u\n", &Zlmin, &Zlmax, &dlZ, &nlz );
	fscanf( data, "%lf %lf %lf %u\n", &rsmin, &rsmax, &drs, &nrs );

	dlti = 1.0 / dlt;
	dldi = 1.0 / dld;
	dlZi = 1.0 / dlZ;
	drsi = 1.0 / drs;

	for ( irs = 0; irs < nrs; irs++ ) {
		for ( ilz = 0; ilz < nlz; ilz++ ) {
			for ( ild = 0; ild < nld; ild++ ) {

				ct_crit = 0.0;
				for ( ilt = 0; ilt < nlt; ilt++ ) {
					fscanf( data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &d[0],
						&d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8],
						&cdum, &hdum );

					if ( d[0] >= 3.2 && ct_crit == 0.0 ) {
						ct_crit = cdum;
					}

					cdum = max( cdum, smallrate );
					cdum = max( cdum, ct_crit );
					hdum = max( hdum, smallrate );

					coolcl[irs][ilz][ild][ilt] = (cdum-hdum) * 1e23 * AL_SD;
					f_ion[irs][ilz][ild][ilt] = d[5] / pow( 10.0, d[1] );
				}
			}
		}
	}

	fclose(data);

#endif /* CLOUDY_COOLING */

#ifdef SD93_COOLING 
	#error "Sutherland & Dopita cooling curves not implemented yet!"
#endif

}

void set_cooling_redshift( double aexp ) {
	double rs, rs1, rs2, ac, bc;
	int irs, irs1, irs2, ilz, ild, ilt;

	rs = max( 1.0 / aexp - 1.0, 0.0 );

	/* find redshift bin */
	irs = (int)((rs - rsmin)*drsi);
	irs1 = max(irs,0);
	irs1 = min(irs1, nrs-1);
	irs2 = min(irs+1,nrs-1);
	irs2 = max(irs2,0);

	if ( irs1 == irs2 ) {
		/* just copy over specific redshift */
		for ( ilz = 0; ilz < nlz; ilz++ ) {
			for ( ild = 0; ild < nld; ild++ ) {
				for ( ilt = 0; ilt < nlt; ilt++ ) {
					ccl_rs[ilz][ild][ilt] = coolcl[irs1][ilz][ild][ilt];
				}
			}
		}
	} else {
		/* need to interpolate to specific redshift */
		rs1 = rsmin + drs*(double)irs1;
		rs2 = rsmin + drs*(double)irs2;

		for ( ilz = 0; ilz < nlz; ilz++ ) {
			for ( ild = 0; ild < nld; ild++ ) {
				for ( ilt = 0; ilt < nlt; ilt++ ) {
					ac = ( coolcl[irs2][ilz][ild][ilt] - coolcl[irs1][ilz][ild][ilt] ) /
						(rs2 - rs1);
					bc = coolcl[irs1][ilz][ild][ilt] - ac * rs1;

					ccl_rs[ilz][ild][ilt] = ac*rs + bc;
				}
			}
		}
	}
}

double cooling_rate( double rhogl, double T_g, double Z_met ) {
	double Tlog;
	int it1, it2;
	int id1, id2;
	int iz1, iz2;
	double dd,td,zd;
	double d1,d2,d3,t1,t2,t3;

#ifdef CLOUDY_COOLING
	/* compute temperature bin */
	Tlog = log10(T_g) + 4.0;
	it1 = (int)((Tlog - tlmin)*dlti);
	it2 = it1 + 1;

	it1 = max(it1,0);
	it1 = min(it1,nlt-1);
	it2 = max(it2,0);
	it2 = min(it2,nlt-1);

	/* compute density bin */
	id1 = (int)((rhogl - dlmin)*dldi);
	id2 = id1 + 1;

	id1 = max(id1,0);
	id1 = min(id1,nld-1);
	id2 = max(id2,0);
	id2 = min(id2,nld-1);

	/* compute metallicity bin */
#ifdef METALCOOLING
	iz1 = (int)((Z_met - Zlmin)*dlZi);
	iz2 = iz1 + 1;
#else
	iz1 = 0;
	iz2 = 0;
#endif /* METALCOOLING */

	iz1 = max(iz1,0);
	iz1 = min(iz1,nlz-1);
	iz2 = max(iz2,0);
	iz2 = min(iz2,nlz-1);

	/* set up interpolation variables */
	td = tlmin + dlt * (double)(it1+1);
	dd = dlmin + dld * (double)(id1+1);
	zd = Zlmin + dlZ * (double)(iz1+1);
	t1 = (td - Tlog) * dlti;
	d1 = 1.0 - t1;
	t2 = (dd - rhogl) * dldi;
	d2 = 1.0 - t2;
	t3 = (zd - Z_met) * dlZi;
	d3 = 1.0 - t3;

/* initial metallicity fails this interpolation
	cart_assert( t1 >= 0.0 && t1 <= 1.0 );
	cart_assert( t2 >= 0.0 && t2 <= 1.0 );
	cart_assert( t3 >= 0.0 && t3 <= 1.0 );
*/

	return	t1*t2*t3 * ccl_rs[iz1][id1][it1] +
		d1*t2*t3 * ccl_rs[iz1][id1][it2] +
		t1*d2*t3 * ccl_rs[iz1][id2][it1] +
		d1*d2*t3 * ccl_rs[iz1][id2][it2] +
		t1*t2*d3 * ccl_rs[iz2][id1][it1] +
		d1*t2*d3 * ccl_rs[iz2][id1][it2] +
		t1*d2*d3 * ccl_rs[iz2][id2][it1] +
		d1*d2*d3 * ccl_rs[iz2][id2][it2];
#endif /* CLOUDY_COOLING */	
}

#endif /* COOLING */


