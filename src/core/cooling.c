#include "config.h"
#if defined(COOLING) && !defined(RADIATIVE_TRANSFER)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cooling.h"
#include "io.h"
#include "units.h"
#include "timing.h"

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

cooling_t coolcl[nrsmax][nlzmax][nldmax][nltmax];
cooling_t ccl_rs[nlzmax][nldmax][nltmax];
double fion[nrsmax][nlzmax][nldmax][nltmax];
double fion_rs[nlzmax][nldmax][nltmax];

char cloudy_table_filename[CONTROL_PARAMETER_STRING_LENGTH] = "clcool.dat";

void config_init_cooling() {
	control_parameter_add(control_parameter_string,cloudy_table_filename,"cloudy-table-file","Full pathname to file which contains cloudy cooling rates (default clcool.dat).");
}

void config_verify_cooling() {
	FILE *f;
	f = fopen( cloudy_table_filename, "r" );
	if ( f == NULL ) {
		cart_error("Unable to locate cloudy cooling table file %s", cloudy_table_filename );
	}
	fclose(f);	
}

void init_cooling() {
	FILE *data;
	int irs, ilz, ild, ilt;
	double d[9];
	double cdum, hdum, ct_crit;

	data = fopen( cloudy_table_filename, "r");
	if ( data == NULL ) {
		cart_error("Unable to open cloudy cooling table file %s", cloudy_table_filename );
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

#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
					coolcl[irs][ilz][ild][ilt] = cdum - hdum;
#else
					coolcl[irs][ilz][ild][ilt].Cooling = max(0.0,cdum);
					coolcl[irs][ilz][ild][ilt].Heating = max(0.0,hdum);
#endif
					fion[irs][ilz][ild][ilt] = d[5] / pow( 10.0, d[1] );
				}
			}
		}
	}

	fclose(data);
}

void set_cooling_redshift( double auni ) {
	double rs, rs1, rs2, ac, bc;
	int irs, irs1, irs2, ilz, ild, ilt;

	start_time( WORK_TIMER );

	rs = max( 1.0 / auni - 1.0, 0.0 );

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
#pragma omp parallel for default(none), private(ilt), shared(irs1,irs2,ilz,ild,ccl_rs,coolcl,nlt)
				for ( ilt = 0; ilt < nlt; ilt++ ) {
					ccl_rs[ilz][ild][ilt] = coolcl[irs1][ilz][ild][ilt];
				}
			}
		}
		for ( ilz = 0; ilz < nlz; ilz++ ) {
			for ( ild = 0; ild < nld; ild++ ) {
#pragma omp parallel for default(none), private(ilt), shared(irs1,irs2,ilz,ild,fion_rs,fion,nlt)
				for ( ilt = 0; ilt < nlt; ilt++ ) {
					fion_rs[ilz][ild][ilt] = fion[irs1][ilz][ild][ilt];
				}
			}
		}
	} else {
		/* need to interpolate to specific redshift */
		rs1 = rsmin + drs*(double)irs1;
		rs2 = rsmin + drs*(double)irs2;

		for ( ilz = 0; ilz < nlz; ilz++ ) {
			for ( ild = 0; ild < nld; ild++ ) {
#pragma omp parallel for default(none), private(ilt,ac,bc), shared(irs1,irs2,ilz,ild,ccl_rs,coolcl,nlt,rs,rs1,rs2)
				for ( ilt = 0; ilt < nlt; ilt++ ) {
					ac = (rs2-rs)/(rs2-rs1);
					bc = (rs-rs1)/(rs2-rs1);
#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
					ccl_rs[ilz][ild][ilt] = ac*coolcl[irs1][ilz][ild][ilt] + bc*coolcl[irs2][ilz][ild][ilt];
#else
					ccl_rs[ilz][ild][ilt].Cooling = ac*coolcl[irs1][ilz][ild][ilt].Cooling + bc*coolcl[irs2][ilz][ild][ilt].Cooling;
					ccl_rs[ilz][ild][ilt].Heating = ac*coolcl[irs1][ilz][ild][ilt].Heating + bc*coolcl[irs2][ilz][ild][ilt].Heating;
#endif
				}
			}
		}
		for ( ilz = 0; ilz < nlz; ilz++ ) {
			for ( ild = 0; ild < nld; ild++ ) {
#pragma omp parallel for default(none), private(ilt,ac,bc), shared(irs1,irs2,ilz,ild,fion_rs,fion,nlt,rs,rs1,rs2)
				for ( ilt = 0; ilt < nlt; ilt++ ) {
					ac = (rs2-rs)/(rs2-rs1);
					bc = (rs-rs1)/(rs2-rs1);
					fion_rs[ilz][ild][ilt] = ac*fion[irs1][ilz][ild][ilt] + bc*fion[irs2][ilz][ild][ilt];
				}
			}
		}
	}

	end_time( WORK_TIMER );
}

cooling_t cooling_rate( double nHlog, double T_g, double Zlog ) {
	double Tlog;
	int it1, it2;
	int id1, id2;
	int iz1, iz2;
	double dd,td,zd;
	double d1,d2,d3,t1,t2,t3;
#ifndef OLDSTYLE_COOLING_EXPLICIT_SOLVER
	cooling_t ret;
#endif

	/* compute temperature bin */
	Tlog = log10(T_g);
	it1 = (int)((Tlog - tlmin)*dlti);
	it2 = it1 + 1;

	it1 = max(it1,0);
	it1 = min(it1,nlt-1);
	it2 = max(it2,0);
	it2 = min(it2,nlt-1);

	/* compute density bin */
	id1 = (int)((nHlog - dlmin)*dldi);
	id2 = id1 + 1;

	id1 = max(id1,0);
	id1 = min(id1,nld-1);
	id2 = max(id2,0);
	id2 = min(id2,nld-1);

	/* compute metallicity bin */
	iz1 = (int)((Zlog - Zlmin)*dlZi);
	iz2 = iz1 + 1;

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
	t2 = (dd - nHlog) * dldi;
	d2 = 1.0 - t2;
	t3 = (zd - Zlog) * dlZi;
	d3 = 1.0 - t3;

/* initial metallicity fails this interpolation
	cart_assert( t1 >= 0.0 && t1 <= 1.0 );
	cart_assert( t2 >= 0.0 && t2 <= 1.0 );
	cart_assert( t3 >= 0.0 && t3 <= 1.0 );
*/

#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER

	return	t1*t2*t3 * ccl_rs[iz1][id1][it1] +
		d1*t2*t3 * ccl_rs[iz1][id1][it2] +
		t1*d2*t3 * ccl_rs[iz1][id2][it1] +
		d1*d2*t3 * ccl_rs[iz1][id2][it2] +
		t1*t2*d3 * ccl_rs[iz2][id1][it1] +
		d1*t2*d3 * ccl_rs[iz2][id1][it2] +
		t1*d2*d3 * ccl_rs[iz2][id2][it1] +
		d1*d2*d3 * ccl_rs[iz2][id2][it2];

#else  /* OLDSTYLE_COOLING_EXPLICIT_SOLVER */

	ret.Cooling =   t1*t2*t3 * ccl_rs[iz1][id1][it1].Cooling +
			d1*t2*t3 * ccl_rs[iz1][id1][it2].Cooling +
			t1*d2*t3 * ccl_rs[iz1][id2][it1].Cooling +
			d1*d2*t3 * ccl_rs[iz1][id2][it2].Cooling +
			t1*t2*d3 * ccl_rs[iz2][id1][it1].Cooling +
			d1*t2*d3 * ccl_rs[iz2][id1][it2].Cooling +
			t1*d2*d3 * ccl_rs[iz2][id2][it1].Cooling +
			d1*d2*d3 * ccl_rs[iz2][id2][it2].Cooling;
	ret.Heating =   t1*t2*t3 * ccl_rs[iz1][id1][it1].Heating +
			d1*t2*t3 * ccl_rs[iz1][id1][it2].Heating +
			t1*d2*t3 * ccl_rs[iz1][id2][it1].Heating +
			d1*d2*t3 * ccl_rs[iz1][id2][it2].Heating +
			t1*t2*d3 * ccl_rs[iz2][id1][it1].Heating +
			d1*t2*d3 * ccl_rs[iz2][id1][it2].Heating +
			t1*d2*d3 * ccl_rs[iz2][id2][it1].Heating +
			d1*d2*d3 * ccl_rs[iz2][id2][it2].Heating;
	return ret;

#endif /* OLDSTYLE_COOLING_EXPLICIT_SOLVER */
}

double cooling_fion( double nHlog, double T_g, double Zlog ) {
	double Tlog;
	int it1, it2;
	int id1, id2;
	int iz1, iz2;
	int irs1, irs2;
	double dd,td,zd;
	double d1,d2,d3,t1,t2,t3;

	/* compute temperature bin */
	Tlog = log10(T_g);
	it1 = (int)((Tlog - tlmin)*dlti);
	it2 = it1 + 1;

	it1 = max(it1,0);
	it1 = min(it1,nlt-1);
	it2 = max(it2,0);
	it2 = min(it2,nlt-1);

	/* compute density bin */
	id1 = (int)((nHlog - dlmin)*dldi);
	id2 = id1 + 1;

	id1 = max(id1,0);
	id1 = min(id1,nld-1);
	id2 = max(id2,0);
	id2 = min(id2,nld-1);

	/* compute metallicity bin */
	iz1 = (int)((Zlog - Zlmin)*dlZi);
	iz2 = iz1 + 1;

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
	t2 = (dd - nHlog) * dldi;
	d2 = 1.0 - t2;
	t3 = (zd - Zlog) * dlZi;
	d3 = 1.0 - t3;

	return	t1*t2*t3 * fion_rs[iz1][id1][it1] +
		d1*t2*t3 * fion_rs[iz1][id1][it2] +
		t1*d2*t3 * fion_rs[iz1][id2][it1] +
		d1*d2*t3 * fion_rs[iz1][id2][it2] +
		t1*t2*d3 * fion_rs[iz2][id1][it1] +
		d1*t2*d3 * fion_rs[iz2][id1][it2] +
		t1*d2*d3 * fion_rs[iz2][id2][it1] +
		d1*d2*d3 * fion_rs[iz2][id2][it2];
}

#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
void test_cooling() {
	double a, rhogl, T_g, Z_met, cr;
	double Tl_min, Tl_max, d_Tl;
	int nt;
	int it;
	FILE *output;
	char filename[128];

	for ( a = 0.019608; a < 1.0; a += 0.01 ) {
		set_cooling_redshift(a);

		sprintf(filename, "%s/coolingrate_%6.4f.dat", output_directory, a );
		output = fopen( filename, "w" );

		for ( rhogl = dlmin; rhogl < dlmax; rhogl += 0.1*dld ) {
			Tl_min = 1.0;
			Tl_max = 9.0;
			d_Tl = 0.05;
			nt = (int)((Tl_max - Tl_min)/d_Tl);
			for ( it = 0; it < nt; it++ ) {
				T_g = pow(10.0, Tl_min + (float)(it)*d_Tl);
				for ( Z_met = Zlmin; Z_met < Zlmax; Z_met += 0.1*dlZ ) {
					cr = cooling_rate( rhogl, T_g, Z_met );

					fprintf(output, "%e %e %e %24.12e\n", rhogl, T_g, Z_met, cr );
				}
			}
		}

		fclose(output);

		cart_error("done outputting table");
	}
}
#endif /* OLDSTYLE_COOLING_EXPLICIT_SOLVER */

#endif /* COOLING && !RADIATIVE_TRANSFER */


