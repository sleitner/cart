#include "config.h"
#ifdef HYDRO

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "tree.h"

#define gammin  	(1.01)
#define gammax  	(10.0)
#define eps			(1e-6)
#define maxit		(50)
#define diffusion	(0.1)
#define dviscmax	(0.1)
#define drhomax		(0.2)
#define small_R     1.0e-20


extern int smooth_density_gradients;


void riemann( double stl[5], double str[5], double sta[4] ) {
	double p_l, p_r, p0, bgam_l, rgam_l, gam_l, xl2, xl3;
	double al, bl, cl, q_l, bgam_r, rgam_r, gam_r;
	double xr2, xr3, ar, br, cr, q_r, ul, ur;
	double p1, ul_0, ur_0, p_0, p_1, xxl, w2l, ul1, xxr, w2r, ur1;
	double p2, devi, dev, u, rho_s, u_s, p_s, bgam_s, gam_s;
	double a_s, b_s, c_s, w_s, rho, gam, xx4, xx5;
	double a2, a3, fs, indd, ind_r;

	int iter;

	/* initial guess for secant iteration */
	p_l     = stl[2];
	p_r     = str[2];
	p0      = 0.5 * ( p_l + p_r );
	bgam_l  = stl[3];
	rgam_l  = 1.0 / bgam_l;
	gam_l   = stl[4];
	xl2     = 2.0 * gam_l * rgam_l - 1.0;
	cl      = xl2 * p_l;
	xl3     = 0.5 * ( gam_l - 1.0 ) * xl2;
	bl      = xl3 * p_l;
	al      = gam_l - xl3;
	q_l     = sqrt( p_l * stl[0] * ((1.0+bgam_l) * p0/p_l + bgam_l - 1.0 ));
	
	bgam_r  = str[3];
	rgam_r  = 1.0 / bgam_r;
	gam_r   = str[4];
	xr2	= 2.0 * gam_r * rgam_r - 1.0;
	cr	= xr2 * p_r;
	xr3	= 0.5 * ( gam_r - 1.0 ) * xr2;
	br	= xr3 * p_r;
	ar	= gam_r - xr3;
	q_r	= sqrt( p_r * str[0] * ((1.0+bgam_r) * p0/p_r + bgam_r - 1.0 ));

	ul	= stl[1];
	ur	= str[1];
	p1	= (((ul-ur)*q_l*0.707106781 + p_l)*q_r + p_r*q_l) / (q_l + q_r);
	ul_0	= ul;
	ur_0	= ur;
	p_0	= p0;
	p_1	= max( small_R, p1 );

	/* Riemann solver - secant iterations for pressure */

	/* first iteration */
	xxl = ( al * p_1 + bl ) / ( p_1 + cl );
	w2l = 1.0/sqrt(max(small_R, xxl * stl[0] * (p_1 + stl[2])));
	ul1 = stl[1] + ( stl[2] - p_1 ) * w2l;
	xxr = ( ar * p_1 + br ) / ( p_1 + cr );
	w2r = 1.0/sqrt(max(small_R, xxr * str[0] * (p_1 + str[2])));
	ur1 = str[1] + ( p_1 - str[2] ) * w2r;
	p2 = max( small_R, 1.0000001 * p_1 - ( ur1 - ul1 ) 
			* fabs( p_1 - p_0 )
			/ ( fabs( ur1 - ur_0 )
			   +fabs( ul1 - ul_0 )
			   +small_R ) );
	p_0 = p_1;
	p_1 = p2;
	ul_0 = ul1;
	ur_0 = ur1;
	devi = fabs( p2 - p_0 ) / ( p2 + p_0 );
	dev = max( 0.0, devi );

	iter = 1;
	while ( iter <= maxit && dev > eps ) {
		xxl = ( al * p_1 + bl ) / ( p_1 + cl );
		w2l = 1.0/sqrt(max(small_R, xxl * stl[0] * (p_1 + stl[2])));
		ul1 = stl[1] + ( stl[2] - p_1 ) * w2l;
		xxr = ( ar * p_1 + br ) / ( p_1 + cr );
		w2r = 1.0/sqrt(max(small_R, xxr * str[0] * (p_1 + str[2])));
		ur1 = str[1] + ( p_1 - str[2] ) * w2r;
		p2 = max( small_R, 1.0000001 * p_1 - ( ur1 - ul1 )
			* fabs( p_1 - p_0 )
			/ ( fabs( ur1 - ur_0 )
			+fabs( ul1 - ul_0 )
			+small_R ) );
		dev = fabs( p2 - p_1 ) / ( p2 + p_1 );
		p_0 = p_1;
                p_1 = p2;
                ul_0 = ul1;
                ur_0 = ur1;
		iter++;
	}

	if ( dev > eps ) {
		cart_error("Riemann solver did not converge!\ndev = %e, p_0 = %e, p_1 = %e\nLeft: %e %e %e %e %e\nRight: %e %e %e %e %e",
			dev, p_0, p_1,
			stl[0], stl[1], stl[2], stl[3], stl[4],
			str[0], str[1], str[2], str[3], str[4] );
	}

	u	= 0.5 * ( ul_0 + ur_0 );
	ind_r = floor( 0.9 - sign( 0.5, u ) );
	cart_assert( ind_r == 0 || ind_r == 1.0 );

	if ( u >= 0.0 ) {
		rho_s = stl[0];
		u_s = stl[1];
		p_s = stl[2];
		bgam_s = stl[3];
		gam_s = stl[4];
		a_s = al;
		b_s = bl;
		c_s = cl;
		indd = -1.0;
	} else {
		rho_s = str[0];
		u_s = str[1];
		p_s = str[2];
		bgam_s = str[3];
		gam_s = str[4];
		a_s = ar;
		b_s = br;
		c_s = cr;
		indd = 1.0;
	}

	w_s	= ( a_s * p_1 + b_s ) / ( p_1 + c_s );
	w_s	= max( small_R, w_s * rho_s * ( p_1 + p_s ) );
	rho	= max( small_R, rho_s / ( 1.0 - rho_s * ( p_1 - p_s ) / w_s ) );
	gam	= gam_s + 2.0 * ( gam_s - 1.0 ) * ( 1.0 - gam_s / bgam_s ) 
			/ ( p_1 + p_s ) * ( p_1 - p_s );

	if ( (p_1 - p_s) > 0 ) {
		a2 = indd * u_s + sqrt( w_s ) / rho_s;
		a3 = 1.0e-10;
	} else {
		xx4 = indd * u_s + sqrt( bgam_s * p_s/rho_s );
		xx5 = indd * u + sqrt( 0.5 * ( stl[3]+str[3] ) * p_1 / rho );
		a2 = xx5;
		a3 = xx4 - xx5;
	}

	fs = -a2/a3;
	if ( fs < 0.0 ) fs = 0.0;
	if ( fs > 1.0 ) fs = 1.0;

	sta[0] = rho	+ fs * ( rho_s - rho );
	sta[1] = u	    + fs * ( u_s - u );
	sta[2] = p_1	+ fs * ( p_s - p_1 );
	sta[3] = gam	+ fs * ( gam_s - gam );
}       

#ifdef GRAVITY_IN_RIEMANN
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double g[2], double c[2], double f[num_hydro_vars-1] ) 
#else 
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double c[2], double f[num_hydro_vars-1] )
#endif
{
	int j;

	double dv0, dv1, dv2, dv0a, dv1a, dv2a;
	double dv11, dv20, dlq0, dlq1, dv00, dv01;
	double rhor_l, u_l, a_l, c_l, c2_l, cp_l, cm_l, x_l, rhow_l;
	double pw_l, uw_l, vw_l, ww_l, gamw_l, b_l, x0_l, rho0_l;
	double v0_l, w0_l, gam0_l, xx1, xx2, b0_l, v_l, w_l, p_l, gam_l;
	double rhor_r, u_r, a_r, c_r, c2_r, cp_r, cm_r, x_r, rhow_r;
	double pw_r, uw_r, vw_r, ww_r, gamw_r, b_r, x0_r, rho0_r;
	double v0_r, w0_r, gam0_r, b0_r, v_r, w_r, p_r, gam_r;
	double rho, vu, pre, gam, xup_r, xup_l, vv, vw, fmass, predtx;
	double pvu, vudtx;
	double fu, fv, fw, fe;
	double fl, fr;
	double dv[6][2];
	double stl[5], str[5], sta[4];

	/* slopes */
	for ( j = 0; j < 6; j++ ) {
		dv0 = v[j][1] - v[j][0];
		dv1 = v[j][2] - v[j][1];
		dv2 = v[j][3] - v[j][2];

		dv0a = 2.0*fabs(dv0);
		dv1a = 2.0*fabs(dv1);
		dv2a = 2.0*fabs(dv2);

		dv11 = v[j][2] - v[j][0];
		dv20 = v[j][3] - v[j][1];

		dlq0 = ( dv1*dv0 < 0.0 ) ? 0.0 : min( dv0a, dv1a );
		dlq1 = ( dv1*dv2 < 0.0 ) ? 0.0 : min( dv1a, dv2a );

		dv00 = sign( min( 0.5*fabs(dv11), dlq0 ), dv11 ) * c[0];
		dv01 = sign( min( 0.5*fabs(dv20), dlq1 ), dv20 ) * c[1];
				  
		dv[j][0] = sign( min( min( fabs(dv00), dv1a ), dv0a ), dv11 );
		dv[j][1] = sign( min( min( fabs(dv01), dv2a ), dv1a ), dv20 );
	}

	/* left states */
	rhor_l	= 1.0/v[0][1];
	u_l	= v[2][1];
	a_l	= sqrt( v[6][1] * v[1][1] * rhor_l );
	c_l	= a_l * v[0][1];
	c2_l	= c_l*c_l;
	cp_l	= u_l + a_l;
	cm_l	= u_l - a_l;
	x_l	= 0.5 * (1.0 - dtx * max( (float)cp_l, 0.0 ));
	rhow_l	= max( small_R, (float)(v[0][1] + x_l*dv[0][0]));
	pw_l	= max( small_R, (float)(v[1][1] + x_l*dv[1][0]));
	uw_l	= v[2][1] + x_l * dv[2][0];
	vw_l	= v[3][1] + x_l * dv[3][0];
	ww_l	= v[4][1] + x_l * dv[4][0];
	gamw_l	= v[5][1] + x_l * dv[5][0];
	b_l	= ( cm_l < 0.0 ) ? 0.0 : dtx2 * rhor_l * ( dv[1][0] / c_l - dv[2][0] );

	x0_l	= 0.5 - dtx2 * v[2][1];
	rho0_l	= v[0][1] + x0_l * dv[0][0];
	v0_l 	= v[3][1] + x0_l * dv[3][0];
	w0_l	= v[4][1] + x0_l * dv[4][0];
	gam0_l	= v[5][1] + x0_l * dv[5][0];
	xx1 	= ( v[2][1] > 0 ) ? 1.0 : 0.0;
	xx2 	= 1.0 - xx1;
	b0_l	= xx1 * a_l * dtx2 * ( 	  dv[0][0] / ( rhow_l * rho0_l ) 
					- dv[1][0] / c2_l );
	v_l	= xx1*v0_l + xx2*vw_l;
	w_l	= xx1*w0_l + xx2*ww_l;

	stl[0]	= max( small_R,(float)(rhow_l / ( 1.0 - ( b0_l + b_l ) * rhow_l ) ) );

#ifdef GRAVITY_IN_RIEMANN
	stl[1]	= uw_l - b_l * c_l + g[0];
#else
	stl[1]	= uw_l - b_l * c_l;
#endif

	p_l	= max( small_R,(float)(pw_l + b_l * c2_l) );
	stl[2]	= p_l;
	stl[3]	= v[6][1];
	gam_l	= xx1 * ( gam0_l + 2.0 * ( 1.0 - v[5][1] / v[6][1] )
					*( v[5][1] - 1.0 )
					*( p_l - v[1][1] )
					/( p_l + v[1][1] ) )
			+ xx2 * gamw_l;
	stl[4]	= max( gammin, min( gammax, (float)(gam_l) ) );

	/* right states */
	rhor_r  = 1.0/v[0][2];
	u_r     = v[2][2];
	a_r     = sqrt( v[6][2] * v[1][2] * rhor_r );
	c_r     = a_r * v[0][2];
	c2_r    = c_r*c_r;
	cp_r    = u_r + a_r;
	cm_r    = u_r - a_r;
	x_r     = -(float)0.5*(1.0 + dtx * min( (float)cm_r, 0.0 ));
	rhow_r  = max( small_R, (float)(v[0][2] + x_r*dv[0][1]));
	pw_r    = max( small_R, (float)(v[1][2] + x_r*dv[1][1]));
	uw_r    = v[2][2] + x_r * dv[2][1];
	vw_r    = v[3][2] + x_r * dv[3][1];
	ww_r    = v[4][2] + x_r * dv[4][1];
	gamw_r  = v[5][2] + x_r * dv[5][1];
	b_r	= ( cp_r >= 0.0 ) ? 0.0 : dtx2 * rhor_r * ( dv[1][1] / c_r + dv[2][1] );

	x0_r    = -(float)0.5*(1.0 + dtx * v[2][2]);
	rho0_r  = v[0][2] + x0_r * dv[0][1];
	v0_r    = v[3][2] + x0_r * dv[3][1];
	w0_r    = v[4][2] + x0_r * dv[4][1];
	gam0_r  = v[5][2] + x0_r * dv[5][1];
	xx1     = ( v[2][2] < 0 ) ? 1.0 : 0.0;
	xx2     = 1.0 - xx1;
	b0_r    = xx1 * a_r * dtx2 * (    dv[1][1] / c2_r
					- dv[0][1] / ( rhow_r * rho0_r ) );
	v_r     = xx1*v0_r + xx2*vw_r;
	w_r     = xx1*w0_r + xx2*ww_r;

	str[0]  = max( small_R,(float)(rhow_r / ( 1.0 - ( b0_r + b_r ) * rhow_r ) ) );

#ifdef GRAVITY_IN_RIEMANN
	str[1]	= uw_r - b_r * c_r + g[1];
#else
	str[1]	= uw_r - b_r * c_r;
#endif

	p_r     = max( small_R,(float)(pw_r + b_r * c2_r) );
	str[2]  = p_r;
	str[3]  = v[6][2];
	gam_r   = xx1 * ( gam0_r + 2.0 * ( 1.0 - v[5][2] / v[6][2] )
					*( v[5][2] - 1.0 )
					*( p_r - v[1][2] )
					/( p_r + v[1][2] ) )
			+ xx2 * gamw_r;
	str[4]  = max( gammin, min( gammax, (float)(gam_r) ) );

	/* call to riemann function */
	riemann( stl, str, sta );

	/* compute fluxes of hydro variables */
	rho	= sta[0];
	vu	= sta[1];
	pre	= sta[2];
	gam	= sta[3];

	if ( vu < 0.0 ) {
		vv 	= v_r;
		vw 	= w_r;
	        xup_r   = 1.0;
		xup_l   = 0.0;
	} else {
		vv 	= v_l;
		vw	= w_l;
		xup_r	= 0.0;
		xup_l	= 1.0;
	}

	vudtx	= vu * dtx;
	fmass	= rho*vudtx;
	predtx	= pre * dtx;
	fu	= fmass * vu + predtx;
	fv	= fmass * vv;
	fw	= fmass * vw;
	pvu	= vu * predtx / ( gam - 1.0 );
	fe	= pvu * gam + 0.5 * fmass * ( vu*vu + vv*vv + vw*vw );

	f[0]	= fmass;
	f[1]	= fu;
	f[2]	= fv;
	f[3]	= fw;
	f[4]	= fe;
	f[5]	= pvu;
	f[6]	= vu;

#ifdef ELECTRON_ION_NONEQUILIBRIUM 
	dv0 = v[7][1] - v[7][0];
	dv1 = v[7][2] - v[7][1];
	dv2 = v[7][3] - v[7][2];
	dv11 = v[7][2] - v[7][0];
	dv20 = v[7][3] - v[7][1];

	dlq0 = ( dv1*dv0 < 0.0 ) ? 0.0 : 2.0*min( fabs(dv0), fabs(dv1) );
	dlq1 = ( dv1*dv2 < 0.0 ) ? 0.0 : 2.0*min( fabs(dv1), fabs(dv2) );

	fl = v[7][1] + ( 0.5 - dtx2 * vu ) * sign( min( 0.5 * fabs(dv11), dlq0 ), dv11 ) * c[0];
	fr = v[7][2] - ( 0.5 + dtx2 * vu ) * sign( min( 0.5 * fabs(dv20), dlq1 ), dv20 ) * c[1];

	f[7] = vudtx * ( fl * xup_l + fr * xup_r );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	/* compute fluxes of advected species */
	for ( j = num_hydro_vars-num_chem_species-1; j < num_hydro_vars - 1; j++ ) {
		dv0 = v[j][1] - v[j][0];
		dv1 = v[j][2] - v[j][1];
		dv2 = v[j][3] - v[j][2];
		dv11 = v[j][2] - v[j][0];
		dv20 = v[j][3] - v[j][1];

		dlq0 = ( dv1*dv0 < 0.0 ) ? 0.0 : 2.0*min( fabs(dv0), fabs(dv1) );
		dlq1 = ( dv1*dv2 < 0.0 ) ? 0.0 : 2.0*min( fabs(dv1), fabs(dv2) );

		fl = v[j][1] + ( 0.5 - dtx2 * vu ) * sign( min( 0.5 * fabs(dv11), dlq0 ), dv11 ) * c[0];
		fr = v[j][2] - ( 0.5 + dtx2 * vu ) * sign( min( 0.5 * fabs(dv20), dlq1 ), dv20 ) * c[1];

		/* advected variables in v[j] must be in dimensionless units (var/density) */
		f[j] = fmass * ( fl * xup_l + fr * xup_r );
	}
}

void lapidus( double dtx2, int L1, int R1, int sweep_direction, int j3, int j4, int j5, 
		double v[num_hydro_vars-1][4], double f[num_hydro_vars-1] ) {
	int i, j;
	double diffk;
	double diff;
	double gvisc;
	double dvisc;
	int neighborsL1[num_neighbors];
	int neighborsR1[num_neighbors];
	int iCh00, iCh01, iCh10, iCh11;
	double v00, v01, v10, v11;
	double xx;
	int j0, j1;

	diffk = dtx2*diffusion;

	gvisc = 2.0 * ( v[2][1] - v[2][2] );

	/* Compute neighbors */
	cell_all_neighbors( L1, neighborsL1 );
	cell_all_neighbors( R1, neighborsR1 );

	/* velocity difference in directions orthogonal to sweep direction */
	for ( i = 0; i < nDim; i++ ) {
		j0 = 2*i;
		j1 = j0+1;

		if ( j0 != sweep_direction && j1 != sweep_direction ) {
			iCh00 = neighborsL1[j0];
			iCh01 = neighborsL1[j1];
			iCh10 = neighborsR1[j0];
			iCh11 = neighborsR1[j1];

			v00 = cell_momentum(iCh00,i) / cell_gas_density(iCh00);
			v01 = cell_momentum(iCh01,i) / cell_gas_density(iCh01);
			v10 = cell_momentum(iCh10,i) / cell_gas_density(iCh10);
			v11 = cell_momentum(iCh11,i) / cell_gas_density(iCh11);

			gvisc += v10 - v11 + v00 - v01;
		}
	}
		
	diff = diffk * max( 0.0, gvisc );

	if(smooth_density_gradients)
	  {
	    xx = drhomax * max( v[0][1], v[0][2] ) / ( min( v[0][1], v[0][2] ) );
	    dvisc = max( 0.0, dviscmax * ( xx - 1.0 ) / ( xx + 1.0 ) );
	    diff = max( diff, dvisc );
	  }

	f[0] += diff * ( cell_gas_density(L1) - cell_gas_density(R1) );
	f[1] += diff * ( cell_momentum(L1,j3) - cell_momentum(R1,j3) );
	f[2] += diff * ( cell_momentum(L1,j4) - cell_momentum(R1,j4) );
	f[3] += diff * ( cell_momentum(L1,j5) - cell_momentum(R1,j5) );
	f[4] += diff * ( cell_gas_energy(L1) - cell_gas_energy(R1) );
	f[5] += diff * ( cell_gas_internal_energy(L1) - cell_gas_internal_energy(R1) );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	f[7] += diff * ( cell_electron_internal_energy(L1) - cell_electron_internal_energy(R1) );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_chem_species; j++ ) {
		f[j+num_hydro_vars-num_chem_species-1] += diff * ( cell_advected_variable(L1,j) - cell_advected_variable(R1,j) );
	}
}

#endif /* HYDRO */
