#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "io.h"
#include "sfc.h"
#include "units.h"
#include "constants.h"
#include "auxiliary.h"

int num_root_cells;
float aexp, aexpv;
float astep;
double vcorr;

int *particle_count;
float *density;
float *momentum;

void particle_callback( particle_struct *particle ) {
	int i, j;
	int ix, ix1, iy, iy1, iz, iz1;
	double dx0, dx1, dy0, dy1, dz0, dz1;
	double d00, d01, d10, d11;
	double mass;
	particle_struct particle2;
	int dx, dy, dz;

	if ( particle->specie != 0 ) {
		particle2.mass = particle->mass / 8.;
		particle2.specie = 0;

		particle2.v[0] = particle->v[0];
		particle2.v[1] = particle->v[1];
		particle2.v[2] = particle->v[2];

		for ( dx = -1; dx <= 1; dx += 2 ) {
			particle2.x[0] = particle->x[0] + (double)dx*0.5;

			for ( dy = -1; dy <= 1; dy += 2 ) {
				particle2.x[1] = particle->x[1] + (double)dy*0.5;

				for ( dz = -1; dz <= 1; dz += 2 ) {
					particle2.x[2] = particle->x[2] + (double)dz*0.5;
					particle_callback( &particle2 );
				}
			}
		}

		return;
	}

	ix = (int)particle->x[0];
	iy = (int)particle->x[1];
	iz = (int)particle->x[2];

	particle_count[iz+num_grid*(iy+num_grid*ix)]++;

	/* compute momentum and density */
	mass = particle->mass * (Omegab0 / Omega0) / ( 1.0 - Omegab0/Omega0 ); 

	ix = ((int)(particle->x[0]-0.5) + num_grid ) % num_grid;
	ix1 = (ix+1) % num_grid;
	iy = ((int)(particle->x[1]-0.5) + num_grid ) % num_grid;
	iy1 = (iy+1) % num_grid;
	iz = ((int)(particle->x[2]-0.5) + num_grid ) % num_grid;
	iz1 = (iz+1) % num_grid;

	dx1 = particle->x[0] - floor(particle->x[0]) - 0.5;
	if ( dx1 < 0.0 ) dx1 += 1.0;

	dy1 = particle->x[1] - floor(particle->x[1]) - 0.5;
	if ( dy1 < 0.0 ) dy1 += 1.0;

	dz1 = particle->x[2] - floor(particle->x[2]) - 0.5;
	if ( dz1 < 0.0 ) dz1 += 1.0;

	dx0 = 1.0 - dx1;
	dy0 = 1.0 - dy1;
	dz0 = 1.0 - dz1;

	dx0 *= mass;
	dx1 *= mass;

	d00 = dx0*dy0;
	d01 = dx0*dy1;
	d10 = dx1*dy0;
	d11 = dx1*dy1;

	density[iz+num_grid*(iy+num_grid*ix)] += d00*dz0;
	density[iz1+num_grid*(iy+num_grid*ix)] += d00*dz1;
	density[iz+num_grid*(iy1+num_grid*ix)] += d01*dz0;
	density[iz1+num_grid*(iy1+num_grid*ix)] += d01*dz1;
	density[iz+num_grid*(iy+num_grid*ix1)] += d10*dz0;
	density[iz1+num_grid*(iy+num_grid*ix1)] += d10*dz1;
	density[iz+num_grid*(iy1+num_grid*ix1)] += d11*dz0;
	density[iz1+num_grid*(iy1+num_grid*ix1)] += d11*dz1;

	for ( j = 0; j < nDim; j++ ) {
		momentum[j*num_root_cells+iz+num_grid*(iy+num_grid*ix)] += d00*dz0*vcorr*particle->v[j];
		momentum[j*num_root_cells+iz1+num_grid*(iy+num_grid*ix)] += d00*dz1*vcorr*particle->v[j];
		momentum[j*num_root_cells+iz+num_grid*(iy1+num_grid*ix)] += d01*dz0*vcorr*particle->v[j];
		momentum[j*num_root_cells+iz1+num_grid*(iy1+num_grid*ix)] += d01*dz1*vcorr*particle->v[j];
		momentum[j*num_root_cells+iz+num_grid*(iy+num_grid*ix1)] += d10*dz0*vcorr*particle->v[j];
		momentum[j*num_root_cells+iz1+num_grid*(iy+num_grid*ix1)] += d10*dz1*vcorr*particle->v[j];
		momentum[j*num_root_cells+iz+num_grid*(iy1+num_grid*ix1)] += d11*dz0*vcorr*particle->v[j];
		momentum[j*num_root_cells+iz1+num_grid*(iy1+num_grid*ix1)] += d11*dz1*vcorr*particle->v[j];
	}
}

int main ( int argc, char *argv[]) {
	int i, j, k;
	FILE *output;
	int size;
	float boxh;
	double TinitK, Tinit, a_th;
	double min_density, max_density;
	double kinetic_energy;
	float *internal, *energy;
	particle_header header;
	int endian, nbody_flag;
	double total_mass;
	char filename[256];
	int max_count, min_count;
	int zero_count;
	
	if ( argc != 7 ) {
		cart_error("usage nbody_to_ics Lbox Omegab0 num_grid PMcrd.DAT PMcrs.DAT output_directory");
	}

	Lbox = atof(argv[1]);
	Omegab0 = atof(argv[2]);
	num_grid = atoi(argv[3]);

	num_root_cells = num_grid*num_grid*num_grid;

	read_particle_header( argv[4], &header, &endian, &nbody_flag );

	aexpn = header.aexpn;
	Omega0 = header.Om0;
	OmegaL0 = header.Oml0;
	hubble = header.hubble;

	init_units();

	aexp = header.aexpn;
	aexpv = header.aexpn - 0.5*header.astep;

	vcorr = pow( (aexpn/aexpv), 1.5 ) * sqrt( Omega0 + OmegaL0*aexpn*aexpn*aexpn ) / 
		sqrt( Omega0 + OmegaL0 * aexpv*aexpv*aexpv );

	cart_debug("aexpv = %f", aexpv );
	cart_debug("aexpn = %f", aexpn );
	cart_debug("vcorr = %e", vcorr );

	particle_count = cart_alloc( num_root_cells * sizeof(int) );
	density = cart_alloc( num_root_cells * sizeof(float) );
	momentum = cart_alloc( 3*num_root_cells * sizeof(float) );
	internal = cart_alloc( num_root_cells * sizeof(float) );
	energy = cart_alloc( num_root_cells * sizeof(float) );

	read_particles( argv[4], argv[5], particle_callback );

	a_th = 1.0 / (1e3*pow(Omegab0*hubble*hubble, 0.4));

	if ( aexpn < a_th ) {
		TinitK = T_CMB0 / aexpn;
	} else {
		TinitK = T_CMB0 / a_th * (a_th/aexpn)*(a_th/aexpn);
	}

	/* second gamma -1 factor to convert from T -> e */
	Tinit = TinitK * aexpn*aexpn / ( T0 * ( gamma - 1.0 ) * ( gamma - 1.0 ));

	total_mass = 0.0;
	max_density = -1e20;
	min_density = 1e20;

	max_count = -1e6;
	min_count = 1e6;
	zero_count = 0;

	/* set energy variables */
	for ( i = 0; i < num_root_cells; i++ ) {
		internal[i] = Tinit * density[i];

		kinetic_energy = 0.0;
		for ( j = 0; j < nDim; j++ ) {
			kinetic_energy += momentum[num_root_cells*j+i]*momentum[num_root_cells*j+i];
		}
		kinetic_energy *= 0.5/density[i];

		energy[i] = kinetic_energy + internal[i];	

		total_mass += density[i];

		min_density = min( min_density, density[i] );
		max_density = max( max_density, density[i] );

		max_count = min( max_count, particle_count[i] );
		min_count = max( min_count, particle_count[i] );

		if ( particle_count[i] == 0 ) {
			zero_count++;
		}
	}

	cart_debug("particle count = %u %u", min_count, max_count );
	cart_debug("zero count = %u", zero_count );

	cart_debug("total_gas_mass = %f, min = %e, max = %e", total_mass, min_density, max_density );

	/* reduce particle masses by Omegab factor */
	for ( i = 0; i < num_particle_species; i++ ) {
		total_mass += particle_species_num[i]*particle_species_mass[i];
	}

	cart_debug("total_mass = %f", total_mass );

	/* now write out initial conditions file */
	//sprintf( filename,  "%s/PMcrd.DAT", output_directory );
	//sprintf( filename2, "%s/PMcrs0.DAT", output_directory );
	//write_particles( filename, filename2, NULL, NULL );

	sprintf( filename, "%s/tr_ic.dat", argv[6] );
	output = fopen( filename, "w" );

	boxh = Lbox;

	size = sizeof(float);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( &boxh, sizeof(float), 1, output ); 
	fwrite( &size, sizeof(int), 1, output );

	size = 2*sizeof(float);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( &aexp, sizeof(float), 1, output );
	fwrite( &astep, sizeof(float), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	size = sizeof(int);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( &num_root_cells, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	size = num_root_cells*sizeof(float);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( density, sizeof(float), num_root_cells, output );
	fwrite( &size, sizeof(int), 1, output );

	fwrite( &size, sizeof(int), 1, output );
	fwrite( momentum, sizeof(float), num_root_cells, output );
	fwrite( &size, sizeof(int), 1, output );

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &momentum[num_root_cells], sizeof(float), num_root_cells, output );
	fwrite( &size, sizeof(int), 1, output );

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &momentum[2*num_root_cells], sizeof(float), num_root_cells, output );
	fwrite( &size, sizeof(int), 1, output );

	fwrite( &size, sizeof(int), 1, output );
	fwrite( energy, sizeof(float), num_root_cells, output );
	fwrite( &size, sizeof(int), 1, output );

	fwrite( &size, sizeof(int), 1, output );
	fwrite( internal, sizeof(float), num_root_cells, output );
	fwrite( &size, sizeof(int), 1, output );	

	fclose(output);

	return 0;
}
