#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "io.h"
#include "units.h"
#include "sfc.h"
#include "auxiliary.h"

char output_directory[256];
char logfile_directory[256];
char jobname[256];
int num_output_files = 1;

double particle_species_mass[10];
int particle_species_indices[11];
int particle_species_num[10];
int num_particle_species;

const float cell_delta[8][3] = {
        { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 },
        { 1.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 },
        { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }
};

void reorder( char *buffer, int size ) {
        int i;
        char tmp;
                                                                                
        for ( i = 0; i < (size/2); i++ ) {
                tmp = buffer[i];
                buffer[i] = buffer[size - i - 1];
                buffer[size - i - 1] = tmp;
        }
}

int compare_indices( const void *a, const void *b ) {
        return *(int *)a - *(int *)b;
}

void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag ) {
	int i;
	FILE *input;
	nbody_particle_header nbody_header;
	char desc[46];
	int size;

	*nbody_flag = 0;
	*endian = 0;

	/* read file header */
	input = fopen( header_filename, "r");
	if ( input == NULL ) {
		cart_error("Unable to open particle file %s", header_filename );
	}

	fread( &size, sizeof(int), 1, input );
	fread( desc, sizeof(char), 45, input );
	desc[45] = '\0';

	cart_debug( "Particle header file: %s", desc );

	if ( size != sizeof(particle_header)+45 ) {
		if ( size == sizeof(nbody_particle_header)+45 ) {
			*nbody_flag = 1;
		} else {
			reorder( (char *)&size, sizeof(int) );

			if ( size != sizeof(particle_header)+45 ) {
				if ( size == sizeof(nbody_particle_header)+45 ) {
					*endian = 1;
					*nbody_flag = 1;
				} else {
					cart_error("Size mismatch in reading particle file header %s (%u vs %u)",
							header_filename, size, sizeof(particle_header)+45 );
				}
			} else {
				*endian = 1;
			}
		}
	}

	if ( *nbody_flag ) {
		cart_debug("USING OLD NBODY FILE FORMAT (hope that's what you wanted...)");

		fread( &nbody_header, sizeof(nbody_particle_header), 1, input );

		header->aexpn = nbody_header.aexpn;
		header->aexp0 = nbody_header.aexp0;
		header->amplt = nbody_header.amplt;
		header->astep = nbody_header.astep;
		header->istep = nbody_header.istep;
		header->partw = nbody_header.partw;
		header->tintg = nbody_header.tintg;
		header->ekin = nbody_header.ekin;
		header->ekin1 = nbody_header.ekin1;
		header->ekin2 = nbody_header.ekin2;
		header->au0 = nbody_header.au0;
		header->aeu0 = nbody_header.aeu0;
		header->Nrow = nbody_header.Nrow;
		header->Ngrid = nbody_header.Ngrid;
		header->Nspecies = nbody_header.Nspecies;
		header->Nseed = nbody_header.Nseed;
		header->Om0 = nbody_header.Om0;
		header->Oml0 = nbody_header.Oml0;
		header->hubble = nbody_header.hubble;
		header->Wp5 = nbody_header.Wp5;
		header->Ocurv = nbody_header.Ocurv;

		for ( i = 0; i < 10; i++ ) {
			header->mass[i] = nbody_header.mass[i];
			header->num[i] = nbody_header.num[i];
		}

	} else {
		fread( header, sizeof(particle_header), 1, input );
	}

	fread( &size, sizeof(int), 1, input );
	fclose(input);

	if ( *endian ) {
		reorder( (char *)&header->aexpn, sizeof(float) );
		reorder( (char *)&header->aexp0, sizeof(float) );
		reorder( (char *)&header->amplt, sizeof(float) );
		reorder( (char *)&header->astep, sizeof(float) );
		reorder( (char *)&header->istep, sizeof(int) );
		reorder( (char *)&header->partw, sizeof(float) );
		reorder( (char *)&header->tintg, sizeof(float) );
		reorder( (char *)&header->ekin, sizeof(float) );
		reorder( (char *)&header->ekin1, sizeof(float) );
		reorder( (char *)&header->ekin2, sizeof(float) );
		reorder( (char *)&header->au0, sizeof(float) );
		reorder( (char *)&header->aeu0, sizeof(float) );
		reorder( (char *)&header->Nrow, sizeof(int) );
		reorder( (char *)&header->Ngrid, sizeof(int) );
		reorder( (char *)&header->Nspecies, sizeof(int) );
		reorder( (char *)&header->Nseed, sizeof(int) );
		reorder( (char *)&header->Om0, sizeof(float) );
		reorder( (char *)&header->Oml0, sizeof(float) );
		reorder( (char *)&header->hubble, sizeof(float) );
		reorder( (char *)&header->Wp5, sizeof(float) );
		reorder( (char *)&header->Ocurv, sizeof(float) );
		reorder( (char *)&header->Omb0, sizeof(float) );

		for ( i = 0; i < header->Nspecies; i++ ) {
			reorder( (char *)&header->mass[i], sizeof(float) );
			reorder( (char *)&header->num[i], sizeof(int) );
		}
	}
}

void read_particles( char *header_filename, char *data_filename, void callback( particle_struct * ) ) {
	int i, j, k;
	char desc[47];
	int proc;
	int num_parts;
	particle_struct particle;
	particle_float *input_page, *x, *y, *z, *vx, *vy, *vz;
	double dt;
	double r;
	int ipart;
	int ihalo;
	int sfc;
	int num_read;
	int coords[nDim];
	int num_parts_in_page, num_parts_per_page;
	int num_parts_per_proc_page;
	int num_pages, index;
	int current_id, current_type, local_count;
	int size, endian, dt_endian, stellar_endian;
	FILE *input;
	FILE *timestep_input;
	particle_header header;
	int nbody_flag;
	int grid_change_flag;
	float rfact, vfact;
	float grid_shift;

	long num_particles_total;
	double total_stellar_mass, total_stellar_initial_mass;

	read_particle_header( header_filename, &header, &endian, &nbody_flag );

	if ( header.Ngrid != num_grid ) {
		rfact = (float)num_grid / (float)header.Ngrid;
		grid_change_flag = 1;
		cart_debug("rfact = %f", rfact );
	} else {
		cart_debug("header.Ngrid = %d, num_grid = %d", header.Ngrid, num_grid );

		grid_change_flag = 0;
	}

	if ( nbody_flag ) {
		vfact = 2.0/sqrt(Omega0);
		grid_shift = 1.5;
	} else {
		vfact = 1.0;
		grid_shift = 1.0;
	}
	
	num_particle_species = header.Nspecies;
	
	particle_species_indices[0] = 0;
	for ( i = 0; i < num_particle_species; i++ ) {
		if ( grid_change_flag ) {
			particle_species_mass[i] = header.mass[i]*rfact*rfact*rfact;
		} else {
			particle_species_mass[i] = header.mass[i];
		}
		particle_species_indices[i+1] = header.num[i];
		particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
	}

	num_particles_total = particle_species_indices[num_particle_species];
	num_parts_per_page = header.Nrow*header.Nrow;
	num_pages = (num_particles_total-1) / num_parts_per_page + 1;

	input_page = cart_alloc( 2*nDim*num_parts_per_page*sizeof(particle_float) );

	x = input_page;
	y = &input_page[num_parts_per_page];
	z = &input_page[2*num_parts_per_page];
	vx = &input_page[3*num_parts_per_page];
	vy = &input_page[4*num_parts_per_page];
	vz = &input_page[5*num_parts_per_page];

	/* start loading actual particle data */
	input = fopen( data_filename, "r" );
	if ( input == NULL ) {
		cart_error( "Unable to open particle file %s for reading!", data_filename );
	}

	current_id = 0;
	current_type = 0;

	for ( i = 0; i < num_pages; i++ ) {
		cart_debug("page %u/%u", i, num_pages );

		if ( i == num_pages - 1 ) {
			num_parts_in_page = num_particles_total - num_parts_per_page*(num_pages-1);
		} else {
			num_parts_in_page = num_parts_per_page;
		}

		num_read = fread( input_page, sizeof(particle_float), 2*nDim*num_parts_per_page, input );
		if ( num_read != 2*nDim*num_parts_per_page ) {
			cart_error("Error reading from particle file %s: insufficient data", data_filename );
		}

		if ( endian ) {
			for ( j = 0; j < num_parts_in_page; j++ ) {
				reorder( (char *)&x[j], sizeof(particle_float) );
				reorder( (char *)&y[j], sizeof(particle_float) );
				reorder( (char *)&z[j], sizeof(particle_float) );
				reorder( (char *)&vx[j], sizeof(particle_float) );
				reorder( (char *)&vy[j], sizeof(particle_float) );
				reorder( (char *)&vz[j], sizeof(particle_float) );
			}
		}

		for ( j = 0; j < num_parts_in_page; j++ ) {
			/* convert to our coordinates 0->num_grid */
			x[j] -= grid_shift;
			y[j] -= grid_shift;
			z[j] -= grid_shift;

			if ( grid_change_flag ) {
				x[j] *= rfact;
				y[j] *= rfact;
				z[j] *= rfact;
			}

			if ( nbody_flag || grid_change_flag ) {
				vx[j] *= vfact;
				vy[j] *= vfact;
				vz[j] *= vfact;
			}

			/* enforce periodic boundary conditions */
			if ( x[j] < 0.0 ) {
				x[j] += (double)num_grid;
			} else if ( x[j] >= (double)num_grid ) {
				x[j] -= (double)num_grid;
			}

			if ( y[j] < 0.0 ) {
				y[j] += (double)num_grid;
			} else if ( y[j] >= (double)num_grid ) {
				y[j] -= (double)num_grid;
			}

			if ( z[j] < 0.0 ) {
				z[j] += (double)num_grid;
			} else if ( z[j] >= (double)num_grid ) {
				z[j] -= (double)num_grid;
			}

			if ( current_id >= particle_species_indices[current_type+1] ) {
				current_type++;
			}

			particle.x[0] = x[j];
			particle.x[1] = y[j];
			particle.x[2] = z[j];
			particle.v[0] = vx[j];
			particle.v[1] = vy[j];
			particle.v[2] = vz[j];
			particle.mass = particle_species_mass[current_type];
			particle.specie = current_type;

			callback( &particle );

			current_id++;
		}
	}
	
	cart_free( input_page );
	fclose(input);
}
