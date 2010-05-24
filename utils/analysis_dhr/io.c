#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "io.h"
#include "units.h"
#include "sfc.h"
#include "analysis.h"
#include "auxiliary.h"

char output_directory[256];
char logfile_directory[256];
char jobname[256];
int num_output_files = 1;

const double cell_delta[8][nDim] = {
                { -0.5, -0.5, -0.5 }, {  0.5, -0.5, -0.5 },
                { -0.5,  0.5, -0.5 }, {  0.5,  0.5, -0.5 },
                { -0.5, -0.5,  0.5 }, {  0.5, -0.5,  0.5 },
                { -0.5,  0.5,  0.5 }, {  0.5,  0.5,  0.5 }
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

#ifdef PARTICLES

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

void read_particles( char *header_filename, char *data_filename, char *stellar_filename,
			halo_list *halos, void callback( halo_struct *, particle_struct * ) ) {
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
	float rfact, vfact;
	float grid_shift;

	long num_particles_total;
	double particle_species_mass[10];
	int particle_species_indices[11];
	int particle_species_num[10];
	int num_particle_species;
	double total_stellar_mass, total_stellar_initial_mass;

#ifdef STARFORM
	FILE *stellar_input;
	int num_stars;
	double st, sa;
	int first_star, first_star_index, num_stars_to_read;
	long seek_amount, var_first;
	float *pw, *pw0, *tbirth, *zstII, *zstIa;
#ifdef ENRICH
#ifdef ENRICH_SNIa
	#define num_star_variables	5
#else
	#define num_star_variables	4
#endif /* ENRICH_SNIa */
#else
	#define num_star_variables	3
#endif /* ENRICH */
#endif /* STARFORM */

	read_particle_header( header_filename, &header, &endian, &nbody_flag );

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
		particle_species_mass[i] = header.mass[i];
		particle_species_indices[i+1] = header.num[i];
		particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
	}

	num_particles_total = particle_species_indices[num_particle_species];
	num_parts_per_page = header.Nrow*header.Nrow;
	num_pages = (num_particles_total-1) / num_parts_per_page + 1;

#ifdef STARFORM
	if ( stellar_filename == NULL ) {
		/* ICs don't include stars, initialize species to 0 */
		if ( num_particle_species+1 > MAX_PARTICLE_SPECIES ) {
			cart_error("header.Nspecies > MAX_PARTICLE_SPECIES.  Increase and rerun.");
		}
		particle_species_indices[num_particle_species+1] = particle_species_indices[num_particle_species];
		particle_species_num[num_particle_species] = 0;
		particle_species_mass[num_particle_species] = 0.0;
		num_particle_species++;

		total_stellar_mass = 0.0;
		total_stellar_initial_mass = 0.0;
	} else {
		stellar_input = fopen( stellar_filename, "r" );
		if ( stellar_input == NULL ) {
			cart_error("Unable to open file %s for reading.", stellar_filename );
		}

		/* read in header */
		fread( &size, sizeof(int), 1, stellar_input );
		if ( size != 2*sizeof(double) ) {
			reorder( (char *)&size, sizeof(int) );

			if ( size != 2*sizeof(double) ) {
				cart_error("Error reading from %s.", stellar_filename );
			}

			stellar_endian = 1;
		} else {
			stellar_endian = 0;
		}

		fread( &st, sizeof(double), 1, stellar_input );
		fread( &sa, sizeof(double), 1, stellar_input );
		fread( &size, sizeof(int), 1, stellar_input );

		fread( &size, sizeof(int), 1, stellar_input );
		fread( &num_stars, sizeof(int), 1, stellar_input );
		fread( &size, sizeof(int), 1, stellar_input );

		if ( stellar_endian ) {
			reorder( (char *)&num_stars, sizeof(int) );
		}

		if ( num_stars != particle_species_num[num_particle_species-1] ) {
			cart_error("num_stars in %s doesn't match last particle specie.", stellar_filename );
		}

		fread( &size, sizeof(int), 1, stellar_input );
		fread( &total_stellar_mass, sizeof(double), 1, stellar_input );
		fread( &total_stellar_initial_mass, sizeof(double), 1, stellar_input );
		fread( &size, sizeof(int), 1, stellar_input );

		if ( stellar_endian ) {
			reorder( (char *)&total_stellar_mass, sizeof(double) );
			reorder( (char *)&total_stellar_initial_mass, sizeof(double) );
		}

		pw = cart_alloc( num_parts_per_page * sizeof(float) );
		pw0 = cart_alloc( num_parts_per_page * sizeof(float) );
		tbirth = cart_alloc( num_parts_per_page * sizeof(float) );

#ifdef ENRICH
		zstII = cart_alloc( num_parts_per_page * sizeof(float) );
#ifdef ENRICH_SNIa
		zstIa = cart_alloc( num_parts_per_page * sizeof(float) );
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
	}
#endif /* STARFORM */

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

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
			/* have we reached star particles yet? */
			if ( current_id+num_parts_in_page-1 >= particle_species_indices[num_particle_species-1] ) {
				/* load page values for stars */
				first_star = max( current_id - particle_species_indices[num_particle_species-1], 0 );
				cart_assert( first_star >= 0 && first_star < num_stars );

				first_star_index = first_star + particle_species_indices[num_particle_species-1] - current_id;
				cart_assert( first_star_index >= 0 && first_star_index < num_parts_per_page );

				num_stars_to_read = min( num_parts_in_page - first_star_index, 
						particle_species_indices[num_particle_species] - first_star );

				cart_assert( num_stars_to_read >= 0 && num_stars_to_read <= num_parts_in_page );

				var_first = 2*sizeof(int) + 2*sizeof(double) +		/* t, a */
					+ 2*sizeof(int) + sizeof(int) +			/* num_stars */
					+ 2*sizeof(int) + 2*sizeof(double) +		/* ws_old, ws_oldi */
					+ sizeof(int) + first_star*sizeof(float);	/* skip over first pw */

				seek_amount = (num_stars - num_stars_to_read)*sizeof(float) + 2*sizeof(int);

				/* weights */
				if ( fseek( stellar_input, var_first, SEEK_SET ) ) {
					cart_error("Error seeking to %u in file %s",
							num_stars*sizeof(float) + sizeof(int), stellar_filename );
				}

				num_read = fread( &pw[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
				if ( num_read != num_stars_to_read ) {
					cart_error("Error reading from stellar file %s", stellar_filename );
				}

				if ( stellar_endian ) {
					for ( j = 0; j < num_stars_to_read; j++ ) {
						reorder( (char *)&pw[first_star_index+j], sizeof(float) );
					}
				}

				/* original weights */
				if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
					cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
				}

				num_read = fread( &pw0[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
				if ( num_read != num_stars_to_read ) {
					cart_error("Error reading from stellar file %s", stellar_filename );
				}

				if ( stellar_endian ) {
					for ( j = 0; j < num_stars_to_read; j++ ) {
						reorder( (char *)&pw0[first_star_index+j], sizeof(float) );
					}
				}

				/* birth times */
				if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
					cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
				}

				num_read = fread( &tbirth[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
				if ( num_read != num_stars_to_read ) {
					cart_error("Error reading from stellar file %s", stellar_filename );
				}

				if ( stellar_endian ) {
					for ( j = 0; j < num_stars_to_read; j++ ) {
						reorder( (char *)&tbirth[first_star_index+j], sizeof(float) );
					}
				}

#ifdef ENRICH
				/* metallicity */
				if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
					cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
				}

				num_read = fread( &zstII[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
				if ( num_read != num_stars_to_read ) {
					cart_error("Error reading from stellar file %s", stellar_filename );
				}

				if ( stellar_endian ) {
					for ( j = 0; j < num_stars_to_read; j++ ) {
						reorder( (char *)&zstII[first_star_index+j], sizeof(float) );
					}
				}

#endif /* ENRICH */
#ifdef ENRICH_SNIa
				if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
					cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
				}

				num_read = fread( &zstIa[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
				if ( num_read != num_stars_to_read ) {
					cart_error("Error reading from stellar file %s", stellar_filename );
				}

				if ( stellar_endian ) {
					for ( j = 0; j < num_stars_to_read; j++ ) {
						reorder( (char *)&zstIa[first_star_index+j], sizeof(float) );
					}
				}
#endif /* ENRICH_SNIa */
			}
		}
#endif /* STARFORM */

		for ( j = 0; j < num_parts_in_page; j++ ) {
			/* convert to our coordinates 0->num_grid */
			x[j] -= grid_shift;
			y[j] -= grid_shift;
			z[j] -= grid_shift;

			if ( nbody_flag ) {
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

			for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
				/* check if particle matches */
				r = compute_distance_periodic( halos->list[ihalo].pos, particle.x );

				if ( r < halos->list[ihalo].analysis_radius ) {
					particle.v[0] = vx[j];
					particle.v[1] = vy[j];
					particle.v[2] = vz[j];

#ifdef STARFORM	
					if ( current_type == num_particle_species-1 ) {
						particle.is_star = 1;
	
						particle.mass = pw[j];
						particle.initial_mass = pw0[j];
						particle.star_tbirth = tbirth[j];
						particle.star_metallicity_II = zstII[j];
						particle.star_metallicity_Ia = zstIa[j];
						
					} else {
						particle.is_star = 0;
						particle.mass = particle_species_mass[current_type];
					}
#else
					particle.is_star = 0;
					particle.mass = particle_species_mass[current_type];
#endif

					callback( &halos->list[ihalo], &particle );
				}
			}

			current_id++;
		}
	}
	
#ifdef STARFORM
	if ( stellar_filename != NULL ) {
#ifdef ENRICH_SNIa
		cart_free( zstIa );
#endif /* ENRICH_SNIa */
#ifdef ENRICH
		cart_free( zstII );
#endif /* ENRICH */

		cart_free( tbirth );
		cart_free( pw0 );
		cart_free( pw );
	}
#endif /* STARFORM */
	cart_free( input_page );
	fclose(input);
}

typedef struct {
	int id;
	particle_float x[nDim];
	particle_float v[nDim];
} particle_indexed_struct;

void read_indexed_particles( char *header_filename, char *data_filename, char *stellar_filename,
			halo_list *halos, halo_list *subhalos, void callback( halo_struct *, particle_struct * ) ) {
	int i, j, k;
	char desc[47];
	int proc;
	int num_parts;
	particle_struct particle;
	particle_indexed_struct *input_particles;
	int isub;
	halo_list *restricted_subhalos;
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
	int size, endian;
	FILE *input;
	particle_header header;
	int nbody_flag;

	long num_particles_total;
	double particle_species_mass[10];
	int particle_species_indices[11];
	int particle_species_num[10];
	int num_particle_species;
	double total_stellar_mass, total_stellar_initial_mass;

	int current_type;
	int count;
	int num_particles;
	int num_star_particles;
	int num_to_read;
	int coords2[nDim];
	int *sfc_list;
	int dx;
	long *root_index;

#ifdef STARFORM
	FILE *stellar_input;
	float *star_vars;
	long *stellar_root_index;
#endif /* STARFORM */

	read_particle_header( header_filename, &header, &endian, &nbody_flag );
	init_sfc();

	num_particle_species = header.Nspecies;

	particle_species_indices[0] = 0;
	for ( i = 0; i < num_particle_species; i++ ) {
		particle_species_mass[i] = header.mass[i];
		particle_species_indices[i+1] = header.num[i];
		particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
	}


	/* start loading actual particle data */
	input = fopen( data_filename, "r" );
	if ( input == NULL ) {
		cart_error( "Unable to open particle file %s for reading!", data_filename );
	}

	root_index = cart_alloc( (long)num_grid*(long)num_grid*(long)num_grid*sizeof(long) );
	fread( root_index, sizeof(long), (long)num_grid*(long)num_grid*(long)num_grid, input );

#ifdef STARFORM
	stellar_input = fopen( stellar_filename, "r" );
	if ( stellar_input == NULL ) {
		cart_error("Unable to open stellar file %s for reading!", stellar_filename );
	}

	stellar_root_index = cart_alloc( (long)num_grid*(long)num_grid*(long)num_grid* sizeof(long) );
	fread( stellar_root_index, sizeof(long), (long)num_grid*(long)num_grid*(long)num_grid, stellar_input );
#endif /* STARFORM */

	restricted_subhalos = cart_alloc( sizeof(halo_list) );

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		restricted_subhalos->num_halos = 0;
		for ( isub = 0; isub < subhalos->num_halos; isub++ ) {
			if ( halos->list[ihalo].id != subhalos->list[isub].id ) {
				r = compute_distance_periodic( halos->list[ihalo].pos, subhalos->list[isub].pos );

				if ( r < rbinmax + subhalo_exclude_radius*subhalos->list[isub].rhalo ) {
					restricted_subhalos->num_halos++;
				}
			}
		}

		restricted_subhalos->list = cart_alloc( restricted_subhalos->num_halos * sizeof(halo_struct) );

		i = 0;
		for ( isub = 0; isub < subhalos->num_halos; isub++ ) {
			if ( halos->list[ihalo].id != subhalos->list[isub].id ) {
				r = compute_distance_periodic( halos->list[ihalo].pos, subhalos->list[isub].pos );

				if ( r < rbinmax + subhalo_exclude_radius*subhalos->list[isub].rhalo ) {
					restricted_subhalos->list[i++] = subhalos->list[isub];
				}
			}
		}

		/* create list of sfc indices */
		for ( i = 0; i < nDim; i++ ) {
			coords[i] = (int)(halos->list[ihalo].pos[i]);
		}

		dx = (int)((double)num_grid*halos->list[ihalo].analysis_radius/Lbox) + 1;
		count = 0;

		sfc_list = cart_alloc((2*dx+1)*(2*dx+1)*(2*dx+1) * sizeof(int) );

		for ( i = coords[0]-dx; i <= coords[0]+dx; i++ ) {
			if ( i < 0 ) {
				coords2[0] = i+num_grid;
			} else if ( i >= num_grid ) {
				coords2[0] = i-num_grid;
			} else {
				coords2[0] = i;
			}

			for ( j = coords[1]-dx; j <= coords[1]+dx; j++ ) {
				if ( j < 0 ) {
					coords2[1] = j+num_grid;
				} else if ( j >= num_grid ) {
					coords2[1] = j-num_grid;
				} else {
					coords2[1] = j;
				}

				for ( k = coords[2]-dx; k <= coords[2]+dx; k++ ) {
					if ( k < 0 ) {
						coords2[2] = k+num_grid;
					} else if ( k >= num_grid ) {
						coords2[2] = k-num_grid;
					} else {
						coords2[2] = k;
					} 

					index = sfc_index(coords2);
					cart_assert( index >= 0 && index < max_sfc_index );
					sfc_list[count++] = index;
				}
			}
		}

		/* sort list */
		qsort( sfc_list, count, sizeof(int), compare_indices );

		cart_debug("selected %u root trees for halo %d...", count, halos->list[ihalo].id  );

		input_particles = cart_alloc( 1024*1024*sizeof(particle_indexed_struct) );

#ifdef STARFORM
		star_vars = cart_alloc( 5*1024*1024*sizeof(float) );
#endif /* STARFORM */

		for ( sfc = 0; sfc < count; sfc++ ) {
			cart_assert( sfc_list[sfc] >= 0 && sfc_list[sfc] < max_sfc_index );

			fseek( input, root_index[sfc_list[sfc]], SEEK_SET );

			fread( &num_particles, sizeof(int), 1, input );


#ifdef STARFORM
			fseek( stellar_input, stellar_root_index[sfc_list[sfc]], SEEK_SET );
			fread( &num_star_particles, sizeof(int), 1, stellar_input );
			cart_assert( num_star_particles <= num_particles );
#else
			num_star_particles = 0;
#endif /* STARFORM */

			/* read non-star particles */
			current_type = 0;
			particle.is_star = 0;

			while ( num_particles > num_star_particles ) {
				num_to_read = min( 1024*1024, min( num_particles, num_particles-num_star_particles ) );
				fread( input_particles, sizeof(particle_indexed_struct), 
						num_to_read, input );

				for ( ipart = 0; ipart < num_to_read; ipart++ ) {
					particle.x[0] = input_particles[ipart].x[0];
					particle.x[1] = input_particles[ipart].x[1];
					particle.x[2] = input_particles[ipart].x[2];
					particle.v[0] = input_particles[ipart].v[0];
					particle.v[1] = input_particles[ipart].v[1];
					particle.v[2] = input_particles[ipart].v[2];

					if ( input_particles[ipart].id >= particle_species_indices[current_type+1] ) {
						current_type++;
					}

					particle.mass = particle_species_mass[current_type];

					particle.subhalo_flag = 0;
					for ( isub = 0; isub < restricted_subhalos->num_halos; isub++ ) {
						r = compute_distance_periodic( particle.x, restricted_subhalos->list[isub].pos );

						if ( r < subhalo_exclude_radius*restricted_subhalos->list[isub].rhalo) {
							particle.subhalo_flag = 1;
							break;
						}
					}

					callback( &halos->list[ihalo], &particle );
				}

				num_particles -= num_to_read;
			}

#ifdef STARFORM
			cart_assert( num_particles == num_star_particles );

			particle.is_star = 1;
			while ( num_star_particles > 0 ) {
				num_to_read = min( 1024*1024, num_star_particles );
                                fread( input_particles, sizeof(particle_indexed_struct), 
                                                num_to_read, input );
				fread( star_vars, 5*sizeof(float), num_to_read, stellar_input );

				for ( ipart = 0; ipart < num_to_read; ipart++ ) {
					particle.x[0] = input_particles[ipart].x[0];
					particle.x[1] = input_particles[ipart].x[1];
					particle.x[2] = input_particles[ipart].x[2];
					particle.v[0] = input_particles[ipart].v[0];
					particle.v[1] = input_particles[ipart].v[1];
					particle.v[2] = input_particles[ipart].v[2];

					particle.mass = star_vars[5*ipart];
					particle.initial_mass = star_vars[5*ipart+1];
					particle.star_tbirth = star_vars[5*ipart+2];
					particle.star_metallicity_II = star_vars[5*ipart+3];
					particle.star_metallicity_Ia = star_vars[5*ipart+4];

					particle.subhalo_flag = 0;
					for ( isub = 0; isub < restricted_subhalos->num_halos; isub++ ) {
						r = compute_distance_periodic( particle.x, restricted_subhalos->list[isub].pos );

						if ( r < subhalo_exclude_radius*restricted_subhalos->list[isub].rhalo) {
							particle.subhalo_flag = 1;
							break;
						}
                                        }

					callback( &halos->list[ihalo], &particle );
				}
	
				num_star_particles -= num_to_read;
			}
#endif /* STARFORM */

		}

		cart_free( input_particles );

#ifdef STARFORM
		cart_free( star_vars );
#endif /* STARFORM */

		cart_free( sfc_list );
	}

	cart_free( restricted_subhalos );

	cart_debug("done with halos");

	cart_free( root_index );

#ifdef STARFORM
	cart_free( stellar_root_index );
	fclose( stellar_input );
#endif /* STARFORM */

	cart_debug("done");	
	fclose(input);
}

#endif /* PARTICLES */

void read_indexed_grid( char *filename, halo_list *halos, halo_list *subhalos, void callback( halo_struct *, cell_struct * ) ) {
	int i, j, k, m, n;
	FILE *input;
	int size;
	int endian;
	int step;
	int sfc;
	int ihalo;
	double aexpn;
	double t, dt;
	float adum, ainit;
	float Lbox, Om0, Oml0, Omb0, h;
	char job[256];
	int nextras;
	float extra[10];
	char lextra[10][256];
	int maxlevel, minlevel;
	double tl[10], dtl[10], tl_old[10], dtl_old[10];
	int level_sweep_dir[10];
	int icell, ncell0;
	int count;
	halo_list *restricted_subhalos;
	float refinement_volume_min[nDim], refinement_volume_max[nDim];
	float star_formation_volume_min[nDim], star_formation_volume_max[nDim];
	long lsize;
	int level;
	long level_count, next_level_count;
	int num_in_page;
	int local_file_root_cells;
	long total_cells;
	int index;
	int num_cells;
	int num_file_vars;
	int current_root_index;
	int *refined;
	float *vars;

	cell_struct cell;
	int density_index;
	int energy_index;
	int pressure_index;
	int internal_energy_index;
	int momentum_x_index;
	int momentum_y_index;
	int momentum_z_index;
	int metallicity_II_index;
	int metallicity_Ia_index;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	int electron_internal_energy_index;
#endif 	

	long *root_index;
	int *root_level_count;
	int *root_next_level_count;
	int cells_per_level[10];
	int *cell_refined;
	int cell_refined_to_read;
	int num_read;
	char varname[32];
	int dim;
	float *pos_level;
	float *pos_next_level;
	int coords[nDim];
	int coords2[nDim];
	int *sfc_list;
	int dx;

	int isub;
	double r;
	
	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s for reading!", filename );
	}

	fread(&size, sizeof(int), 1, input );
	endian = 0;
	if ( size != 256*sizeof(char) ) {
		reorder( (char *)&size, sizeof(int) );
		if ( size != 256*sizeof(char) ) {
			cart_error("Error: file %s is corrupted", filename );
		} else {
			endian = 1;
			cart_debug("Reordering bytes (file endianness is opposite program)");
		}
	}

	fread(job, sizeof(char), 256, input );
	fread(&size, sizeof(int), 1, input );

	/* istep, t, dt, adum, ainit */
	fread( &size, sizeof(int), 1, input );
	fread( &step, sizeof(int), 1, input );
	fread( &t, sizeof(double), 1, input );
	fread( &dt, sizeof(double), 1, input );
	fread( &adum, sizeof(float), 1, input );
	fread( &ainit, sizeof(float), 1, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&step, sizeof(int) );
		reorder( (char *)&t, sizeof(double) );
		reorder( (char *)&dt, sizeof(double) );
		reorder( (char *)&adum, sizeof(float) );
		reorder( (char *)&ainit, sizeof(float) );
	}

	/* boxh, Om0, Oml0, Omb0, hubble */
	fread( &size, sizeof(int), 1, input );
	fread( &Lbox, sizeof(float), 1, input );
	fread( &Om0, sizeof(float), 1, input );
	fread( &Oml0, sizeof(float), 1, input );
	fread( &Omb0, sizeof(float), 1, input );
	fread( &h, sizeof(float), 1, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&Lbox, sizeof(float) );
		reorder( (char *)&Om0, sizeof(float) );
		reorder( (char *)&Oml0, sizeof(float) );
		reorder( (char *)&Omb0, sizeof(float) );
		reorder( (char *)&h, sizeof(float) );
	}

	size = 5*sizeof(float);

	/* nextra (no evidence extras are used...) extra lextra */
	fread( &size, sizeof(int), 1, input );
	fread( &nextras, sizeof(int), 1, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&nextras, sizeof(int) );
	}

	/* extra */
	fread( &size, sizeof(int), 1, input );
	fread( extra, sizeof(float), nextras, input );
	fread( &size, sizeof(int), 1, input );

	/* lextra */
	fread( &size, sizeof(int), 1, input );
	fread( lextra, 256*sizeof(char), nextras, input );
	fread( &size, sizeof(int), 1, input );

	/* Minlevel, MaxLevelNow */
	fread( &size, sizeof(int), 1, input );
	fread( &minlevel, sizeof(int), 1, input );
	fread( &maxlevel, sizeof(int), 1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&minlevel, sizeof(int) );
		reorder( (char *)&maxlevel, sizeof(int) );
	}

	/* tl */
	fread( &size, sizeof(int), 1, input );
	fread( &tl, sizeof(double), maxlevel-minlevel+1, input );
	fread( &size, sizeof(int), 1, input);

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&tl[i], sizeof(double) );
		}
	}

	/* dtl */
	fread( &size, sizeof(int), 1, input );
	fread( &dtl, sizeof(double), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&dtl[i], sizeof(double) );
		}
	}

	/* tl_old */
	fread( &size, sizeof(int), 1, input );
	fread( &tl_old, sizeof(double), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&tl_old[i], sizeof(double) );
		}
	}

	/* dtl_old */
	fread( &size, sizeof(int), 1, input );
	fread( &dtl_old, sizeof(double), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&dtl_old[i], sizeof(double) );
		}
	}

	/* iSO */
	fread( &size, sizeof(int), 1, input );
	fread( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&level_sweep_dir[i], sizeof(int) );
		}
	}

	/* sfc ordering used */
	fread( &size, sizeof(int), 1, input );
	fread( &sfc_order, sizeof(int), 1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&sfc_order, sizeof(int) );
	}

	/* refinement volume */
	fread( &size, sizeof(int), 1, input );
	fread( refinement_volume_min, sizeof(float), nDim, input );
	fread( refinement_volume_max, sizeof(float), nDim, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = 0; i < nDim; i++ ) {
			reorder( (char *)&refinement_volume_min[i], sizeof(float) );
			reorder( (char *)&refinement_volume_max[i], sizeof(float) );
		}
	}

	/* star formation volume (when included) */
	fread( &size, sizeof(int), 1, input );

	if ( size == 6*sizeof(float) ) {
		fread( star_formation_volume_min, sizeof(float), nDim, input );
		fread( star_formation_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&star_formation_volume_min, sizeof(float) );
				reorder( (char *)&star_formation_volume_max, sizeof(float) );
			}
		}

		/* ncell0 */
		fread( &size, sizeof(int), 1, input );
	}

	/* read in variable strings */
	fread( &num_file_vars, sizeof(int), 1, input );
	fread( &size, sizeof(int), 1, input );

	cart_debug("num_file_vars = %d", num_file_vars );

	for ( i = 0; i < num_file_vars; i++ ) {
		fread( &size, sizeof(int), 1, input );
		fread( varname, sizeof(char), 32, input );
		fread( &size, sizeof(int), 1, input );

		cart_debug("file_var[%d] = %s", i, varname );

		if ( strcmp( varname, "hydro_gas_density" ) == 0 ) {
			density_index = i;
		} else if ( strcmp( varname, "hydro_gas_energy" ) == 0 ) {
			energy_index = i;
		} else if ( strcmp( varname, "hydro_momentum_x" ) == 0 ) {
			momentum_x_index = i;
		} else if ( strcmp( varname, "hydro_momentum_y" ) == 0 ) {
			momentum_y_index = i;
		} else if ( strcmp( varname, "hydro_momentum_z" ) == 0 ) {
			momentum_z_index = i;
		} else if ( strcmp( varname, "hydro_gas_pressure" ) == 0 ) {
			pressure_index = i;
		} else if ( strcmp( varname, "hydro_gas_internal_energy" ) == 0 ) {
			internal_energy_index = i;
		} else if ( strcmp( varname, "hydro_metallicity_II" ) == 0 ) {
			metallicity_II_index = i;
		} else if ( strcmp( varname, "hydro_metallicity_Ia" ) == 0 ) {
			metallicity_Ia_index = i;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
		} else if ( strcmp( varname, "hydro_electron_internal_energy" ) == 0 ) {
			electron_internal_energy_index = i;
#endif
		}
	}

	fread( &size, sizeof(int), 1, input );
        fread( &ncell0, sizeof(int), 1, input);
        fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&ncell0, sizeof(int) );
	}

	init_sfc();

	root_index = cart_alloc( ncell0 * sizeof(long) );

	fread( &size, sizeof(int), 1, input );
	fread( root_index, sizeof(long), ncell0, input );
	fread( &size, sizeof(int), 1, input );

	restricted_subhalos = (halo_list *)cart_alloc( sizeof(halo_list) );

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {

		restricted_subhalos->num_halos = 0;
		for ( isub = 0; isub < subhalos->num_halos; isub++ ) {
			r = compute_distance_periodic( halos->list[ihalo].pos, subhalos->list[isub].pos );
		
			if ( r < rbinmax + subhalo_exclude_radius*subhalos->list[isub].rhalo ) {
				restricted_subhalos->num_halos++;
			}	
		}

		restricted_subhalos->list = cart_alloc( restricted_subhalos->num_halos * sizeof(halo_struct) );

		i = 0;
		for ( isub = 0; isub < subhalos->num_halos; isub++ ) {
			r = compute_distance_periodic( halos->list[ihalo].pos, subhalos->list[isub].pos );

			if ( r < rbinmax + subhalo_exclude_radius*subhalos->list[isub].rhalo ) {
				restricted_subhalos->list[i++] = subhalos->list[isub];
			}
		}

		/* create list of sfc indices */
		for ( i = 0; i < nDim; i++ ) {
			coords[i] = (int)(halos->list[ihalo].pos[i]);
		}

		dx = (int)((double)num_grid*halos->list[ihalo].analysis_radius/Lbox) + 1;
		count = 0;

		sfc_list = cart_alloc((2*dx+1)*(2*dx+1)*(2*dx+1) * sizeof(int) );

		for ( i = coords[0]-dx; i <= coords[0]+dx; i++ ) {
			if ( i < 0 ) {
				coords2[0] = i+num_grid;
			} else if ( i >= num_grid ) {
				coords2[0] = i-num_grid;
			} else {
				coords2[0] = i;
			}

			for ( j = coords[1]-dx; j <= coords[1]+dx; j++ ) {
				if ( j < 0 ) {
					coords2[1] = j+num_grid;
				} else if ( j >= num_grid ) {
					coords2[1] = j-num_grid;
				} else {
					coords2[1] = j;
				}

				for ( k = coords[2]-dx; k <= coords[2]+dx; k++ ) {
					if ( k < 0 ) {
						coords2[2] = k+num_grid;
					} else if ( k >= num_grid ) {
						coords2[2] = k-num_grid;
					} else {
						coords2[2] = k;
					} 

					index = sfc_index(coords2);
					sfc_list[count++] = index;
				}
			}
		}

		/* sort list */
		qsort( sfc_list, count, sizeof(int), compare_indices );

		cart_debug("selected %u root trees for halo %d...", count, halos->list[ihalo].id  );

		for ( sfc = 0; sfc < count; sfc++ ) {
			fseek( input, root_index[sfc_list[sfc]], SEEK_SET );

			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
				if ( size != (maxlevel-minlevel+1)*sizeof(int)) {
					cart_error("Error, file is corrupt!");
				}
			}

			fread( cells_per_level, sizeof(int), maxlevel-minlevel+1, input );
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( level = minlevel; level <= maxlevel; level++ ) {
					reorder( (char *)&cells_per_level[level], sizeof(int) );
				}
			}

			total_cells = 0;
			for ( level = 0; level <= maxlevel; level++ ) {
				total_cells += cells_per_level[level];
			}

			cell_refined_to_read = total_cells-cells_per_level[maxlevel];

			cell_refined = cart_alloc( cell_refined_to_read*sizeof(int) );
			vars = cart_alloc( num_file_vars*total_cells*sizeof(float) );

			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );

				if ( size != cell_refined_to_read*sizeof(int)) {
					cart_error("Error: file is corrupt reading cell_refined!");
				}
			}
			fread( cell_refined, sizeof(int), cell_refined_to_read, input );
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( i = 0; i < cell_refined_to_read; i++ ) {
					reorder( (char *)&cell_refined[i], sizeof(int) );
				}
			}

			fread( &size, sizeof(int), 1, input );
			fread( vars, sizeof(float), num_file_vars*total_cells, input );
			fread( &size, sizeof(int), 1, input );

			pos_level = cart_alloc( nDim*sizeof(float) );

			sfc_coords( sfc_list[sfc], coords );

			for ( dim = 0; dim < nDim; dim++ ) {
				pos_level[dim] = (double)coords[dim]+0.5;
			}

			icell = 0;

			for ( level = minlevel; level <= maxlevel; level++ ) {
				if ( level < maxlevel ) {
					pos_next_level = cart_alloc( nDim*cells_per_level[level+1]*sizeof(float) );
					next_level_count = 0;
				} else {
					pos_next_level = NULL;
				}

				for ( j = 0; j < cells_per_level[level]; j++ ) {
					if ( level == maxlevel || !cell_refined[icell+j] ) {
						cell.level = level;
						for ( k = 0; k < nDim; k++ ) {
							cell.pos[k] = pos_level[nDim*j+k];
						}

						cell.subhalo_flag = 0;
						for ( isub = 0; isub < restricted_subhalos->num_halos; isub++ ) {
							r = compute_distance_periodic( cell.pos, restricted_subhalos->list[isub].pos );

							if ( r < restricted_subhalos->list[isub].rhalo) {
								cell.subhalo_flag = 1;
								break;
							}
						}

						cell.gas_density = vars[density_index*total_cells + icell + j];
						cell.gas_energy = vars[energy_index*total_cells + icell + j];
						cell.gas_pressure = vars[pressure_index*total_cells + icell + j];
						cell.gas_internal_energy = vars[internal_energy_index*total_cells + icell + j];
						cell.momentum[0] = vars[momentum_x_index*total_cells + icell + j];
						cell.momentum[1] = vars[momentum_y_index*total_cells + icell + j];
						cell.momentum[2] = vars[momentum_z_index*total_cells + icell + j];

#ifdef ENRICH
						cell.metallicity_II = vars[metallicity_II_index*total_cells + icell + j];
#ifdef ENRICH_SNIa
						cell.metallicity_Ia = vars[metallicity_Ia_index*total_cells + icell + j];
#endif
#endif

#ifdef ELECTRON_ION_NONEQUILIBRIUM
						cell.electron_internal_energy = vars[electron_internal_energy_index*total_cells + icell + j];
#endif

						cell.subhalo_flag = 0;
						for ( isub = 0; isub < subhalos->num_halos; isub++ ) {
							r = compute_distance_periodic( cell.pos, subhalos->list[isub].pos );

							if ( r < subhalos->list[isub].rhalo) {
								cell.subhalo_flag = 1;
								break;
							}
						}

						/*
						r = compute_distance_periodic( cell.pos, halos->list[ihalo].pos );
						if ( r > 3.0/r0 && cell.gas_density > 100.0 * Omegab0 / Omega0 ) {
							cell.subhalo_flag = 1;
						}
						*/

						/* callback here */
						callback( &halos->list[ihalo], &cell );

					} else if ( cells_per_level[level+1] > 0 ) {
						for ( k = 0; k < 8; k++ ) {
							for ( dim = 0; dim < nDim; dim++ ) {
								pos_next_level[nDim*next_level_count+nDim*k+dim] =
									pos_level[nDim*j+dim] + 0.5*cell_size[level]*cell_delta[k][dim];
							}
						}

						next_level_count += 8;
					}
				}

				icell += cells_per_level[level];

				cart_free(pos_level);
				pos_level = pos_next_level;
			}

			cart_free( cell_refined );
			cart_free( vars );
		}

		cart_free( sfc_list );
		cart_free( restricted_subhalos->list );
	}

	cart_free(restricted_subhalos);

	cart_free( root_index );
	fclose(input);
}

