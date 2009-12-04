#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include "sfc.h"

#define min(x,y)        (((x) < (y)) ? (x): (y))
#define max(x,y)        (((x) > (y)) ? (x): (y))

#define MAX_PARTICLE_SPECIES	10
typedef float particle_float;

typedef struct {
	float aexpn;
	float aexp0;
	float amplt;
	float astep;
	int   istep;
	float partw;
	float tintg;
	float ekin;
	float ekin1;
	float ekin2;
	float au0;
	float aeu0;
	int   Nrow;
	int   Ngrid;
	int   Nspecies;
	int   Nseed;
	float Om0;
	float Oml0;
	float hubble;
	float Wp5;
	float Ocurv;
	float Omb0;  
	float mass[10];
	unsigned int   num[10];
	float fill[80];
} particle_header;

typedef struct {
	float aexpn;
	float aexp0;
	float amplt;
	float astep;
	int   istep;
	float partw;
	float tintg;
	float ekin;
	float ekin1;
	float ekin2;
	float au0;
	float aeu0;
	int   Nrow;
	int   Ngrid;
	int   Nspecies;
	int   Nseed;
	float Om0;
	float Oml0;
	float hubble;
	float Wp5;
	float Ocurv;
	float mass[10];
	unsigned int   num[10];
	float fill[80];
} nbody_particle_header;

typedef struct {
	int id;
	particle_float x[3];
	particle_float v[3];
} particle_struct;

int *pid;
int *pindex;

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
	int a1 = *(int *)a;
	int b1 = *(int *)b;

	if ( pindex[a1] == pindex[b1] ) {
		return pid[a1] - pid[b1];
	} else {
		return pindex[a1] - pindex[b1];
	}
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
		fprintf(stderr,"Unable to open particle file %s\n", header_filename );
	}

	fread( &size, sizeof(int), 1, input );
	fread( desc, sizeof(char), 45, input );
	desc[45] = '\0';

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
					fprintf(stderr, "Size mismatch in reading particle file header %s (%u vs %lu)\n",
							header_filename, size, sizeof(particle_header)+45 );
				}
			} else {
				*endian = 1;
			}
		}
	}

	if ( *nbody_flag ) {
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

int main ( int argc, char *argv[]) {
	int i, j, k;
	int coords[nDim];
	long num_particles_total;
	int num_parts_in_page, num_parts_per_page;
        int num_parts_per_proc_page;
	int num_particle_species;
	float particle_species_mass[MAX_PARTICLE_SPECIES];
	int particle_species_num[MAX_PARTICLE_SPECIES];
	int particle_species_indices[MAX_PARTICLE_SPECIES+1];
        int num_pages, index;
	int current_id, current_type;
	particle_struct *particles;
	particle_header header;
	int nbody_flag;
	int endian;
	int count;
	particle_float *input_page, *x, *y, *z, *vx, *vy, *vz;
	long *file_index;
	int *order;
	FILE *input, *output;
	int *root_particle_count;
        int num_read;
        int size;
        float vfact, grid_shift;

	FILE *stellar_input, *stellar_output;
	int stellar_endian;
	int *root_star_count, *stellar_root_index;
	long *star_file_index;
	int *root_index;
        int num_stars;
        double st, sa;
	double total_stellar_mass, total_stellar_initial_mass;
        float *pw, *pw0, *tbirth,*zstII, *zstIa;
	float *star_vars;
	//#ifdef ENRICH
	//#ifdef ENRICH_SNIa
	//        #define num_star_variables      5
	//#else
	//        #define num_star_variables      4
	//#endif /* ENRICH_SNIa */
	//#else
	//        #define num_star_variables      3
	//#endif /* ENRICH */
	int enrich_flag;
	int starform_flag;
	int num_star_variables;

	if ( argc == 4){
	  fprintf(stderr,"no starformation/enrichment\n");
	  enrich_flag = 0;
	  starform_flag = 0;
	  
	}else if( argc == 7 ) {
	  fprintf(stderr,"starformation\n");
	  starform_flag = 1;

	  enrich_flag = atoi( argv[6] );
	  fprintf(stderr,"enrichment by %d(0=NONE,1=SNII,2=SNII+SNIA)\n",enrich_flag);
	  if ( enrich_flag > 2 || enrich_flag < 0){ fprintf(stderr,"bad enrich flag! %d",enrich_flag);exit(-1);}


	}else{
		fprintf(stderr,"Usage: ./cart_create_indexed_particles input_dph input_dxv  output_dxv_indexed input_dst output_dst_indexed enrich(0/1/2)\n");
		fprintf(stderr,"enrichment by %d(0=NONE,1=SNII,2=SNII+SNIA)\n",enrich_flag);
		exit(-1);
	}
	
	num_star_variables = 3+enrich_flag;
	
	read_particle_header( argv[1], &header, &endian, &nbody_flag );
	if ( nbody_flag ) {
		vfact = 2.0/sqrt(header.Om0);
		grid_shift = 1.5;
	} else {
		vfact = 1.0;
		grid_shift = 1.0;
	}

	num_grid = header.Ngrid;
	num_particle_species = header.Nspecies;
	for ( i = 0; i < num_particle_species; i++ ) {
		particle_species_mass[i] = header.mass[i];
		particle_species_indices[i+1] = header.num[i];
		particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
	}

	num_particles_total = particle_species_indices[num_particle_species];
	num_parts_per_page = header.Nrow*header.Nrow;
#ifdef ADHOC_PAGESIZE 
	num_parts_per_page = 1024*1024;
#endif
	num_pages = (num_particles_total-1) / num_parts_per_page + 1;

	init_sfc();



	root_particle_count = malloc( num_grid*num_grid*num_grid*sizeof(int) );
	for ( i = 0; i < num_grid*num_grid*num_grid; i++ ) {
		root_particle_count[i] = 0;
	}

	if ( starform_flag == 1 ) {
		root_star_count = malloc( num_grid*num_grid*num_grid*sizeof(int) );
		for ( i = 0; i < num_grid*num_grid*num_grid; i++ ) {
			root_star_count[i] = 0;
		}
	}
	
	input_page = malloc( 2*nDim*num_parts_per_page*sizeof(particle_float) );
        pindex = malloc( num_particles_total*sizeof(int) );
	root_index = malloc( num_particles_total*sizeof(int) );
	particles = malloc( num_particles_total*sizeof(particle_struct) );
	if ( input_page == NULL || pindex == NULL || particles == NULL || root_index == NULL ) {
		fprintf(stderr,"Ran out of memory allocating buffers!\n" );
		exit(1);
	}

	x = input_page;
	y = &input_page[num_parts_per_page];
	z = &input_page[2*num_parts_per_page];
	vx = &input_page[3*num_parts_per_page];
	vy = &input_page[4*num_parts_per_page];
	vz = &input_page[5*num_parts_per_page];

	/* first scan through particle file and count particles in sfc cells */	
	input = fopen( argv[2], "r" );
	if ( input == NULL ) {
		fprintf( stderr, "Unable to open particle file %s for reading!\n", argv[2] );
		exit(1);
	}

	current_id = 0;
	current_type = 0;

	for ( i = 0; i < num_pages; i++ ) {
               if ( i % 10 == 0 ) {
			printf("page %d/%d\n", i, num_pages );
		}

		if ( i == num_pages - 1 ) {
			num_parts_in_page = num_particles_total - num_parts_per_page*(num_pages-1);
		} else {
			num_parts_in_page = num_parts_per_page;
		}

		num_read = fread( input_page, sizeof(particle_float), 2*nDim*num_parts_per_page, input );
		if ( num_read != 2*nDim*num_parts_per_page ) {
			fprintf(stderr,"Error reading from particle file %s: insufficient data!\n", argv[2] );
			exit(1);
		}

		if ( endian ) {
			for ( j = 0; j < num_parts_in_page; j++ ) {
				reorder( (char *)&x[j], sizeof(particle_float) );
				reorder( (char *)&y[j], sizeof(particle_float) );
				reorder( (char *)&z[j], sizeof(particle_float) );
			}
		}

		for ( j = 0; j < num_parts_in_page; j++ ) {
			/* convert to our coordinates 0->num_grid */
			x[j] -= grid_shift;
			y[j] -= grid_shift;
			z[j] -= grid_shift;

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

			coords[0] = (int)(x[j]);
			coords[1] = (int)(y[j]);
			coords[2] = (int)(z[j]);

			index = sfc_index( coords );	
		
			pindex[current_id] = index;	
			root_particle_count[index]++;

			if ( current_id >= particle_species_indices[current_type+1] ) {
				current_type++;
			}

			if ( starform_flag == 1 && current_type == num_particle_species-1 ) {
				root_star_count[index]++;
			}

			current_id++;
		}
	}

	printf("done with first sweep...\n");

	root_index[0] = 0;
	for ( i = 1; i < num_grid*num_grid*num_grid; i++ ) {
		root_index[i] = root_index[i-1] + root_particle_count[i-1];
	}

	rewind(input);

	current_id = current_type = 0;

	for ( i = 0; i < num_pages; i++ ) {
	        if ( i % 10 == 0 ) {
			printf("Page %d/%d\n", i, num_pages );
		}

		if ( i == num_pages - 1 ) {
			num_parts_in_page = num_particles_total - num_parts_per_page*(num_pages-1);
		} else {
			num_parts_in_page = num_parts_per_page;
		}

		num_read = fread( input_page, sizeof(particle_float), 2*nDim*num_parts_per_page, input );
		if ( num_read != 2*nDim*num_parts_per_page ) {
			fprintf(stderr,"Error reading from particle file %s: insufficient data!\n", argv[2] );
			exit(1);
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

			k = root_index[pindex[current_id]];
			root_index[pindex[current_id]]++;

			particles[k].id = current_id;
			particles[k].x[0] = x[j];
			particles[k].x[1] = y[j];
			particles[k].x[2] = z[j];
			particles[k].v[0] = vx[j];
			particles[k].v[1] = vy[j];
			particles[k].v[2] = vz[j];

			current_id++;
		}
	}

	fclose(input);

	/* now that we know how many particles start writing output file */
	output = fopen( argv[3], "w" );
	if ( output == NULL ) {
		fprintf(stderr, "Unable to open %s for writing!\n", argv[3] );
		exit(1);
	}

	file_index = malloc( (long)num_grid*(long)num_grid*(long)num_grid*sizeof(long) );
	if ( file_index == NULL ) {
		fprintf(stderr, "Ran out of memory allocating file_index!\n" );
		exit(1);
	}

	file_index[0] = (long)num_grid*(long)num_grid*(long)num_grid*sizeof(long);
	for ( i = 1; i < num_grid*num_grid*num_grid; i++ ) {
		file_index[i] = file_index[i-1] + sizeof(int) +
			sizeof(particle_struct)*root_particle_count[i-1];
	}

	fwrite( file_index, sizeof(long), num_grid*num_grid*num_grid, output );

	printf("done writing particle index...\n");

	count = 0;
	for ( index = 0; index < num_grid*num_grid*num_grid; index++ ) {
		if ( (index % 1024*1024) == 0 ) {
		        printf("writing grid index %d...\n", index );
		}
		fwrite( &root_particle_count[index], sizeof(int), 1, output );
		fwrite( &particles[count], sizeof(particle_struct), root_particle_count[index], output );
		count += root_particle_count[index];
	}

	fclose(output);
	free( particles );
	free( root_index );

	if ( starform_flag == 1 ) {
		printf("reading/writing star variables...\n");

		stellar_input = fopen( argv[4], "r" );
		if ( stellar_input == NULL ) {
			fprintf(stderr, "Unable to open file %s for reading.\n", argv[4] );
		}

		/* read in header */
		fread( &size, sizeof(int), 1, stellar_input );
		if ( size != 2*sizeof(double) ) {
			reorder( (char *)&size, sizeof(int) );

			if ( size != 2*sizeof(double) ) {
				fprintf(stderr,"Error reading from %s.\n", argv[4] );
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
			fprintf(stderr,"num_stars in %s doesn't match last particle specie.\n", argv[4] );
			exit(1);
		}
		fprintf(stderr,"num_stars=%d\n",num_stars);

		fread( &size, sizeof(int), 1, stellar_input );
		fread( &total_stellar_mass, sizeof(double), 1, stellar_input );
		fread( &total_stellar_initial_mass, sizeof(double), 1, stellar_input );
		fread( &size, sizeof(int), 1, stellar_input );

		pw = malloc( num_stars * sizeof(float) );
		pw0 = malloc( num_stars * sizeof(float) );
		tbirth = malloc( num_stars * sizeof(float) );
		if(enrich_flag > 0){ //#ifdef ENRICH
		  zstII = malloc( num_stars * sizeof(float) );
		  if(enrich_flag == 2){//#ifdef ENRICH_SNIa
		    zstIa = malloc( num_stars * sizeof(float) );
		  }//#endif
		}//#endif
		
		

		fread( &size, sizeof(int), 1, stellar_input );
		fread( pw, sizeof(float), num_stars, stellar_input );
		fread( &size, sizeof(int), 1, stellar_input );

		fread( &size, sizeof(int), 1, stellar_input );
		fread( pw0, sizeof(float), num_stars, stellar_input );
		fread( &size, sizeof(int), 1, stellar_input );

		fread( &size, sizeof(int), 1, stellar_input );
		fread( tbirth, sizeof(float), num_stars, stellar_input );
		fread( &size, sizeof(int), 1, stellar_input );

		if(enrich_flag > 0){//#ifdef ENRICH
		  fread( &size, sizeof(int), 1, stellar_input );
		  fread( zstII, sizeof(float), num_stars, stellar_input );
		  fread( &size, sizeof(int), 1, stellar_input );
		  if(enrich_flag == 2){//#ifdef ENRICH_SNIa
		    fread( &size, sizeof(int), 1, stellar_input );
		    fread( zstIa, sizeof(float), num_stars, stellar_input );
		    fread( &size, sizeof(int), 1, stellar_input );
		  }//#endif
		}//#endif
		

		if ( stellar_endian ) {
			for ( i = 0; i < num_stars; i++ ) {
				reorder( (char *)&pw[i], sizeof(float) );
				reorder( (char *)&pw0[i], sizeof(float) );
				reorder( (char *)&tbirth[i], sizeof(float) );
				if(enrich_flag > 0){//#ifdef ENRICH
				  reorder( (char *)&zstII[i], sizeof(float) );
				  if(enrich_flag == 2){//#ifdef ENRICH_SNIa
				    reorder( (char *)&zstIa[i], sizeof(float) );
				  }//#endif
				}//#endif
			}
		}

		stellar_root_index = malloc( num_grid*num_grid*num_grid*sizeof(int) );

		stellar_root_index[0] = 0;
		for ( i = 1; i < num_grid*num_grid*num_grid; i++ ) {
		  stellar_root_index[i] = stellar_root_index[i-1] + root_star_count[i-1];
		  //fprintf(stderr,"grid_index=%d stellar_root_index=%d root_star_count %d \n",i,stellar_root_index[i], root_star_count[i-1] );
		  if(root_star_count[i-1]>0){
		    fprintf(stderr,"grid#=%d root_star_count %d...\n",i-1, root_star_count[i-1] );
		  }
		}

		star_vars = malloc( num_star_variables*num_stars*sizeof(float) );

		for ( i = 0; i < num_stars; i++ ) {
			index = pindex[i+particle_species_indices[num_particle_species-1]];
			j = stellar_root_index[index];

			star_vars[num_star_variables*j] = pw[i];
			star_vars[num_star_variables*j+1] = pw0[i];
			star_vars[num_star_variables*j+2] = tbirth[i];
			if(enrich_flag > 0){//#ifdef ENRICH
			  star_vars[num_star_variables*j+3] = zstII[i];
			  if(enrich_flag == 2){//#ifdef ENRICH_SNIa
			    star_vars[num_star_variables*j+4] = zstIa[i];
			  }//#endif
			}//#endif

			stellar_root_index[index]++;
		}


		//--------------------------==================--------------------

		stellar_output = fopen( argv[5], "w");       
		if ( stellar_output == NULL ) {
			fprintf(stderr,"Unable to open %s for writing!\n", argv[5] );
			exit(1);
		}

		star_file_index = malloc( (long)num_grid*(long)num_grid*(long)num_grid*sizeof(long) );
		if ( star_file_index == NULL ) {
			fprintf(stderr, "Ran out of memory allocating star_file_index!\n" );
			exit(1);
		}

		star_file_index[0] = (long)num_grid*(long)num_grid*(long)num_grid*sizeof(long);
		for ( i = 1; i < num_grid*num_grid*num_grid; i++ ) {
			star_file_index[i] = star_file_index[i-1] + sizeof(int) +
				num_star_variables*sizeof(float)*root_star_count[i-1];
			}

		fwrite( star_file_index, sizeof(long), num_grid*num_grid*num_grid, stellar_output );

		count = 0;
		for ( index = 0; index < num_grid*num_grid*num_grid; index++ ) {
			if ( index % 1024 == 0 ) {
			        printf("writing grid index %d...%d\n", index, stellar_root_index[index] );
			}
			fwrite( &root_star_count[index], sizeof(int), 1, stellar_output );
			fwrite( &star_vars[count], num_star_variables*sizeof(float), root_star_count[index], stellar_output );
			count += num_star_variables*root_star_count[index];
		}

		free( stellar_root_index );
		free( star_file_index );
		free( pw );
		free( pw0 );
		free( tbirth );
		
		if(enrich_flag > 0){//#ifdef ENRICH
		  free( zstII );
		  if(enrich_flag == 2){  //#ifdef ENRICH_SNIa
		    free( zstIa );
		  }//#endif
		}//#endif
		free( star_vars );

		fclose(stellar_output);
	}

	return 0;
}
