#include "config.h"

#ifdef COSMOLOGY

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "times.h"
#include "units.h"

#include "halos.h"

void load_halo_finder_epochs( char *filename, int *num_epochs, float **epoch ) {
	FILE *input;
	float aexpn;
	float *list;
	int num;
	
	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s", filename );
	}

	/* count number of lines */
	num = 0;
	while ( fscanf( input, "%f", &aexpn ) != EOF ) {
		num++;
	}

	rewind(input);

	list = cart_alloc(float, num );

	num = 0;
	while ( fscanf( input, "%f", &aexpn ) != EOF ) {
		list[num] = aexpn;
		num++;
	}		

	qsort( list, num, sizeof(float), compare_floats );

	*epoch = list;
	*num_epochs = num;
}

halo_list *load_halo_finder_catalog( const char *filename, int nmem_min, float mvir_min, float vmax_min, float rvir_min, int max_num_halos) {
	int i;
	FILE *input;
	char line[256];
	int id;
	int pid;
	float px, py, pz, vx, vy, vz;
	float rvir, rhalo, mvir;
	float vmax, rmax, rs;
	int np, coords[nDim];
	halo_list *halos;
	halo *h;
	float a, OmM, OmL, OmB, h100; 
	int cart_hfind_format_flag = 0;
	
	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s", filename );
	}

	/*
	//  Skip the job name
	 */
	if(fgets(line,1024,input) == NULL) cart_error("Error in reading halo catalog %s", filename );

	/* check for # character, denoting Fortran hfind vs cart hfind */
	if ( line[0] == '#' ) {
		cart_hfind_format_flag = 1;
	}

	/*
	//  Check the scale factor
	 */
	if ( cart_hfind_format_flag ) {
		if(fscanf(input,"# step = %*u, auni = %f, abox = %*f\n", &a) != 1) cart_error("Error in reading halo catalog %s", filename );
	} else {
		if(fscanf(input," A=%f A0=%*f Ampl=%*f Step=%*f\n",&a) != 1) cart_error("Error in reading halo catalog %s", filename );
	}
	if(fabs(a-auni[min_level]) > 1.0e-3) {
		cart_debug("Scalar factor in HLIST file (%f) is different from the current value (%f)",a,auni[min_level]);
	}

	if ( cart_hfind_format_flag ) {
		if(fscanf(input,"# Cosmology: OmM = %f, OmL = %f, OmB = %f, h = %f, DeltaDC = %*f\n",
				&OmM,&OmL,&OmB,&h100) != 4) cart_error("Error in reading halo catalog %s", filename );
	} else {
		if(fgets(line,1024,input) == NULL) cart_error("Error in reading halo catalog %s", filename );
		if(fscanf(input," Nrow=%*d Ngrid=%*d  Omega_0=%f OmLam_0=%f  Omegab_0=%f Hubble=%f\n",&OmM,&OmL,&OmB,&h100) != 4) cart_error("Error in reading halo catalog %s", filename );
	}

	if(fabs(OmM/cosmology->OmegaM-1.0) > 1.0e-2) {
		cart_debug("OmegaM in HLIST file (%f) is different from the current value (%f)",OmM,cosmology->OmegaM);
	}
	if(fabs(OmL/cosmology->OmegaL-1.0) > 1.0e-2) {
		cart_debug("OmegaL in HLIST file (%f) is different from the current value (%f)",OmL,cosmology->OmegaL);
	}
	if(fabs(OmB/cosmology->OmegaB-1.0) > 1.0e-2) {
		cart_debug("OmegaB in HLIST file (%f) is different from the current value (%f)",OmB,cosmology->OmegaB);
	}
	if(fabs(h100/cosmology->h-1.0) > 1.0e-2) {
		cart_debug("Hubble in HLIST file (%f) is different from the current value (%f)",h100,cosmology->h);
	}

	/*
	//  Skip the rest
	 */
	if ( cart_hfind_format_flag ) {
		do {
			if(fgets(line,1024,input) == NULL) cart_error("Error in reading halo catalog %s", filename );
		} while ( line[0] == '#' && line[1] != '#' );
	} else {
		for(i=4; i<17; i++) {
			if(fgets(line,1024,input) == NULL) cart_error("Error in reading halo catalog %s", filename );
		}
	}

	/* count halos */
	halos = halo_list_alloc(100);

	do {
		if ( cart_hfind_format_flag ) {
			if ( fscanf( input, "%u %e %e %e %e %e %e %e %e %u %e %e", &id, &px, &py, &pz,
					&vx, &vy, &vz, &rvir, &mvir, &np, &vmax, &rmax ) != 12 ) {
				break;
			}
		
			rhalo = rvir;
			rs = 0.0;
			pid = 0;
		} else {
			if ( fscanf( input, "%u %e %e %e %e %e %e %e %e %e %u %e %e %e %u",
					&id, &px, &py, &pz, &vx, &vy, &vz, &rvir, &rhalo,
					&mvir, &np, &vmax, &rmax, &rs, &pid ) != 15 ) {
				break;
			}
		}

		if(np>=nmem_min && mvir>=mvir_min && vmax>=vmax_min && rvir>=rvir_min) {
			h = halo_list_add_halo(halos);

			h->id = id;

			/* convert position to code units */
			h->pos[0] = px/units->length_in_chimps;
			h->pos[1] = py/units->length_in_chimps;
			h->pos[2] = pz/units->length_in_chimps;

			h->vel[0] = constants->kms*vx/units->velocity;
			h->vel[1] = constants->kms*vy/units->velocity;
			h->vel[2] = constants->kms*vz/units->velocity;

			h->rvir = rvir/( 1e3*units->length_in_chimps );
			h->rhalo = rhalo/( 1e3*units->length_in_chimps );
			h->mvir = constants->Msun/cosmology->h*mvir/units->mass;
			h->vmax = constants->kms*vmax/units->velocity;
			h->np = np;

			coords[0] = (int)h->pos[0];
			coords[1] = (int)h->pos[1];
			coords[2] = (int)h->pos[2];

			if ( num_procs > 1 ) {
				h->proc = processor_owner( sfc_index( coords ) );
			} else {
				h->proc = local_proc_id;
			}
		}
	} while ( halos->num_halos < max_num_halos );

	fclose(input);

	cart_debug("Done loading halo list; read %d halos.",halos->num_halos);

	return halos;
}

float *binding_energy;

int compare_binding_energy( const void *a, const void *b ) {
	int p1 = *(int *)a;
	int p2 = *(int *)b;

	if ( binding_energy[p1] < binding_energy[p2] ) {
		return -1;
	} else {
		return 1;
	}
}

void load_halo_particle_mapping( char *filename, halo_list *halos ) {
	int i, j;
	FILE *input;
	halo *h;
	int nh, np, ih, nhp;
	int size;
	float aexpn;
	int endian = 0;

	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s", filename );
	}

	fread( &size, sizeof(int), 1, input );

	if ( size != sizeof(float) ) {
		reorder( (char *)&size, sizeof(int) );
	
		if ( size != sizeof(float) ) {
			cart_error("Error reading from file %s\n", filename );
		}

		endian = 1;
	}

	fread( &aexpn, sizeof(float), 1, input );
	if ( endian ) {
		reorder( (char *)&aexpn, sizeof(float) );
	}
	if(fabs(aexpn-auni[min_level]) > 1.0e-3) {
		cart_debug("Scalar factor in halo particle file %s (%f) is different from the current value (%f)",
			filename, aexpn,auni[min_level]);
    }
	fread( &size, sizeof(int), 1, input );

	fread( &size, sizeof(int), 1, input );
	fread( &nh, sizeof(int), 1, input );
	fread( &np, sizeof(int), 1, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&nh, sizeof(int) );
		reorder( (char *)&np, sizeof(int) );
	}

	if ( nh != halos->num_halos ) {
		cart_error("Error: number of halos in %s (%u) don't match provided halo_list (%u)", 
			filename, nh, halos->num_halos );
	}

	for ( i = 0; i < nh; i++ ) {
		fread( &size, sizeof(int), 1, input );
		fread( &ih, sizeof(int), 1, input );
		fread( &nhp, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&ih, sizeof(int) );
			reorder( (char *)&nhp, sizeof(int) );
		}

		h = find_halo_by_id( halos, ih );

		if ( h != NULL ) {
			h->particles = cart_alloc(int, nhp );
			h->binding_order = cart_alloc(int, nhp );

			binding_energy = cart_alloc(float, nhp );

			fread( h->particles, sizeof(int), nhp, input );
			fread( binding_energy, sizeof(float), nhp, input );
			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				for ( j = 0; j < nhp; j++ ) {
					reorder( (char *)&h->particles[j], sizeof(int) );
					reorder( (char *)&binding_energy[j], sizeof(float) );
				}
			}

			for ( j = 0; j < nhp; j++ ) {
				/* convert to CART indexes DHR - hfind now outputs in 0 based indices
				   halos->list[ih].particles[j] -= 1;
				 */
				h->binding_order[j] = j;
			}

			/* sort particles by binding energy */
			qsort( h->binding_order, nhp, sizeof(int), compare_binding_energy );

			cart_free( binding_energy );

			for ( j = 0; j < nhp; j++ ) {
				h->binding_order[j] = h->particles[h->binding_order[j]];
			}

			/* sort by particle index for faster searching */
			qsort( h->particles, nhp, sizeof(int), compare_ints );

			/* force consistency between particle mapping and np */
			h->np = nhp;
		}	
	}
	
	fclose(input);
}

void destroy_halo_list( halo_list *halos ) {
	int ihalo;
	if ( halos->list != NULL ) {
		for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
			if ( halos->list[ihalo].particles != NULL ) {
				cart_free( halos->list[ihalo].particles );
			}
		}

		cart_free( halos->list );
	}

	if ( halos->map != NULL ) cart_free( halos->map );

	cart_free( halos );
}

halo* find_halo_by_id(halo_list *halos, int id)
{
  int ih;

  cart_assert(halos);
  /* DHR - would replace with binary search, but maybe we'll sort list by something else? */
  for(ih=0; ih<halos->num_halos; ih++)
    {
      if(halos->list[ih].id == id) return &halos->list[ih];
    }
  return NULL;
}

halo_list *halo_list_alloc( int size ) {
	halo_list *halos = cart_alloc(halo_list, 1);

	halos->list = cart_alloc(halo, size);
	halos->size = size;
	halos->num_halos = 0;
	halos->map = NULL;

	return halos;
}

halo *halo_list_add_halo( halo_list *halos ) {
	halo *h;
	halo *tmp;

	if(halos->num_halos == halos->size) {
		halos->size *= 2;
		tmp = cart_alloc(halo, halos->size );
		memcpy(tmp,halos->list,halos->num_halos*sizeof(halo));
		cart_free(halos->list);
		halos->list = tmp;
	}

	h = &halos->list[halos->num_halos];
	halos->num_halos++;
	memset( (void *)h, 0, sizeof(halo) );
	h->id = -1;
	h->np = -1;
	h->particles = NULL;
	h->binding_order = NULL;
	h->flag = 0;
	return h;
}	

#endif /* COSMOLOGY */
