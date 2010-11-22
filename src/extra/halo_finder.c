#include "config.h"

#ifdef COSMOLOGY

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cooling.h"
#include "cosmology.h"
#include "iterators.h"
#include "hydro.h"
#include "io.h"
#include "parallel.h"
#include "particle.h"
#include "rt_utilities.h"
#include "sfc.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "halo_finder.h"

#ifdef ANALYZE_XRAY
#include "xrays.h"
#endif


#define rbinmax		(5e3/cosmology->h)
#define rbinmin		(1.0/cosmology->h)
#define max_bins	100

#define rbinvirmax	4.0
#define rbinvirmin	0.01
#define num_vir_bins	80
#define virial_radius_index	2	/* r500, 0 = rvir */

#define	points_per_cell 1

#define Tcold		(1e5)
#define new_star_age	(0.1)

/* Chandra parameters for X-ray spectroscopic fitting (Vikhlinin 2006) */
#define alpha_xray	(0.875)
#define beta_xray	(1.0)
#define delta_xray_1	(0.19)
#define delta_xray_2	(0.25)

int num_radii = 6;
char *radii_label[] = { "rout",
			"r2500c",
			"r500c",
			"r200c",
			"r200m",
			"r180m"	};

float radii_delta[] = {	0.0,
			2500.0,
			500.0,
			200.0,
			200.0,
			180.0 };

char radii_units[] = {	'r',
			'c',
			'c',
			'c',
			'm',
			'm' };

int num_halo_leafs;
int *leaf_index;

int build_leaf_index( int icell, int level ) {
	cart_assert( icell >= 0 && icell < num_cells );

        if ( cell_is_leaf(icell) ) {
		leaf_index[num_halo_leafs++] = icell;
	}

	return 0;
}


void recv_int_bins( int *output_buffer, int num_bins, int proc, int tag ) {
	int bin;
	int recv_buffer[max_bins];

	MPI_Recv( recv_buffer, num_bins, MPI_INT, proc, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

	for ( bin = 0; bin < num_bins; bin++ ) {
		output_buffer[bin] += recv_buffer[bin];
	}
}

void recv_float_bins( float *output_buffer, int num_bins, int proc, int tag ) {
	int bin;
	float recv_buffer[max_bins];

	MPI_Recv( recv_buffer, num_bins, MPI_FLOAT, proc, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

	for ( bin = 0; bin < num_bins; bin++ ) {
		output_buffer[bin] += recv_buffer[bin];
	}
}

double log_interpolate( double *binned_var, int bin, double rlout, double rri, double rll ) {
	double aM1, aM2, ah, bh;

	if ( binned_var[bin-1] > 0.0 ) {
		aM1 = log10( binned_var[bin-1] );
	} else {
		aM1 = -15.0;
	}
	if ( binned_var[bin] > 0.0 ) {
		aM2 = log10( binned_var[bin] );
	} else {
		aM2 = -15.0;
	}
	ah = ( aM2 - aM1 ) * rri;
	bh = aM1 - ah * rll;
	return pow( 10.0, ah*rlout + bh );
}

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

halo_list *load_halo_finder_catalog( const char *filename, int nmem_min, float mvir_min, float vmax_min, float rvir_min)
{
  int i, size;
  FILE *input;
  char line[256];
  int id;
  int pid;
  float px, py, pz, vx, vy, vz;
  float rvir, rhalo, mvir;
  float vmax, rmax, rs;
  int np, coords[nDim];
  halo_list *halos;
  halo *tmp;
  float a, OmM, OmL, OmB, h100, r0; 

  halos = cart_alloc(halo_list, 1 );

  input = fopen( filename, "r" );
  if ( input == NULL )
    {
      cart_error("Unable to open %s", filename );
    }

  /*
  //  Skip the job name
  */
  if(fgets(line,1024,input) == NULL) cart_error("Error in reading halo catalog %s", filename );
      
  /*
  //  Check the scale factor
  */
  if(fscanf(input," A=%f A0=%*f Ampl=%*f Step=%*f\n",&a) != 1) cart_error("Error in reading halo catalog %s", filename );
  if(fabs(a-auni[min_level]) > 1.0e-3)
    {
      cart_debug("Scalar factor in HLIST file (%f) is different from the current value (%f)",a,auni[min_level]);
    }

  if(fgets(line,1024,input) == NULL) cart_error("Error in reading halo catalog %s", filename );

  if(fscanf(input," Nrow=%*d Ngrid=%*d  Omega_0=%f OmLam_0=%f  Omegab_0=%f Hubble=%f\n",&OmM,&OmL,&OmB,&h100) != 4) cart_error("Error in reading halo catalog %s", filename );
  if(fabs(OmM/cosmology->OmegaM-1.0) > 1.0e-2)
    {
      cart_debug("OmegaM in HLIST file (%f) is different from the current value (%f)",OmM,cosmology->OmegaM);
    }
  if(fabs(OmL/cosmology->OmegaL-1.0) > 1.0e-2)
    {
      cart_debug("OmegaL in HLIST file (%f) is different from the current value (%f)",OmL,cosmology->OmegaL);
    }
  if(fabs(OmB/cosmology->OmegaB-1.0) > 1.0e-2)
    {
      cart_debug("OmegaB in HLIST file (%f) is different from the current value (%f)",OmB,cosmology->OmegaB);
    }
  if(fabs(h100/cosmology->h-1.0) > 1.0e-2)
    {
      cart_debug("Hubble in HLIST file (%f) is different from the current value (%f)",h100,cosmology->h);
    }

  /*
  //  Skip the rest
  */
  for(i=4; i<17; i++)
    {
      if(fgets(line,1024,input) == NULL) cart_error("Error in reading halo catalog %s", filename );
    }

  r0 = units->length_in_chimps;

  /* count halos */
  halos->num_halos = 0;
  size = 100;
  halos->list = cart_alloc(halo, size );

  while (fscanf( input, "%u %e %e %e %e %e %e %e %e %e %u %e %e %e %u",
		 &id, &px, &py, &pz, &vx, &vy, &vz, &rvir, &rhalo,
		 &mvir, &np, &vmax, &rmax, &rs, &pid ) == 15 )
    {
      if(np>=nmem_min && mvir>=mvir_min && vmax>=vmax_min && rvir>=rvir_min)
	{
	  if(halos->num_halos == size)
	    {
	      size *= 2;
	      tmp = cart_alloc(halo, size );
	      memcpy(tmp,halos->list,halos->num_halos*sizeof(halo));
	      cart_free(halos->list);
	      halos->list = tmp;
	    }

	  halos->list[halos->num_halos].id = id;
	
	  /* convert position to code units */
	  halos->list[halos->num_halos].pos[0] = px/r0;
	  halos->list[halos->num_halos].pos[1] = py/r0;
	  halos->list[halos->num_halos].pos[2] = pz/r0;

	  halos->list[halos->num_halos].vel[0] = vx/units->velocity;
	  halos->list[halos->num_halos].vel[1] = vy/units->velocity;
	  halos->list[halos->num_halos].vel[2] = vz/units->velocity;

	  halos->list[halos->num_halos].rvir = rvir/( 1e3*r0 );
	  halos->list[halos->num_halos].rhalo = rhalo/( 1e3*r0 );
	  halos->list[halos->num_halos].mvir = constants->Msun/cosmology->h*mvir/units->mass;
	  halos->list[halos->num_halos].vmax = constants->kms*vmax/units->velocity;
	  halos->list[halos->num_halos].np = np;

	  coords[0] = (int)halos->list[halos->num_halos].pos[0];
	  coords[1] = (int)halos->list[halos->num_halos].pos[1];
	  coords[2] = (int)halos->list[halos->num_halos].pos[2];

	  if ( num_procs > 1 )
	    {
	      halos->list[halos->num_halos].proc = processor_owner( sfc_index( coords ) );
	    }
	  else
	    {
	      halos->list[halos->num_halos].proc = local_proc_id;
	    }

	  halos->list[halos->num_halos].particles = NULL;
	  halos->list[halos->num_halos].binding_order = NULL;

	  halos->num_halos++;
	}
    }

  fclose(input);

  cart_debug("Done loading halo list; read %d halos.",halos->num_halos);

  halos->map = -1;

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
		cart_error("Error: number of halos in %s don't match provided halo_list", filename );
	}

	for ( i = 0; i < nh; i++ ) {
		fread( &size, sizeof(int), 1, input );
		cart_debug("size = %d", size );
		fread( &ih, sizeof(int), 1, input );
		fread( &nhp, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&ih, sizeof(int) );
			reorder( (char *)&nhp, sizeof(int) );
		}

		ih--;
		cart_assert( ih >= 0 && ih < halos->num_halos );

		cart_debug("ih = %d, nhp = %d", ih, nhp ); fflush(stdout);
		halos->list[ih].particles = cart_alloc(int, nhp );
		halos->list[ih].binding_order = cart_alloc(int, nhp );

		binding_energy = cart_alloc(float, nhp );

		fread( halos->list[ih].particles, sizeof(int), nhp, input );
		fread( binding_energy, sizeof(float), nhp, input );

		if ( endian ) {
			for ( j = 0; j < nhp; j++ ) {
				reorder( (char *)&halos->list[ih].particles[j], sizeof(int) );
				reorder( (char *)&binding_energy[j], sizeof(float) );
			}
		}

		for ( j = 0; j < nhp; j++ ) {
			/* convert to CART indexes */
			halos->list[ih].particles[j] -= 1;
			halos->list[ih].binding_order[j] = j;
		}

		cart_debug("about to sort"); fflush(stdout);

		/* sort particles by binding energy */
		qsort( halos->list[ih].binding_order, nhp, sizeof(int), compare_binding_energy );

		cart_free( binding_energy );

		for ( j = 0; j < nhp; j++ ) {
			halos->list[ih].binding_order[j] = halos->list[ih].particles[halos->list[ih].binding_order[j]];
		}

		/* sort by particle index for faster searching */
		qsort( halos->list[ih].particles, nhp, sizeof(int), compare_ints );

		halos->list[ih].np = nhp;

		fread( &size, sizeof(int), 1, input );
		cart_debug("size = %d", size ); fflush(stdout);
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

	cart_free( halos );
}

#ifdef STARFORM_OLD

void crude_stellar_mass_fractions( halo_list *halos ) {
	int i, j, k;
	int ipart, ihalo;
	double r, halo_pos[nDim];
	double stellar_mass[max_halos], total_stellar_mass[max_halos];
	double dm_mass[max_halos], total_dm_mass[max_halos];
	double overdensity;
	FILE *output;
	char filename[256];
	int index, icell;
	int coords[nDim];
	float rhalomax;

	if ( local_proc_id == MASTER_NODE ) {
		sprintf(filename, "%s/stellar_fractions_a%5.3f.dat", output_directory, auni[min_level] );
		output = fopen(filename, "w");
	}

	for ( ihalo = 0; ihalo < num_halos; ihalo++ ) {
		halo_pos[0] = halos->list[ihalo].pos[0];
		halo_pos[1] = halos->list[ihalo].pos[1];
		halo_pos[2] = halos->list[ihalo].pos[2];

		stellar_mass[ihalo] = dm_mass[ihalo] = 0.0;

		rhalomax = hrvir[ihalo];

		for ( i = -(int)(2*rhalomax+1.5); i <= (int)(2*rhalomax+1.5); i++ ) {
			coords[0] = ((int)(halo_pos[0]) + i ) % num_grid;
			for ( j = -(int)(2*rhalomax+1.5); j <= (int)(2*rhalomax+1.5); j++ ) {
				coords[1] = ((int)(halo_pos[1]) + j ) % num_grid;
				for ( k = -(int)(2*rhalomax+1.5); k <= (int)(2*rhalomax+1.5); k++ ) {
					coords[2] = ((int)(halo_pos[2]) + k ) % num_grid;
					icell = root_cell_location( sfc_index( coords ) );
					if ( index != NULL_OCT ) {
						ipart = cell_particle_list[icell];
						while ( ipart != NULL_PARTICLE ) {
							r = compute_distance_periodic( halo_pos, particle_x[ipart] );
							if ( r < rhalomax ) {
								if ( particle_is_star(ipart) ) {
									stellar_mass[ihalo] += particle_mass[ipart];
								} else {
									dm_mass[ihalo] += particle_mass(ipart);
								}
							}

							ipart = particle_list_next[ipart];
						}
					}
				}	
			}
		}

		if ( local_proc_id == MASTER_NODE ) {
			cart_debug("halo %u", halos->list[ihalo].id);
		} 
	}


	MPI_Reduce( stellar_mass, total_stellar_mass, num_halos, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( dm_mass, total_dm_mass, num_halos, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

	if ( local_proc_id == MASTER_NODE ) {
		for ( ihalo = 0; ihalo < num_halos; ihalo++ ) {
			overdensity = total_dm_mass[ihalo] / ( 4.*M_PI/3. * hrvir[ihalo]*hrvir[ihalo]*hrvir[ihalo] );
			total_stellar_mass[ihalo] *= aM0*hubble;
			total_dm_mass[ihalo] *= aM0*hubble;
			fprintf(output,"%u %e %e %e %e %e\n", halos->list[ihalo].id, 
				total_stellar_mass[ihalo], total_dm_mass[ihalo], aM0*hubble*halos->list[ihalo].mvir, 
				overdensity, hrvir[ihalo]*r0 );
		}

		fclose(output);
	}
}

#endif /* STARFORM */

void compute_halo_properties( char *analysis_directory, int halo_section, halo_list *halos ) {
	int i, j, k;
	int coords[nDim];
	int ihalo, icell, ipart;
	int irmax;
	int ix, iy, iz;
	int index;
	int level;
	int leaf;
	int found;
	int bin, proc;
	int throw_point;
	double r, rrl, rll, rri, rout, rlout;
	double halo_pos[nDim];
	double point_pos[nDim];
	double cell_pos[nDim];
	int max_bin, min_bin;
	int num_bins;
	int processor_mask[MAX_PROCS];
	double rlmin, rlmax, drl, dlout;
	double rlvirmin, rlvirmax, drlvir;
	float cell_vx, cell_vy, cell_vz;
	double rr[max_bins], rl[max_bins], bin_volume[max_bins], bin_volume_cumulative[max_bins];
	double rvir_bin_volume[max_bins];
	double rmid[max_bins];
	double rrvir[num_vir_bins], rlvir[num_vir_bins], rmidvir[num_vir_bins];
	double inverse_mass, vmean, cell_mass;
	double Tcell, Tcell_kev, Tfact, vfact, rfact;
	double Pcell, Pfact, Scell, Sfact;
	double dtfact, tcool, dEfact, dEcell;
	double Ycell, szfact;
	double rhogi, rhogl;
	double Zdum, Zldum;
	double tnewstar;
#ifdef COOLING
	cooling_t coolrate;
#endif
	int irvir, irflag;
	double rvir, rdout, rvdout, rmass;
	double dbi1, dbi2, dlbi1, dlbi2;
	double ah, bh, aM1, aM2;
	double aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_dark, aM_baryons, aM_total;
	double Ysz, Ysz_total, Tx;
	double total_dark_mass;
	double total_star_mass;
	double total_new_star_mass;
	double total_gas_mass;
	double total_cold_gas_mass;
	double avg_gas_metallicity_II;
	double avg_gas_metallicity_Ia;
	double avg_star_metallicity_II;
	double avg_star_metallicity_Ia;
	double avg_new_star_metallicity_II;
	double avg_new_star_metallicity_Ia;
	double avg_star_age;
	double avg_gas_mass;
	double avg_new_star_mass;
	double avg_star_mass;
	double xz, omega, a3, E2;
	double virial_overdensity;
	double overdensity;
	double vmax, rmax;
	double rcirc, vcirc_gas, vcirc_dm, vcirc_stars, vcirc_total;
	double age_star;
	double mass_fraction;
	double volume_fraction;
	double fact_nH;
	float mass;
	double Fline, Fcont, Tline, Tcont, x, f_line;
	double Tcont1, Tcont2, avgE;
	double xray_Tcont1_cell, xray_Tcont2_cell;
	double xray_cT, xray_lambda, xray_fT, xray_w;
	double xray_Tcont1, xray_Tcont2;
	double xray_Fcont_cell, xray_Fline_cell, xray_avgE_cell;

	double radii_overdensity[num_radii];
	int iflag_r[num_radii];
	double delta_r[num_radii];
	double aM_gas_r[num_radii];
	double aM_cold_gas_r[num_radii];
	double aM_stars_r[num_radii];
	double aM_new_stars_r[num_radii];
	double aM_dark_r[num_radii];
	double aM_baryons_r[num_radii];
	double aM_total_r[num_radii];
	double Ysz_r[num_radii];
	double Tx_r[num_radii];
	double avg_gas_metallicity_II_r[num_radii];
	double avg_gas_metallicity_Ia_r[num_radii];
	double avg_star_metallicity_II_r[num_radii];
	double avg_star_metallicity_Ia_r[num_radii];
	double avg_new_star_metallicity_II_r[num_radii];
	double avg_new_star_metallicity_Ia_r[num_radii];
	double avg_star_age_r[num_radii];
	double avg_gas_mass_r[num_radii];
	double avg_new_star_mass_r[num_radii];
	double avg_star_mass_r[num_radii];

	FILE *rlist[num_radii];

	char filename[256];
	char prefix[128];
	char suffix[128];
	FILE *blist;
	FILE *bszlist;
	FILE *btxlist;
	FILE *bmpro;
	FILE *bvpro;
	FILE *bzpro;
	FILE *bgpro;

	FILE *bmprovir;
	FILE *bvprovir;
	FILE *bzprovir;
	FILE *bgprovir;

	float bin_volume_fraction[max_bins];

	int bin_dark_num[max_bins];
	float bin_dark_mass[max_bins];
	float bin_dark_momentum[nDim][max_bins];
	float bin_dark_vrms[max_bins];

	int bin_star_num[max_bins];
	float bin_star_mass[max_bins];
	float bin_star_momentum[nDim][max_bins];
	float bin_star_vrms[max_bins];
	float bin_star_age[max_bins];
	float bin_star_metallicity_II[max_bins];
	float bin_star_metallicity_Ia[max_bins];

	int bin_new_star_num[max_bins];
	float bin_new_star_mass[max_bins];
	float bin_new_star_momentum[nDim][max_bins];
	float bin_new_star_vrms[max_bins];
	float bin_new_star_metallicity_II[max_bins];
	float bin_new_star_metallicity_Ia[max_bins];

	float bin_gas_mass[max_bins];
	float bin_gas_velocity[nDim][max_bins];
	float bin_gas_vrms[max_bins];
	float bin_cold_gas_mass[max_bins];
	float bin_gas_temperature[max_bins];
	float bin_gas_entropy[max_bins];
	float bin_gas_pressure[max_bins];
	float bin_gas_coolingrate[max_bins];
	float bin_gas_tcool[max_bins];
	float bin_gas_metallicity_II[max_bins];
	float bin_gas_metallicity_Ia[max_bins];
	float bin_sz_flux[max_bins];
	float bin_xray_Fcont[max_bins];
	float bin_xray_Fline[max_bins];
	float bin_xray_avgE[max_bins];
	float bin_xray_Tcont1[max_bins];
	float bin_xray_Tcont2[max_bins];

	double bin_total_dark_mass[max_bins];
	double bin_total_star_mass[max_bins];
	double bin_total_new_star_mass[max_bins];
	double bin_total_gas_mass[max_bins];
	double bin_total_cold_gas_mass[max_bins];
	double bin_total_density[max_bins];
	double bin_total_sz_flux[max_bins];

	double bin_total_xray_Fcont[max_bins];
	double bin_total_xray_Fline[max_bins];
	double bin_total_xray_avgE[max_bins];
	double bin_total_xray_Tcont1[max_bins];
	double bin_total_xray_Tcont2[max_bins];

	/* set up conversion constants */
#ifdef CHECK_LEGACY_UNITS
	double r0 = legacy_units->r0;
	double aM0 = legacy_units->M0;
#ifdef HYDRO
	Tfact = legacy_units->T0 * ( constants->gamma - 1.0 ) / ( abox[min_level]*abox[min_level] );
	Pfact = legacy_units->P0/(abox[min_level]*abox[min_level]*abox[min_level]*abox[min_level]*abox[min_level]);
	Sfact = legacy_units->S0;
	szfact = 3.9207e-15 * ( constants->gamma - 1.0 ) * legacy_units->T0 * legacy_units->r0 * cosmology->OmegaM * cosmology->h / 
			(abox[min_level]*abox[min_level]*abox[min_level]*abox[min_level]);
#endif
#else  /* CHECK_LEGACY_UNITS */
#error "This function uses legacy units that have not been converted yet. Define CHECK_LEGACY_UNITS."
#endif /* CHECK_LEGACY_UNITS */

	rfact = 1.0e3 * units->length * constants->pc; /* proper kpc -> code units */
	vfact = units->velocity;

#ifdef STARFORM
	tnewstar = constants->Gyr * new_star_age / units->time ;
#endif

#ifdef COOLING
	dtfact = 1e6 * units->time / constants->yr;
	dEfact = units->time / constants->yr / Pfact;

	fact_nH = log10( 1.12e-5 * cosmology->Omh2 * ( 1.0 - constants->Yp ) / 
			(abox[min_level]*abox[min_level]*abox[min_level]) );
#endif

	/* set up binning */
	irmax = (int)( rbinmax / ( 1000.0 * r0 ) ) + 1;
	rlmin = log10(rbinmin/rfact);
	rlmax = log10(rbinmax/rfact);
	drl = (rlmax - rlmin)/(float)(max_bins-1);
	num_bins = (rlmax - rlmin)/drl + 1;

	cart_debug("rlmin = %e", rlmin );
	cart_debug("rlmax = %e", rlmax );
	cart_debug("drl = %e", drl );
	cart_debug("num_bins = %u", num_bins );
	cart_debug("irmax = %u", irmax );

	cart_assert( num_bins <= max_bins );

	rl[0] = 0.0;
	rr[0] = rbinmin/rfact;
	rmid[0] = 0.5*rbinmin/rfact;
	bin_volume[0] = 4.0*M_PI/3.0 * rr[0]*rr[0]*rr[0];
	bin_volume_cumulative[0] = bin_volume[0];

	for ( i = 1; i < num_bins; i++ ) {
		rl[i] = rr[i-1];
		rr[i] = pow( 10.0, rlmin + (float)i*drl );
		rmid[i] = pow( 10.0, rlmin + (float)(i-0.5)*drl );
		bin_volume[i] = 4.0*M_PI/3.0 * ( rr[i]*rr[i]*rr[i] - rl[i]*rl[i]*rl[i] );
		bin_volume_cumulative[i] = 4.0*M_PI/3.0*rr[i]*rr[i]*rr[i];
        }

	/* set up binning in units of rvir (virial_radius_index radius) */
	rlvirmin = log10(rbinvirmin);
	rlvirmax = log10(rbinvirmax);
	drlvir = (rlvirmax-rlvirmin)/(float)(num_vir_bins-1);

	rlvir[0] = 0.0;
	rmidvir[0] = 0.5*rbinmin;
	rrvir[0] = rbinvirmin;
	rvir_bin_volume[0] = 4.0*M_PI/3.0 * rrvir[0]*rrvir[0]*rrvir[0];

	for ( i = 1; i < num_vir_bins; i++ ) {
		rlvir[i] = rrvir[i-1];
		rrvir[i] = pow( 10.0, rlvirmin + (float)i*drlvir );
		rmidvir[i] = pow( 10.0, rlvirmin + (float)(i-0.5)*drlvir );

		rvir_bin_volume[i] = 4.0*M_PI/3.0 * ( rrvir[i]*rrvir[i]*rrvir[i] - 
					rlvir[i]*rlvir[i]*rlvir[i] );
	}

	/* compute virial overdensity using Bryan & Norman (1998) */
	a3 = abox[min_level] * abox[min_level] * abox[min_level];
	E2 = cosmology->OmegaM / a3 + cosmology->OmegaL;
	omega = cosmology->OmegaM / a3 / E2;
	xz = omega - 1.0;
	virial_overdensity = ( 18.0*M_PI*M_PI + 82.0*xz - 39.0*xz*xz ) / ( 1 + xz );
	dlout = log10(virial_overdensity);

	cart_debug("E2 = %e", E2 );
	cart_debug("omega = %e", omega );
	cart_debug("xz = %e", xz );
	cart_debug("virial_overdensity = %e", virial_overdensity );

	if ( num_procs == 1 ) {
		sprintf( suffix, "_a%6.4f.dat", abox[min_level] );
	} else {
		sprintf( suffix, "_a%6.4f_%04u.dat", abox[min_level], local_proc_id );
	}

	if ( halo_section == -1 ) {
		sprintf( prefix, "%s/h_", analysis_directory );
	} else {
		sprintf( prefix, "%s/h_%06u_", analysis_directory, halo_section );
	}

	/* convert overdensities and open files */
	for ( i = 0; i < num_radii; i++ ) {
		if ( radii_units[i] == 'c' ) {
			radii_overdensity[i] = radii_delta[i] / cosmology->OmegaM;
		} else {
			radii_overdensity[i] = radii_delta[i];
		}

		sprintf( filename, "%sblist_%s%s", prefix, radii_label[i], suffix );
		rlist[i] = fopen( filename, "w" );
		if ( rlist[i] == NULL ) {
			cart_error( "Unable to open %s for writing", filename );
		}
	}

	/* open output files */
	sprintf( filename, "%sblist%s", prefix, suffix );
	blist = fopen( filename, "w" );
	if ( blist == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

	sprintf( filename, "%sbmpro%s", prefix, suffix );
	bmpro = fopen( filename, "w" );
	if ( bmpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

#ifdef HYDRO
	sprintf( filename, "%sszlist%s", prefix, suffix );
	bszlist = fopen( filename, "w" );
	if ( bszlist == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

#ifdef ANALYZE_XRAY
	sprintf( filename, "%stxlist%s", prefix, suffix );
	btxlist = fopen( filename, "w" );
	if ( btxlist == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* ANALYZE_XRAY */

	sprintf( filename, "%sbgpro%s", prefix, suffix );
	bgpro = fopen( filename, "w" );
	if ( bgpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* HYDRO */

#ifdef STARFORM
	sprintf( filename, "%sbzpro%s", prefix, suffix );
	bzpro = fopen( filename, "w" );
	if ( bzpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* STARFORM */

	sprintf( filename, "%shvpro%s", prefix, suffix );
	bvpro = fopen( filename, "w" );
	if ( bvpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

	/* profiles binned in units of virial radius */
	sprintf( filename, "%sbmpro_%s%s", prefix, radii_label[virial_radius_index], suffix );
	bmprovir = fopen( filename, "w" );
	if ( bmprovir == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

#ifdef HYDRO
	sprintf( filename, "%sbgpro_%s%s", prefix, radii_label[virial_radius_index], suffix );
	bgprovir = fopen( filename, "w" );
	if ( bgprovir == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* HYDRO */

#ifdef STARFORM
	sprintf( filename, "%sbzpro_%s%s", prefix, radii_label[virial_radius_index], suffix );
	bzprovir = fopen( filename, "w" );
	if ( bzprovir == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* STARFORM */

	sprintf( filename, "%sbvpro_%s%s", prefix, radii_label[virial_radius_index], suffix );
	bvprovir = fopen( filename, "w" );
	if ( bvprovir == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

	/* put headers on files */
	if ( local_proc_id == MASTER_NODE ) {
		fprintf( blist, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( blist, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( blist, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( blist, "# virial overdensity = %.3f\n", virial_overdensity );

		fprintf( blist, "# Columns:\n");
		fprintf( blist, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( blist, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( blist, "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( blist, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bszlist, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( bszlist, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bszlist, "# Overdensities: vir (%.3f)", virial_overdensity );

		for ( i = 0; i < num_radii; i++ ) {
			fprintf( bszlist, ", %s", radii_label[i] );
		}

		fprintf( bszlist, "\n");

		fprintf( bszlist, "# Columns:\n" );
		fprintf( bszlist, "# id Mg Mdm Mtotal Ysz (for each overdensity)\n" );


#ifdef ANALYZE_XRAY
		fprintf( btxlist, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( btxlist, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( btxlist, "# Overdensities: vir (%.3f)", virial_overdensity );

		for ( i = 0; i < num_radii; i++ ) {
			fprintf( btxlist, ", %s", radii_label[i] );
		}

		fprintf( btxlist, "\n");

		fprintf( btxlist, "# Columns:\n" );
		fprintf( btxlist, "# id Mg Mdm Mtotal Tx (for each overdensity)\n" );
#endif

		/* density profiles */
		fprintf( bmpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( bmpro, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bmpro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( bmpro, "# virial overdensity = %.3f\n", virial_overdensity );

                fprintf( bmpro, "# Header Columns:\n");
                fprintf( bmpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
                fprintf( bmpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bmpro, "# vmax [km/s] rmax [/h kpc]\n" );
                fprintf( bmpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bmpro, "# Profile Columns:\n");
		fprintf( bmpro, "# rmid rr [/h kpc] Mdm Mg Mcg M* M*new (< rr) [/h Msolar]\n");

		fprintf( bmprovir, "# Binning: %u bins, %f to %f\n", num_vir_bins, rbinvirmin, rbinvirmax );
		fprintf( bmprovir, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bmprovir, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( bmprovir, "# virial overdensity = %s, %.3f\n", radii_label[virial_radius_index], radii_delta[virial_radius_index] );

		fprintf( bmprovir, "# Header Columns:\n");
		fprintf( bmprovir, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( bmprovir, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bmprovir, "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( bmprovir, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bmprovir, "# Profile Columns:\n");
		fprintf( bmprovir, "# rmid rr [/h kpc] Mdm Mg Mcg M* M*new (< rr) [/h Msolar]\n");

		/* velocity profiles */
		fprintf( bvpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( bvpro, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bvpro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( bvpro, "# virial overdensity = %.3f\n", virial_overdensity );

		fprintf( bvpro, "# Header Columns:\n");
		fprintf( bvpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( bvpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bvpro, "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( bvpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bvpro, "# Profile Columns:\n");
		fprintf( bvpro, "# rmid rr [/h kpc] Vrms_dm Vrms_gas Vrms_* Vrms_*new [km/s]\n");
		fprintf( bvpro, "# Vc_dm Vc_gas Vc_* Vc_total\n");

		fprintf( bvprovir, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( bvprovir, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bvprovir, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( bvprovir, "# virial overdensity = %s, %.3f\n", radii_label[virial_radius_index], radii_delta[virial_radius_index] );

		fprintf( bvprovir, "# Header Columns:\n");
		fprintf( bvprovir, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( bvprovir, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bvprovir, "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( bvprovir, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bvprovir, "# Profile Columns:\n");
		fprintf( bvprovir, "# rmid rr [/h kpc] Vrms_dm Vrms_gas Vrms_* Vrms_*new [km/s]\n");
		fprintf( bvprovir, "# Vc_dm Vc_gas Vc_* Vc_total\n");

#ifdef HYDRO
		fprintf( bgpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( bgpro, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bgpro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( bgpro, "# virial overdensity = %.3f\n", virial_overdensity );

		fprintf( bgpro, "# Header Columns:\n");
		fprintf( bgpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( bgpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bgpro, "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( bgpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bgpro, "# Profile Columns:\n");
		fprintf( bgpro, "# rmid rr [/h kpc] Mg Mcg [/h Msolar] T [K] P [ergs cm^{-3}] \n");
		fprintf( bgpro, "# S [keV cm^2] dE_cool [ergs s^-1] <tcool> [Myr] <Zg_II> <Zg_Ia>\n");

		fprintf( bgprovir, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( bgprovir, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bgprovir, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( bgprovir, "# virial overdensity = %.3f\n", virial_overdensity );

		fprintf( bgprovir, "# Header Columns:\n");
		fprintf( bgprovir, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( bgprovir, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bgprovir, "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( bgprovir, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bgprovir, "# Profile Columns:\n");
		fprintf( bgprovir, "# rmid rr [/h kpc] Mg Mcg [/h Msolar] T [K] P [ergs cm^{-3}] \n");
		fprintf( bgprovir, "# S [keV cm^2] <tcool> [Myr] dE_cool [ergs s^-1] <Zg_II> <Zg_Ia>\n");
#endif /* HYDRO */

		/* star and metallicity profiles */
#ifdef STARFORM
                fprintf( bzpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
                fprintf( bzpro, "# Monte carlo points per cell = %u\n", points_per_cell );
                fprintf( bzpro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
                fprintf( bzpro, "# virial overdensity = %.3f\n", virial_overdensity );

                fprintf( bzpro, "# Header Columns:\n");
                fprintf( bzpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
                fprintf( bzpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bzpro, "# vmax [km/s] rmax [/h kpc]\n" );
                fprintf( bzpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

                fprintf( bzpro, "# Profile Columns:\n");
                fprintf( bzpro, "# rmid rr [/h kpc] Zg_II Zg_Ia Z*_II Z*_Ia Z*new_II Z*new_Ia <t*> [Gyr]\n");

		fprintf( bzprovir, "# Binning: %u bins, %f to %f\n", num_vir_bins, rbinvirmin, rbinvirmax );
		fprintf( bzprovir, "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( bzprovir, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( bzprovir, "# virial overdensity = %s, %.3f\n", radii_label[virial_radius_index], radii_delta[virial_radius_index] );

		fprintf( bzprovir, "# Header Columns:\n");
		fprintf( bzprovir, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( bzprovir, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
		fprintf( bzprovir, "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( bzprovir, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fprintf( bzprovir, "# Profile Columns:\n");
		fprintf( bzprovir, "# rmid rr [/h kpc] Zg_II Zg_Ia Z*_II Z*_Ia Z*new_II Z*new_Ia <t*> [Gyr]\n");
#endif /* STARFORM */

		for ( i = 0; i < num_radii; i++ ) {
			fprintf( rlist[i], "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
			fprintf( rlist[i], "# Monte carlo points per cell = %u\n", points_per_cell );
			fprintf( rlist[i], "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
			fprintf( rlist[i], "# overdensity = %s, %.3f\n", radii_label[i], radii_delta[i] );

			if ( radii_units == "r" ) {
				fprintf( rlist[i], "# Quantities at fixed radius rout\n" );
			} else if ( radii_units == "c" ) {
				fprintf( rlist[i], "# Quantities at fixed overdensity %f with respect to critical\n", radii_delta[i] );
			} else if ( radii_units == "m" ) {
				fprintf( rlist[i], "# Quantities at fixed matter overdensity %f\n", radii_delta[i] );
			}

			fprintf( rlist[i], "# Columns:\n");
			fprintf( rlist[i], "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
			fprintf( rlist[i], "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rout] Mvir [/h Msolar]\n");
			fprintf( rlist[i], "# vmax [km/s] rmax [/h kpc]\n" );
			fprintf( rlist[i], "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

			fflush(rlist[i]);
		}
	}

	/* process all halos */
	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		if ( halo_section == -1 || halos->list[ihalo].section == halo_section ) {
			ix = (int)(halos->list[ihalo].pos[0]);
			iy = (int)(halos->list[ihalo].pos[1]);
			iz = (int)(halos->list[ihalo].pos[2]);

			for ( i = 0; i < nDim; i++ ) {
				halo_pos[i] = halos->list[ihalo].pos[i];
			}

			num_halo_leafs = 0;

			/* rmass = halos->list[ihalo].rhalo; */
			rmass = rr[num_bins-1];

			if ( halos->list[ihalo].proc == local_proc_id ) {
				found = 1;
	
				/* determine which procs to receive from */
				for ( proc = 0; proc < num_procs; proc++ ) {
					processor_mask[proc] = 0;
				}

				for ( i = -irmax; i <= irmax; i++ ) {
					coords[0] = ( ix + i + num_grid ) % num_grid;

					for ( j = -irmax; j <= irmax; j++ ) {
						coords[1] = ( iy + j + num_grid ) % num_grid;

						for ( k = -irmax; k <= irmax; k++ ) {
							coords[2] = ( iz + k + num_grid ) % num_grid;

							index = sfc_index( coords );

							if ( num_procs == 1 || root_cell_is_local(index) ) {
								num_halo_leafs += cell_num_child_leaves( root_cell_location(index), min_level );
							} else {
								proc = processor_owner( index );
								cart_assert( proc >= 0 && proc < num_procs );
								processor_mask[proc] = 1;
							}
						}
					}
				}
			} else {
				found = 0;
				for ( i = -irmax; i <= irmax; i++ ) {
					coords[0] = ( ix + i + num_grid ) % num_grid;

					for ( j = -irmax; j <= irmax; j++ ) {
						coords[1] = ( iy + j + num_grid ) % num_grid;

						for ( k = -irmax; k <= irmax; k++ ) {
							coords[2] = ( iz + k + num_grid ) % num_grid;

							index = sfc_index( coords );
							if ( root_cell_is_local(index) ) {
								found = 1;
								num_halo_leafs += cell_num_child_leaves( root_cell_location(index), min_level );
							}
						}
					}
				}
			}

			if ( found ) {
				cart_debug("found halo %u, num_halo_leafs = %d", halos->list[ihalo].id, num_halo_leafs );

				/* clear out bins */
				for ( bin = 0; bin < num_bins; bin++ ) {
					bin_dark_num[bin] = 0;
					bin_dark_mass[bin] = 0.0;
					bin_dark_momentum[0][bin] = 0.0;
					bin_dark_momentum[1][bin] = 0.0;
					bin_dark_momentum[2][bin] = 0.0;
					bin_dark_vrms[bin] = 0.0;
					bin_star_num[bin] = 0;
					bin_star_mass[bin] = 0.0;
					bin_star_momentum[0][bin] = 0.0;
					bin_star_momentum[1][bin] = 0.0;
					bin_star_momentum[2][bin] = 0.0;
					bin_star_vrms[bin] = 0.0;
					bin_star_age[bin] = 0.0;
					bin_star_metallicity_II[bin] = 0.0;
					bin_star_metallicity_Ia[bin] = 0.0;
					bin_new_star_num[bin] = 0;
					bin_new_star_mass[bin] = 0.0;
					bin_new_star_momentum[0][bin] = 0.0;
					bin_new_star_momentum[1][bin] = 0.0;
					bin_new_star_momentum[2][bin] = 0.0;
					bin_new_star_vrms[bin] = 0.0;
					bin_new_star_metallicity_II[bin] = 0.0;
					bin_new_star_metallicity_Ia[bin] = 0.0;
					bin_gas_mass[bin] = 0.0;
					bin_gas_velocity[0][bin] = 0.0;
					bin_gas_velocity[1][bin] = 0.0;
					bin_gas_velocity[2][bin] = 0.0;
					bin_gas_vrms[bin] = 0.0;
					bin_cold_gas_mass[bin] = 0.0;
					bin_gas_temperature[bin] = 0.0;
					bin_gas_entropy[bin] = 0.0;
					bin_gas_pressure[bin] = 0.0;
					bin_gas_coolingrate[bin] = 0.0;
					bin_gas_tcool[bin] = 0.0;
					bin_gas_metallicity_II[bin] = 0.0;
					bin_gas_metallicity_Ia[bin] = 0.0;
					bin_sz_flux[bin] = 0.0;
#ifdef ANALYZE_XRAY
					bin_xray_Fcont[bin] = 0.0;
					bin_xray_Fline[bin] = 0.0;
					bin_xray_avgE[bin] = 0.0;
					bin_xray_Tcont1[bin] = 0.0;
					bin_xray_Tcont2[bin] = 0.0;
#endif /* ANALYZE_XRAY */
				}

#ifdef PARTICLES
				/* assign particles to bins */
				for ( ipart = 0; ipart < num_particles; ipart++ ) {
					if ( particle_level[ipart] != FREE_PARTICLE_LEVEL ) {
						r = compute_distance_periodic( halo_pos, particle_x[ipart] );

						if ( r < rr[0] ) {
							bin = 0;
						} else {
							bin = (int)((log10(r) - rlmin)/drl) + 1;
						}

						if ( bin < num_bins ) {
#ifdef STARFORM
							if ( particle_is_star(ipart) ) {
								bin_star_num[bin]++;
								bin_star_mass[bin] += particle_mass[ipart];

								bin_star_momentum[0][bin] += particle_v[ipart][0]*particle_mass[ipart];
								bin_star_momentum[1][bin] += particle_v[ipart][1]*particle_mass[ipart];
								bin_star_momentum[2][bin] += particle_v[ipart][2]*particle_mass[ipart];

								bin_star_vrms[bin] += particle_mass[ipart] * (
										particle_v[ipart][0]*particle_v[ipart][0] +
										particle_v[ipart][1]*particle_v[ipart][1] +
										particle_v[ipart][2]*particle_v[ipart][2] );

								bin_star_age[bin] += star_tbirth[ipart]*particle_mass[ipart];
								bin_star_metallicity_II[bin] += star_metallicity_II[ipart]*particle_mass[ipart];
								bin_star_metallicity_Ia[bin] += star_metallicity_Ia[ipart]*particle_mass[ipart];

								if ( particle_t[ipart] - star_tbirth[ipart] <= tnewstar ) {
									bin_new_star_num[bin]++;
									bin_new_star_mass[bin] += particle_mass[ipart];

									bin_new_star_momentum[0][bin] += particle_v[ipart][0]*particle_mass[ipart];
									bin_new_star_momentum[1][bin] += particle_v[ipart][1]*particle_mass[ipart];
									bin_new_star_momentum[2][bin] += particle_v[ipart][2]*particle_mass[ipart];

									bin_new_star_vrms[bin] += particle_mass[ipart] * (
											particle_v[ipart][0]*particle_v[ipart][0] +
											particle_v[ipart][1]*particle_v[ipart][1] +
											particle_v[ipart][2]*particle_v[ipart][2] );
									bin_new_star_metallicity_II[bin] += star_metallicity_II[ipart]*particle_mass[ipart];
									bin_new_star_metallicity_Ia[bin] += star_metallicity_Ia[ipart]*particle_mass[ipart];
								}
							} else {
								mass = particle_mass[ipart];
								bin_dark_num[bin]++;
								bin_dark_mass[bin] += mass;

								bin_dark_momentum[0][bin] += particle_v[ipart][0]*mass;
								bin_dark_momentum[1][bin] += particle_v[ipart][1]*mass;
								bin_dark_momentum[2][bin] += particle_v[ipart][2]*mass;

								bin_dark_vrms[bin] += mass * (
										particle_v[ipart][0]*particle_v[ipart][0] +
										particle_v[ipart][1]*particle_v[ipart][1] +
										particle_v[ipart][2]*particle_v[ipart][2] );

							}
#else
							bin_dark_num[bin]++;
							mass = particle_mass[ipart];
							bin_dark_mass[bin] += mass;
							bin_dark_momentum[0][bin] += particle_v[ipart][0]*mass;
							bin_dark_momentum[1][bin] += particle_v[ipart][1]*mass;
							bin_dark_momentum[2][bin] += particle_v[ipart][2]*mass;
							bin_dark_vrms[bin] += mass * (
									particle_v[ipart][0]*particle_v[ipart][0] +
									particle_v[ipart][1]*particle_v[ipart][1] +
									particle_v[ipart][2]*particle_v[ipart][2] );
#endif /* STARFORM */
						}
					}
				}

				cart_debug("done assigning particles");
#endif /* PARTICLES */

#ifdef HYDRO
				/* build list of leaf cells within halo */
				leaf_index = cart_alloc(int, num_halo_leafs );

				num_halo_leafs = 0;

				for ( i = -irmax; i <= irmax; i++ ) {
					coords[0] = ( ix + i + num_grid ) % num_grid;

					for ( j = -irmax; j <= irmax; j++ ) {
						coords[1] = ( iy + j + num_grid ) % num_grid;

						for ( k = -irmax; k <= irmax; k++ ) {
							coords[2] = ( iz + k + num_grid ) % num_grid;

							index = sfc_index( coords );

							if ( num_procs == 1 || root_cell_is_local(index) ) {
								tree_traversal( root_cell_location(index), build_leaf_index );
							} 
						}
					}
				}

				for ( leaf = 0; leaf < num_halo_leafs; leaf++ ) {
					icell = leaf_index[leaf];
					cart_assert( icell >= 0 && icell < num_cells );
					cart_assert( cell_is_leaf(icell) );

					level = cell_level(icell);

					rhogi = 1.0 / cell_gas_density(icell);

					/* grab cell properties */
					Tcell = Tfact * cell_gas_internal_energy(icell) * rhogi;
					Tcell_kev = Tcell * 8.6170868e-8;

					Pcell = Pfact * cell_gas_pressure(icell);
					Scell = Sfact * cell_gas_internal_energy(icell)*pow(rhogi,constants->gamma);
					Ycell = szfact * cell_gas_internal_energy(icell)*cell_volume[level];

#if defined(COOLING) && !defined(RADIATIVE_TRANSFER)
					/* take code density -> log10(n_H [cm^-3]) */
					rhogl = log10(cell_gas_density(icell)) + fact_nH;
#ifdef ENRICH
					Zdum = max(1.0e-10,cell_gas_metal_density(icell)/(constants->Zsun*cell_gas_density(icell)));
					Zldum = log10(Zdum);
#else
					Zdum = 0.0;
					Zldum = 0.0;
#endif /* ENRICH */
#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
					dEcell = cooling_rate( rhogl, Tcell*1e-4, Zldum ) * 
						cell_gas_density(icell)*cell_gas_density(icell) *
						abox[level];
#else
					coolrate = cooling_rate( rhogl, Tcell*1e-4, Zldum);
					dEcell = (coolrate.Cooling-coolrate.Heating) * 
						cell_gas_density(icell)*cell_gas_density(icell) *
						abox[level];	
#endif
					tcool = cell_gas_internal_energy(icell) / dEcell / dtfact;
					dEcell /= dEfact;
#endif /* COOLING && !RADIATIVE_TRANSFER */

#ifdef ANALYZE_XRAY
					/* Alexey's weighting (Vikhlinin 2006) */
					xray_calibration( Tcell_kev, &xray_cT, &xray_lambda, &xray_fT );

					/* eq 6 */
					xray_w = xray_cT * cell_gas_density(icell)*cell_gas_density(icell)*pow(Tcell_kev,-alpha_xray );

					/* numerator & denominator of eq 4, computes <T>_cont */
					xray_Tcont1 = xray_w*Tcell_kev*cell_volume[level];
					xray_Tcont2 = xray_w*cell_volume[level];

					/* eq 9 */
					xray_Fcont_cell = xray_cT * cell_gas_density(icell)*cell_gas_density(icell)*cell_volume[level];
				
					/* eq 11 */
					xray_Fline_cell = xray_lambda * Zdum * cell_gas_density(icell)*cell_gas_density(icell)*cell_volume[level];
					xray_avgE_cell = xray_fT * xray_Fline_cell;

/*
					if ( Tcell_kev > 0.1 ) {
						cart_debug("Tcell = %e", Tcell );
						cart_debug("Tcell = %e", Tcell_kev );
						cart_debug("xray_cT = %e", xray_cT );
						cart_debug("xray_lambda = %e", xray_lambda );
						cart_debug("xray_fT = %e", xray_fT );
						cart_debug("xray_w = %e", xray_w );
						cart_debug("xray_Tcont1 = %e", xray_Tcont1 );
						cart_debug("xray_Tcont2 = %e", xray_Tcont2 );
						cart_debug("xray_Fcont_cell = %e", xray_Fcont_cell );
						cart_debug("xray_Fline_cell = %e", xray_Fline_cell );
						cart_debug("xray_avgE_cell = %e", xray_avgE_cell );
					}
*/

#endif /* ANALYZE_XRAY */

					cell_mass = cell_gas_density(icell)*cell_volume[level];
					cell_vx = cell_momentum(icell,0)*rhogi;
					cell_vy = cell_momentum(icell,1)*rhogi;
					cell_vz = cell_momentum(icell,2)*rhogi;

					/* find which bin the center of this cell lies in */
					cell_center_position( icell, cell_pos );

					r = compute_distance_periodic( halos->list[ihalo].pos, cell_pos );

					if ( r < rr[0] ) {
						bin = 0;
					} else {
						bin = (int)((log10(r) - rlmin)/drl) + 1;
					}

#if points_per_cell > 1
					max_bin = min_bin = bin;

					/* now find if corners of cell lie in different bins */
					for ( i = 0; i < num_children; i++ ) {
						for ( j = 0; j < nDim; j++ ) {
							point_pos[j] = cell_pos[j] + cell_delta[i][j]*cell_size[level];
						}

						r = compute_distance_periodic( halo_pos, point_pos );

						if ( r < rr[0] ) {
							bin = 0;
						} else {
							bin = (int)((log10(r) - rlmin)/drl) + 1;
						}

						cart_assert( bin >= num_bins || (r >= rl[bin] && r <= rr[bin]) );

						min_bin = min( min_bin, bin );
						max_bin = max( max_bin, bin );
					}

					if ( min_bin < num_bins ) {
						/* if only in one bin, just add completely */
						if ( max_bin == min_bin ) {
							bin_gas_mass[min_bin] += cell_mass;

							if ( Tcell <= Tcold ) {
								bin_cold_gas_mass[min_bin] += cell_mass;
							}

							bin_gas_pressure[min_bin] += Pcell*cell_mass;
							bin_gas_entropy[min_bin] += Scell*cell_mass;
							bin_gas_temperature[min_bin] += Tcell*cell_mass;
							bin_sz_flux[min_bin] += Ycell;


#ifdef COOLING
							bin_gas_coolingrate[min_bin] += dEcell*cell_volume[level];
								bin_gas_tcool[min_bin] += cell_mass*tcool;
#endif /* COOLING */

#ifdef ENRICH
							bin_gas_metallicity_II[min_bin] += cell_gas_metal_density_II(icell)*rhogi*cell_mass;
							bin_gas_metallicity_Ia[min_bin] += cell_gas_metal_density_Ia(icell)*rhogi*cell_mass;
#endif /* ENRICH */

							bin_gas_velocity[0][min_bin] += cell_vx*cell_mass;
							bin_gas_velocity[1][min_bin] += cell_vy*cell_mass;
							bin_gas_velocity[2][min_bin] += cell_vz*cell_mass;

							bin_gas_vrms[min_bin] += (cell_vx*cell_vx + cell_vy*cell_vy + cell_vz*cell_vz)*cell_mass;

#ifdef ANALYZE_XRAY
							/* X-ray quantities */
							bin_xray_Fcont[min_bin] += xray_Fcont_cell;
							bin_xray_Fline[min_bin] += xray_Fline_cell;
							bin_xray_avgE[min_bin] += xray_avgE_cell;
							bin_xray_Tcont1[min_bin] += xray_Tcont1;
							bin_xray_Tcont2[min_bin] += xray_Tcont2;
#endif /* ANALYZE_XRAY */
						} else {
							/* monte carlo cell properties into bins */
							for ( bin = 0; bin < num_bins; bin++ ) {
								bin_volume_fraction[bin] = 0;
							}

							for ( throw_point = 0; throw_point < points_per_cell; throw_point++ ) {
								for ( j = 0; j < nDim; j++ ) {
									point_pos[j] = (cart_rand()-0.5)*cell_size[level] +
										cell_pos[j];
								}

								r = compute_distance_periodic( halo_pos, point_pos );

								if ( r < rr[0] ) {
									bin = 0;
								} else {
									bin = (int)((log10(r) - rlmin)/drl) + 1;
								}

								if ( bin < num_bins ) {
									bin_volume_fraction[bin]++;
								}
							}

							for ( bin = 0; bin < num_bins; bin++ ) {
								if ( bin_volume_fraction[bin] > 0 ) {
									volume_fraction = (double)bin_volume_fraction[bin] /
										(double)points_per_cell;
									mass_fraction = cell_mass*volume_fraction;

									bin_gas_mass[bin] += mass_fraction;

									if ( Tcell <= Tcold ) {
										bin_cold_gas_mass[bin] += mass_fraction;
									}

									bin_gas_pressure[bin] += Pcell*mass_fraction;
									bin_gas_entropy[bin] += Scell*mass_fraction;
									bin_gas_temperature[bin] += Tcell*mass_fraction;
									bin_sz_flux[bin] += Ycell*volume_fraction;

#ifdef COOLING
									bin_gas_coolingrate[bin] += dEcell*cell_volume[level]*volume_fraction;
									bin_gas_tcool[bin] += tcool*mass_fraction;
#endif /* COOLING */

#ifdef ENRICH
									bin_gas_metallicity_II[bin] += cell_gas_metal_density_II(icell) *
										rhogi * mass_fraction; 
									bin_gas_metallicity_Ia[bin] += cell_gas_metal_density_Ia(icell) *
										rhogi * mass_fraction;

#endif /* ENRICH */

									bin_gas_velocity[0][bin] += cell_vx*mass_fraction;
									bin_gas_velocity[1][bin] += cell_vy*mass_fraction;
									bin_gas_velocity[2][bin] += cell_vz*mass_fraction;

#ifdef ANALYZE_XRAY
									/* X-ray quantities */
									bin_xray_Fcont[bin] += xray_Fcont_cell*volume_fraction;
									bin_xray_Fline[bin] += xray_Fline_cell*volume_fraction;
									bin_xray_avgE[bin] += xray_avgE_cell*volume_fraction;
									bin_xray_Tcont1[bin] += xray_Tcont1*volume_fraction;
									bin_xray_Tcont2[bin] += xray_Tcont2*volume_fraction;
#endif /* ANALYZE_XRAY */
								}
							}
						}
					}
#else
					if ( bin < num_bins ) {
						bin_gas_mass[bin] += cell_mass;

						if ( Tcell <= Tcold ) {
							bin_cold_gas_mass[bin] += cell_mass;
						}

						bin_gas_pressure[bin] += Pcell*cell_mass;
						bin_gas_entropy[bin] += Scell*cell_mass;
						bin_gas_temperature[bin] += Tcell*cell_mass;
						bin_sz_flux[bin] += Ycell;

#ifdef COOLING
						bin_gas_coolingrate[bin] += dEcell*cell_volume[level]*volume_fraction;
						bin_gas_tcool[bin] += tcool*mass_fraction;
#endif /* COOLING */
#ifdef ENRICH
						bin_gas_metallicity_II[bin] += cell_gas_metal_density_II(icell)*rhogi*cell_mass;
						bin_gas_metallicity_Ia[bin] += cell_gas_metal_density_Ia(icell)*rhogi*cell_mass;
#endif /* ENRICH */
						bin_gas_velocity[0][bin] += cell_vx*cell_mass;
						bin_gas_velocity[1][bin] += cell_vy*cell_mass;
						bin_gas_velocity[2][bin] += cell_vz*cell_mass;

						bin_gas_vrms[bin] += (cell_vx*cell_vx + cell_vy*cell_vy + cell_vz*cell_vz)*cell_mass;

#ifdef ANALYZE_XRAY
						/* X-ray quantities */
						bin_xray_Fcont[bin] += xray_Fcont_cell;
						bin_xray_Fline[bin] += xray_Fline_cell;
						bin_xray_avgE[bin] += xray_avgE_cell;
						bin_xray_Tcont1[bin] += xray_Tcont1;
						bin_xray_Tcont2[bin] += xray_Tcont2;
#endif /* ANALYZE_XRAY */
					}
#endif /* points_per_cell > 1 */
				}
#endif /* HYDRO */

				cart_debug("done assigning gas to grid");

				if ( halos->list[ihalo].proc == local_proc_id ) {
					/* receive bin information from other processors */
					for ( proc = 0; proc < num_procs; proc++ ) {
						if ( processor_mask[proc] ) {
							recv_int_bins( bin_dark_num, num_bins, proc, ihalo );
							recv_float_bins( bin_dark_mass, num_bins, proc, ihalo );
							recv_float_bins( bin_dark_momentum[0], num_bins, proc, ihalo );
							recv_float_bins( bin_dark_momentum[1], num_bins, proc, ihalo );
							recv_float_bins( bin_dark_momentum[2], num_bins, proc, ihalo );
							recv_float_bins( bin_dark_vrms, num_bins, proc, ihalo );
							recv_int_bins( bin_star_num, num_bins, proc, ihalo );
							recv_float_bins( bin_star_mass, num_bins, proc, ihalo );
							recv_float_bins( bin_star_momentum[0], num_bins, proc, ihalo );
							recv_float_bins( bin_star_momentum[1], num_bins, proc, ihalo );
							recv_float_bins( bin_star_momentum[2], num_bins, proc, ihalo );
							recv_float_bins( bin_star_vrms, num_bins, proc, ihalo );
							recv_float_bins( bin_star_age, num_bins, proc, ihalo );
							recv_float_bins( bin_star_metallicity_II, num_bins, proc, ihalo );
							recv_float_bins( bin_star_metallicity_Ia, num_bins, proc, ihalo );
							recv_int_bins( bin_new_star_num, num_bins, proc, ihalo );
							recv_float_bins( bin_new_star_mass, num_bins, proc, ihalo );
							recv_float_bins( bin_new_star_momentum[0], num_bins, proc, ihalo );
							recv_float_bins( bin_new_star_momentum[1], num_bins, proc, ihalo );
							recv_float_bins( bin_new_star_momentum[2], num_bins, proc, ihalo );
							recv_float_bins( bin_new_star_vrms, num_bins, proc, ihalo );
							recv_float_bins( bin_new_star_metallicity_II, num_bins, proc, ihalo );
							recv_float_bins( bin_new_star_metallicity_Ia, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_mass, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_velocity[0], num_bins, proc, ihalo );
							recv_float_bins( bin_gas_velocity[1], num_bins, proc, ihalo );
							recv_float_bins( bin_gas_velocity[2], num_bins, proc, ihalo );
							recv_float_bins( bin_gas_vrms, num_bins, proc, ihalo );
							recv_float_bins( bin_cold_gas_mass, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_pressure, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_entropy, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_temperature, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_tcool, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_coolingrate, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_metallicity_II, num_bins, proc, ihalo );
							recv_float_bins( bin_gas_metallicity_Ia, num_bins, proc, ihalo );
							recv_float_bins( bin_sz_flux, num_bins, proc, ihalo );
#ifdef ANALYZE_XRAY
							recv_float_bins( bin_xray_Fcont, num_bins, proc, ihalo );
							recv_float_bins( bin_xray_Fline, num_bins, proc, ihalo );
							recv_float_bins( bin_xray_avgE, num_bins, proc, ihalo );
							recv_float_bins( bin_xray_Tcont1, num_bins, proc, ihalo );
							recv_float_bins( bin_xray_Tcont2, num_bins, proc, ihalo );
#endif /* ANALYZE_XRAY */
						}
					}

					cart_debug("done receiving bin data from other processors, processing...");

					/* process bins */
					for ( bin = 0; bin < num_bins; bin++ ) {
						if ( bin_gas_mass[bin] > 0.0 ) {
							inverse_mass = 1.0 / bin_gas_mass[bin];

							bin_gas_pressure[bin] *= inverse_mass;
							bin_gas_temperature[bin] *= inverse_mass;
							bin_gas_entropy[bin] *= inverse_mass;
							bin_gas_tcool[bin] *= inverse_mass;
							bin_gas_coolingrate[bin] /= bin_volume[bin];
							bin_gas_metallicity_II[bin] *= inverse_mass;
							bin_gas_metallicity_Ia[bin] *= inverse_mass;
							bin_gas_velocity[0][bin] *= inverse_mass;
							bin_gas_velocity[1][bin] *= inverse_mass;
							bin_gas_velocity[2][bin] *= inverse_mass;
							vmean = bin_gas_velocity[0][bin]*bin_gas_velocity[0][bin] +
								bin_gas_velocity[1][bin]*bin_gas_velocity[1][bin] +
								bin_gas_velocity[2][bin]*bin_gas_velocity[2][bin];
							bin_gas_vrms[bin] = sqrt(fabs(bin_gas_vrms[bin]*inverse_mass - vmean))*vfact;
							bin_gas_velocity[0][bin] *= vfact;
							bin_gas_velocity[1][bin] *= vfact;
							bin_gas_velocity[2][bin] *= vfact;
						}	
					}

					/* compute cumulative quantities */
					total_dark_mass = 0.0;
					total_star_mass = 0.0;
					total_new_star_mass = 0.0;
					total_gas_mass = 0.0;
					total_cold_gas_mass = 0.0;

					Ysz_total = 0.0;

					avg_gas_metallicity_II = 0.0;
					avg_gas_metallicity_Ia = 0.0;
					avg_star_metallicity_II = 0.0;
					avg_star_metallicity_Ia = 0.0;
					avg_new_star_metallicity_II = 0.0;
					avg_new_star_metallicity_Ia = 0.0;
					avg_star_age = 0.0;

					avg_new_star_mass = 0.0;
					avg_star_mass = 0.0;
					avg_gas_mass = 0.0;

					irvir = 0;
					irflag = 0;

					rvir = 0.0;
					rout = 0.0;

					for ( i = 0; i < num_radii; i++ ) {
						avg_gas_metallicity_II_r[i] = 0.0;
						avg_gas_metallicity_Ia_r[i] = 0.0;
						avg_star_metallicity_II_r[i] = 0.0;
						avg_star_metallicity_Ia_r[i] = 0.0;
						avg_new_star_metallicity_II_r[i] = 0.0;
						avg_new_star_metallicity_Ia_r[i] = 0.0;
						avg_star_age_r[i] = 0.0;

						avg_new_star_mass_r[i] = 0.0;
						avg_star_mass_r[i] = 0.0;
						avg_gas_mass_r[i] = 0.0;

						iflag_r[i] = 0;
						delta_r[i] = 0.0;
					}

					vmax = 0.0;

					for ( bin = 0; bin < num_bins; bin++ ) {
						cart_assert( bin_dark_mass[bin] >= 0.0 );
						cart_assert( bin_star_mass[bin] >= 0.0 );
						cart_assert( bin_new_star_mass[bin] >= 0.0 );

						total_dark_mass += bin_dark_mass[bin];
						total_star_mass += bin_star_mass[bin];
						total_new_star_mass += bin_new_star_mass[bin];
						total_gas_mass += bin_gas_mass[bin];
						total_cold_gas_mass += bin_cold_gas_mass[bin];

						Ysz_total += bin_sz_flux[bin];

						bin_total_dark_mass[bin] = total_dark_mass;
						bin_total_star_mass[bin] = total_star_mass;
						bin_total_new_star_mass[bin] = total_new_star_mass;
						bin_total_gas_mass[bin] = total_gas_mass;
						bin_total_cold_gas_mass[bin] = total_cold_gas_mass;

						bin_total_sz_flux[bin] = Ysz_total;

#ifdef ANALYZE_XRAY
						bin_total_xray_Fcont[bin] += bin_xray_Fcont[bin];
						bin_total_xray_Fline[bin] += bin_xray_Fline[bin];
						bin_total_xray_avgE[bin] += bin_xray_avgE[bin];
						bin_total_xray_Tcont1[bin] += bin_xray_Tcont1[bin];
						bin_total_xray_Tcont2[bin] += bin_xray_Tcont2[bin];
#endif /* ANALYZE_XRAY */

						if ( irflag == 0 ) {
							avg_gas_metallicity_II += bin_gas_metallicity_II[bin] * bin_gas_mass[bin];
							avg_gas_metallicity_Ia += bin_gas_metallicity_Ia[bin] * bin_gas_mass[bin];
							avg_star_metallicity_II += bin_star_metallicity_II[bin];
							avg_star_metallicity_Ia = bin_star_metallicity_Ia[bin];
							avg_new_star_metallicity_II += bin_new_star_metallicity_II[bin];
							avg_new_star_metallicity_Ia += bin_new_star_metallicity_Ia[bin];
							avg_star_age += bin_star_age[bin];

							avg_gas_mass += bin_gas_mass[bin];
							avg_new_star_mass += bin_new_star_mass[bin];
							avg_star_mass += bin_star_mass[bin];
						}

						/* find maximum circular velocity within min(rvir,rhalo) */
						if ( bin > 0 && irflag == 0 && rr[bin] < halos->list[ihalo].rhalo) {
							if ( sqrt( (	bin_total_dark_mass[bin] +
											bin_total_gas_mass[bin] +
											bin_total_star_mass[bin] ) / rr[bin] ) >= vmax ) {
								rmax = rr[bin];
								vmax = sqrt( ( 	bin_total_dark_mass[bin] +
											bin_total_gas_mass[bin] +
											bin_total_star_mass[bin] ) / rr[bin] );
							}
						}

						if ( bin == 0 ) {
							dbi1 = 0.0;
						} else {
							dbi1 = ( bin_total_dark_mass[bin-1] + bin_total_star_mass[bin-1] + 
									bin_total_gas_mass[bin-1] ) / bin_volume_cumulative[bin-1];
						}

						dbi2 = (bin_total_dark_mass[bin] + bin_total_star_mass[bin] + 
								bin_total_gas_mass[bin]) / bin_volume_cumulative[bin];

						/* compute virial radius */
						if ( bin > 0 && ( dbi1 >= virial_overdensity && dbi2 < virial_overdensity ) ) {
							rrl = log10(rr[bin]);
							rll = log10(rl[bin]);
							rri = 1.0/(rrl-rll);
							dlbi1 = log10(dbi1);
							dlbi2 = log10(dbi2);
							irvir = 1;
							rvir = pow( 10.0, (dlout*(rrl-rll) + rll*dlbi2 - 
										rrl*dlbi1)/(dlbi2-dlbi1));
						}

						if ( bin > 0 && irflag == 0 &&
								( ( rl[bin] <= rmass && rr[bin] > rmass ) || 
								  ( dbi1 >= virial_overdensity && dbi2 < virial_overdensity ) ) ) {

							irflag = 1;
							rrl = log10(rr[bin]);
							rll = log10(rl[bin]);
							rri = 1.0/(rrl-rll);

							if ( dbi1 >= virial_overdensity && dbi2 < virial_overdensity ) {
								dlbi1 = log10(dbi1);
								dlbi2 = log10(dbi2);
								rout = pow( 10.0, (dlout*(rrl-rll) + rll*dlbi2 
											- rrl*dlbi1)/(dlbi2-dlbi1));
							} else {
								rout = rmass;
							}

							rlout = log10(rout);

							/* interpolate total quantities to rout */
							aM_gas = log_interpolate( bin_total_gas_mass, bin, rlout, rri, rll );
							aM_cold_gas = log_interpolate( bin_total_cold_gas_mass, bin, rlout, rri, rll );

#ifdef STARFORM
							aM_stars = log_interpolate( bin_total_star_mass, bin, rlout, rri, rll );
							aM_new_stars = log_interpolate( bin_total_new_star_mass, bin, rlout, rri, rll );
#else
							aM_stars = 0.0;
							aM_new_stars = 0.0;
#endif /* STARFORM */

							aM_dark = log_interpolate( bin_total_dark_mass, bin, rlout, rri, rll );
						
							Ysz = log_interpolate( bin_total_sz_flux, bin, rlout, rri, rll );

#ifdef ANALYZE_XRAY
							Fcont = log_interpolate( bin_total_xray_Fcont, bin, rlout, rri, rll );
							Fline = log_interpolate( bin_total_xray_Fline, bin, rlout, rri, rll );
							avgE = log_interpolate( bin_total_xray_avgE, bin, rlout, rri, rll );
							Tcont1 = log_interpolate( bin_total_xray_Tcont1, bin, rlout, rri, rll );
							Tcont2 = log_interpolate( bin_total_xray_Tcont2, bin, rlout, rri, rll );

							Tcont = Tcont1/Tcont2;

							avgE /= Fline;

							Tline = xray_calibrated_line_temperature(avgE);

							f_line = Fline / ( Fline + Fcont );
							x = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
								exp( -pow(f_line/delta_xray_2,8.0) ) *
								exp( -pow(f_line/delta_xray_2,8.0) );

							Tx = x*Tcont + (1.0-x)*Tline;
#endif /* ANALYZE_XRAY */

							/* convert to output units */
							aM_gas *= aM0 * cosmology->h;
							aM_cold_gas *= aM0 * cosmology->h;
							aM_stars *= aM0 * cosmology->h;
							aM_new_stars *= aM0 * cosmology->h;
							aM_dark *= aM0 * cosmology->h;

							Ysz *= rfact*rfact / 1e6;

						}

						for ( i = 0; i < num_radii; i++ ) {
							if ( iflag_r[i] == 0 ) {
								avg_gas_metallicity_II_r[i] += bin_gas_metallicity_II[bin] * bin_gas_mass[bin];
								avg_gas_metallicity_Ia_r[i] += bin_gas_metallicity_Ia[bin] * bin_gas_mass[bin];
								avg_star_metallicity_II_r[i] += bin_star_metallicity_II[bin];
								avg_star_metallicity_Ia_r[i] = bin_star_metallicity_Ia[bin];
								avg_new_star_metallicity_II_r[i] += bin_new_star_metallicity_II[bin];
								avg_new_star_metallicity_Ia_r[i] += bin_new_star_metallicity_Ia[bin];

								avg_star_age_r[i] += bin_star_age[bin];
								avg_gas_mass_r[i] += bin_gas_mass[bin];
								avg_new_star_mass_r[i] += bin_new_star_mass[bin];
								avg_star_mass_r[i] += bin_star_mass[bin];
							}

							if ( bin > 0 && iflag_r[i] == 0 &&
									( ( i > 0 && dbi1 >= radii_overdensity[i] && dbi2 < radii_overdensity[i] ) ||
									  ( i == 0 && rl[bin] < halos->list[ihalo].rhalo && rr[bin] >= halos->list[ihalo].rhalo ) ) ) {

								iflag_r[i] = 1;
								rrl = log10(rr[bin]);
								rll = log10(rl[bin]);
								rri = 1.0/(rrl-rll);

								if ( i == 0 ) {
									delta_r[i] = halos->list[ihalo].rhalo;
								} else {
									dlbi1 = log10(dbi1);
									dlbi2 = log10(dbi2);
									delta_r[i] = pow( 10.0, (log10(radii_overdensity[i])*(rrl-rll) + rll*dlbi2
												- rrl*dlbi1)/(dlbi2-dlbi1));
								}

								rlout = log10(delta_r[i]);

								/* interpolate total quantities to rlout */
								aM_gas_r[i] = log_interpolate( bin_total_gas_mass, bin, rlout, rri, rll );
								aM_cold_gas_r[i] = log_interpolate( bin_total_cold_gas_mass, bin, rlout, rri, rll );

#ifdef STARFORM
								aM_stars_r[i] = log_interpolate( bin_total_star_mass, bin, rlout, rri, rll );
								aM_new_stars_r[i] = log_interpolate( bin_total_new_star_mass, bin, rlout, rri, rll );
#else
								aM_stars_r[i] = 0.0;
								aM_new_stars_r[i] = 0.0;
#endif /* STARFORM */

								aM_dark_r[i] = log_interpolate( bin_total_dark_mass, bin, rlout, rri, rll );

								Ysz_r[i] = log_interpolate( bin_total_sz_flux, bin, rlout, rri, rll ); 

#ifdef ANALYZE_XRAY
								Fcont = log_interpolate( bin_total_xray_Fcont, bin, rlout, rri, rll );
								Fline = log_interpolate( bin_total_xray_Fline, bin, rlout, rri, rll );
								avgE = log_interpolate( bin_total_xray_avgE, bin, rlout, rri, rll );
								Tcont1 = log_interpolate( bin_total_xray_Tcont1, bin, rlout, rri, rll );
								Tcont2 = log_interpolate( bin_total_xray_Tcont2, bin, rlout, rri, rll );

								Tcont = Tcont1/Tcont2;
								avgE /= Fline;
								Tline = xray_calibrated_line_temperature(avgE);

								f_line = Fline / ( Fline + Fcont );
								x = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
									exp( -pow(f_line/delta_xray_2,8.0) ) *
									exp( -pow(f_line/delta_xray_2,8.0) );

								Tx_r[i] = x*Tcont + (1.0-x)*Tline;
#endif /* ANALYZE_XRAY */

								aM_gas_r[i] *= aM0 * cosmology->h;
								aM_cold_gas_r[i] *= aM0 * cosmology->h;
								aM_stars_r[i] *= aM0 * cosmology->h;
								aM_new_stars_r[i] *= aM0 * cosmology->h;
								aM_dark_r[i] *= aM0 * cosmology->h;

								Ysz *= rfact*rfact/1e6;
							}
						}
					}

					if ( irvir == 0 ) {
						/* never reached virial radius */
						rout = rbinmax/rfact;
						rvir = rbinmax/rfact; 
						aM_gas = total_gas_mass * aM0 * cosmology->h;
						aM_cold_gas = total_cold_gas_mass * aM0 * cosmology->h;
						aM_stars = total_star_mass * aM0 * cosmology->h;
						aM_new_stars = total_new_star_mass * aM0 * cosmology->h;
						aM_dark = total_dark_mass * aM0 * cosmology->h;

						Ysz = Ysz_total * rfact*rfact/1e6;

#ifdef ANALYZE_XRAY
						Fcont = bin_total_xray_Fcont[max_bins-1];
						Fline = bin_total_xray_Fline[max_bins-1];
						avgE = bin_total_xray_avgE[max_bins-1];
						Tcont1 = bin_total_xray_Tcont1[max_bins-1];
						Tcont2 = bin_total_xray_Tcont2[max_bins-1];

						Tcont = Tcont1/Tcont2;
						avgE /= Fline;
						Tline = xray_calibrated_line_temperature(avgE);

						f_line = Fline / ( Fline + Fcont );
						x = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
							exp( -pow(f_line/delta_xray_2,8.0) ) *
							exp( -pow(f_line/delta_xray_2,8.0) );

						Tx_r[i] = x*Tcont + (1.0-x)*Tline;
#endif /* ANALYZE_XRAY */
					}

					/* compute average quantities */
					if ( avg_gas_mass > 0.0 ) {
						avg_gas_metallicity_II /= avg_gas_mass*constants->Zsun; 
						avg_gas_metallicity_Ia /= avg_gas_mass*constants->Zsun;
					}

					if ( avg_star_mass > 0.0 ) {
						avg_star_metallicity_II /= avg_star_mass*constants->Zsun;
						avg_star_metallicity_Ia /= avg_star_mass*constants->Zsun;
						avg_star_age = tphys_from_tcode( avg_star_age / avg_star_mass );
					}

					if ( avg_new_star_mass > 0.0 ) {
						avg_new_star_metallicity_II /= avg_new_star_mass*constants->Zsun;
						avg_new_star_metallicity_Ia /= avg_new_star_mass*constants->Zsun;
					}

					aM_baryons = aM_stars + aM_gas;
					aM_total = aM_baryons + aM_dark;

					rdout = rout * rfact * cosmology->h;
					rvdout = rvir * rfact * cosmology->h;

					rmax *= rfact * cosmology->h;
					vmax *= 2.07498e-3 * sqrt(aM0/rfact);

					/* output */
					fprintf( blist, "%u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n", 
							halos->list[ihalo].id, rdout, rvdout, aM_gas, 
							aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
							aM_dark, aM_total, aM0*cosmology->h*halos->list[ihalo].mvir, vmax, rmax, 
							avg_gas_metallicity_II, avg_gas_metallicity_Ia,
							avg_star_metallicity_II, avg_star_metallicity_Ia, 
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );
					fflush(blist);

					for ( i = 0; i < num_radii; i++ ) {
						/* compute average quantities */
						if ( avg_gas_mass_r[i] > 0.0 ) {
							avg_gas_metallicity_II_r[i] /= avg_gas_mass_r[i]*constants->Zsun;
							avg_gas_metallicity_Ia_r[i] /= avg_gas_mass_r[i]*constants->Zsun;
						}

						if ( avg_star_mass_r[i] > 0.0 ) {
							avg_star_metallicity_II_r[i] /= avg_star_mass_r[i]*constants->Zsun;
							avg_star_metallicity_Ia_r[i] /= avg_star_mass_r[i]*constants->Zsun;
							avg_star_age_r[i] = tphys_from_tcode( avg_star_age_r[i] / avg_star_mass_r[i] );
						}

						if ( avg_new_star_mass_r[i] > 0.0 ) {
							avg_new_star_metallicity_II_r[i] /= avg_new_star_mass_r[i]*constants->Zsun;
							avg_new_star_metallicity_Ia_r[i] /= avg_new_star_mass_r[i]*constants->Zsun;
						}

						aM_baryons_r[i] = aM_stars_r[i] + aM_gas_r[i];
						aM_total_r[i] = aM_baryons_r[i] + aM_dark_r[i];

						fprintf( rlist[i], "%u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
								halos->list[ihalo].id, 
								delta_r[i]*rfact*cosmology->h, rvdout,
								aM_gas_r[i], aM_cold_gas_r[i], aM_stars_r[i], aM_new_stars_r[i], 
								aM_baryons_r[i], aM_dark_r[i], aM_total_r[i], aM_total, 
								vmax, rmax, avg_gas_metallicity_II_r[i], 
								avg_gas_metallicity_Ia_r[i], avg_star_metallicity_II_r[i], 
								avg_star_metallicity_Ia_r[i], avg_new_star_metallicity_II_r[i], 
								avg_new_star_metallicity_Ia_r[i], avg_star_age_r[i] );

						fflush(rlist[i]);
					}

#ifdef HYDRO
					fprintf( bszlist, "%u %e %e %e %e", halos->list[ihalo].id, aM_gas, aM_dark, aM_total, Ysz );
					for ( i = 0; i < num_radii; i++ ) {
						fprintf( bszlist, " %e %e %e %e", aM_gas_r[i], aM_dark_r[i], aM_total_r[i], Ysz_r[i] );
					}
					fprintf( bszlist, "\n" );
#endif /* HYDRO */

#ifdef ANALYZE_XRAY
					fprintf( btxlist, "%u %e %e %e %e", halos->list[ihalo].id, aM_gas, aM_dark, aM_total, Tx );
					for ( i = 0; i < num_radii; i++ ) {
						fprintf( btxlist, " %e %e %e %e", aM_gas_r[i], aM_dark_r[i], aM_total_r[i], Tx_r[i] );
					}
					fprintf( btxlist, "\n" );
#endif /* ANALYZE_XRAY */


					fprintf( bmpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id, rdout, rvdout, 
							aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
							aM_dark, aM_total, aM0*cosmology->h*halos->list[ihalo].mvir, vmax, rmax,
							avg_gas_metallicity_II, avg_gas_metallicity_Ia,
							avg_star_metallicity_II, avg_star_metallicity_Ia, 
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, 
							avg_star_age );

					/* write out profiles */
					for ( bin = 0; bin < num_bins; bin++ ) {
						fprintf( bmpro, "%.3f %.3f %e %e %e %e %e\n",
								rmid[bin]*rfact*cosmology->h, rr[bin]*rfact*cosmology->h,
								bin_total_dark_mass[bin]*aM0*cosmology->h,
								bin_total_gas_mass[bin]*aM0*cosmology->h,
								bin_total_cold_gas_mass[bin]*aM0*cosmology->h,
								bin_total_star_mass[bin]*aM0*cosmology->h,
								bin_total_new_star_mass[bin]*aM0*cosmology->h );
					}

					fflush(bmpro);

					fprintf( bvpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id, 
							rdout, rvdout, aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
							aM_dark, aM_total, aM0*cosmology->h*halos->list[ihalo].mvir, vmax, rmax,
							avg_gas_metallicity_II, avg_gas_metallicity_Ia,
							avg_star_metallicity_II, avg_star_metallicity_Ia, 
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

					for ( bin = 0; bin < num_bins; bin++ ) {
						rcirc = rr[bin]*rfact; /* code units to proper kpc */
						vcirc_dm = 2.07498e-3 * sqrt( aM0*bin_total_dark_mass[bin] / rcirc );
						vcirc_gas = 2.07498e-3 * sqrt( aM0*bin_total_gas_mass[bin] / rcirc );
						vcirc_stars = 2.07498e-3 * sqrt( aM0*bin_total_star_mass[bin] / rcirc );
						vcirc_total = 2.07498e-3 * sqrt( aM0*(bin_total_dark_mass[bin]+
									bin_total_gas_mass[bin]+bin_total_star_mass[bin])/rcirc );

						bin_dark_vrms[bin] = vfact*sqrt(bin_dark_vrms[bin]);
						bin_gas_vrms[bin] = vfact*sqrt(bin_gas_vrms[bin]);
						bin_star_vrms[bin] = vfact*sqrt(bin_star_vrms[bin]);
						bin_new_star_vrms[bin] = vfact*sqrt(bin_new_star_vrms[bin]);

						fprintf( bvpro, "%.3f %.3f %e %e %e %e %.3f %.3f %.3f %.3f\n", 
								rmid[bin]*rfact*cosmology->h, rr[bin]*rfact*cosmology->h,
								bin_dark_vrms[bin], bin_gas_vrms[bin], 
								bin_star_vrms[bin], bin_new_star_vrms[bin],
								vcirc_dm, vcirc_gas, vcirc_stars, vcirc_total );
					}

					fflush(bvpro);

#ifdef HYDRO
					fprintf( bgpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id,
							rdout, rvdout, aM_gas, aM_cold_gas,
							aM_stars, aM_new_stars, aM_baryons,
							aM_dark, aM_total, aM0*cosmology->h*halos->list[ihalo].mvir,
							vmax, rmax, avg_gas_metallicity_II, avg_gas_metallicity_Ia,
							avg_star_metallicity_II, avg_star_metallicity_Ia,
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

					for ( bin = 0; bin < num_bins; bin++ ) {
						fprintf( bgpro, "%.3f %.3f %e %e %e %e %e %.3f %.3f %e %e\n",
								rmid[bin]*rfact*cosmology->h, rr[bin]*rfact*cosmology->h,
								bin_total_gas_mass[bin]*aM0*cosmology->h,
								bin_total_cold_gas_mass[bin]*aM0*cosmology->h,
								bin_gas_temperature[bin],
								bin_gas_pressure[bin],
								bin_gas_entropy[bin],
								bin_gas_coolingrate[bin],
								bin_gas_tcool[bin],
								bin_gas_metallicity_II[bin]/constants->Zsun,
								bin_gas_metallicity_Ia[bin]/constants->Zsun
						       );
					}

					fflush(bgpro);
#endif /* HYDRO */

#ifdef STARFORM
					fprintf( bzpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id, rdout, rvdout, 
							aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
							aM_dark, aM_total, aM0*cosmology->h*halos->list[ihalo].mvir, vmax, rmax,
							avg_gas_metallicity_II, avg_gas_metallicity_Ia,
							avg_star_metallicity_II, avg_star_metallicity_Ia, 
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, 
							avg_star_age );

					for ( bin = 0; bin < num_bins; bin++ ) {
						if ( bin_star_age[bin] > 0 ) {
							age_star = tphys_from_tcode( bin_star_age[bin] );
						} else {
							age_star = -1.0;
						}

						fprintf( bzpro, "%.3f %.3f %e %e %e %e %e %e %.3f\n",
								rmid[bin]*rfact*cosmology->h, rr[bin]*rfact*cosmology->h,
								bin_gas_metallicity_II[bin]/constants->Zsun, bin_gas_metallicity_Ia[bin]/constants->Zsun,
								bin_star_metallicity_II[bin]/constants->Zsun, bin_star_metallicity_Ia[bin]/constants->Zsun,
								bin_new_star_metallicity_II[bin]/constants->Zsun, 
								bin_new_star_metallicity_Ia[bin]/constants->Zsun,
								age_star );
					}

					fflush(bzpro);
#endif /* STARFORM */
				} else {
					/* send bin values to halo owner */
					MPI_Send( bin_dark_num, num_bins, MPI_INT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_mass, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_momentum[0], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_momentum[1], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_momentum[2], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_vrms, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_num, num_bins, MPI_INT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_mass, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_momentum[0], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_momentum[1], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_momentum[2], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_vrms, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_age, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_metallicity_II, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_metallicity_Ia, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_num, num_bins, MPI_INT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_mass, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_momentum[0], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_momentum[1], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_momentum[2], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_vrms, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_metallicity_II, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_metallicity_Ia, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_mass, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_velocity[0], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_velocity[1], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_velocity[2], num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_vrms, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_cold_gas_mass, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_pressure, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_entropy, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_temperature, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_tcool, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_coolingrate, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_metallicity_II, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_metallicity_Ia, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_sz_flux, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
#ifdef ANALYZE_XRAY
					MPI_Send( bin_xray_Fcont, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_xray_Fline, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_xray_avgE, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_xray_Tcont1, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_xray_Tcont2, num_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
#endif /* ANALYZE_XRAY */
				}

				cart_debug("done processing bins or sending info to parent");

				/* set up binning in units of the virial radius (any virial_radius_index),
				 *  or the tidal radius if it's smaller */
				if ( halos->list[ihalo].proc == local_proc_id ) {
					for ( proc = 0; proc < num_procs; proc++ ) {
						if ( processor_mask[proc] ) {
							MPI_Send( delta_r, num_radii, MPI_FLOAT, proc, ihalo,
								MPI_COMM_WORLD );
						}
					}
				} else {
					cart_assert( num_procs > 1 );

					MPI_Recv( delta_r, num_radii, MPI_FLOAT, halos->list[ihalo].proc, 
							ihalo, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}

				if ( delta_r[virial_radius_index] < 1e-6 ) {
					rvir = delta_r[0];
				} else {
					rvir = min( delta_r[0], delta_r[virial_radius_index] );
				}
				rmax /= (rvir*cosmology->h*rfact);

				/* clear out bins again to bin in units of virial radius (virial_radius_index) */
				for ( bin = 0; bin < num_bins; bin++ ) {
					bin_dark_num[bin] = 0;
					bin_dark_mass[bin] = 0.0;
					bin_dark_momentum[0][bin] = 0.0;
					bin_dark_momentum[1][bin] = 0.0;
					bin_dark_momentum[2][bin] = 0.0;
					bin_dark_vrms[bin] = 0.0;
					bin_star_num[bin] = 0;
					bin_star_mass[bin] = 0.0;
					bin_star_momentum[0][bin] = 0.0;
					bin_star_momentum[1][bin] = 0.0;
					bin_star_momentum[2][bin] = 0.0;
					bin_star_vrms[bin] = 0.0;
					bin_star_age[bin] = 0.0;
					bin_star_metallicity_II[bin] = 0.0;
					bin_star_metallicity_Ia[bin] = 0.0;
					bin_new_star_num[bin] = 0;
					bin_new_star_mass[bin] = 0.0;
					bin_new_star_momentum[0][bin] = 0.0;
					bin_new_star_momentum[1][bin] = 0.0;
					bin_new_star_momentum[2][bin] = 0.0;
					bin_new_star_vrms[bin] = 0.0;
					bin_new_star_metallicity_II[bin] = 0.0;
					bin_new_star_metallicity_Ia[bin] = 0.0;
					bin_gas_mass[bin] = 0.0;
					bin_gas_velocity[0][bin] = 0.0;
					bin_gas_velocity[1][bin] = 0.0;
					bin_gas_velocity[2][bin] = 0.0;
					bin_gas_vrms[bin] = 0.0;
					bin_cold_gas_mass[bin] = 0.0;
					bin_gas_temperature[bin] = 0.0;
					bin_gas_entropy[bin] = 0.0;
					bin_gas_pressure[bin] = 0.0;
					bin_gas_tcool[bin] = 0.0;
					bin_gas_coolingrate[bin] = 0.0;
					bin_gas_metallicity_II[bin] = 0.0;
					bin_gas_metallicity_Ia[bin] = 0.0;
				}

#ifdef PARTICLES
				cart_debug("assigning particles to virial bins");

				/* assign particles to bins */
				for ( ipart = 0; ipart < num_particles; ipart++ ) {
					if ( particle_level[ipart] != FREE_PARTICLE_LEVEL ) {
						r = compute_distance_periodic( halo_pos, particle_x[ipart] ) / rvir;

						if ( r < rrvir[0] ) {
							bin = 0;
						} else {
							bin = (int)((log10(r) - rlvirmin)/drlvir) + 1;
						}

						cart_assert( bin >= 0 );
						if ( bin < num_vir_bins ) {
#ifdef STARFORM
							if ( particle_is_star(ipart) ) {
								bin_star_num[bin]++;
								bin_star_mass[bin] += particle_mass[ipart];

								bin_star_momentum[0][bin] += particle_v[ipart][0]*particle_mass[ipart];
								bin_star_momentum[1][bin] += particle_v[ipart][1]*particle_mass[ipart];
								bin_star_momentum[2][bin] += particle_v[ipart][2]*particle_mass[ipart];

								bin_star_vrms[bin] += particle_mass[ipart] * (
										particle_v[ipart][0]*particle_v[ipart][0] +
										particle_v[ipart][1]*particle_v[ipart][1] +
										particle_v[ipart][2]*particle_v[ipart][2] );

								bin_star_age[bin] += star_tbirth[ipart]*particle_mass[ipart];
								bin_star_metallicity_II[bin] += star_metallicity_II[ipart]*particle_mass[ipart];
								bin_star_metallicity_Ia[bin] += star_metallicity_Ia[ipart]*particle_mass[ipart];

								if ( particle_t[ipart] - star_tbirth[ipart] <= tnewstar ) {
									bin_new_star_num[bin]++;
									bin_new_star_mass[bin] += particle_mass[ipart];

									bin_new_star_momentum[0][bin] += particle_v[ipart][0]*particle_mass[ipart];
									bin_new_star_momentum[1][bin] += particle_v[ipart][1]*particle_mass[ipart];
									bin_new_star_momentum[2][bin] += particle_v[ipart][2]*particle_mass[ipart];

									bin_new_star_vrms[bin] += particle_mass[ipart] * (
											particle_v[ipart][0]*particle_v[ipart][0] +
											particle_v[ipart][1]*particle_v[ipart][1] +
											particle_v[ipart][2]*particle_v[ipart][2] );
									bin_new_star_metallicity_II[bin] += star_metallicity_II[ipart]*particle_mass[ipart];
									bin_new_star_metallicity_Ia[bin] += star_metallicity_Ia[ipart]*particle_mass[ipart];
								}
							} else {
								bin_dark_num[bin]++;
								mass = particle_mass[ipart];
								bin_dark_mass[bin] += mass;

								bin_dark_momentum[0][bin] += particle_v[ipart][0]*mass;
								bin_dark_momentum[1][bin] += particle_v[ipart][1]*mass;
								bin_dark_momentum[2][bin] += particle_v[ipart][2]*mass;

								bin_dark_vrms[bin] += mass * (
										particle_v[ipart][0]*particle_v[ipart][0] +
										particle_v[ipart][1]*particle_v[ipart][1] +
										particle_v[ipart][2]*particle_v[ipart][2] );

							}
#else
							bin_dark_num[bin]++;
							bin_dark_mass[bin] += mass;
							bin_dark_momentum[0][bin] += particle_v[ipart][0]*mass;
							bin_dark_momentum[1][bin] += particle_v[ipart][1]*mass;
							bin_dark_momentum[2][bin] += particle_v[ipart][2]*mass;
							bin_dark_vrms[bin] += mass * (
									particle_v[ipart][0]*particle_v[ipart][0] +
									particle_v[ipart][1]*particle_v[ipart][1] +
									particle_v[ipart][2]*particle_v[ipart][2] );
#endif /* STARFORM */
						}
					}
				}

				cart_debug("done assigning particles to virial bins"); fflush(stdout);
#endif /* PARTICLES */

#ifdef HYDRO
				for ( leaf = 0; leaf < num_halo_leafs; leaf++ ) {
					icell = leaf_index[leaf];

					cart_assert( icell >= 0 && icell < num_cells );

					level = cell_level(icell);

					cart_assert( level >= min_level && level <= max_level );
					cart_assert( cell_gas_density(icell) != 0.0 );

					rhogi = 1.0 / cell_gas_density(icell);

					/* grab cell properties */
					Tcell = Tfact * cell_gas_internal_energy(icell) * rhogi;
					Pcell = Pfact * cell_gas_pressure(icell);
					Scell = Sfact * cell_gas_internal_energy(icell)*pow(rhogi,constants->gamma);

#if defined(COOLING) && !defined(RADIATIVE_TRANSFER)
					/* take code density -> log10(n_H [cm^-3]) */
					rhogl = log10(cell_gas_density(icell)) + fact_nH;
#ifdef ENRICH
					Zdum = max(1.0e-10,cell_gas_metal_density(icell)/(constants->Zsun*cell_gas_density(icell)));
					Zldum = log10( Zdum );
#else
					Zdum = 0.0;
					Zldum = 0.0;
#endif /* ENRICH */
#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
					dEcell = cooling_rate( rhogl, Tcell*1e-4, Zldum ) *
						cell_gas_density(icell)*cell_gas_density(icell) *
						abox[level];
#else
					coolrate = cooling_rate(rhogl, Tcell*1e-4, Zldum);
					dEcell = (coolrate.Cooling-coolrate.Heating) *
						cell_gas_density(icell)*cell_gas_density(icell) *
						abox[level];
					
#endif

					tcool = cell_gas_internal_energy(icell) / dEcell / dtfact;
					dEcell /= dEfact;
#endif /* COOLING && !RADIATIVE_TRANSFER */

					cell_mass = cell_gas_density(icell)*cell_volume[level];
					cell_vx = cell_momentum(icell,0)*rhogi;
					cell_vy = cell_momentum(icell,1)*rhogi;
					cell_vz = cell_momentum(icell,2)*rhogi;

					/* find which bin the center of this cell lies in */
					cell_center_position( icell, cell_pos );

					cart_assert( rvir > 0.0 );
					r = compute_distance_periodic( halos->list[ihalo].pos, cell_pos ) / rvir;

					if ( r < rrvir[0] ) {
						bin = 0;
					} else {
						bin = (int)((log10(r) - rlvirmin)/drlvir) + 1;
					}

#if (points_per_cell > 1)
					max_bin = min_bin = bin;

					/* now find if corners of cell lie in different bins */
					for ( i = 0; i < num_children; i++ ) {
						for ( j = 0; j < nDim; j++ ) {
							point_pos[j] = cell_pos[j] + cell_delta[i][j]*cell_size[level];
						}

						r = compute_distance_periodic( halo_pos, point_pos ) / rvir;

						if ( r < rrvir[0] ) {
							bin = 0;
						} else {
							bin = (int)((log10(r) - rlvirmin)/drlvir) + 1;
						}

						cart_assert( bin >= num_vir_bins || (r >= rlvir[bin] && r <= rrvir[bin]) );

						min_bin = min( min_bin, bin );
						max_bin = max( max_bin, bin );
					}

					if ( min_bin < num_vir_bins ) {
						/* if only in one bin, just add completely */
						if ( max_bin == min_bin ) {
							bin_gas_mass[min_bin] += cell_mass;

							if ( Tcell <= Tcold ) {
								bin_cold_gas_mass[min_bin] += cell_mass;
							}

							bin_gas_pressure[min_bin] += Pcell*cell_mass;
							bin_gas_entropy[min_bin] += Scell*cell_mass;
							bin_gas_temperature[min_bin] += Tcell*cell_mass;

#ifdef COOLING
							bin_gas_coolingrate[min_bin] += dEcell*cell_volume[level];
							bin_gas_tcool[min_bin] += tcool*cell_mass;
#endif /* COOLING */

#ifdef ENRICH
							bin_gas_metallicity_II[min_bin] += cell_gas_metal_density_II(icell)*rhogi*cell_mass;
							bin_gas_metallicity_Ia[min_bin] += cell_gas_metal_density_Ia(icell)*rhogi*cell_mass;
#endif /* ENRICH */

							bin_gas_velocity[0][min_bin] += cell_vx*cell_mass;
							bin_gas_velocity[1][min_bin] += cell_vy*cell_mass;
							bin_gas_velocity[2][min_bin] += cell_vz*cell_mass;

							bin_gas_vrms[min_bin] += (cell_vx*cell_vx + cell_vy*cell_vy + cell_vz*cell_vz)*cell_mass;
						} else {
							/* monte carlo cell properties into bins */
							for ( bin = 0; bin < num_vir_bins; bin++ ) {
								bin_volume_fraction[bin] = 0;
							}

							for ( throw_point = 0; throw_point < points_per_cell; throw_point++ ) {
								for ( j = 0; j < nDim; j++ ) {
									point_pos[j] = cell_pos[j] + (cart_rand()-0.5)*cell_size[level];
								}

								r = compute_distance_periodic( halo_pos, point_pos ) / rvir;

								if ( r < rrvir[0] ) {
									bin = 0;
								} else {
									bin = (int)((log10(r) - rlvirmin)/drlvir) + 1;
								}

								if ( bin < num_vir_bins ) {
									bin_volume_fraction[bin]++;
								}
							}

							for ( bin = min_bin; bin <= min( max_bin, num_vir_bins-1 ); bin++ ) {
								cart_assert( bin >= 0 && bin < num_vir_bins );

								if ( bin_volume_fraction[bin] > 0 ) {
									volume_fraction = (double)bin_volume_fraction[bin] /
										(double)points_per_cell;
									mass_fraction = cell_mass*volume_fraction;

									bin_gas_mass[bin] += mass_fraction;

									if ( Tcell <= Tcold ) {
										bin_cold_gas_mass[bin] += mass_fraction;
									}

									bin_gas_pressure[bin] += Pcell*mass_fraction;
									bin_gas_entropy[bin] += Scell*mass_fraction;
									bin_gas_temperature[bin] += Tcell*mass_fraction;

#ifdef COOLING
									bin_gas_coolingrate[min_bin] += dEcell*volume_fraction*cell_volume[level];
									bin_gas_tcool[min_bin] += tcool*mass_fraction;
#endif /* COOLING */

#ifdef ENRICH
									bin_gas_metallicity_II[bin] += cell_gas_metal_density_II(icell)*rhogi*mass_fraction;
									bin_gas_metallicity_Ia[bin] += cell_gas_metal_density_Ia(icell)*rhogi*mass_fraction;
#endif /* ENRICH */

									bin_gas_velocity[0][bin] += cell_vx*mass_fraction;
									bin_gas_velocity[1][bin] += cell_vy*mass_fraction;
									bin_gas_velocity[2][bin] += cell_vz*mass_fraction;
								}
							}
						}
					}
#else 
					cart_assert( bin >= 0 );
					if ( bin < num_vir_bins ) {
						bin_gas_mass[bin] += cell_mass;

						if ( Tcell <= Tcold ) {
							bin_cold_gas_mass[bin] += cell_mass;
						}

						bin_gas_pressure[bin] += Pcell*cell_mass;
						bin_gas_entropy[bin] += Scell*cell_mass;
						bin_gas_temperature[bin] += Tcell*cell_mass;
		
#ifdef COOLING
						bin_gas_coolingrate[bin] += dEcell*cell_volume[level];
						bin_gas_tcool[bin] += tcool*cell_mass;
#endif /* COOLING */

#ifdef ENRICH
						bin_gas_metallicity_II[bin] += cell_gas_metal_density_II(icell)*rhogi*cell_mass;
						bin_gas_metallicity_Ia[bin] += cell_gas_metal_density_Ia(icell)*rhogi*cell_mass;
#endif /* ENRICH */
						bin_gas_velocity[0][bin] += cell_vx*cell_mass;
						bin_gas_velocity[1][bin] += cell_vy*cell_mass;
						bin_gas_velocity[2][bin] += cell_vz*cell_mass;

						bin_gas_vrms[bin] += (cell_vx*cell_vx + cell_vy*cell_vy + cell_vz*cell_vz)*cell_mass;
					}
#endif /* points per bin > 1 */
				}

				cart_debug("done assigning gas to bins in virial radius"); fflush(stdout);
#endif /* HYDRO */

				if ( halos->list[ihalo].proc == local_proc_id ) {
					/* receive bin information from other processors */
					for ( proc = 0; proc < num_procs; proc++ ) {
						if ( processor_mask[proc] ) {
							recv_int_bins( bin_dark_num, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_dark_mass, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_dark_momentum[0], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_dark_momentum[1], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_dark_momentum[2], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_dark_vrms, num_vir_bins, proc, ihalo );
							recv_int_bins( bin_star_num, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_mass, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_momentum[0], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_momentum[1], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_momentum[2], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_vrms, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_age, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_metallicity_II, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_star_metallicity_Ia, num_vir_bins, proc, ihalo );
							recv_int_bins( bin_new_star_num, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_new_star_mass, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_new_star_momentum[0], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_new_star_momentum[1], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_new_star_momentum[2], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_new_star_vrms, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_new_star_metallicity_II, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_new_star_metallicity_Ia, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_mass, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_velocity[0], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_velocity[1], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_velocity[2], num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_vrms, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_cold_gas_mass, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_pressure, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_entropy, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_temperature, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_tcool, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_coolingrate, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_metallicity_II, num_vir_bins, proc, ihalo );
							recv_float_bins( bin_gas_metallicity_Ia, num_vir_bins, proc, ihalo );
						}
					}

					/* process bins */
					for ( bin = 0; bin < num_vir_bins; bin++ ) {
						if ( bin_gas_mass[bin] > 0.0 ) {
							inverse_mass = 1.0 / bin_gas_mass[bin];

							bin_gas_pressure[bin] *= inverse_mass;
							bin_gas_entropy[bin] *= inverse_mass;
							bin_gas_temperature[bin] *= inverse_mass;
							bin_gas_coolingrate[bin] /= rvir_bin_volume[bin]*rvir*rvir*rvir;
							bin_gas_tcool[bin] *= inverse_mass;
							bin_gas_metallicity_II[bin] *= inverse_mass;
							bin_gas_metallicity_Ia[bin] *= inverse_mass;
							bin_gas_velocity[0][bin] *= inverse_mass;
							bin_gas_velocity[1][bin] *= inverse_mass;
							bin_gas_velocity[2][bin] *= inverse_mass;
							vmean = bin_gas_velocity[0][bin]*bin_gas_velocity[0][bin] +
								bin_gas_velocity[1][bin]*bin_gas_velocity[1][bin] +
								bin_gas_velocity[2][bin]*bin_gas_velocity[2][bin];
							bin_gas_vrms[bin] = sqrt(fabs(bin_gas_vrms[bin]*inverse_mass - vmean))*vfact;
							bin_gas_velocity[0][bin] *= vfact;
							bin_gas_velocity[1][bin] *= vfact;
							bin_gas_velocity[2][bin] *= vfact;
						}
					}

					/* compute cumulative quantities */
					total_dark_mass = 0.0;
					total_star_mass = 0.0;
					total_new_star_mass = 0.0;
					total_gas_mass = 0.0;
					total_cold_gas_mass = 0.0;

					avg_gas_metallicity_II = 0.0;
					avg_gas_metallicity_Ia = 0.0;
					avg_star_metallicity_II = 0.0;
					avg_star_metallicity_Ia = 0.0;
					avg_new_star_metallicity_II = 0.0;
					avg_new_star_metallicity_Ia = 0.0;
					avg_star_age = 0.0;

					avg_new_star_mass = 0.0;
					avg_star_mass = 0.0;
					avg_gas_mass = 0.0;

					for ( bin = 0; bin < num_vir_bins; bin++ ) {
						cart_assert( bin_dark_mass[bin] >= 0.0 );
						cart_assert( bin_star_mass[bin] >= 0.0 );
						cart_assert( bin_new_star_mass[bin] >= 0.0 );

						total_dark_mass += bin_dark_mass[bin];
						total_star_mass += bin_star_mass[bin];
						total_new_star_mass += bin_new_star_mass[bin];
						total_gas_mass += bin_gas_mass[bin];
						total_cold_gas_mass += bin_cold_gas_mass[bin];

						bin_total_dark_mass[bin] = total_dark_mass;
						bin_total_star_mass[bin] = total_star_mass;
						bin_total_new_star_mass[bin] = total_new_star_mass;
						bin_total_gas_mass[bin] = total_gas_mass;
						bin_total_cold_gas_mass[bin] = total_cold_gas_mass;

						avg_gas_metallicity_II += bin_gas_metallicity_II[bin] * bin_gas_mass[bin];
						avg_gas_metallicity_Ia += bin_gas_metallicity_Ia[bin] * bin_gas_mass[bin];
						avg_star_metallicity_II += bin_star_metallicity_II[bin];
						avg_star_metallicity_Ia = bin_star_metallicity_Ia[bin];
						avg_new_star_metallicity_II += bin_new_star_metallicity_II[bin];
						avg_new_star_metallicity_Ia += bin_new_star_metallicity_Ia[bin];
						avg_star_age += bin_star_age[bin];

						avg_gas_mass += bin_gas_mass[bin];
						avg_new_star_mass += bin_new_star_mass[bin];
						avg_star_mass += bin_star_mass[bin];
					}

					/* never reached virial radius */
					aM_gas = total_gas_mass * aM0 * cosmology->h;
					aM_cold_gas = total_cold_gas_mass * aM0 * cosmology->h;
					aM_stars = total_star_mass * aM0 * cosmology->h;
					aM_new_stars = total_new_star_mass * aM0 * cosmology->h;
					aM_dark = total_dark_mass * aM0 * cosmology->h;

					/* compute average quantities */
					if ( avg_gas_mass > 0.0 ) {
						avg_gas_metallicity_II /= avg_gas_mass*constants->Zsun;
						avg_gas_metallicity_Ia /= avg_gas_mass*constants->Zsun;
					}

					if ( avg_star_mass > 0.0 ) {
						avg_star_metallicity_II /= avg_star_mass*constants->Zsun;
						avg_star_metallicity_Ia /= avg_star_mass*constants->Zsun;
						avg_star_age = tphys_from_tcode( avg_star_age / avg_star_mass );
					}

					if ( avg_new_star_mass > 0.0 ) {
						avg_new_star_metallicity_II /= avg_new_star_mass*constants->Zsun;
						avg_new_star_metallicity_Ia /= avg_new_star_mass*constants->Zsun;
					}

					aM_baryons = aM_stars + aM_gas;
					aM_total = aM_baryons + aM_dark;

					rdout = rvir * rfact * cosmology->h;
					rvdout = rvir * rfact * cosmology->h;

					fprintf( bmprovir, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id, rdout, rvdout, aM_gas, aM_cold_gas, 
							aM_stars, aM_new_stars, aM_baryons, aM_dark, aM_total, 
							aM0*cosmology->h*halos->list[ihalo].mvir, vmax, rmax, avg_gas_metallicity_II, 
							avg_gas_metallicity_Ia, avg_star_metallicity_II, avg_star_metallicity_Ia,
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

					/* write out profiles */
					for ( bin = 0; bin < num_vir_bins; bin++ ) {
						fprintf( bmprovir, "%.3f %.3f %e %e %e %e %e\n",
								rmidvir[bin], rrvir[bin],
								bin_total_dark_mass[bin]*aM0*cosmology->h,
								bin_total_gas_mass[bin]*aM0*cosmology->h,
								bin_total_cold_gas_mass[bin]*aM0*cosmology->h,
								bin_total_star_mass[bin]*aM0*cosmology->h,
								bin_total_new_star_mass[bin]*aM0*cosmology->h );
					}

					fflush(bmprovir);

					fprintf( bvprovir, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id, rdout, rvdout, aM_gas, aM_cold_gas, 
							aM_stars, aM_new_stars, aM_baryons, aM_dark, aM_total, 
							aM0*cosmology->h*halos->list[ihalo].mvir, vmax, rmax, avg_gas_metallicity_II, 
							avg_gas_metallicity_Ia, avg_star_metallicity_II, avg_star_metallicity_Ia,
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

					for ( bin = 0; bin < num_vir_bins; bin++ ) {
						rcirc = rrvir[bin]*rvir*rfact; /* code units to proper kpc */
						vcirc_dm = 2.07498e-3 * sqrt( aM0*bin_total_dark_mass[bin] / rcirc );
						vcirc_gas = 2.07498e-3 * sqrt( aM0*bin_total_gas_mass[bin] / rcirc );
						vcirc_stars = 2.07498e-3 * sqrt( aM0*bin_total_star_mass[bin] / rcirc );
						vcirc_total = 2.07498e-3 * sqrt( aM0*(bin_total_dark_mass[bin]+
									bin_total_gas_mass[bin]+bin_total_star_mass[bin])/rcirc );

						bin_dark_vrms[bin] = vfact*sqrt(bin_dark_vrms[bin]);
						bin_gas_vrms[bin] = vfact*sqrt(bin_gas_vrms[bin]);
						bin_star_vrms[bin] = vfact*sqrt(bin_star_vrms[bin]);
						bin_new_star_vrms[bin] = vfact*sqrt(bin_new_star_vrms[bin]);

						fprintf( bvprovir, "%.3f %.3f %e %e %e %e %.3f %.3f %.3f %.3f\n",
								rmidvir[bin], rrvir[bin],
								bin_dark_vrms[bin], bin_gas_vrms[bin],
								bin_star_vrms[bin], bin_new_star_vrms[bin],
								vcirc_dm, vcirc_gas, vcirc_stars, vcirc_total );
					}

					fflush(bvprovir);

#ifdef HYDRO
					fprintf( bgprovir, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id,
							rdout, rvdout, aM_gas, aM_cold_gas,
							aM_stars, aM_new_stars, aM_baryons,
							aM_dark, aM_total, aM0*cosmology->h*halos->list[ihalo].mvir,
							vmax, rmax, avg_gas_metallicity_II, avg_gas_metallicity_Ia,
							avg_star_metallicity_II, avg_star_metallicity_Ia,
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

					for ( bin = 0; bin < num_vir_bins; bin++ ) {
						fprintf( bgprovir, "%.3f %.3f %e %e %e %e %e %.3f %.3f %e %e\n",
								rmidvir[bin], rrvir[bin],
								bin_total_gas_mass[bin]*aM0*cosmology->h,
								bin_total_cold_gas_mass[bin]*aM0*cosmology->h,
								bin_gas_temperature[bin],
								bin_gas_pressure[bin],
								bin_gas_entropy[bin],
								bin_gas_coolingrate[bin],	
								bin_gas_tcool[bin],
								bin_gas_metallicity_II[bin]/constants->Zsun,
								bin_gas_metallicity_Ia[bin]/constants->Zsun
						       );
					}

					fflush(bgprovir);
#endif /* HYDRO */

#ifdef STARFORM
					fprintf( bzprovir, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
							halos->list[ihalo].id, 
							rdout, rvdout, aM_gas, aM_cold_gas, 
							aM_stars, aM_new_stars, aM_baryons,
							aM_dark, aM_total, aM0*cosmology->h*halos->list[ihalo].mvir, 
							vmax, rmax, avg_gas_metallicity_II, avg_gas_metallicity_Ia,
							avg_star_metallicity_II, avg_star_metallicity_Ia,
							avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

					for ( bin = 0; bin < num_vir_bins; bin++ ) {
						if ( bin_star_age[bin] > 0 ) {
							age_star = tphys_from_tcode( bin_star_age[bin] );
						} else {
							age_star = -1.0;
						}

						fprintf( bzprovir, "%.3f %.3f %e %e %e %e %e %e %.3f\n",
								rmidvir[bin], rrvir[bin],
								bin_gas_metallicity_II[bin]/constants->Zsun, 
								bin_gas_metallicity_Ia[bin]/constants->Zsun,
								bin_star_metallicity_II[bin]/constants->Zsun, 
								bin_star_metallicity_Ia[bin]/constants->Zsun,
								bin_new_star_metallicity_II[bin]/constants->Zsun, 
								bin_new_star_metallicity_Ia[bin]/constants->Zsun,
								age_star );
					}

					fflush(bzprovir);
#endif /* STARFORM */
				} else {
					/* send bin values to halo owner */
					MPI_Send( bin_dark_num, num_vir_bins, MPI_INT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_mass, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_momentum[0], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_momentum[1], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_momentum[2], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_dark_vrms, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_num, num_vir_bins, MPI_INT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_mass, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_momentum[0], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_momentum[1], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_momentum[2], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_vrms, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_age, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_metallicity_II, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_star_metallicity_Ia, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_num, num_vir_bins, MPI_INT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_mass, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_momentum[0], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_momentum[1], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_momentum[2], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_vrms, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_metallicity_II, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_new_star_metallicity_Ia, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_mass, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_velocity[0], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_velocity[1], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_velocity[2], num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_vrms, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_cold_gas_mass, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_pressure, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_entropy, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_temperature, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_tcool, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_coolingrate, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_metallicity_II, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
					MPI_Send( bin_gas_metallicity_Ia, num_vir_bins, MPI_FLOAT, halos->list[ihalo].proc, ihalo, MPI_COMM_WORLD );
				}

#ifdef HYDRO
				cart_free( leaf_index );
#endif
			}
		}
	}

	/* finish up */
	fclose( blist );
	fclose( bmpro );
	fclose( bvpro );
	fclose( bmprovir );
	fclose( bvprovir );

#ifdef HYDRO
	fclose( bszlist );
#ifdef ANALYZE_XRAY
	fclose( btxlist );
#endif /* ANALYZE_XRAY */
	fclose( bgpro );
	fclose( bgprovir );
#endif /* HYDRO */

#ifdef STARFORM
	fclose( bzpro );
	fclose( bzprovir );
#endif /* STARFORM */

	for ( i = 0; i < num_radii; i++ ) {
		fclose( rlist[i] );
	}
}

#endif /* COSMOLOGY */


int halo_level( const halo *h, MPI_Comm local_comm )
{
  int llev, glev, cell;

  cart_assert(h);

  cell = cell_find_position((double *)h->pos);

  if(cell < 0)
    {
      llev = min_level - 1;
    }
  else
    {
      llev = cell_level(cell);
    }

  MPI_Allreduce(&llev,&glev,1,MPI_INT,MPI_MAX,local_comm);

  return glev;
}


halo* find_halo_by_id(halo_list *halos, int id)
{
  int ih;

  cart_assert(halos);
  for(ih=0; ih<halos->num_halos; ih++)
    {
      if(halos->list[ih].id == id) return &halos->list[ih];
    }
  return NULL;
}


void dump_region_around_halo(const char *filename, const halo *h, float size)
{
  int i, j, n, nbuf;
  int *ids;
  FILE *f;

  cart_assert(h != NULL);

  for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
    {
      if(compute_distance_periodic(particle_x[j],(double *)h->pos) < h->rvir*size)
        {
          n++;
        }
    }

  nbuf = n + 1;  /* in case n is zero */
  ids = cart_alloc(int,nbuf);
  
  for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
    {
      /*
      //  If particle is within the size*Rvir, save its id
      */
      if(compute_distance_periodic(particle_x[j],(double *)h->pos) < h->rvir*size)
        {
          ids[n++] = particle_id[j];
        }
    }

  /*
  //  Write from the master node
  */
  if(local_proc_id == MASTER_NODE)
    {

      cart_debug("Dumping region for halo %d",h->id);

      f = fopen(filename,"w");
      cart_assert(f != NULL);

      for(j=0; j<n; j++)
	{
	  fprintf(f,"%d\n",ids[j]);
	}

      for(i=1; i<num_procs; i++)
        {
          MPI_Recv(&n,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          if(n > nbuf)
            {
              cart_free(ids);
	      nbuf = n;
              ids = cart_alloc(int,nbuf);
            }
          MPI_Recv(ids,n,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  for(j=0; j<n; j++)
	    {
	      fprintf(f,"%d\n",ids[j]);
	    }
        }

      fclose(f);
    }
  else
    {
      MPI_Send(&n,1,MPI_INT,MASTER_NODE,0,MPI_COMM_WORLD);
      MPI_Send(ids,n,MPI_INT,MASTER_NODE,0,MPI_COMM_WORLD);
   }

  cart_free(ids);
}


/*
//  Set cell_var(c,var) with the halo id for each halo, or 0 if belongs to 
//  none; a cell belongs to a halo if it is inside its size_factor*Rtrunc, 
//  and satellites are accounted for properly.
*/
void map_halos(int var, int resolution_level, halo_list *halos, float size_factor)
{
  int j, ih, iold, *halo_levels;
  MESH_RUN_DECLARE(level,cell);
  double pos[3], dx, r2, r2old, r2Cut;

  cart_assert(halos != NULL);
  cart_assert(var>=0 && var<num_vars);

  if(halos->map == var) return; /* Already mapped */

  halos->map = var;

  halo_levels = cart_alloc(int,halos->num_halos);
  for(ih=0; ih<halos->num_halos; ih++) halo_levels[ih] = halo_level(&halos->list[ih],MPI_COMM_WORLD);

  /*
  //  Loop over levels first to avoid selecting cells multiple times
  */
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  cart_debug("Mapping level %d...",level);

  /*
  //  Zero map array
  */
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,var)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell)) cell_var(cell,var) = 0.0;
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  for(ih=0; ih<halos->num_halos; ih++) if(halo_levels[ih] >= resolution_level)
    {
      r2Cut = pow(size_factor*halos->list[ih].rhalo,2.0);

      /*
      //  Map halo indices(+1), not ids, first (to simplify inter-comparison)
      */
#pragma omp parallel for default(none), private(_Index,cell,j,dx,r2,iold,r2old,pos), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,var,r2Cut,halos,ih)
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
      if(cell_is_leaf(cell))
	{
	  cell_center_position(cell,pos);
	  for(j=0, r2=0.0; j<nDim; j++)
	    {
	      dx = pos[j] - halos->list[ih].pos[j];
	      if(dx < -0.5*num_grid) dx += num_grid;
	      if(dx >  0.5*num_grid) dx -= num_grid;
	      r2 += dx*dx;
	    }
	  
	  if(r2 < r2Cut)
	    {
	      iold = (int)(0.5+cell_var(cell,var));
	      if(iold == 0)
		{
		  cell_var(cell,var) = ih + 1;
		}
	      else
		{
		  cart_assert(iold>=1 && iold<=halos->num_halos);
		  iold--;
		  for(j=0, r2old=0.0; j<nDim; j++)
		    {
		      dx = pos[j] - halos->list[iold].pos[j];
		      if(dx < -0.5*num_grid) dx += num_grid;
		      if(dx >  0.5*num_grid) dx -= num_grid;
		      r2old += dx*dx;
		    }
		  if(r2old > r2)
		    {
		      /*
		      //  This cells belongs to a satellite
		      */
		      cell_var(cell,var) = ih + 1;
		    }
		}
	    }
	}

      MESH_RUN_OVER_CELLS_OF_LEVEL_END;

    }

  MESH_RUN_OVER_LEVELS_END;

  cart_free(halo_levels);
}

