#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "index_hash.h"
#include "io.h"
#include "io_artio.h"
#include "io_cart.h" 
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "rt_transfer.h"
#include "sfc.h"
#include "skiplist.h"
#include "starformation.h"
#include "timing.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#ifdef RADIATIVE_TRANSFER
#include "frt/frt_c.h"
#endif /* RADIATIVE_TRANSFER */

#include "../tools/artio/artio.h"

extern int step;

DECLARE_LEVEL_ARRAY(double,dtl);
DECLARE_LEVEL_ARRAY(double,tl_old);
DECLARE_LEVEL_ARRAY(double,dtl_old);
DECLARE_LEVEL_ARRAY(int,level_sweep_dir);
DECLARE_LEVEL_ARRAY(int,time_refinement_factor);
DECLARE_LEVEL_ARRAY(int,time_refinement_factor_old);

#ifdef COSMOLOGY
DECLARE_LEVEL_ARRAY(double,abox_old);
#endif

extern double auni_init;

int num_artio_grid_files = 1;
int artio_grid_allocation_strategy = ARTIO_ALLOC_EQUAL_SFC;

#ifdef PARTICLES
int num_artio_particle_files = 1;
int artio_particle_allocation_strategy = ARTIO_ALLOC_EQUAL_SFC;
#endif /* PARTICLES */

void artio_restart_load_balance( artio_file handle );

void control_parameter_set_allocation_strategy(const char *value, void *ptr, int ind) {
	if ( strcmp( value, "ARTIO_ALLOC_EQUAL_SFC" ) == 0 ) {
		*(int *)ptr = ARTIO_ALLOC_EQUAL_SFC;
	} else if ( strcmp( value, "ARTIO_ALLOC_EQUAL_PROC" ) == 0 ) {
		*(int *)ptr = ARTIO_ALLOC_EQUAL_PROC;
	} else {
		cart_error("Invalid allocation strategy %s", value );
	}
}

void control_parameter_list_allocation_strategy(FILE *stream, const void *ptr) {
	switch( *(int *)ptr ) {
		case ARTIO_ALLOC_EQUAL_SFC :
			fprintf(stream,"ARTIO_ALLOC_EQUAL_SFC");
			break;
		case ARTIO_ALLOC_EQUAL_PROC :
			fprintf(stream,"ARTIO_ALLOC_EQUAL_PROC");
			break;
		default :
			fprintf(stream,"UNKNOWN");
	}
}

void config_init_io_artio() {
	ControlParameterOps control_parameter_allocation_strategy = { 
		control_parameter_set_allocation_strategy, 
		control_parameter_list_allocation_strategy };

	num_artio_grid_files = num_procs;
	control_parameter_add2(control_parameter_int,&num_artio_grid_files,"io:num-grid-files","num_grid_files","Number of output grid files. Defaults to the number of MPI tasks.");
	control_parameter_add2(control_parameter_allocation_strategy,&artio_grid_allocation_strategy,"io:grid-file-allocation-strategy","grid_allocation_strategy","Determine how root cells are divided amongst the output files.  Supported options: ARTIO_ALLOC_EQUAL_SFC, ARTIO_ALLOC_EQUAL_PROC");

#ifdef PARTICLES
	num_artio_particle_files = num_procs;
	control_parameter_add2(control_parameter_int,&num_artio_particle_files,"io:num-particle-files","num_particle_files","Number of output particle files.  Defaults to the number of MPI tasks.");
	control_parameter_add2(control_parameter_allocation_strategy,&artio_particle_allocation_strategy,"io:particle-file-allocation-strategy","particle_allocation_strategy","Determine how root cells are divided amongst the output files.  Supported options: ARTIO_ALLOC_EQUAL_SFC, ARTIO_ALLOC_EQUAL_PROC");
#endif /* PARTICLES */
}

void config_verify_io_artio() {
	VERIFY(io:num-grid-files, num_artio_grid_files > 0 && num_artio_grid_files < max_sfc_index );
	VERIFY(io:grid-file-allocation-strategy, artio_grid_allocation_strategy == ARTIO_ALLOC_EQUAL_SFC || artio_grid_allocation_strategy == ARTIO_ALLOC_EQUAL_PROC );
#ifdef PARTICLES 
	VERIFY(io:num-particle-files, num_artio_particle_files > 0 && num_artio_particle_files < max_sfc_index );
	VERIFY(io:particle-file-allocation-strategy, artio_particle_allocation_strategy == ARTIO_ALLOC_EQUAL_SFC || artio_particle_allocation_strategy == ARTIO_ALLOC_EQUAL_PROC );
#endif /* PARTICLES */
}

void read_artio_grid( artio_file handle, int file_max_level );
void write_artio_grid( artio_file handle, int num_file_vars, int *var_indices );

#ifdef PARTICLES
void read_artio_particles( artio_file handle );
void write_artio_particles( artio_file handle, int *root_tree_particle_list,
            int *num_particles_per_species_per_root_tree );
#endif /* PARTICLES */

void define_file_variables(int *num_variables, char *variable_labels[num_vars],
		int variable_indices[num_vars]) {
	int j;

#define add_variable(label,index)    \
		variable_labels[*num_variables] = cart_alloc(char, strlen(label)+1); \
		strcpy( variable_labels[*num_variables], label ); \
		variable_indices[*num_variables] = index; \
		(*num_variables)++;

	*num_variables = 0;

#ifdef HYDRO
	add_variable( "HVAR_GAS_DENSITY", HVAR_GAS_DENSITY );
	add_variable( "HVAR_GAS_ENERGY", HVAR_GAS_ENERGY );
	add_variable( "HVAR_PRESSURE", HVAR_PRESSURE );
	add_variable( "HVAR_GAMMA", HVAR_GAMMA );
	add_variable( "HVAR_INTERNAL_ENERGY", HVAR_INTERNAL_ENERGY );
	add_variable( "HVAR_MOMENTUM_X", HVAR_MOMENTUM+0 );
	add_variable( "HVAR_MOMENTUM_Y", HVAR_MOMENTUM+1 );
	add_variable( "HVAR_MOMENTUM_Z", HVAR_MOMENTUM+2 );
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	add_variable( "HVAR_ELECTRON_INTERNAL_ENERGY", HVAR_ELECTRON_INTERNAL_ENERGY );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#ifdef RADIATIVE_TRANSFER
	add_variable( "RT_HVAR_HI", RT_HVAR_OFFSET+0 );
	add_variable( "RT_HVAR_HII", RT_HVAR_OFFSET+1 );
	add_variable( "RT_HVAR_HeI", RT_HVAR_OFFSET+2 );
	add_variable( "RT_HVAR_HeII", RT_HVAR_OFFSET+3 );
	add_variable( "RT_HVAR_HeIII", RT_HVAR_OFFSET+4 );
	add_variable( "RT_HVAR_H2", RT_HVAR_OFFSET+5 );
#endif /* RADIATIVE_TRANSFER */
#ifdef ENRICHMENT
	add_variable( "HVAR_METAL_DENSITY_II", HVAR_METAL_DENSITY_II );
#ifdef ENRICHMENT_SNIa
	add_variable( "HVAR_METAL_DENSITY_Ia", HVAR_METAL_DENSITY_Ia );
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef BLASTWAVE_FEEDBACK
	add_variable( "HVAR_BLASTWAVE_TIME", HVAR_BLASTWAVE_TIME );
#endif /* BLASTWAVE_FEEDBACK */
#endif /* HYDRO */

#ifdef GRAVITY
	add_variable( "VAR_POTENTIAL", VAR_POTENTIAL );
#ifdef HYDRO 
	add_variable( "VAR_POTENTIAL_HYDRO", VAR_POTENTIAL_HYDRO );
#endif /* HYDRO */
#endif /* GRAVITY */

#ifdef RADIATIVE_TRANSFER
	for( j = 0; j < rt_num_disk_vars; j++ ) { 
		variable_labels[*num_variables] = cart_alloc(char, 64);
		snprintf(variable_labels[*num_variables], 64, "RT_DISK_VAR_%d", j );
		variable_indices[*num_variables] = rt_disk_offset + j;
		(*num_variables)++;
	}
#endif /* RADIATIVE_TRANSFER */

#undef add_variable
}

#ifdef PARTICLES
void define_particle_variables( int *num_species, char *species_labels[MAX_PARTICLE_SPECIES],
		int num_primary_variables[MAX_PARTICLE_SPECIES], char **primary_variable_labels[MAX_PARTICLE_SPECIES],
		int num_secondary_variables[MAX_PARTICLE_SPECIES], char **secondary_variable_labels[MAX_PARTICLE_SPECIES] ) {

	int i, j;
	int num_nbody_species;
	const char dim[] = { 'X', 'Y', 'Z' };

	*num_species = num_particle_species;

#ifdef STAR_FORMATION
	num_nbody_species = num_particle_species - 1;
#else
	num_nbody_species = num_particle_species;
#endif /* STAR_FORMATION */

	for ( i = 0; i < num_nbody_species; i++ ) {
		species_labels[i] = cart_alloc(char, 64);
		snprintf(species_labels[i], 64, "N-BODY");
		num_primary_variables[i] = 2*nDim+1;
		num_secondary_variables[i] = 0;
	}

#ifdef STAR_FORMATION 
	species_labels[num_nbody_species] = cart_alloc(char, 64);
	snprintf(species_labels[num_nbody_species], 64, "STAR");
	num_primary_variables[num_nbody_species] = 2*nDim+1;

#ifdef ENRICHMENT
#ifdef ENRICHMENT_SNIa
	num_secondary_variables[num_nbody_species] = 5;
#else
	num_secondary_variables[num_nbody_species] = 4;
#endif /* ENRICHMENT_SNIa */
#else
	num_secondary_variables[num_nbody_species] = 3;
#endif /* ENRICHMENT */

	secondary_variable_labels[num_nbody_species] = cart_alloc( char *, num_secondary_variables[num_nbody_species] );

	secondary_variable_labels[num_nbody_species][0] = cart_alloc( char, 64 );
	snprintf( secondary_variable_labels[num_nbody_species][0], 64, "BIRTH_TIME" );

	secondary_variable_labels[num_nbody_species][1] = cart_alloc( char, 64 );
	snprintf( secondary_variable_labels[num_nbody_species][1], 64, "INITIAL_MASS" );

	secondary_variable_labels[num_nbody_species][2] = cart_alloc( char, 64 );
	snprintf( secondary_variable_labels[num_nbody_species][2], 64, "MASS" );	

#ifdef ENRICHMENT
	secondary_variable_labels[num_nbody_species][3] = cart_alloc( char, 64 );
	snprintf( secondary_variable_labels[num_nbody_species][3], 64, "METALLICITY_SNII" );    
#ifdef ENRICHMENT_SNIa
	secondary_variable_labels[num_nbody_species][4] = cart_alloc( char, 64 );
    snprintf( secondary_variable_labels[num_nbody_species][4], 64, "METALLICITY_SNIa" );
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */

#endif /* STAR_FORMATION */

	/* copy primary variable names */
	for ( i = 0; i < *num_species; i++ ) {
		primary_variable_labels[i] = cart_alloc( char *, num_primary_variables[i] );
		for ( j = 0; j < nDim; j++ ) {
			primary_variable_labels[i][j] = cart_alloc( char, 64 );
			snprintf( primary_variable_labels[i][j], 64, "POSITION_%c", dim[j] );
		}

		for ( j = 0; j < nDim; j++ ) {
            primary_variable_labels[i][nDim+j] = cart_alloc( char, 64 );
            snprintf( primary_variable_labels[i][nDim+j], 64, "VELOCITY_%c", dim[j] );
        }

		primary_variable_labels[i][2*nDim] = cart_alloc( char, 64 );
		snprintf( primary_variable_labels[i][2*nDim], 64, "TIMESTEP" );
	}
}
#endif /* PARTICLES */

void pack_cell_vars(int icell, int num_pack_vars, int *var_indices,
		float *variables) {
	int i;

	for (i = 0; i < num_pack_vars; i++) {
		variables[i] = cell_var(icell, var_indices[i]);
	}
}

void unpack_cell_vars(int icell, int num_pack_vars, int *sim_var_indices,
		int *file_var_indices, float *variables) {
	int i;

	for (i = 0; i < num_pack_vars; i++) {
		cell_var(icell, sim_var_indices[i]) = variables[file_var_indices[i]];
	}
}

#ifdef PARTICLES
int compare_particle_species_id( const void *a, const void *b ) {
	int id1 = particle_id[*(int *)a];
	int id2 = particle_id[*(int *)b];
	int species1 = particle_species(id1);
	int species2 = particle_species(id2);

	return ( species1 == species2 ) ? id1 - id2 : species1 - species2;
}
#endif /* PARTICLES */

void create_artio_filename( int filename_flag, char *label, char *filename ) {
	switch(filename_flag) {
		case WRITE_SAVE:
			sprintf( filename, "%s/%s_%s", output_directory, jobname, label );
			break;
		case WRITE_BACKUP:
			sprintf( filename, "%s/%s_2", output_directory, jobname );
			break;
		case WRITE_GENERIC:
			sprintf( filename, "%s/%s", output_directory, jobname );
			break;
	}
}

void write_artio_restart_worker( char *filename, int fileset_write_options );
#ifdef HYDRO_TRACERS
char filename_tracers[256];
#endif /* HYDRO_TRACERS */

void write_artio_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag ) { 
	/* ignores tracer_filename_flag for now... */
	char label[256];
	char filename[256];
	int fileset_options = 0;

#ifdef COSMOLOGY
	sprintf(label,"a%06.4f",auni[min_level]);
#else
	sprintf(label,"%06d",step);
#endif /* COSMOLOGY */

#ifdef HYDRO_TRACERS
	/* DHR - hack to allow writing cart format tracer files for now */
	if ( tracer_filename_flag != NO_WRITE ) {
		fileset_options |= WRITE_TRACERS;

		switch(tracer_filename_flag) {
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
			case WRITE_SAVE:
				sprintf( filename_tracers, "%s/%s_%s.dtr", output_directory, jobname, label );
				break;
			case WRITE_BACKUP:
				sprintf( filename_tracers, "%s/%s_2.dtr", output_directory, jobname );
				break;
			case WRITE_GENERIC:
				sprintf( filename_tracers, "%s/%s.dtr", output_directory, jobname );
				break;
#else  /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
			case WRITE_SAVE:
				sprintf( filename_tracers, "%s/tracers_%s.dat", output_directory, label );
				break;
			case WRITE_BACKUP:
				sprintf( filename_tracers, "%s/tracers_2.dat", output_directory );
				break;
			case WRITE_GENERIC:
				sprintf( filename_tracers, "%s/tracers.dat", output_directory );
				break;
#endif /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
			default:
				cart_error("Invalid value for tracer_filename_flag in write_artio_restart!");
		}
	}
#endif /* HYDRO_TRACERS */	

	if ( particle_filename_flag == NO_WRITE ) {
		if ( grid_filename_flag != NO_WRITE ) {
			create_artio_filename(grid_filename_flag, label, filename);
			write_artio_restart_worker( filename, fileset_options | WRITE_GRID );
		} else if ( fileset_options ) {
			create_artio_filename(grid_filename_flag, label, filename);
			write_artio_restart_worker( filename, fileset_options );
		}
	} else {
		if ( particle_filename_flag == grid_filename_flag ) {
			create_artio_filename(particle_filename_flag, label, filename);
			write_artio_restart_worker( filename, fileset_options | WRITE_GRID | WRITE_PARTICLES );
		} else {
			if ( grid_filename_flag == WRITE_BACKUP || grid_filename_flag == WRITE_GENERIC ) {
				create_artio_filename(grid_filename_flag, label, filename);
				write_artio_restart_worker( filename, fileset_options | WRITE_GRID | WRITE_PARTICLES );

				create_artio_filename(particle_filename_flag, label, filename);
				write_artio_restart_worker( filename, WRITE_PARTICLES );
			} else if ( grid_filename_flag == WRITE_SAVE ) {
				/* implies particle_filename_flag == BACKUP || GENERIC */
				create_artio_filename(particle_filename_flag, label, filename);
				write_artio_restart_worker( filename, fileset_options | WRITE_GRID | WRITE_PARTICLES );

				create_artio_filename(grid_filename_flag, label, filename);
				write_artio_restart_worker( filename, WRITE_GRID );
			} else {
				/* implies particle_filename_flag == WRITE_SAVE */
				create_artio_filename(particle_filename_flag, label, filename);
				write_artio_restart_worker( filename, fileset_options | WRITE_PARTICLES );
			}
		}
	}
}

void write_artio_restart_worker( char *filename, int fileset_write_options ) {
	int i, j;
	int sfc;
	int icell;
	int iroot, ipart;
	int level;
	long total;
	int num_level_cells, *level_cells;
	artio_file handle;
	int num_file_vars;
	int var_indices[num_vars];
	char *var_labels[num_vars];
	int num_levels;
	int *num_levels_per_root_tree;
	int *num_octs_per_root_tree;
	int64_t local_proc_sfc_index[MAX_PROCS+1];
	FILE *restart;
	char restart_filename[256];
#ifdef RADIATIVE_TRANSFER
	frt_intg n;
	frt_real *data;
	float *buffer;
#endif /* RADIATIVE_TRANSFER */
	struct artio_context_struct con = { mpi.comm.run };

#ifdef PARTICLES
	int num_species;
	char *species_labels[MAX_PARTICLE_SPECIES];
	int num_primary_variables[MAX_PARTICLE_SPECIES];
	char **primary_variable_labels[MAX_PARTICLE_SPECIES];
	int num_secondary_variables[MAX_PARTICLE_SPECIES];
	char **secondary_variable_labels[MAX_PARTICLE_SPECIES];

	int *root_tree_particle_list;
	int *order;
	int *num_particles_per_species_per_root_tree;
#endif

	if ( fileset_write_options & WRITE_GRID || fileset_write_options & WRITE_PARTICLES ) {
		cart_debug("creating parallel file %s jobname %s", filename, jobname );
		handle = artio_fileset_create( filename, num_root_cells,
				proc_sfc_index[local_proc_id],
				proc_sfc_index[local_proc_id+1]-1, &con );
		if ( handle == NULL ) {
			cart_error("Unable to create parallel file %s", filename );
		}

		/* write header variables */
		num_levels = max_level_now_global(mpi.comm.run) - min_level + 1;

		artio_parameter_set_string( handle, "jobname", (char *)jobname );

		artio_parameter_set_int( handle, "sfc", SFC );
		artio_parameter_set_int( handle, "max_refinement_level", num_levels-1 );

		/* timestepping variables */
		artio_parameter_set_int( handle, "step", step );
		artio_parameter_set_double_array( handle, "tl", num_levels, tl );
		artio_parameter_set_double_array( handle, "tl_old", num_levels, tl_old );
		artio_parameter_set_double_array( handle, "dtl", num_levels, dtl );
		artio_parameter_set_double_array( handle, "dtl_old", num_levels, dtl_old );

		artio_parameter_set_int_array( handle, "time_refinement_factor", num_levels, time_refinement_factor );
		artio_parameter_set_int_array( handle, "time_refinement_factor_old", num_levels, time_refinement_factor_old );

#ifdef HYDRO	
		artio_parameter_set_int_array( handle, "hydro_sweep_direction", num_levels, level_sweep_dir );
#endif

		/* unit parameters */
		artio_parameter_set_double( handle, "box_size", box_size );
#ifdef COSMOLOGY 
		artio_parameter_set_double( handle, "OmegaM", cosmology->OmegaM );
		artio_parameter_set_double( handle, "OmegaL", cosmology->OmegaL );
		artio_parameter_set_double( handle, "OmegaB", cosmology->OmegaB );
		artio_parameter_set_double( handle, "hubble", cosmology->h );
		artio_parameter_set_double( handle, "DeltaDC", cosmology->DeltaDC );

		artio_parameter_set_double( handle, "auni_init", auni_init );
		artio_parameter_set_double_array( handle, "abox", num_levels, abox );
		artio_parameter_set_double_array( handle, "auni", num_levels, auni );

		artio_parameter_set_double( handle, "Hbox", Hubble(abox[min_level]) );
#endif /* COSMOLOGY */

		/* write unit conversion factors independent of cosmology flag */
		artio_parameter_set_double( handle, "mass_unit", primary_units->mass );
		artio_parameter_set_double( handle, "time_unit", primary_units->time );
		artio_parameter_set_double( handle, "length_unit", primary_units->length );

#ifdef PARTICLES
		/* energy conservation variables */	
		artio_parameter_set_double( handle, "energy:tintg", tintg );
		artio_parameter_set_double( handle, "energy:ekin", ekin );
		artio_parameter_set_double( handle, "energy:ekin1", ekin1 );
		artio_parameter_set_double( handle, "energy:ekin2", ekin2 );
		artio_parameter_set_double( handle, "energy:au0", au0 );
		artio_parameter_set_double( handle, "energy:aeu0", aeu0 );
		artio_parameter_set_double( handle, "energy:ap0", ap0 );
#endif /* PARTICLES */

		/* refinement boundaries */
		artio_parameter_set_float_array( handle, "refinement_volume_min", 
				nDim, refinement_volume_min );
		artio_parameter_set_float_array( handle, "refinement_volume_max", 
				nDim, refinement_volume_max );

#ifdef STAR_FORMATION
		artio_parameter_set_float_array( handle, "star_formation_volume_min", 
				nDim, star_formation_volume_min );
		artio_parameter_set_float_array( handle, "star_formation_volume_max", 
				nDim, star_formation_volume_max );

		artio_parameter_set_double( handle, "total_stellar_mass", total_stellar_mass );
		artio_parameter_set_double( handle, "total_stellar_initial_mass", total_stellar_initial_mass );
#endif /* STAR_FORMATION */

#ifdef RADIATIVE_TRANSFER
		n = 0;
		frtCall(packradiationfields)(&n,0);
		if(n < 1)
		{
			cart_error("Unable to pack Radiation Field data.");
		}

		data = cart_alloc(frt_real, n );
		frtCall(packradiationfields)(&n,data);

		if(sizeof(float) == sizeof(frt_real))
		{
			buffer = (float *)data;
		}
		else
		{
			buffer = cart_alloc(float,n);
			for(i=0; i<n; i++) buffer[i] = data[i];
		}

		artio_parameter_set_float_array( handle, "radiation_background", n, buffer );

		if(sizeof(float) == sizeof(frt_real))
		{
		}
		else
		{
			cart_free(buffer);
		}

		cart_free(data);

#ifdef RT_SINGLE_SOURCE
		artio_parameter_set_int( handle, "rt:single_source:level", rtSingleSourceLevel );
		artio_parameter_set_float( handle, "rt:single_source:value", rtSingleSourceValue );
		artio_parameter_set_double_array( handle, "rt:single_source:pos", nDim, rtSingleSourcePos );
#endif
#endif /* RADIATIVE_TRANSFER */

		/* load balance parameters */
		artio_parameter_set_int( handle, "num_octs_per_mpi_task", num_octs );
#ifdef PARTICLES
		artio_parameter_set_int( handle, "num_particles_per_mpi_task", num_particles );
#ifdef STAR_FORMATION
		artio_parameter_set_int( handle, "num_star_particles_per_mpi_task", num_star_particles );
#endif /* STAR_FORMATION */
#endif /* PARTICLES */

		for ( i = 0; i < num_procs+1; i++ ) {
			local_proc_sfc_index[i] = proc_sfc_index[i];
		}

		artio_parameter_set_long_array( handle, "mpi_task_sfc_index", num_procs+1, local_proc_sfc_index );

		if ( fileset_write_options & WRITE_GRID ) {
			start_time( GAS_WRITE_IO_TIMER );

			/* build list of variables to write */
			define_file_variables(&num_file_vars, var_labels, var_indices);

			/* collect tree information */
			num_levels_per_root_tree = cart_alloc( int, num_cells_per_level[min_level] );
			num_octs_per_root_tree = cart_alloc( int, num_cells_per_level[min_level] );

			for ( sfc = proc_sfc_index[local_proc_id]; sfc < proc_sfc_index[local_proc_id+1]; sfc++ ) {
				icell = root_cell_location(sfc);
				num_levels_per_root_tree[sfc-proc_sfc_index[local_proc_id]] = tree_max_level( icell );
				num_octs_per_root_tree[sfc-proc_sfc_index[local_proc_id]] = tree_oct_count( icell );
			}

			if ( artio_fileset_add_grid( handle,
					num_artio_grid_files, artio_grid_allocation_strategy,
					num_file_vars, var_labels,
					num_levels_per_root_tree,
					num_octs_per_root_tree ) != ARTIO_SUCCESS ) {
				cart_error("Error adding grid to fileset %s", filename );
			}

			for (j=0;j<num_file_vars;j++) {
				cart_free(var_labels[j]);
			}

			cart_free( num_levels_per_root_tree );
			cart_free( num_octs_per_root_tree );

			write_artio_grid( handle, num_file_vars, var_indices );

			end_time( GAS_WRITE_IO_TIMER );
		}

#ifdef PARTICLES
		if ( fileset_write_options & WRITE_PARTICLES ) {
			start_time( PARTICLE_WRITE_IO_TIMER );

			define_particle_variables( &num_species, species_labels,
					num_primary_variables, primary_variable_labels,
					num_secondary_variables, secondary_variable_labels );

			/* redo particle linked list */
			root_tree_particle_list = cart_alloc( int, num_cells_per_level[min_level] );

			for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
				ipart = cell_particle_list[i];
				if ( ipart != NULL_PARTICLE ) {
					while (particle_list_next[ipart] != NULL_PARTICLE ) {
						ipart = particle_list_next[ipart];
					}
				}

				root_tree_particle_list[i] = ipart;
			}

			for ( level = min_level+1; level <= max_level; level++ ) {
				select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

				for ( i = 0; i < num_level_cells; i++ ) {
					icell = level_cells[i];

					if ( cell_particle_list[icell] != NULL_PARTICLE ) {
						iroot = cell_parent_root_cell(icell);
						ipart = cell_particle_list[icell];
						particle_list_prev[ipart] = root_tree_particle_list[iroot];

						/* find tail of list */
						while ( particle_list_next[ipart] != NULL_PARTICLE ) {
							ipart = particle_list_next[ipart];
						}

						root_tree_particle_list[iroot] = ipart;
					}
				}

				cart_free( level_cells );
			}

			num_particles_per_species_per_root_tree = cart_alloc( int, 
					num_cells_per_level[min_level]*num_particle_species );

			for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {		
				/* count particle species in root cell */
				for ( j = 0; j < num_particle_species; j++ ) {
					num_particles_per_species_per_root_tree[num_particle_species*i+j] = 0;
				}

				ipart = root_tree_particle_list[i];
				while ( ipart != NULL_PARTICLE ) {
					num_particles_per_species_per_root_tree[
						num_particle_species*i+particle_species(particle_id[ipart])]++;
					ipart = particle_list_prev[ipart];
				}

				total = 0;
				for ( j = 0; j < num_particle_species; j++ ) {
					total += num_particles_per_species_per_root_tree[num_particle_species*i+j];
				}

				/* sort linked list by species, id */
				order = cart_alloc( int, total );

				ipart = root_tree_particle_list[i];
				j = 0;
				while ( ipart != NULL_PARTICLE ) {
					order[j++] = ipart;
					ipart = particle_list_prev[ipart];
				}

				qsort( order, total, sizeof(int), compare_particle_species_id );

				if ( total > 0 ) {
					root_tree_particle_list[i] = order[0];
					for ( j = 0; j < total-1; j++ ) {
						particle_list_prev[order[j]] = order[j+1];
					}
					particle_list_prev[order[total-1]] = NULL_PARTICLE;
				}

				cart_free( order );	
			}

			if ( artio_fileset_add_particles( handle,
					num_artio_particle_files, 
					artio_particle_allocation_strategy,
					num_species, 
					species_labels,
					num_primary_variables,
					num_secondary_variables,
					primary_variable_labels,
					secondary_variable_labels,
					num_particles_per_species_per_root_tree ) != ARTIO_SUCCESS ) {
				cart_error("Error adding particles to fileset %s", filename );
			}

			for ( i = 0; i < num_species; i++ ) {
				for ( j = 0; j < num_primary_variables[i]; j++ ) {
					cart_free( primary_variable_labels[i][j] );
				}
				cart_free( primary_variable_labels[i] );

				for ( j = 0; j < num_secondary_variables[i]; j++ ) {
					cart_free( secondary_variable_labels[i][j] );
				}
				if ( num_secondary_variables[i] > 0 ) {
					cart_free( secondary_variable_labels[i] );
				}
			}

			artio_parameter_set_int_array( handle, "particle_species_num", 
					num_particle_species, particle_species_num );
			artio_parameter_set_float_array( handle, "particle_species_mass", 
					num_particle_species, particle_species_mass );

			write_artio_particles( handle, root_tree_particle_list, 
					num_particles_per_species_per_root_tree );
			cart_free( num_particles_per_species_per_root_tree );
			cart_free( root_tree_particle_list );

			/* return particle_list_prev to doublely-linked-list state */
			rebuild_particle_list();

			end_time( PARTICLE_WRITE_IO_TIMER );
		}
#endif /* PARTICLES */

		if ( artio_fileset_close(handle) != ARTIO_SUCCESS ) {
			cart_error("Error finalizing fileset %s, be careful, it may be corrupted!",
				filename );
		}
	}

#ifdef HYDRO_TRACERS
	if ( fileset_write_options & WRITE_TRACERS ) {
		start_time( PARTICLE_WRITE_IO_TIMER );
		write_cart_hydro_tracers( filename_tracers );                                                                                                  
		end_time( PARTICLE_WRITE_IO_TIMER );
	}
#endif /* HYDRO_TRACERS */

	if ( local_proc_id == MASTER_NODE &&
			fileset_write_options & WRITE_GRID 
#ifdef PARTICLES
			&& fileset_write_options & WRITE_PARTICLES
#endif /* PARTICLES */
#ifdef HYDRO_TRACERS 
			&& fileset_write_options & WRITE_TRACERS 
#endif /* HYDRO_TRACERS */
	) {
		/* write out restart file */
		sprintf( restart_filename, "%s/restart.dat", output_directory );
		restart = fopen( restart_filename, "w" );
		if ( restart == NULL ) {
			cart_error("Unable to open restart file %s for writing!", restart_filename );
		}
		fprintf( restart, "%s\n", filename );
#ifdef HYDRO_TRACERS
		fprintf( restart, "%s\n", filename_tracers );
#endif /* HYDRO_TRACERS */
        fclose(restart);
	}
}

void write_artio_grid( artio_file handle, int num_file_vars, int *var_indices ) {
	int i, j;
	int64_t sfc;
	int icell, ioct;
	int refined[num_children];
	int num_level_octs, num_next_level_octs;
	int *level_octs, *next_level_octs;
	int num_octs_per_level[max_level - min_level + 1];
	int level;
	float variables[num_vars * num_children];                                                          
	
	for ( sfc = proc_sfc_index[local_proc_id]; sfc < proc_sfc_index[local_proc_id+1]; sfc++) {
		icell = root_cell_location(sfc);
		pack_cell_vars( icell, num_file_vars, var_indices, variables );
		root_tree_level_oct_count( icell, num_octs_per_level );

		level = 1;
		while ( level <= max_level && num_octs_per_level[level-1] > 0 ) level++;

		if ( artio_grid_write_root_cell_begin(handle,sfc,variables, level-1, 
				num_octs_per_level) != ARTIO_SUCCESS ) {
			cart_error("Error writing grid root cell at sfc %ld", sfc );
		}

		level = min_level+1;

		if ( cell_is_refined(icell) ) {
			num_next_level_octs = 1;
			next_level_octs = cart_alloc( int, num_next_level_octs );
			next_level_octs[0] = cell_child_oct[icell];
		} else {
			num_next_level_octs = 0;
		}

		while ( num_next_level_octs > 0 ) {
			num_level_octs = num_next_level_octs;
			level_octs = next_level_octs;

			num_next_level_octs = 0;
			if ( level < max_level && num_octs_per_level[level] > 0 ) {
				next_level_octs = cart_alloc( int, num_children*num_level_octs );
			}

			if ( artio_grid_write_level_begin(handle,level) != ARTIO_SUCCESS ) {
				cart_error("Error writing level %u in sfc %ld", level, sfc );
			}

			for ( i = 0; i < num_level_octs; i++ ) {
				ioct = level_octs[i];
				for ( j = 0; j < num_children; j++ ) {
					icell = oct_child( ioct, j );
					pack_cell_vars( icell, num_file_vars, var_indices, &variables[num_file_vars*j] );
					if ( cell_is_refined(icell) ) {
						next_level_octs[num_next_level_octs++] = cell_child_oct[icell];
						refined[j] = 1;
					} else {
						refined[j] = 0;
					}
				}

				if ( artio_grid_write_oct(handle, variables, refined ) != ARTIO_SUCCESS ) {
					cart_error("Error writing oct %d on level %d, sfc = %ld", 
						ioct, level, sfc );
				}
			}

			cart_free( level_octs );
			artio_grid_write_level_end(handle);
			level++;
		}

		artio_grid_write_root_cell_end(handle);
	}
}

#ifdef PARTICLES
void write_artio_particles( artio_file handle, int *root_tree_particle_list, 
			int *num_particles_per_species_per_root_tree ) {
	int i, j, k, m;
	int sfc;
	int ipart;
	int64_t id;
	int *num_particles_per_species;
	double primary_variables[2*nDim+1];
	float secondary_variables[5];
	int subspecie = 0;

	for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
		sfc = root_cell_sfc_index(i);

		num_particles_per_species = &num_particles_per_species_per_root_tree[num_particle_species*i];

		if ( artio_particle_write_root_cell_begin(handle,sfc,
				num_particles_per_species) != ARTIO_SUCCESS ) {
			cart_error("Error writing particle root cell sfc %ld", sfc );
		}

		/* this code assumes two types of particles: N-body, Stars */
		ipart = root_tree_particle_list[i];
		for ( j = 0; j < num_particle_species; j++ ) {
			if ( artio_particle_write_species_begin( handle, j ) != ARTIO_SUCCESS ) {
				cart_error("Error writing particle species %d, sfc %ld", j, sfc );
			}

			for ( k = 0; k < num_particles_per_species[j]; k++ ) {
				id = particle_id[ipart];
				cart_assert( particle_species(particle_id[ipart]) == j );

				for ( m = 0; m < nDim; m++ ) {
					primary_variables[m] = particle_x[ipart][m];
				}
				for ( m = 0; m < nDim; m++ ) {
					primary_variables[nDim+m] = particle_v[ipart][m];
				}
				primary_variables[2*nDim] = particle_dt[ipart];
#ifdef STAR_FORMATION
				if ( j == num_particle_species - 1 ) {
					cart_assert( particle_is_star(ipart) );

#ifdef STAR_PARTICLE_TYPES
					subspecie = star_particle_type[ipart];
#endif /* STAR_PARTICLE_TYPES */

					secondary_variables[0] = star_tbirth[ipart];
					secondary_variables[1] = star_initial_mass[ipart];
					secondary_variables[2] = particle_mass[ipart];
#ifdef ENRICHMENT
					secondary_variables[3] = star_metallicity_II[ipart];
#ifdef ENRICHMENT_SNIa
					secondary_variables[4] = star_metallicity_Ia[ipart];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef STAR_PARTICLE_TYPES
				} else {
					subspecie = 0;
#endif /* STAR_PARTICLE_TYPES */
				}
#endif /* STAR_FORMATION */

				if ( artio_particle_write_particle( handle, id, subspecie, 
						primary_variables, secondary_variables ) != ARTIO_SUCCESS ) {
					cart_error("Error writing particle id %ld, species %d, sfc %ld",
						id, j, sfc );
				}
				ipart = particle_list_prev[ipart];
			}

			artio_particle_write_species_end(handle);
		}

		artio_particle_write_root_cell_end(handle);
	}
}
#endif /* PARTICLES */

void artio_restart_load_balance( artio_file handle ) {
	int i;
	int *constrained_quantities;
	float *cell_work;
	int page, file_max_level, num_oct_levels;
	int64_t sfc, end_sfc;
	float *variables;
	int *num_octs_per_level;
	int num_file_variables;
	int num_file_species;
	int *num_particles_per_species;
	int grid_flag = 0;
	int particle_flag = 0;
	
	if ( num_procs == 1 ) {
		proc_sfc_index[0] = 0;
		proc_sfc_index[1] = num_root_cells;
		init_tree();
		return;
	}

	if ( local_proc_id == MASTER_NODE ) {
		/* do load balancing */
		constrained_quantities = cart_alloc(int, num_constraints*num_root_cells );
		cell_work = cart_alloc(float, num_root_cells );

		for ( i = 0; i < num_root_cells; i++ ) {
			cell_work[i] = 0.0;
		}

		for ( i = 0; i < num_constraints*num_root_cells; i++ ) {
			constrained_quantities[i] = 0;
		}

		/* load grid information */
		if ( artio_parameter_get_int(handle, "num_grid_variables", &num_file_variables) == ARTIO_SUCCESS ) {
			grid_flag = 1;
			artio_parameter_get_int( handle, "max_refinement_level", &file_max_level );
			variables = cart_alloc(float, num_file_variables);
			num_octs_per_level = cart_alloc( int, file_max_level );

			for ( page = 0, sfc = 0; page < num_grid; page++ ) {
				end_sfc = MIN( sfc + num_grid*num_grid, num_root_cells ) - 1;
				if ( artio_grid_cache_sfc_range(handle, sfc, end_sfc ) != ARTIO_SUCCESS ) {
					cart_error("Error caching grid sfc range %ld to %ld", sfc, end_sfc );
				}

				for ( ; sfc <= end_sfc; sfc++ ) {
					if ( artio_grid_read_root_cell_begin(handle, sfc, variables, 
								&num_oct_levels, num_octs_per_level ) != ARTIO_SUCCESS ) {
						cart_error("Error reading grid root cell %ld", sfc );
					}

					/* contribution from root cell */
					constrained_quantities[num_constraints*sfc]++;
					cell_work[sfc] += cost_per_cell;

					for ( i = 0; i < num_oct_levels; i++ ) {
						constrained_quantities[num_constraints*sfc] += num_children*num_octs_per_level[i];
						cell_work[sfc] += cost_per_cell*(float)(2<<i)*num_children*num_octs_per_level[i];
					}

					artio_grid_read_root_cell_end(handle);
				}	
			}

			cart_free( num_octs_per_level );
			cart_free( variables );
		}

#ifdef PARTICLES
		if ( artio_parameter_get_int( handle, "num_particle_species", &num_file_species) == ARTIO_SUCCESS ) {
			particle_flag = 1;
			num_particles_per_species = cart_alloc(int, num_file_species);

			for ( page = 0, sfc = 0; page < num_grid; page++ ) {
				end_sfc = MIN( sfc + num_grid*num_grid, num_root_cells ) - 1;
				if ( artio_particle_cache_sfc_range(handle, sfc, end_sfc ) != ARTIO_SUCCESS ) {
					cart_error("Error caching particle sfc range %ld to %ld", sfc, end_sfc );
				}

				for ( ; sfc <= end_sfc; sfc++ ) {
					if ( artio_particle_read_root_cell_begin( handle, sfc, 
								num_particles_per_species ) != ARTIO_SUCCESS ) {
						cart_error("Error reading particle root cell %ld", sfc );
					}

					for ( i = 0; i < num_file_species; i++ ) {
						constrained_quantities[num_constraints*sfc+1] += num_particles_per_species[i];
						cell_work[sfc] += cost_per_particle*num_particles_per_species[i];
					}

					artio_particle_read_root_cell_end(handle);
				}   
			}

			cart_free( num_particles_per_species );
		}
#endif /* PARTICLES */

		/* ensure that some work has been estiamated */
		if ( !grid_flag && !particle_flag ) {
			cart_error("Unable to load balance based on an artio file without either grid or particle data!");
		}

		load_balance_entire_volume( cell_work, constrained_quantities, proc_sfc_index );

		cart_free( cell_work );
		cart_free( constrained_quantities );
	}

	/* let all other processors know what their new workload is */
	MPI_Bcast( proc_sfc_index, num_procs+1, MPI_INT, MASTER_NODE, mpi.comm.run );
	init_tree();
}

void read_artio_restart( const char *label ) {
	int i;
	artio_file handle;
	int sfc_order;
	int64_t num_file_root_cells;
	int level;
	int file_max_level;
	int num_levels;
	char filename[256];
	int64_t mpi_task_sfc_index[MAX_PROCS+1];
	int num_grid_files, num_particle_files;
	int num_file_procs, num_file_octs, num_file_particles, num_file_star_particles;
	FILE *restart;
	char str[CONTROL_PARAMETER_STRING_LENGTH];
#ifdef RADIATIVE_TRANSFER
	frt_intg n;
	frt_real *data;
	float *buffer;
#endif /* RADIATIVE_TRANSFER */
	struct artio_context_struct con = { mpi.comm.run };

#ifdef COSMOLOGY
	double OmM0, OmB0, OmL0, h100, DelDC;
#else
	double mass_unit, time_unit, length_unit;
#endif

	if ( label == NULL ) {
		/* try to load from restart.dat */
		sprintf( filename, "%s/restart.dat", output_directory );
		restart = fopen( filename, "r" );

		if ( restart == NULL ) {
			cart_debug("Unable to locate restart.dat, trying default filename!");
			sprintf( filename, "%s/%s", output_directory, jobname );
#ifdef HYDRO_TRACERS
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
			sprintf( filename_tracers, "%s/%s.dtr", output_directory, jobname );
#else
			sprintf( filename_tracers, "%s/tracers.dat", output_directory );
#endif
#endif /* HYDRO_TRACERS */
		} else {
			fscanf( restart, "%s\n", filename );
#ifdef HYDRO_TRACERS
			fscanf( restart, "%s\n", filename_tracers );
#endif /* HYDRO_TRACERS */
			fclose(restart);
		}
	} else {
		sprintf( filename, "%s/%s_%s", output_directory, jobname, label );
#ifdef HYDRO_TRACERS
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
        sprintf( filename_tracers, "%s/%s_%s.dtr", output_directory, jobname, label );
#else
        sprintf( filename_tracers, "%s/tracers_%s.dat", output_directory, label );
#endif
#endif /* HYDRO_TRACERS */
	}	

	handle = artio_fileset_open(filename, ARTIO_OPEN_HEADER, &con);
	if( handle == NULL ) {
		cart_error("Unable to open ARTIO fileset %s",filename);
	}

	if ( artio_parameter_get_long(handle, "num_root_cells", &num_file_root_cells) != ARTIO_SUCCESS ) {
		cart_error("Unable to load number of root cells from artio header file!");
	}

	if (num_file_root_cells != num_root_cells) {
		cart_error( "Number of root cells in file %s header does not match compiled value: %d vs %d",
				filename, num_file_root_cells, num_root_cells);
	}

    if ( artio_parameter_get_int( handle, "num_grid_files", &num_grid_files ) != ARTIO_SUCCESS ||
			artio_fileset_open_grid(handle) != ARTIO_SUCCESS ) {
        cart_debug("Warning: the fileset does not contain grid data, the code cannot be restarted from a complete snapshot");
    }

#ifdef PARTICLES
    if ( artio_parameter_get_int( handle, "num_particle_files", &num_particle_files ) != ARTIO_SUCCESS ||
			artio_fileset_open_particles(handle) != ARTIO_SUCCESS ) {
        cart_debug("Warning: fileset does not contain particles, the code cannot be restarted from a complete snapshot");
    }
#endif /* PARTICLES */

	/* try to load balance */
	if ( artio_parameter_get_array_length( handle, "mpi_task_sfc_index", &num_file_procs ) != ARTIO_SUCCESS
			|| num_file_procs != num_procs+1 
			|| artio_parameter_get_int( handle, "num_octs_per_mpi_task", &num_file_octs ) != ARTIO_SUCCESS 
			|| num_file_octs > num_octs 
#ifdef PARTICLES
			|| artio_parameter_get_int( handle, "num_particles_per_mpi_task", &num_file_particles ) != ARTIO_SUCCESS 
			|| num_file_particles > num_particles
#ifdef STAR_FORMATION
			|| artio_parameter_get_int( handle, "num_star_particles_per_mpi_task", &num_file_star_particles ) != ARTIO_SUCCESS 
			|| num_file_star_particles > num_star_particles 
#endif /* STAR_FORMATION */
#endif /* PARTICLES */
	) {
		artio_restart_load_balance( handle );
	} else {
		artio_parameter_get_long_array( handle, "mpi_task_sfc_index", num_procs+1, mpi_task_sfc_index );
		for ( i = 0; i < num_procs+1; i++ ) {
			proc_sfc_index[i] = mpi_task_sfc_index[i];
		}
		init_tree();
	}

	/* load all simulation parameters here */
	artio_parameter_get_string( handle, "jobname", str, CONTROL_PARAMETER_STRING_LENGTH);
	set_jobname( str );
	cart_debug("jobname: %s", jobname );

	if ( artio_parameter_get_int( handle, "sfc", &sfc_order ) != ARTIO_SUCCESS || sfc_order != SFC ) {
		cart_error("Grid fileset has undefined or different sfc indexing than compiled code!");
	}

	/* check maximum level */
	if ( artio_parameter_get_int( handle, "max_refinement_level", &file_max_level ) != ARTIO_SUCCESS ||
			file_max_level > max_level ) {
		cart_error("Grid fileset has undefined refinement levels or contains more levels than compiled code!");
	}
	num_levels = file_max_level+1;
    
#ifdef HYDRO    
	artio_parameter_get_int_array( handle, "hydro_sweep_direction", num_levels, level_sweep_dir );
	for ( level = file_max_level+1; level < max_level; level++ ) {
		level_sweep_dir[level] = 0;
	}
#endif
    
	/* unit parameters */
	if ( artio_parameter_get_double( handle, "box_size", &box_size ) != ARTIO_SUCCESS ) {
		cart_error("Unable to read box size from artio header file!");
	}

#ifdef COSMOLOGY 
	if ( artio_parameter_get_double( handle, "auni_init", &auni_init ) != ARTIO_SUCCESS ||
			artio_parameter_get_double( handle, "OmegaM", &OmM0 ) != ARTIO_SUCCESS ||
			artio_parameter_get_double( handle, "OmegaL", &OmL0 ) != ARTIO_SUCCESS ||
			artio_parameter_get_double( handle, "OmegaB", &OmB0 ) != ARTIO_SUCCESS ||
			artio_parameter_get_double( handle, "hubble", &h100 ) != ARTIO_SUCCESS ||
			artio_parameter_get_double( handle, "DeltaDC", &DelDC ) ) {
		cart_error("Unable to read cosmology information from artio header file!");
	}

	cosmology_set(OmegaM,OmM0);
	cosmology_set(OmegaL,OmL0);
	cosmology_set(OmegaB,OmB0);
	cosmology_set(h,h100);
	cosmology_set(DeltaDC,DelDC);
		
#else 
	if ( artio_parameter_get_double( handle, "mass_unit", &mass_unit ) != ARTIO_SUCCESS ||
			artio_parameter_get_double( handle, "time_unit", &time_unit ) != ARTIO_SUCCESS ||
			artio_parameter_get_double( handle, "length_unit", &length_unit ) != ARTIO_SUCCESS ) {
		cart_error("Unable to read unit information from artio header file!");
	}

	units_set( mass_unit, time_unit, length_unit );
#endif /* COSMOLOGY */

	/* timestepping variables */
	artio_parameter_get_int( handle, "step", &step );
	artio_parameter_get_double_array( handle, "tl", num_levels, tl );
	artio_parameter_get_double_array( handle, "tl_old", num_levels, tl_old );
	artio_parameter_get_double_array( handle, "dtl", num_levels, dtl );
	artio_parameter_get_double_array( handle, "dtl_old", num_levels, dtl_old );

	artio_parameter_get_int_array( handle, "time_refinement_factor", num_levels, time_refinement_factor );
	artio_parameter_get_int_array( handle, "time_refinement_factor_old", num_levels, time_refinement_factor_old );                                              

#ifdef COSMOLOGY
	artio_parameter_get_double( handle, "auni_init", &auni_init );
	artio_parameter_get_double_array( handle, "abox", num_levels, abox );
	artio_parameter_get_double_array( handle, "auni", num_levels, auni );

	for ( level = min_level; level <= file_max_level; level++ ) {
		abox_old[level] = abox_from_tcode(tl_old[level]);
	}
#endif /* COSMOLOGY */

	for ( level = file_max_level+1; level < max_level; level++ ) {
		tl[level] = tl[file_max_level];
		tl_old[level] = tl_old[file_max_level];
		dtl[level] = dtl[level-1];
		dtl_old[level] = dtl[level-1];

#ifdef COSMOLOGY
		abox[level] = abox[file_max_level];
		auni[level] = abox[file_max_level];
		abox_old[level] = abox_old[file_max_level];
#endif /* COSMOLOGY */
	}

#ifdef PARTICLES
	/* energy conservation variables */ 
	artio_parameter_get_double( handle, "energy:tintg", &tintg );
	artio_parameter_get_double( handle, "energy:ekin", &ekin );
	artio_parameter_get_double( handle, "energy:ekin1", &ekin1 );
	artio_parameter_get_double( handle, "energy:ekin2", &ekin2 );
	artio_parameter_get_double( handle, "energy:au0", &au0 );
	artio_parameter_get_double( handle, "energy:aeu0", &aeu0 );
	artio_parameter_get_double( handle, "energy:ap0", &ap0 );
#endif /* PARTICLES */

    /* refinement boundaries */
	artio_parameter_get_float_array( handle, "refinement_volume_min", 
			nDim, refinement_volume_min );
	artio_parameter_get_float_array( handle, "refinement_volume_max", 
			nDim, refinement_volume_max );

#ifdef STAR_FORMATION
	artio_parameter_get_float_array( handle, "star_formation_volume_min", 
			nDim, star_formation_volume_min );
	artio_parameter_get_float_array( handle, "star_formation_volume_max", 
			nDim, star_formation_volume_max );

	artio_parameter_get_double( handle, "total_stellar_mass", &total_stellar_mass );
	artio_parameter_get_double( handle, "total_stellar_initial_mass", &total_stellar_initial_mass );
#endif /* STAR_FORMATION */

#ifdef RADIATIVE_TRANSFER
	n = 0;
	frtCall(unpackradiationfields)(&n,0);
	if(n < 1)
	  {
	    cart_error("Unable to unpack Radiation Field data.");
	  }

	data = cart_alloc(frt_real, n );

	if(sizeof(float) == sizeof(frt_real))
	  {
	    buffer = (float *)data;
	  }
	else
	  {
	    buffer = cart_alloc(float,n);
	  }

	artio_parameter_get_float_array( handle, "radiation_background", n, buffer );

	if(sizeof(float) == sizeof(frt_real))
	  {
	  }
	else
	  {
	    for(i=0; i<n; i++) data[i] = buffer[i];
	    cart_free(buffer);
	  }

	frtCall(unpackradiationfields)(&n,data);
	cart_free(data);

#ifdef RT_SINGLE_SOURCE
	artio_parameter_get_int( handle, "rt:single_source:level", &rtSingleSourceLevel );
	artio_parameter_get_float( handle, "rt:single_source:value", &rtSingleSourceValue );
	artio_parameter_get_double_array( handle, "rt:single_source:pos", nDim, rtSingleSourcePos );
#endif
#endif /* RADIATIVE_TRANSFER */

	if ( artio_parameter_get_int( handle, "num_grid_files", &num_grid_files ) == ARTIO_SUCCESS ) {
		read_artio_grid(handle, file_max_level);
		cart_debug("done reading grid");
	}

#ifdef PARTICLES
	if ( artio_parameter_get_int( handle, "num_particle_files", &num_particle_files ) == ARTIO_SUCCESS ) {
		read_artio_particles(handle);
		cart_debug("num_local_particles = %u", num_local_particles );
		build_particle_list();
	}
#endif /* PARTICLES */

	artio_fileset_close(handle);

#ifdef HYDRO_TRACERS
	start_time( PARTICLE_READ_IO_TIMER );
	read_cart_hydro_tracers( filename_tracers );
	end_time( PARTICLE_READ_IO_TIMER );
#endif /* HYDRO_TRACERS */

}

void read_artio_grid( artio_file handle, int file_max_level ) {
	int i, j;
	int64_t sfc;
	int icell, ioct;
	int *oct_order, *next_level_order;
	int next_level_octs;

	int num_sim_variables, num_file_variables;
	int file_var_indices[num_vars];
	int sim_var_indices[num_vars];
	char *sim_var_labels[num_vars];

	float * root_variables;
	int num_tree_levels;
	int *num_octs_per_level;
	int level;
	float * oct_variables;
	int oct_refined[num_children];
	char ** file_variables;

	start_time( GAS_READ_IO_TIMER );

	if ( artio_parameter_get_int(handle, "num_grid_variables", &num_file_variables) != ARTIO_SUCCESS ) {
		cart_error("Unable to locate grid data in artio fileset");
	}

	/* load list of variables the code expects */
	define_file_variables(&num_sim_variables, sim_var_labels, sim_var_indices);

	if (num_file_variables < num_sim_variables) {
		cart_error("Not enough variables in file header!" );
	} else if (num_file_variables > num_sim_variables) {
		cart_debug(
				"WARNING: file %s contains more variables than code expects (%d vs %d), hope you know what you're doing...",
				num_file_variables, num_sim_variables);
	}

	root_variables = cart_alloc(float, num_file_variables);
	oct_variables = cart_alloc(float, 8 * num_file_variables);
	file_variables = cart_alloc(char *, num_file_variables);
	for(i=0; i<num_file_variables; i++) file_variables[i] = cart_alloc(char, 256);

	if ( artio_parameter_get_string_array(handle, "grid_variable_labels", num_file_variables, file_variables, 256 ) != ARTIO_SUCCESS ) {
		cart_error("Unable to read grid variable labels from artio header file, which may be corrupt.");
	}

	/* match expected variables with variables in file */
	for (i = 0; i < num_sim_variables; i++) {
		for (j = 0; j < num_file_variables; j++) {
			if (!strcmp(sim_var_labels[i], file_variables[j])) {
				file_var_indices[i] = j;
				break;
			}
		}

		if (j == num_file_variables) {
			cart_error("Unable to locate expected variable %s in header",sim_var_labels[i]);
		}

		cart_free(sim_var_labels[i]);
	}

	for(i=0; i<num_file_variables; i++) cart_free( file_variables[i] );
	cart_free( file_variables );

	num_octs_per_level = cart_alloc( int, file_max_level );

	/* cache file offsets */
	if ( artio_grid_cache_sfc_range( handle, proc_sfc_index[local_proc_id], 
				proc_sfc_index[local_proc_id+1]-1 ) != ARTIO_SUCCESS ) {
		cart_error("Error caching grid sfc range %ld to %ld", 
				proc_sfc_index[local_proc_id],
				proc_sfc_index[local_proc_id+1]-1 );
	}

	/* load each space-filling-curve index in turn */
	for (sfc = proc_sfc_index[local_proc_id]; 
			sfc < proc_sfc_index[local_proc_id + 1]; sfc++) {

		if ( artio_grid_read_root_cell_begin(handle, sfc, root_variables,
					&num_tree_levels, num_octs_per_level) != ARTIO_SUCCESS ) {
			cart_error( "Error reading grid root cell sfc %ld", sfc );
		}
	
		icell = root_cell_location(sfc);
		unpack_cell_vars(icell, num_sim_variables, sim_var_indices,
				file_var_indices, root_variables);
		
		if (num_tree_levels > 0) {
			if ( split_cell(icell) ) {
				cart_error("Unable to split root cell, sfc = %ld", sfc );
			}

			next_level_order = cart_alloc(int, 1);
			next_level_order[0] = cell_child_oct[icell];
			next_level_octs = 1;
		}

		for (level = 0; level < num_tree_levels; level++) {
			if ( artio_grid_read_level_begin(handle, level+1) != ARTIO_SUCCESS ) {
				cart_error("Error reading level %d, sfc %ld", level+1, sfc );
			}

			oct_order = next_level_order;
			cart_assert(num_octs_per_level[level] == next_level_octs);

			if (level < num_tree_levels-1) {
				next_level_order = cart_alloc(int, num_octs_per_level[level+1] );
			} else {
				next_level_order = NULL;
			}
			next_level_octs = 0;

			for (i = 0; i < num_octs_per_level[level]; i++) {
				ioct = oct_order[i];
				if ( artio_grid_read_oct(handle, oct_variables, oct_refined) != ARTIO_SUCCESS ) {
					cart_error("Error reading oct %d on level %d, sfc %ld", ioct, 
						level+1, sfc );
				}

				for (j = 0; j < num_children; j++) {
					icell = oct_child(ioct, j);

					unpack_cell_vars(icell, num_file_variables, sim_var_indices, 
							file_var_indices, &oct_variables[num_file_variables*j]);

					if (oct_refined[j]) {
						cart_assert( level < num_tree_levels );
						cart_assert( next_level_octs < num_octs_per_level[level+1] );

						if ( split_cell(icell) ) {
							cart_error("Unable to split cell, sfc = %ld, level = %d, icell = %d",
								sfc, level+1, icell );
						}

						next_level_order[next_level_octs++]
							= cell_child_oct[icell];
					}
				}
			}

			cart_free(oct_order);
			artio_grid_read_level_end(handle);
		}

		artio_grid_read_root_cell_end(handle);
	}

	cart_free( num_octs_per_level );
	cart_free(root_variables);
	cart_free (oct_variables);

	end_time( GAS_READ_IO_TIMER );
}

#ifdef PARTICLES
void read_artio_particles( artio_file handle ) {
	int i, j;
	int sfc;
	int64_t pid;
	int id, ipart;
	int species;
	int subspecies;
	long num_particles_local = 0L;
	int num_species;

	double primary_variables[7];
	float secondary_variables[5];

	int num_particles_per_species[MAX_PARTICLE_SPECIES];

	start_time( PARTICLE_READ_IO_TIMER );

	artio_parameter_get_int( handle, "num_particle_species", &num_species);

	if ( num_species > MAX_PARTICLE_SPECIES ) {
		cart_error("Ran out of particle species!");
	}

	num_particle_species = num_species;

	if ( artio_parameter_get_int_array( handle, "particle_species_num", num_species, particle_species_num ) != ARTIO_SUCCESS ||
			artio_parameter_get_float_array( handle, "particle_species_mass", num_species, particle_species_mass ) != ARTIO_SUCCESS ) {
		cart_error("Unable to load particle species data from artio header file!");
	}

	for ( i = 0; i < num_species; i++ ) {
		cart_debug("particle species %d: %u particles, %e mass", i, particle_species_num[i], particle_species_mass[i] );
	}

	particle_species_indices[0] = 0;                                                                                     
	for ( i = 0; i < num_particle_species; i++ ) {
		particle_species_indices[i+1] = particle_species_indices[i] + particle_species_num[i];
	}

	/* cache file offsets */
	if ( artio_particle_cache_sfc_range( handle, proc_sfc_index[local_proc_id], 
				proc_sfc_index[local_proc_id+1]-1 ) != ARTIO_SUCCESS ) {
		cart_error("Error caching particle sfc range %ld to %ld",
			proc_sfc_index[local_proc_id],
			proc_sfc_index[local_proc_id+1]-1 );
	}

	/* load each space-filling-curve index in turn */
	for (sfc = proc_sfc_index[local_proc_id]; 
			sfc < proc_sfc_index[local_proc_id + 1]; sfc++) {

		if ( artio_particle_read_root_cell_begin(handle, sfc, num_particles_per_species ) != 
				ARTIO_SUCCESS ) {
			cart_error("Error reading particle root cell %ld", sfc );
		}

		for ( species = 0; species < num_species; species++ ) {
			if ( artio_particle_read_species_begin(handle, species) != ARTIO_SUCCESS ) {
				cart_error("Error reading particle specie %d, sfc = %ld", species, sfc );
			}

			num_particles_local += num_particles_per_species[species];

			for ( i = 0; i < num_particles_per_species[species]; i++ ) {
				if ( artio_particle_read_particle(handle, &pid, &subspecies,
							primary_variables, secondary_variables ) != ARTIO_SUCCESS ) {
					cart_error("Error reading particle %d, species %d, sfc %ld",
						i, species, sfc );
				}

				/* unpack variables */
				if ( pid >= INT_MAX ) {
					cart_error("Error storing particle id %ld in 32-bit integer!", pid );
				}
				id = (int)pid;
				ipart = particle_alloc( id );

				for ( j = 0; j < nDim; j++ ) {
					particle_x[ipart][j] = primary_variables[j];
				}

				for ( j = 0; j < nDim; j++ ) {
					particle_v[ipart][j] = primary_variables[nDim+j];
				}

				particle_dt[ipart] = primary_variables[2*nDim];
				particle_t[ipart] = tl[min_level];

#ifdef STAR_FORMATION
				if ( species == num_species - 1 ) {
#ifdef STAR_PARTICLE_TYPES
					star_particle_type[ipart] = subspecies;
#endif /* STAR_PARTICLE_TYPES */

					star_tbirth[ipart] = secondary_variables[0];
					star_initial_mass[ipart] = secondary_variables[1];
					particle_mass[ipart] = secondary_variables[2];

#ifdef ENRICHMENT
					star_metallicity_II[ipart] = secondary_variables[3];
#ifdef ENRICHMENT_SNIa
					star_metallicity_Ia[ipart] = secondary_variables[4];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
				} else {
					particle_mass[ipart] = particle_species_mass[particle_species(id)];
				}
#else
				particle_mass[ipart] = particle_species_mass[particle_species(id)];
#endif /* STAR_FORMATION */
			}

			artio_particle_read_species_end(handle);
		}

		artio_particle_read_root_cell_end(handle);
	}

/*  DHR - due to deleted particles, this code is incorrect, i.e. the sum of the local active particles
    does not equal num_particles_total
	MPI_Allreduce( &num_particles_local, &num_particles_total, 1, MPI_LONG, MPI_SUM, mpi.comm.run );
*/
	num_particles_total = particle_species_indices[num_particle_species];

	end_time( PARTICLE_READ_IO_TIMER );
}
#endif /* PARTICLES */
