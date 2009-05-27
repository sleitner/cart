#include "defs.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "units.h"
#include "timing.h"
#include "timestep.h"
#include "io.h"
#include "load_balance.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "auxiliary.h"
#include "starformation.h"
#include "hydro.h"
#include "parallel.h"


void read_config( char *filename ) {
	FILE *configfile;
	char line[1024];
	char *tag, *value, *p;
	char *output;
	int i;
	int indicator, refinement_level_min, refinement_level_max;

	/* set number = 0 here so we can allow for multiple lines */
	num_outputs = 0;

	configfile = fopen( filename, "r" );
	if ( configfile == NULL ) {
		cart_error("Unable to load configuration file %s", filename );
	}

	while ( fgets( line, 1024, configfile ) ) {
		/* check for different states */
		if ( line[0] != '#' ) {
			/* translate all tabs/carriage return to spaces */
			for ( i = 0; i < strlen(line); i++ ) {
				if ( line[i] == '\t' || line[i] == '\n' ) {
					line[i] = ' ';
				}
			}

			/* skip over leading spaces */
			p = line;
			while ( *p == ' ' ) {
				*p++;
			}

			tag = p;

			/* convert tag to lower case */
			while ( *p != ' ' ) {
				*p = (char)tolower( (int)*p );
				*p++;
			}

			/* null-terminate tag string */
			*p = '\0';
			*p++;

			/* skip spaces between tag and output */
			while ( *p == ' ' ) {
				*p++;
			}

			value = p;

			/* remove trailing spaces */
			i = strlen(value) - 1;
			while ( value[i] == ' ' ) {
				value[i] = '\0';
				i--;
			}

			if ( strcmp( tag, "omega0" ) == 0 || strcmp( tag, "omegam" ) == 0 ) {
			        cosmology_set(OmegaM,atof( value ));
			} else if ( strcmp( tag, "omegal0" ) == 0 || strcmp( tag, "omegal" ) == 0 ) {
			        cosmology_set(OmegaL,atof( value ));
			} else if ( strcmp( tag, "omegab0" ) == 0 || strcmp( tag, "omegab" ) == 0 ) {
			        cosmology_set(OmegaB,atof( value ));
			} else if ( strcmp( tag, "hubble" ) == 0 || strcmp( tag, "h" ) == 0 ) {
			        cosmology_set(h,atof( value ));
			} else if ( strcmp( tag, "deltadc" ) == 0 ) {
			        cosmology_set(DeltaDC,atof( value ));
			} else if ( strcmp( tag, "lbox" ) == 0 ) {
				Lbox = atof( value );
			} else if ( strcmp( tag, "output_directory" ) == 0 ) {
				strcpy( output_directory, value );
			} else if ( strcmp( tag, "logfile_directory" ) == 0 ) {
				strcpy( logfile_directory, value );
			} else if ( strcmp( tag, "jobname" ) == 0 ) {
				strcpy( jobname, value );
			} else if ( strcmp( tag, "requeue_command" ) == 0 ) {
				strcpy( requeue_command, value );
			} else if( strcmp( tag, "stopfile" ) == 0 ) {
				
			} else if ( strcmp( tag, "a_init" ) == 0 || strcmp( tag, "auni_init" ) == 0 ) {
				auni_init = atof( value );
				#ifdef COSMOLOGY
				t_init = tcode_from_auni( auni_init );
				#endif
			} else if ( strcmp( tag, "a_end" ) == 0 || strcmp( tag, "auni_end" ) == 0 ) {
				auni_end = atof( value );
				#ifdef COSMOLOGY
				t_end = tcode_from_auni( auni_end );
				#endif
			} else if ( strcmp( tag, "t_init" ) == 0 ) {
				t_init = atof( value );
				#ifdef COSMOLOGY
					auni_init = auni_from_tcode(t_init);
				#else
					auni_init = 0.0;
				#endif
			} else if ( strcmp( tag, "t_end" ) == 0 ) {
				t_end = atof( value );
				#ifdef COSMOLOGY
					auni_end = auni_from_tcode(t_end);
				#else
					auni_end = 1.0e35;
				#endif
			} else if ( strcmp( tag, "timelimit" ) == 0 ) {
				timelimit = atof(value);
			} else if ( strcmp( tag, "max_steps" ) == 0 ) {
				max_steps = atoi(value);
			} else if ( strcmp( tag, "output_frequency" ) == 0 ) {
				output_frequency = atoi(value);
			} else if ( strcmp( tag, "restart_frequency" ) == 0 ) {
				restart_frequency = atoi(value);
			} else if ( strcmp( tag, "load_balance_frequency" ) == 0 ) {
				load_balance_frequency = atoi(value);
			} else if ( strcmp( tag, "particle_output_frequency" ) == 0 ) {
				particle_output_frequency = atoi(value);
			} else if ( strcmp( tag, "tracer_output_frequency" ) == 0 ) {
				tracer_output_frequency = atoi(value);
			} else if ( strcmp( tag, "grid_output_frequency" ) == 0 ) {
				grid_output_frequency = atoi(value);
			} else if ( strcmp( tag, "num_output_files" ) == 0 ) {
				num_output_files = atoi(value);
			} else if ( strcmp( tag, "cost_per_cell" ) == 0 ) {
				cost_per_cell = atof(value);
			} else if ( strcmp( tag, "cost_per_particle" ) == 0 ) {
				cost_per_particle = atof(value);
			} else if ( strcmp( tag, "mpi_customization_mode" ) == 0 ) {
				mpi_customization_mode = atoi(value);
#ifdef _OPENMP
			} else if ( strcmp( tag, "omp_num_threads" ) == 0 ) {
			        omp_set_num_threads(atoi(value));
#endif
			} else if ( strcmp( tag, "outputs" ) == 0 ) {
				output = strtok( value, " " );
				while ( output != NULL && num_outputs < MAX_OUTPUTS ) {
					outputs[num_outputs++] = atof( output );
					output = strtok( NULL, " " );
				}
			} else if ( strcmp( tag, "cfl" ) == 0 ) {
				cfl_run = cfl_max = atof( value );
			} else if ( strcmp( tag, "cfl_run" ) == 0 ) {
				cfl_run = atof( value );
				if(cfl_max < cfl_run) cfl_max = cfl_run;
			} else if ( strcmp( tag, "cfl_max" ) == 0 ) {
				cfl_max = atof( value );
				if(cfl_max < cfl_run) cfl_max = cfl_run;
			} else if ( strcmp( tag, "max_cfl_sync_level" ) == 0 ) {
				max_cfl_sync_level = atoi( value );
			} else if ( strcmp( tag, "particle_cfl" ) == 0 ) {
				particle_cfl = atof( value );
			} else if ( strcmp( tag, "max_time_inc" ) == 0 ) {
				max_time_inc = atof( value );
			} else if ( strcmp( tag, "min_time_dec" ) == 0 ) {
				min_time_dec = atof( value );
			} else if ( strcmp( tag, "max_da" ) == 0 ) {
				max_da = atof( value );
			} else if ( strcmp( tag, "max_frac_da" ) == 0 ) {
				max_frac_da = atof( value );
			} else if ( strcmp( tag, "max_dt" ) == 0 ) {
				max_dt = atof( value );
			} else if ( strcmp( tag, "refinement" ) == 0 ) {
				/* which indicator? */
				output = strtok( value, " " );
				indicator = atoi( output );

				cart_assert( indicator >= 0 && indicator < num_refinement_indicators );

				/* weight */
				output = strtok( NULL, " " );
				refinement_indicator_weight[indicator] = atof(output);

				/* levels */
				output = strtok( NULL, " " );
				refinement_level_min = max( min_level, atoi( output ) );

				output = strtok( NULL, " " );
				refinement_level_max = min( max_level-1, atoi( output ) );

				for ( i = refinement_level_min; i <= refinement_level_max; i++ ) {
					use_refinement_indicator[indicator][i] = 1;

					output = strtok( NULL, " " );
					refinement_indicator_threshold[indicator][i] = atof(output);
				}

			} else if ( strcmp( tag, "split_tolerance" ) == 0 ) {
				split_tolerance = atof( value );
			} else if ( strcmp( tag, "wsplit" ) == 0 ) {
				split_tolerance = atof(value);
			} else if ( strcmp( tag, "join_tolerance" ) == 0 ) {
				join_tolerance = atof( value );
			} else if ( strcmp( tag, "wjoin" ) == 0 ) {
				join_tolerance = atof( value );
			} else if ( strcmp( tag, "num_diffusion_steps" ) == 0 ) {
				num_diffusion_steps = atoi( value );
			} else if ( strcmp( tag, "reaction_increment" ) == 0 ) {
				reaction_increment = atof( value );
			} else if ( strcmp( tag, "diffusion_coefficient" ) == 0 ) {
				diffusion_coefficient = atof( value );
			} else if ( strcmp( tag, "momentum_increment" ) == 0 ) {
				momentum_increment = atof( value );
#ifdef HYDRO
			} else if ( strcmp( tag, "gas_temperature_floor" ) == 0 ) {
				gas_temperature_floor = atof( value );
			} else if ( strcmp( tag, "gas_density_floor" ) == 0 ) {
				gas_density_floor = atof( value );
#ifdef PRESSURE_FLOOR
			} else if ( strcmp( tag, "pressure_floor_min_level" ) == 0 ) {
				pressure_floor_min_level = atoi( value );
#endif  /* PRESSURE_FLOOR */
#endif  /* HYDRO */
#ifdef STARFORM
			} else if ( strcmp( tag, "sf_recipe" ) == 0 ) {
				sf_recipe = atoi( value );
			} else if ( strcmp( tag, "alpha_sf" ) == 0 ) {
				alpha_SF = atof( value );
			} else if ( strcmp( tag, "eps_sf" ) == 0 ) {
				eps_SF = atof( value );
			} else if ( strcmp( tag, "dtmin_sf" ) == 0 ) {
				dtmin_SF = atof( value );
			} else if ( strcmp( tag, "tau_sf" ) == 0 ) {
				tau_SF = atof( value );
			} else if ( strcmp( tag, "dm_star_min" ) == 0 ) {
				dm_star_min = atof( value );
			} else if ( strcmp( tag, "rho_sf" ) == 0 ) {
				rho_SF = atof( value );
			} else if ( strcmp( tag, "t_sf" ) == 0 ) {
				T_SF = atof( value );
			} else if ( strcmp( tag, "a_imf" ) == 0 ) {
				a_IMF = atof( value );
			} else if ( strcmp( tag, "am_stl" ) == 0 ) {
				aM_stl = atof( value );
			} else if ( strcmp( tag, "am_stu" ) == 0 ) {
				aM_stu = atof( value );
			} else if ( strcmp( tag, "am_snii" ) == 0 ) {
				aM_SNII = atof( value );
			} else if ( strcmp( tag, "am_snia1" ) == 0 ) {
				aM_SNIa1 = atof( value );
			} else if ( strcmp( tag, "am_snia2" ) == 0 ) {
				aM_SNIa2 = atof( value );
			} else if ( strcmp( tag, "ejm_snia" ) == 0 ) {
				ejM_SNIa = atof( value );
			} else if ( strcmp( tag, "c_snia" ) == 0 ) {
				C_SNIa = atof( value );
			} else if ( strcmp( tag, "t_snia" ) == 0 ) {
				t_SNIa = atof( value );
			} else if ( strcmp( tag, "e_51" ) == 0 ) {
				E_51 = atof( value );
			} else if ( strcmp( tag, "t_fb" ) == 0 ) {
				t_fb = atof( value );
			} else if ( strcmp( tag, "t0_ml" ) == 0 ) {
				T0_ml = atof( value );
			} else if ( strcmp( tag, "c0_ml" ) == 0 ) {
				c0_ml = atof( value );
			} else if ( strcmp( tag, "T_max_feedback" ) == 0 ) {
				T_max_feedback = atof( value );
#endif /* STARFORM */
			} else if ( strcmp( tag, "" ) == 0 ) {
				/* no tag, just ignore */
			} else {
				/* unknown tag */
				cart_debug("Unrecognized tag [%s] in parameter file %s", tag, filename );
			}
		}
	}

	fclose( configfile );
}
