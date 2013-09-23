#include "config.h"
#ifdef REFINEMENT

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "hydro.h"
#include "tree.h"
#include "hydro.h"


int join( int cell ) {
	int i, j;
	int child, neighbor;

	cart_assert( cell_is_refined(cell) );
	
	for ( i = 0; i < num_children; i++ ) {
		child = cell_child( cell, i );

		if ( cell_is_refined(child) ) {
			return -2;
		}

		for ( j = 0; j < nDim; j++ ) {
			neighbor = cell_neighbor( child, external_direction[i][j] );

			if ( neighbor != NULL_OCT && cell_is_refined(neighbor) ) {
				return -2;
			}
		}
	}
	
	return join_cell(cell);
}

#ifdef HYDRO
double cell_internal_energy( int icell ) {
	return cell_gas_energy(icell) - cell_gas_kinetic_energy(icell);
}
#endif /* HYDRO */

int split ( int cell ) {
	int i, j, ind;
	int neighbor, cell_number;
	int result;
	int child_cell;
	int neighbors[nDim];
	float mass;

	double moment;	
#ifdef HYDRO
	double weights[num_hydro_vars-nDim];
#endif /* HYDRO */

	cart_assert( cell >= 0 && cell < num_cells );

#ifdef HYDRO
	for ( i = 0; i < num_hydro_vars-nDim; i++ ) {
		weights[i] = 0.0;
	}
#endif /* HYDRO */

	/* check +/- refinement criterion */
	if ( cell_level(cell) > min_level ) {
		cell_number = cell_child_number(cell);
		for ( i = 0; i < nDim; i++ ) {
			neighbor = cell_neighbor( cell, external_direction[cell_number][i] );
			if ( neighbor != NULL_OCT && 
					( cell_is_leaf(neighbor) && cell_level(neighbor) < cell_level(cell) ) ) {
				return -4;
			}
		}
	}
	
	result = split_cell(cell);
	if ( result == 0 ) {
		cart_assert( cell_is_refined(cell) );

		for ( i = 0; i < num_children; i++ ) {
			child_cell = cell_child( cell, i );
			cart_assert( child_cell >= 0 && child_cell < num_cells );

			/* find the neighbors used in interpolation */
			cell_interpolation_neighbors( cell, i, neighbors );

#ifdef HYDRO
			/* interpolate momentum */
			for ( j = 0; j < nDim; j++ ) {
				cell_momentum(child_cell,j) = cell_interpolate_with_neighbors( cell, HVAR_MOMENTUM+j, neighbors );
			}

			/* sum momentum */
			moment = 0.0;
			for ( j = 0; j < nDim; j++ ) {
				moment += cell_momentum(child_cell,j)*cell_momentum(child_cell,j);
			}
			weights[0] += sqrt(moment);

			/* interpolate density,  pressure and internal energy */
			cell_gas_density(child_cell) = cell_interpolate_with_neighbors( cell, HVAR_GAS_DENSITY, neighbors );
			weights[1] += cell_gas_density(child_cell);

			cell_gas_pressure(child_cell) = cell_interpolate_with_neighbors( cell, HVAR_PRESSURE, neighbors );
			weights[2] += cell_gas_pressure(child_cell);

			cell_gas_internal_energy(child_cell) = cell_interpolate_with_neighbors( cell, HVAR_INTERNAL_ENERGY, neighbors );
			weights[3] += cell_gas_internal_energy(child_cell);

			/* gamma is just copied */
			cell_gas_gamma(child_cell) = cell_gas_gamma(cell);

			/* interpolate potential and add to kinetic to get total energy */

			/*
			   ASK DOUG WHY WE DO IT THAT WAY
			   */
			cell_gas_energy(child_cell) = cell_gas_kinetic_energy(child_cell) + cell_gas_internal_energy(child_cell);
			//				cell_interpolate_function_with_neighbors( cell, cell_internal_energy, neighbors );
			weights[4] += cell_gas_energy(child_cell);

#ifdef	ELECTRON_ION_NONEQUILIBRIUM
			cell_electron_internal_energy(child_cell) = cell_interpolate_with_neighbors( cell, HVAR_ELECTRON_INTERNAL_ENERGY, neighbors );
			weights[5] += cell_electron_internal_energy(child_cell);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

			for ( j = 0; j < num_extra_energy_variables; j++ ) {
				cell_extra_energy_variables(child_cell,j) = 
					cell_interpolate_with_neighbors( cell, HVAR_EXTRA_ENERGY_VARIABLES+j, neighbors );
				weights[j+5+num_electronion_noneq_vars] +=  cell_extra_energy_variables(child_cell,j);
			}

#ifdef EXTRA_PRESSURE_SOURCE
                        /* no smearing -- maintain gradient and \dot{p}=P*A */
			cell_extra_pressure_source(child_cell) = cell_extra_pressure_source(cell);
#endif

			for ( j = 0; j < num_chem_species; j++ ) {
				cell_advected_variable(child_cell,j) = cell_interpolate_with_neighbors( cell, HVAR_ADVECTED_VARIABLES+j, neighbors );
				weights[num_hydro_vars-num_chem_species-nDim+j] += cell_advected_variable(child_cell,j);
			}
#endif /* HYDRO */

#ifdef GRAVITY
			cell_total_mass(child_cell) = cell_interpolate_with_neighbors( cell, VAR_TOTAL_MASS, neighbors );
			cell_potential(child_cell) = cell_interpolate_with_neighbors( cell, VAR_POTENTIAL, neighbors );
#ifdef PARTICLES
			cell_first_species_mass(child_cell) = 0.0;
#endif /* PARTICLES */

#ifdef HYDRO
			cell_potential_hydro(child_cell) = cell_interpolate_with_neighbors( cell, VAR_POTENTIAL_HYDRO, neighbors );
#endif /* HYDRO */
#endif /* GRAVITY */

#ifdef RADIATIVE_TRANSFER
			for ( j = 0; j < rt_num_vars; j++ ) {
				cell_var(child_cell,rt_grav_vars_offset+j) = cell_interpolate_with_neighbors( cell, rt_grav_vars_offset+j, neighbors );
			}
#endif /* RADIATIVE_TRANSFER */
		}

#ifdef HYDRO
		/* calculate scalings to conserve quantities */
		moment = 0.0;
		for ( j = 0; j < nDim; j++ ) {
			moment += cell_momentum(cell,j)*cell_momentum(cell,j);
		}
		weights[0] = (weights[0] == 0.0) ? 0.0: (double)num_children * sqrt(moment) / weights[0];
		weights[1] = (weights[1] == 0.0) ? 0.0: (double)num_children * (double)cell_gas_density(cell) / weights[1];
		weights[2] = (weights[2] == 0.0) ? 0.0: (double)num_children * (double)cell_gas_pressure(cell) / weights[2]; 
		weights[3] = (weights[3] == 0.0) ? 0.0: (double)num_children * (double)cell_gas_internal_energy(cell) / weights[3];
		weights[4] = (weights[4] == 0.0) ? 0.0: (double)num_children * (double)cell_gas_energy(cell) / weights[4];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		weights[5] = (weights[5] == 0.0) ? 0.0: (double)num_children * (double)cell_electron_internal_energy(cell) / weights[5];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
                for ( j = 0; j < num_extra_energy_variables; j++ ) {
                    ind = j+5+num_electronion_noneq_vars;
                    weights[ind] = (weights[ind] == 0.0) ? 0.0: 
                        (double)num_children *
                        (double)cell_extra_energy_variables(cell,j) / 
                        weights[ind];
                }

		for ( j = 0; j < num_chem_species; j++ ) {
			weights[num_hydro_vars-num_chem_species-nDim+j] = 
				(weights[num_hydro_vars-num_chem_species-nDim+j] == 0.0) ? 0.0: 
				(double)num_children * (double)cell_advected_variable(cell,j) / 
				weights[num_hydro_vars-num_chem_species-nDim+j];
		}

		/* enforce conservation laws in children */
		mass = 0.0;
		for ( i = 0; i < num_children; i++ ) {
			child_cell = cell_child(cell, i);
			cart_assert( child_cell >= 0 && child_cell < num_cells );

			for ( j = 0; j < nDim; j++ ) {
				cell_momentum(child_cell,j) *= weights[0];
			}

			cell_gas_density(child_cell) *= weights[1];
			cell_gas_pressure(child_cell) *= weights[2];
			cell_gas_internal_energy(child_cell) *= weights[3];
			cell_gas_energy(child_cell) *= weights[4];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			cell_electron_internal_energy(child_cell) *= weights[5];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

			for ( j = 0; j < num_extra_energy_variables; j++ ) {
				cell_extra_energy_variables(child_cell,j) *= weights[j+5+num_electronion_noneq_vars];
			}

			for ( j = 0; j < num_chem_species; j++ ) {
				cell_advected_variable(child_cell,j) *= weights[num_hydro_vars-num_chem_species-nDim+j];
			}

			mass += cell_gas_density(child_cell) * cell_volume[ cell_level(child_cell) ];
		}

		if ( mass>0.0 && fabs( mass - cell_gas_density(cell)*cell_volume[cell_level(cell)] )/mass > 1e-6 ) {
			cart_error("Error in mass conservation in split_cell: %e %e\n", mass,
					cell_gas_density(cell)*cell_volume[cell_level(cell)] );
		}

		/* ensure new values correspond to variable limits */
		hydro_magic_one_cell(cell);

#endif /* HYDRO */

	} else {
		cart_error("ERROR: split_cell(%u) = %d\n", cell, result );
	}

	return result;
}

#endif /* REFINEMENT */

