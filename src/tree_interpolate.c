#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "auxiliary.h"

void cell_interpolation_neighbors( int cell, int child_number, int neighbor[nDim] ) {
	int i;

	for ( i = 0; i < nDim; i++ ) {
		neighbor[i] = cell_neighbor(cell, pyramid_vertices[child_number][i]);
		cart_assert( neighbor[i] != NULL_OCT );
	}
}

double cell_interpolate_with_neighbors( int icell, int variable, int neighbor[nDim] ) {
	int i;
	double value;

	const double weight =
		#if nDim == 3
			0.25
		#else
			#error Unsupported number of dimensions for weight in cell_interpolate_with_neighbors
		#endif
	;
	
	value = cell_var(icell, variable);
	for ( i = 0; i < nDim; i++ ) {
		value += cell_var(neighbor[i],variable);
	}

	return weight*value;
}

double cell_interpolate( int icell, int child_number, int variable ) {
	int neighbors[nDim];

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( child_number >= 0 && child_number < num_children );
	
	cell_interpolation_neighbors(icell,child_number,neighbors);
	return cell_interpolate_with_neighbors(icell,variable,neighbors);
}

double cell_interpolate_function_with_neighbors( int icell, double function(int), int neighbor[nDim] ) {
	int i;
	double value;

	const double weight =
		#if nDim == 3
			0.25
		#else
			#error Unsupported number of dimensions for weight in cell_interpolate_with_neighbors
		#endif
	;

	value = function(icell);
	for ( i = 0; i < nDim; i++ ) {
		value += function(neighbor[i]);
	}

	return weight*value;
}

double cell_interpolate_function( int icell, int child_number, double function(int) ) {
	int neighbors[nDim];

	cell_interpolation_neighbors(icell,child_number, neighbors);
	return cell_interpolate_function_with_neighbors( icell, function, neighbors );
}
