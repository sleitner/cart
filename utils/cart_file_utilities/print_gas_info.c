#include <stdlib.h>
#include <stdio.h>

void reorder( char *buffer, int size ) {
	int i;
	char tmp;

	for ( i = 0; i < (size/2); i++ ) {
		tmp = buffer[i];
		buffer[i] = buffer[size - i - 1];
		buffer[size - i - 1] = tmp;
	}
}

int main ( int argc, char *argv[] ) {
	FILE *input;
	char job[256];
	int size;
	int endian;
	int step;
	double t, dt;
	float adum, ainit;
	float boxh, Om0, Oml0, Omb0, h;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int maxlevel, minlevel;
	double tl[10], dtl[10], tl_old[10], dtl_old[10];
	int level_sweep_dir[10];
	int ncell0;
	float refinement_volume_min[3], refinement_volume_max[3];
	float star_formation_volume_min[3], star_formation_volume_max[3];
	int sfc_order;
	int i, j;
	long long num_cells;
	int page_count;
	int counts[1024];

	/* open file and read header */
	input = fopen( argv[1],"r");
	if ( input == NULL ) {
		printf( "Unable to open file %s for input!\n", argv[1]);
	}

	fread(&size, sizeof(int), 1, input );
	endian = 0;
	if ( size != 256 ) {
		reorder( (char *)&size, sizeof(int) );
		if ( size != 256 ) {
			printf("Error: file %s is corrupted\n", argv[1] );
		} else {
			endian = 1;
			printf("Reordering bytes (file endianness is opposite program)\n");
		}
	}

	fread(&job, sizeof(char), 256, input );
	fread(&size, sizeof(int), 1, input );

	printf("job: %s\n", job );

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

	printf("step = %u\n", step );
	printf("t = %e\n", t );
	printf("dt = %e\n", dt );
	printf("adum = %e\n", adum );
	printf("ainit = %e\n", ainit );

	/* boxh, Om0, Oml0, Omb0, hubble */
	fread( &size, sizeof(int), 1, input );
	fread( &boxh, sizeof(float), 1, input );
	fread( &Om0, sizeof(float), 1, input );
	fread( &Oml0, sizeof(float), 1, input );
	fread( &Omb0, sizeof(float), 1, input );
	fread( &h, sizeof(float), 1, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&boxh, sizeof(float) );
		reorder( (char *)&Om0, sizeof(float) );
		reorder( (char *)&Oml0, sizeof(float) );
		reorder( (char *)&Omb0, sizeof(float) );
		reorder( (char *)&h, sizeof(float) );
	}

	printf("boxh = %f\n", boxh );
	printf("Om0 = %f\n", Om0 );
	printf("Oml0 = %f\n", Oml0 );
	printf("Omb0 = %e\n", Omb0 );
	printf("h = %f\n", h );

	/* nextra (no evidence extras are used...) extra lextra */
	fread( &size, sizeof(int), 1, input );
	fread( &nextras, sizeof(int), 1, input );
	fread( &size, sizeof(int), 1, input );

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

	printf("minlevel = %u\n", minlevel );
	printf("maxlevel = %u\n", maxlevel );

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

	printf("tl dtl tl_old dtl_old\n");
	for ( i = minlevel; i <= maxlevel; i++ ) {
		printf("%f %f %f %f\n", tl[i], dtl[i], tl_old[i], dtl_old[i] );
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

	printf("sfc_order = %u\n", sfc_order );

	/* refinement volume */
	fread( &size, sizeof(int), 1, input );
	fread( refinement_volume_min, sizeof(float), 3, input );
	fread( refinement_volume_max, sizeof(float), 3, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = 0; i < 3; i++ ) {
			reorder( (char *)&refinement_volume_min[i], sizeof(float) );
			reorder( (char *)&refinement_volume_max[i], sizeof(float) );
		}
	}

	for ( i = 0; i < 3; i++ ) {
		printf("refinement_volume[%u] = %f to %f\n", i, refinement_volume_min[i], refinement_volume_max[i] );
	}

	/* star formation volume */
	fread( &size, sizeof(int), 1, input );

	if ( size == 6*sizeof(float) ) {
		printf("reading star formation information...\n");

		fread( &star_formation_volume_min[i], sizeof(float), 3, input );
		fread( &star_formation_volume_max[i], sizeof(float), 3, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < 3; i++ ) {
				reorder( (char *)&star_formation_volume_min[i], sizeof(float) );
				reorder( (char *)&star_formation_volume_max[i], sizeof(float) );
			}
		}

		/* ncell0 */
		fread( &size, sizeof(int), 1, input );
	}

	fread( &ncell0, sizeof(int), 1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&ncell0, sizeof(int) );
	}

	printf("num_root_cells = %d\n", ncell0 );

	/* figure out how many total cells in this file */
	fread( &size, sizeof(int), 1, input );
	
	if ( endian ) {
		reorder( (char *)&size, sizeof(int) );
	}
	ncell0 = size / sizeof(int);

	printf("num_root_cells in file = %d\n", ncell0 );

	num_cells = 0;
	while ( ncell0 > 0 ) {
		page_count = (ncell0 < 1024) ? ncell0 : 1024;
		fread( counts, sizeof(int), page_count, input );

		if ( endian ) {
			for ( j = 0; j < page_count; j++ ) {
				reorder( (char *)&counts[j], sizeof(int) );
			}
		}

		for ( j = 0; j < page_count; j++ ) {
			num_cells += counts[j];
		}

		ncell0 -= page_count;
	}

	printf("total cells = %ld\n", num_cells );
	printf("total octs = %ld\n", (num_cells-ncell0)/8 );

	return 0;
}
