#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <fitsio.h>

int main ( int argc, char *argv[] ) {
	FILE *input;
	int num_image_pixels;
	long pixels[2];
	double *image;
	double Lbox, width, depth, center[3];
	int status;
	fitsfile *fptr;

	if ( argc != 3 ) {
		fprintf(stderr,"Usage: convert_to_fits input.dat output.fits\n");
		exit(1);
	}

	/* read in image in data format */
	input = fopen( argv[1], "r");
	if ( input == NULL ) {
		fprintf(stderr, "Error opening file %s\n", argv[1] );
		exit(1);
	}

	fread( &num_image_pixels, sizeof(int), 1, input );
	fread( &Lbox, sizeof(double), 1, input );
	fread( center, sizeof(double), 3, input );
	fread( &width, sizeof(double), 1, input );
	fread( &depth, sizeof(double), 1, input );
	
	pixels[0] = pixels[1] = num_image_pixels;

	status = 0;
	fits_create_file( &fptr, argv[2], &status );

	image = malloc( num_image_pixels*num_image_pixels*sizeof(double) );
	fread( image, sizeof(double), num_image_pixels*num_image_pixels, input );

	fits_create_img(fptr, DOUBLE_IMG, 2, pixels, &status);
	fits_update_key( fptr, TDOUBLE, "Lbox", &Lbox, "Box size",&status);
	fits_update_key( fptr, TDOUBLE, "width", &width, "Slice width (Mpc/h comoving)",&status);
	fits_update_key( fptr, TDOUBLE, "depth", &depth, "Slice depth (Mpc/h comoving)",&status);
	fits_update_key( fptr, TDOUBLE, "x", &center[0], "Image Center (Mpc/h comoving)",&status);
	fits_update_key( fptr, TDOUBLE, "y", &center[1], "Image Center (Mpc/h comoving)",&status);
	fits_update_key( fptr, TDOUBLE, "z", &center[2], "Image Center (Mpc/h comoving)",&status);
	
	fits_write_img(fptr, TDOUBLE, 1, num_image_pixels*num_image_pixels, image, &status);

	fclose(input);

	fits_close_file(fptr, &status);            
	fits_report_error(stderr, status);  

	free( image );

	return status;
}
