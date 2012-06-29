#include "config.h"

#include <stdio.h>

#include "auxiliary.h"
#include "times.h"

#include "power.h"
extern int step;

const char *current_working_directory = ".";

int main_analysis(int argc, char **argv) {
	char label[10];
	char filename[256];
	const char *output;
	const char *str;

#ifdef COSMOLOGY
	sprintf(label,"a%6.4f",abox[min_level]);
#else
	sprintf(label,"%05u",step);
#endif /* COSMOLOGY */

	if ( !is_option_present("dark","d",0) &&
			!is_option_present("gas","g",0) &&
			!is_option_present("baryons","b",0) &&
			!is_option_present("stars","s",0) &&	
			!is_option_present("total","t",0) ) {

        cart_error("Usage: cart_power <config file> <aexp> (one or more of -d -g -b -s -t) [options].\n"
            "Valid options are:\n"
            "  -output, --power-directory=<dir>              The directory into which the output power spectra will be placed.\n"
			"  -fold, --power-foldings=<int>                 The number of foldings used to extend power spectra beyond nyquist-freq.\n"
			"  -mesh, --power-mesh=<int>                     The size fft to use to compute power spectrum (power-of-two, defaults to num_grid).\n"
			"  Matter components:\n"
			"  -d, --dark                                    Dark matter.\n"
			"  -g, --gas                                     Gas.\n"
			"  -b, --baryons                                 Total baryons.\n"
			"  -s, --stars                                   Stars.\n"
            "  -t, --total                                   Total matter.\n" );
    }

	output = extract_option1("power-directory","output",NULL);
	if ( output == NULL ) {
		output = current_working_directory;
	}

	str = extract_option1("power-foldings","fold",NULL);
	if ( str != NULL ) {
		if ( sscanf( str, "%d", &num_power_foldings ) != 1 ) {
			cart_error("Error parsing parameter 'power-foldings'");
		}
	}

	str = extract_option1("power-mesh","mesh",NULL);
	if ( str != NULL ) {
		if ( sscanf( str, "%d", &power_mesh_size ) != 1 ) {
			cart_error("Error parsing parameter 'power-mesh'");
		}
	}

	if ( is_option_present("dark","d",0) ) {
		sprintf(filename,"%s/power_dark_%s.dat",output,label);
		compute_power_spectrum(filename, POWER_TYPE_DARK );
	}
	if ( is_option_present("gas","g",0) ) {
        sprintf(filename,"%s/power_gas_%s.dat",output,label);
        compute_power_spectrum(filename, POWER_TYPE_GAS );
    }
	if ( is_option_present("baryons","b",0) ) {
        sprintf(filename,"%s/power_baryons_%s.dat",output,label);
        compute_power_spectrum(filename, POWER_TYPE_BARYONS );
    }
	if ( is_option_present("stars","s",0) ) {
        sprintf(filename,"%s/power_stars_%s.dat",output,label);
        compute_power_spectrum(filename, POWER_TYPE_STARS );
    }
	if ( is_option_present("total","t",0) ) {
        sprintf(filename,"%s/power_total_%s.dat",output,label);
        compute_power_spectrum(filename, POWER_TYPE_TOTAL );
    }

	return 0;
}
