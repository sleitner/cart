#include "config.h"

#include <stdio.h>
#include <string.h>

#include "auxiliary.h"

#include "halos.h"
#include "gic_tools.h"

int mask_mode = GIC_MASK_MODE_PLAIN;
int mask_level = 1;
float mask_rvir_size = 5.0;
int mask_buffer_width = 2;

int main_analysis(int num_options, char **options) {
	int i;
	halo_list *halos;
	halo *h;
	char *tmp;
	const char *str;
	char *str2;
	FILE *input;
	int hid;
	const char *p;

	if ( !is_option_present("halo-catalog","hlist",1) || !is_option_present("halos","h",1) ||
			!is_option_present("output","o",1) ) {

		cart_error("Usage: gic_mask <config file> <aexp> -halo_catalog=<hlist> -halos=<file or , sep. ids> -output=<mask name> [options].\n"
			"Valid options are:\n"
			"  -hlist, --halo_catalog=<file>                 The halo catalog file generated either by hfind or cart_hfind.\n"
			"  -h, --halos=<option>                          The selected halos from halo catalog. Either a file with one halo id\n"
			"                                                per line or a comma separated list of halo ids.\n"
			"  -mode, --mask-mode=<plain,chull>              Controls which selection algorithm to apply in lagrangian coordinates.\n"
			"  -level, --mask-level=<int>                    Sets the highest resolution level of the mask.\n"
			"  -width, --mask-buffer-width=<int>             Sets the size of the lower-resolution particle buffer.\n"
			"  -radius, --mask-virial-radius-factor=<float>  Sets the size of the sphere, in units of a halo's virial radius,\n"
			"                                                within which particles are selected to define the high-resolution region\n"
            "                                                in lagrangian coordinates.\n" );
	}

	str = extract_option1("halo-catalog","hlist",NULL);
	if ( str != NULL ) {
		halos = load_halo_finder_catalog( str, 0, 0.0, 0.0, 0.0, 1e6 );
	}

	str = extract_option1("halos","h",NULL);
	if ( str != NULL ) {
		/* first check if file */
		input = fopen( str, "r" );
		if ( input != NULL ) {
			while ( fscanf( input, "%d", &hid ) == 1 ) {
				h = find_halo_by_id( halos, hid );
				if ( !h ) {
					cart_error("Halo %d not found in catalog!", hid );
				}
				h->flag = 1;
			}
			fclose(input);
		} else {
			/* assume it's a comma separated list of ids */
			p = str;
			while (*p) {
				if ( *p != ',' && ( *p < '0' || *p > '9' ) ) {	
					cart_error("Invalid value for parameter 'halos'");
				}
				p++;
			}

			i = 0;

			str2 = cart_alloc(char,strlen(str)+1);
			strcpy(str2,str); 
			p = strtok(str2, ",");
			for(; p != NULL ;){
				if ( sscanf( p, "%d", &hid ) != 1 ) {
					cart_error("Error parsing parameter 'halos'");
				}
				h = find_halo_by_id( halos, hid );
				if ( !h ) {
					cart_error("Halo %d not found in catalog!", hid );
				}
				h->flag = 1;	
				p = strtok(NULL, ",");
			}
			cart_free(str2);
		}
	}

	str = extract_option1("mask-mode","mode",NULL);
	if ( str != NULL ) {
		if ( strcmp(str,"plain") == 0 ) {
			mask_mode = GIC_MASK_MODE_PLAIN;
		} else if ( strcmp(str,"chull") == 0 ) {
			mask_mode = GIC_MASK_MODE_CHULL;
		} else {
			cart_error("mask-mode must be either 'plain' or 'chull'");
		}
	}

	str = extract_option1("mask-level","level",NULL);
	if ( str != NULL ) {
		if ( sscanf( str, "%d", &mask_level ) != 1 ) {
			cart_error("Error parsing parameter 'mask-level'");
		}
		/* do some error checking on mask level... */
	}

	str = extract_option1("mask-buffer-width","width",NULL);
	if ( str != NULL ) {
		if ( sscanf( str, "%d", &mask_buffer_width ) != 1 ) {
			cart_error("Error parsing parameter 'mask-buffer-width'");
		}
	}

	str = extract_option1("mask-virial-radius-factor","radius",NULL);
	if ( str != NULL ) {
		if ( sscanf( str, "%f", &mask_rvir_size ) != 1 ) {
			cart_error("Error parsing parameter 'mask-virial-radius-factor'");	
		}
	}

	//die_on_unknown_options();

	str = extract_option1("output","o",NULL);
	gicMakeMask( str, halos, mask_rvir_size, mask_mode, mask_level, mask_buffer_width);

	return 0;
}
