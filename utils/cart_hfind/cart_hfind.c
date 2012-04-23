#include "config.h"

#include <stdio.h>

#include "auxiliary.h"

#include "halos.h"
#include "halo_finder.h"

int main_analysis(int argc, char **argv) {
	halo_list *list;

	list = find_halos();
	write_halo_list( list );
	
	return 0;
}
