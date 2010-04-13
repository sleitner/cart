/*
 * cartio.h
 *
 *  Created on: Feb 21, 2010
 *      Author: Yongen YU
 */

#ifndef CARTIO_H_
#define CARTIO_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "cartio_mpi.h"


#define num_root_cells (1<<21)


typedef struct cartio_file_struct
{
	MPI_FFile  ffh;
	int num_variables;
} cartio_file;


void log_debug(char * fmt,...);

/*
 * Description	Open the file
 *
 * filename		The name of the file
 * mode			(read,write and append)
 * num_files	The number of files (1 for the first version)
 * handle		The handle of the file
 */
int cartio_open(char * filename, int mode, int num_files,cartio_file * handle);

/*
 * Description	Close the file
 */
void cartio_close(cartio_file * handle);


/*
 * Description	Create the global SFC index ->offset table to store the offset of each root level cell
 */

void cartio_sfc_index_offset(cartio_file * handle, int * num_levels_per_root_tree,int * num_octs_per_root_tree);


/*
 * Description	Move the File Handle to the very position for each node according to the global SFC index offset table
 */
void root_cell_location(cartio_file * handle,int sfc);


/*
 * Description		Output the variables of the root level cell and the hierarchy of the Oct tree correlated with this root level cell
 *
 * handle			The File handle
 * variables		The variables of the root level cell
 * level			The depth of the Oct tree correlated to the root level cell
 * num_level_octs	The array store the number of Oct nodes each level
 */
void cartio_root_cell_begin(cartio_file * handle, float * variables,int level, int * num_octs_per_level);

/*
 * Description		Do something at the end of writing the root level cell
 */
void cartio_root_cell_end(cartio_file * handle);


/*
 * Description		Do something at the beginning of each level
 */
void cartio_level_begin(cartio_file * handle,int level, int num_level_octs);

/*
 * Description		Do something at the end of each level
 */
void cartio_level_end(cartio_file * handle);

/*
 *Description		Output the data of a special oct tree node to the file
 *
 * Handle			The handle of the file
 * variables 		The array recording the variables of the eight cells belonging to this Octree node.
 */
void cartio_oct_begin(cartio_file * handle, float * variables, int refined[8]);

/*
 * Do something..
 */
void cartio_oct_end(cartio_file * handle);



#endif /* CARTIO_H_ */
