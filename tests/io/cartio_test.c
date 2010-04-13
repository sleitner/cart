/*
 * main.c
 *
 *  Created on: Feb 22, 2010
 *      Author: yongen yu
 */

#include "cartio.h"


#define MAX_LEVEL	6
#define MAX_PROCS	512
#define NUM_VAR		2

int num_procs;
int  proc_sfc_index[MAX_PROCS+1];
char * FILE_NAME = "cartio.grid";

cartio_file handle;

int my_rank;
int num_levels_per_root_tree[num_root_cells];
int num_octs_per_root_tree[num_root_cells];

int num_octs_per_level[]={1,2,4,8,16,32};

float  variables[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
int refined[8]={1,0,0,0,1,0,0,0};

extern void log_debug(char * fmt,...);


int cell_is_refined(int sfc)
{
	return 1;
}

void inital_test_data()
{
	int i;
	int index;

	for( i=0, index = 0; i<num_procs; i++, index +=num_root_cells/num_procs)
	{
		proc_sfc_index[i]=index;
	}
	proc_sfc_index[num_procs]=num_root_cells;

	for(int j=0 ;j<num_root_cells; j++)
	{
		num_levels_per_root_tree[j]=6;
		num_octs_per_root_tree[j]=63;
	}

}

int write_grid()
{
	int sfc;
	MPI_Offset c;
	log_debug("%d: Begin to write grid\n", my_rank);
	cartio_open(FILE_NAME,MPI_MODE_CREATE|MPI_MODE_RDWR,1,&handle);
	handle.num_variables=NUM_VAR;

	//create the global sfc index table
	if(0==my_rank)
	{
		log_debug("Create the global SFC index table\n");
		cartio_sfc_index_offset(&handle, num_levels_per_root_tree,num_octs_per_root_tree);
	}

	// Barrier. All the processes block here until the whole global sfc index table is recorded in the file.
	MPI_Barrier(MPI_COMM_WORLD);

	sfc = proc_sfc_index[my_rank];
	log_debug("%d SFC index %d\n",my_rank,proc_sfc_index[my_rank]);
	root_cell_location(&handle,sfc);

	MPI_File_get_position(handle.ffh.fh, &c);
	log_debug("%d: Position before write octs %lld\n",my_rank,c);

	for(;sfc<proc_sfc_index[my_rank+1];sfc++)
	{
		int level=MAX_LEVEL;
		cartio_root_cell_begin(&handle,variables,level,num_octs_per_level);
		if(cell_is_refined(sfc))
		{
			for( level = 1; level<=MAX_LEVEL;level++)
			{
				int i;
				cartio_level_begin(&handle,level,num_octs_per_level[level]);
				for(i=0;i<num_octs_per_level[level];i++)
				{
					cartio_oct_begin(&handle,variables, refined);
					cartio_oct_end(&handle);
				}
				cartio_level_end(&handle);
			}
		}
		cartio_root_cell_end(&handle);
	}
	cartio_close(&handle);
	return 1;
}

int cartio_oct_read(cartio_file * handle,int sfc)
{

	MPI_Status status;
	MPI_Offset offset;
	float var[2];
	float variables_of_oct[16];
	int level, oct_per_level[6], refiend[8];
	MPI_Offset c;
	root_cell_location(handle,sfc);
	MPI_File_get_position(handle->ffh.fh, &c);
	log_debug("position before read root cell %lld\n",c);
	MPI_File_read(handle->ffh.fh, &var, handle->num_variables, MPI_FLOAT, &status);
	MPI_File_get_position(handle->ffh.fh, &c);

	MPI_File_read(handle->ffh.fh, &level, 1, MPI_INT, &status);
	log_debug("%d: The depth of the Oct tree is %d\n",my_rank,level);

	MPI_File_read(handle->ffh.fh, &oct_per_level, level, MPI_INT, &status);
	log_debug("%d: The octs number of level: 1 %d; 2 %d;  3 %d; 4 %d\n",my_rank,oct_per_level[0],oct_per_level[1],oct_per_level[2],oct_per_level[3]);

	MPI_File_read(handle->ffh.fh, &variables_of_oct, 8*handle->num_variables, MPI_FLOAT, &status);
	log_debug("%d: variable of the OCT: %f %f %f %f %f %f ...\n",my_rank,variables_of_oct[0],variables_of_oct[1],variables_of_oct[2],variables_of_oct[3],variables_of_oct[4],variables_of_oct[5]);
	MPI_File_read(handle->ffh.fh, &refiend, 8, MPI_INT, &status);
	log_debug("%d: refinement: %d %d %d %d %d %d %d %d\n",my_rank,refiend[0],refiend[1],refiend[2],refiend[3],refiend[4],refiend[5],refiend[6],refiend[7]);
	return 0;
}

int read_grid()
{

	cartio_file hd_r;
	log_debug("Begin to read file\n");
	cartio_open(FILE_NAME,MPI_MODE_RDONLY,1,&hd_r);
	hd_r.num_variables=NUM_VAR;
	log_debug("%d SFC index %d\n",my_rank,proc_sfc_index[my_rank]);
	cartio_oct_read(&hd_r,proc_sfc_index[my_rank]);
	cartio_close(&hd_r);
	return 0;
}

int main(argc,argv)
int argc; char * argv[];
{
	int  pool_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


	log_debug("initial data\n");

	inital_test_data();

	if(0==my_rank)
	{
		log_debug("the number of root cells %d\n",num_root_cells);
	}

	log_debug("try to write data\n");

	write_grid();

	MPI_Barrier(MPI_COMM_WORLD);

	read_grid();

	MPI_Finalize();
	exit (0);
}

