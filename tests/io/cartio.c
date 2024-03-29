/*
 * cartio.c
 *
 *  Created on: Feb 21, 2010
 *  Author: Yongen Yu
 */

#include "cartio.h"
#include <limits.h>


inline void log_debug(char * fmt,...)
{
#ifdef DEBUG
	va_list list_p;
	va_start(list_p,fmt);
	vprintf(fmt,list_p);
	va_end(list_p);
#endif
}

/*
 * Open
 */
int cartio_open(char * filename, int mode, int num_files, cartio_file * handle )
{
	int status;
	log_debug("Begin to open file\n", filename);
	status = MPI_File_open(mpi.comm.run,filename, mode ,MPI_INFO_NULL,&(handle->ffh.fh));
	if (status != MPI_SUCCESS) {
	    char error_string[BUFSIZ];
	    int length_of_error_string, error_class;
	    MPI_Error_class(status, &error_class);
	    MPI_Error_string(error_class, error_string, &length_of_error_string);
	    printf("%s\n", error_string);
	    MPI_Error_string(status, error_string, &length_of_error_string);
	    printf("%s\n", error_string);
	    MPI_Abort(mpi.comm.run, status);
	  }
	log_debug("File %s has been open\n", filename);
	handle->ffh.p=0;
	memset(handle->ffh.data,0,MAX_SIZE);
	return 0;
}

void cartio_close(cartio_file * handle)
{
	MPI_Status status;
	MPI_File_fclose(&(handle->ffh));
}

void cartio_sfc_index_offset(cartio_file * handle, int * num_levels_per_root_tree,int * num_octs_per_root_tree)
{
	MPI_Status status;
	MPI_Offset offset_table[num_root_cells];
	MPI_Offset cur=0;
	int inc=0;
	offset_table[0]=cur;

	for(int i=1;i<num_root_cells;i++)
	{
		inc = sizeof(float)*handle->num_variables;
		inc+= sizeof(int)*(1+num_levels_per_root_tree[i-1]);
		inc+= num_octs_per_root_tree[i-1] *(sizeof(float)*handle->num_variables*8+sizeof(int)*8);
		cur+= inc;
		offset_table[i]=cur;
	}

	MPI_File_write(handle->ffh.fh,offset_table,num_root_cells,MPI_LONG_LONG_INT,&status);

}


/*
 * File Handle seeks according to the global sfc index
 */
void root_cell_location(cartio_file * handle,int sfc)
{
	MPI_Status status;
	MPI_Offset offset,c,size_of_sfc_table;
	size_of_sfc_table = sizeof(MPI_Offset)*num_root_cells;
	MPI_File_seek(handle->ffh.fh,0,MPI_SEEK_SET);
	MPI_File_seek(handle->ffh.fh,sizeof(MPI_Offset)*sfc,MPI_SEEK_SET);
	MPI_File_read(handle->ffh.fh, &offset, 1, MPI_LONG_LONG_INT, &status);

	MPI_File_seek(handle->ffh.fh,0,MPI_SEEK_SET);
	MPI_File_seek(handle->ffh.fh,(MPI_Offset)(size_of_sfc_table+offset),MPI_SEEK_SET);

	MPI_File_get_position(handle->ffh.fh, &c);
}


void cartio_root_cell_begin(cartio_file * handle, float * variables,int level, int * num_octs_per_level)
{
	MPI_Status status;
	MPIO_Request req;
	int size;

	MPI_File_fwrite(&(handle->ffh),variables,handle->num_variables,MPI_FLOAT,&status);

	MPI_File_fwrite(&(handle->ffh),&level,1,MPI_INT,&status);


	MPI_File_fwrite(&(handle->ffh),num_octs_per_level,level,MPI_INT,&status);
}

void cartio_root_cell_end(cartio_file * handle)
{
	//
}



void cartio_level_begin(cartio_file * handle,int level, int num_level_octs)
{
	//
}

void cartio_level_end(cartio_file * handle)
{
	//
}

void cartio_oct_begin(cartio_file * handle, float * variables, int cellrefined[8])
{
	//
	MPI_Status status;
	MPIO_Request req;
	int size;
	MPI_File_fwrite(&(handle->ffh),variables,8*handle->num_variables,MPI_INT,&status);

	MPI_File_fwrite(&(handle->ffh),cellrefined,8,MPI_INT,&status);
}

void cartio_oct_end(cartio_file * handle)
{
	//
}



