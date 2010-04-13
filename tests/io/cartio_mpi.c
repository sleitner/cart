/*
 * cartio_mpi.c
 *
 *  Created on: Apr 9, 2010
 *      Author: Yongen Yu
 *
 */

#include "cartio_mpi.h"
#include <stdlib.h>


extern void log_debug(char * fmt,...);

int MPI_File_fwrite(MPI_FFile * handle,void *buf,int count,MPI_Datatype type, MPI_Status *status)
{
	int pos = handle->p;
	int size = sizeof(type)*count;
	memcpy(&(handle->data[pos]),buf,size);
	handle->p+=size;

	if(handle->p>MAX_SIZE)
	{
		log_debug("Buffer overflow\n");
	}
	else if(handle->p>(MAX_SIZE*9/10))
	{
		MPI_File_write(handle->fh,handle->data,handle->p,MPI_BYTE,status);
		handle->p=0;
	}
	return 0;
}

int MPI_File_fflush(MPI_FFile * handle)
{
	MPI_Status  status;
	MPI_File_write(handle->fh,handle->data,handle->p,MPI_BYTE,&status);
	handle->p=0;
	return 0;
}

int MPI_File_fread(MPI_FFile * handle, void *buf, int count,MPI_Datatype type, MPI_Status *status)
{

	return 0;
}

int MPI_File_fclose(MPI_FFile * handle)
{
	MPI_File_fflush(handle);
	MPI_File_close(&(handle->fh));
	return 0;
}
