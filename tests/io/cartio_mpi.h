/*
 * cartio_mpi.h
 *
 *  Created on: Apr 9, 2010
 *      Author: eric-desktop
 */

#ifndef CARTIO_MPI_H_
#define CARTIO_MPI_H_

#include <mpi.h>
#include <string.h>
#include <stdio.h>

#define MAX_SIZE (1<<22)

typedef struct MPI_F
{
	MPI_File fh;
	char data[MAX_SIZE];
	int p;
}MPI_FFile;


int MPI_File_fwrite(MPI_FFile  * handle, void *buf,int count,MPI_Datatype type,MPI_Status *status);

int MPI_File_fflush(MPI_FFile * handle);

int MPI_File_fread(MPI_FFile * handle, void *buf, int count,MPI_Datatype type, MPI_Status *status);

int MPI_File_fclose(MPI_FFile * handle);


#endif /* CARTIO_MPI_H_ */
