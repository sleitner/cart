#include "defs.h"
#ifdef RADIATIVE_TRANSFER

#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "tree.h"

#include "rt_c2f.h"
#include "rt_tree.h"
#include "rt_utilities.h"
#include "rt_transfer.h"


void f2c_wrapper(frtpackradiationfields)(f2c_intg *n, f2c_real *data);
void f2c_wrapper(frtunpackradiationfields)(f2c_intg *n, f2c_real *data);


int rtWriteRFHelper(FILE *f, f2c_intg n, f2c_real *data, int fortran_style)
{
  int size;

  size = n*sizeof(f2c_real);
  if(fortran_style && fwrite(&size,sizeof(int),1,f)!=1) return 1;
  if(fwrite(data,sizeof(f2c_real),n,f) != n) return 1;
  if(fortran_style && fwrite(&size,sizeof(int),1,f)!=1) return 1;

#ifdef RT_SINGLE_SOURCE
  size = sizeof(float) + nDim*sizeof(double);
  if(fortran_style && fwrite(&size,sizeof(int),1,f)!=1) return 1;
  if(fwrite(&rtSingleSourceVal,sizeof(float),1,f) != 1) return 1;
  if(fwrite(rtSingleSourcePos,sizeof(double),nDim,f) != nDim) return 1;
  if(fortran_style && fwrite(&size,sizeof(int),1,f)!=1) return 1;
#endif

  return 0;
}


int rtReadRFHelper(FILE *f, f2c_intg n, f2c_real *data, int fortran_style)
{
  int size;

  size = n*sizeof(f2c_real);
  if(fortran_style && fread(&size,sizeof(int),1,f)!=1) return 1;
  if(fread(data,sizeof(f2c_real),n,f) != n) return 1;
  if(fortran_style && fread(&size,sizeof(int),1,f)!=1) return 1;

#ifdef RT_SINGLE_SOURCE
  size = sizeof(float) + nDim*sizeof(double);
  if(fortran_style && fread(&size,sizeof(int),1,f)!=1) return 1;
  if(fread(&rtSingleSourceVal,sizeof(float),1,f) != 1) return 1;
  if(fread(rtSingleSourcePos,sizeof(double),nDim,f) != nDim) return 1;
  if(fortran_style && fread(&size,sizeof(int),1,f)!=1) return 1;
#endif

  return 0;
}


void rtWriteRadiationFieldData(const char *fileroot, int fortran_style)
{
  FILE *f;
  f2c_intg n;
  f2c_real *data;
  char *filename;

  if(local_proc_id == MASTER_NODE)
    {
      n = 0;
      f2c_wrapper(frtpackradiationfields)(&n,0);
      
      if(n < 1)
	{
	  cart_error("Unable to pack Radiation Field data.");
	}

      filename = cart_alloc(char, (strlen(fileroot)+2) );
      strcpy(filename,fileroot);
      strcat(filename,"rf");

      f = fopen(filename,"w");
      if(f == 0)
	{
	  cart_error("Unable to open file %s for writing.",filename);
	}

      data = cart_alloc(f2c_real, n );
      f2c_wrapper(frtpackradiationfields)(&n,data);

      if(rtWriteRFHelper(f,n,data,fortran_style))
	{
	  cart_error("Error in writing file %s.",filename);
	}

      cart_free(filename);
      cart_free(data);
      fclose(f);
    }
}


void rtReadRadiationFieldData(const char *fileroot, int fortran_style)
{
  FILE *f;
  f2c_intg n;
  f2c_real *data;
  char *filename;

  n = 0;
  f2c_wrapper(frtunpackradiationfields)(&n,0);
  
  if(n < 1)
    {
      cart_error("Unable to unpack Radiation Field data.");
    }

  data = cart_alloc(f2c_real, n );

  if(local_proc_id == MASTER_NODE)
    {
      filename = cart_alloc(char, (strlen(fileroot)+2) );
      strcpy(filename,fileroot);
      strcat(filename,"rf");

      f = fopen(filename,"r");
      if(f == 0)
	{
	  cart_error("Unable to open file %s for reading.",filename);
	}

      if(rtReadRFHelper(f,n,data,fortran_style))
	{
	  cart_error("Error in reading file %s (endiness-awareness is not implemented).",filename);
	}
    }

  if(sizeof(f2c_real) == sizeof(float))
    {
      MPI_Bcast(data,n,MPI_FLOAT,MASTER_NODE,MPI_COMM_WORLD);
    }
  else
    {
      cart_error("Fortran type REAL and C type FLOAT are expected to be the same.");
    }
  f2c_wrapper(frtunpackradiationfields)(&n,data);
  
  cart_free(data);
 
  if(local_proc_id == MASTER_NODE)
    {
      cart_free(filename);
      fclose(f);
    }

#ifdef RT_SINGLE_SOURCE
  MPI_Bcast(&rtSingleSourceVal,1,MPI_FLOAT,MASTER_NODE,MPI_COMM_WORLD);
  MPI_Bcast(rtSingleSourcePos,nDim,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
#endif
}

#endif  // RADIATIVE_TRANSFER

