/*
//  Helper functions for reading GIC files
*/
#include <stdio.h>
#include <string.h>


#include "gic_reader.h"


#define GIC_INTG        GIC_INTEGER  /* just showing off */
#define GIC_INT8        long long    /* ignore non-std g77 for now */


/*
//  Helper macros
*/
#define SWAP(val)   gicSwapBytes((char *)&(val),sizeof(val))

#define SET(type,val)   val = *((type *)((wrong_order ? gicSwapBytes(buffer,sizeof(type)) : buffer ))); buffer += sizeof(type)


/*
//  Helper functions
*/
char* gicSwapBytes(char *buffer, int size)
{
  int i;
  char tmp;
                                                                               
  for(i=0; i<size/2; i++)
    {
      tmp = buffer[i];
      buffer[i] = buffer[size - i - 1];
      buffer[size - i - 1] = tmp;
    }
  return buffer;
}


int gicReadRecordHelper(FILE *f, long size, void* buffer, int *wrong_order)
{
  GIC_RECORD s1, s2;

  if(buffer == NULL) return -1;
  if(f == NULL) return -2;

  /*
  //  Detect endiness
  */
  if(fread(&s1,sizeof(GIC_RECORD),1,f) != 1) return 1;
  if(s1 == size)
    {
      *wrong_order = 0;
    }
  else
    {
      *wrong_order = 1;
      SWAP(s1);
    }
  if(s1 != size) return 1;

  if(fread(buffer,size,1,f) != 1) return 1;

  if(fread(&s2,sizeof(GIC_RECORD),1,f) != 1) return 1;
  if(*wrong_order)
    {
      SWAP(s2);
    }
  if(s2 != size) return 1;

  return 0;
}


/*
//  Interface functions
*/
int gicReadManifest(struct gicFile *f, struct gicManifest *manifest)
{
  int i, ret, wrong_order;
  long InternalBuffer[GIC_MANIFEST_SIZE/sizeof(long)+1];
  char *buffer = (char *)InternalBuffer;

  if(f==NULL || manifest == NULL) return -1;

  ret = gicReadRecordHelper(f->File,GIC_MANIFEST_SIZE,buffer,&wrong_order);
  if(ret != 0) return ret;

  f->WrongOrder = wrong_order;
  f->Nrec = 0;

  memcpy(manifest->name,buffer,256);
  i = 255;
  while(i>=0 && manifest->name[i] == ' ') i--;
  manifest->name[i+1] = 0;

  buffer += 256;

  SET(GIC_REAL,manifest->OmegaB);
  SET(GIC_REAL,manifest->OmegaX);
  SET(GIC_REAL,manifest->OmegaL);
  SET(GIC_REAL,manifest->OmegaN);
  SET(GIC_REAL,manifest->h100);
  SET(GIC_REAL,manifest->dx);
  SET(GIC_REAL,manifest->ns);
  SET(GIC_REAL,manifest->s8);
  SET(GIC_REAL,manifest->kp);

  return 0;
}


int gicReadFileHeader(struct gicFile *f, struct gicFileHeader *header)
{
  int ret, wrong_order;
  long InternalBuffer[GIC_FILEHEADER_SIZE/sizeof(long)+1];
  char *buffer = (char *)InternalBuffer;

  if(f==NULL || header == NULL) return -1;

  ret = gicReadRecordHelper(f->File,GIC_FILEHEADER_SIZE,buffer,&wrong_order);
  if(f->WrongOrder != wrong_order) ret = -3;
  if(ret != 0) return ret;

  SET(GIC_REAL,header->aBegin);
  SET(GIC_REAL,header->DeltaDC);
  SET(GIC_INTG,header->dims[0]);
  SET(GIC_INTG,header->dims[1]);
  SET(GIC_INTG,header->dims[2]);
  SET(GIC_INTG,header->seed);
  SET(GIC_INTG,header->Nrec);
  SET(GIC_INT8,header->Ntot);
  SET(GIC_INTG,header->Lmax);

  f->Nrec = header->Nrec;

  return 0;
}


int gicReadLevelHeader(struct gicFile *f, struct gicLevelHeader *header)
{
  int ret, wrong_order;
  long InternalBuffer[GIC_LEVELHEADER_SIZE/sizeof(long)+1];
  char *buffer = (char *)InternalBuffer;

  if(f==NULL || header == NULL) return -1;

  ret = gicReadRecordHelper(f->File,GIC_LEVELHEADER_SIZE,buffer,&wrong_order);
  if(f->WrongOrder != wrong_order) ret = -3;
  if(ret != 0) return ret;

  SET(GIC_INTG,header->L);
  SET(GIC_INTG,header->Lmax);
  SET(GIC_INT8,header->Nlev);
  SET(GIC_REAL,header->Mlev);
  SET(GIC_INTG,header->ind);

  return 0;
}


int gicReadFortranRecordReal(struct gicFile *f, GIC_REAL *buffer)
{
  int i, ret, wrong_order;

  if(f==NULL || buffer == NULL) return -1;

  ret = gicReadRecordHelper(f->File,f->Nrec*sizeof(GIC_REAL),buffer,&wrong_order);
  if(f->WrongOrder != wrong_order) ret = -3;
  if(ret != 0) return ret;

  if(f->WrongOrder)
    {
      for(i=0; i<f->Nrec; i++) SWAP(buffer[i]);
    }

  return 0;
}


int gicReadFortranRecordInteger(struct gicFile *f, GIC_INTEGER *buffer)
{
  int i, ret, wrong_order;

  if(f==NULL || buffer == NULL) return -1;

  ret = gicReadRecordHelper(f->File,f->Nrec*sizeof(GIC_INTEGER),buffer,&wrong_order);
  if(f->WrongOrder != wrong_order) ret = -3;
  if(ret != 0) return ret;

  if(f->WrongOrder)
    {
      for(i=0; i<f->Nrec; i++) SWAP(buffer[i]);
    }

  return 0;
}


int gicSkipFortranRecordReal(struct gicFile *f)
{
  if(f == NULL) return -1;

  if(fseek(f->File,f->Nrec*sizeof(GIC_REAL)+2*sizeof(GIC_RECORD),SEEK_CUR) != 0) return 2;

  return 0;
}


int gicSkipFortranRecordInteger(struct gicFile *f)
{
  if(f == NULL) return -1;

  if(fseek(f->File,f->Nrec*sizeof(GIC_INTEGER)+2*sizeof(GIC_RECORD),SEEK_CUR) != 0) return 2;

  return 0;
}


