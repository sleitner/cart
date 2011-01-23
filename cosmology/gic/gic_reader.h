#ifndef __GIC_READER_H__
#define __GIC_READER_H__


/*
//  Fortran types (will depend on a platform)
//  Currently ported to: 32/64-bit Linux
*/
#define GIC_REAL        float
#define GIC_INTEGER     int
#define GIC_RECORD      int          /* ignore non-std g77 for now */

#define GIC_MANIFEST_SIZE   (256+1*sizeof(GIC_INTEGER))
#define GIC_FILEHEADER_SIZE   (8+112*sizeof(GIC_REAL)+8*sizeof(GIC_INTEGER))
#define GIC_LEVELHEADER_SIZE   (8+1*sizeof(GIC_REAL)+3*sizeof(GIC_INTEGER))


/*
//  Helper classes
*/
struct gicFile
{
  FILE *File;
  int WrongOrder;
  int Nrec;
};


struct gicManifest
{
  char name[256];
  int id;
};


/*
//  The different order here is for 32-bit/64-bit alignment.
*/
struct gicFileHeader
{
  long Ntot;
  float OmegaB;
  float OmegaX;
  float OmegaL;
  float OmegaN;
  float h100;
  float ns;
  float s8;
  float kp;
  float aBegin;
  float dx;
  float DeltaDC;
  float rmsDC;
  int dims[3];
  int seed[2];
  int Nrec;
  int Lmax;
  int flag;
};


struct gicLevelHeader
{
  long Nlev;
  float Mlev;
  int L;
  int Lmax;
  int ind;
};


int gicReadRecordHelper(FILE *f, long size, void* buffer, int *wrong_order);


int gicReadManifest(struct gicFile *f, struct gicManifest *manifest);
int gicReadFileHeader(struct gicFile *f, struct gicFileHeader *header);
int gicReadLevelHeader(struct gicFile *f, struct gicLevelHeader *header);

int gicReadFortranRecordReal(struct gicFile *f, GIC_REAL* buffer);
int gicReadFortranRecordInteger(struct gicFile *f, GIC_INTEGER *buffer);

int gicSkipFortranRecordReal(struct gicFile *f);
int gicSkipFortranRecordInteger(struct gicFile *f);

#endif  /* __GIC_READER_H__ */
