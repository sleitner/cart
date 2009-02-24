#ifndef __GIC_READER_H__
#define __GIC_READER_H__


/*
//  Fortran types (will depend on a platform)
//  Currently ported to: 32/64-bit Linux
*/
#define GIC_REAL        float
#define GIC_INTEGER     int
#define GIC_RECORD      int          /* ignore non-std g77 for now */

#define GIC_MANIFEST_SIZE   (256+9*sizeof(GIC_REAL))
#define GIC_FILEHEADER_SIZE   (8+2*sizeof(GIC_REAL)+6*sizeof(GIC_INTEGER))
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

/*
//  The different order here is for 32-bit/64-bit alignment.
*/
struct gicManifest
{
  char name[256];
  float OmegaB;
  float OmegaX;
  float OmegaL;
  float OmegaN;
  float h100;
  float dx;
  float ns;
  float s8;
  float kp;
};


struct gicFileHeader
{
  long Ntot;
  float aBegin;
  float DeltaDC;
  int dims[3];
  int seed;
  int Nrec;
  int Lmax;
};


struct gicLevelHeader
{
  long Nlev;
  float Mlev;
  int L;
  int Lmax;
  int ind;
};

  
int gicReadManifest(struct gicFile *f, struct gicManifest *manifest);
int gicReadFileHeader(struct gicFile *f, struct gicFileHeader *header);
int gicReadLevelHeader(struct gicFile *f, struct gicLevelHeader *header);

int gicReadFortranRecordReal(struct gicFile *f, GIC_REAL* buffer);
int gicReadFortranRecordInteger(struct gicFile *f, GIC_INTEGER *buffer);

int gicSkipFortranRecordReal(struct gicFile *f);
int gicSkipFortranRecordInteger(struct gicFile *f);

#endif  /* __GIC_READER_H__ */
