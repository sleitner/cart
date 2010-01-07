#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "kd.h"
#include "smooth.h"


#define BIGNOTQUITEMAXFLOAT		((float)1.0e+37)


int main(int argc,char **argv)
{
	KD kd;
	SMX smx;
	int nBucket,nSmooth,i,j;
	FILE *fp;
	float fPeriod[3];
	float period;

   	nBucket = 16;
	nSmooth = 24;

	if ( argc < 4 || argc > 5 ) {
		fprintf(stderr,"USAGE: smooth PMcrd.DAT PMcrs.DAT [stars.dat] rho_smooth.den\n");
		exit(1);
	}

	kdInit(&kd,nBucket);

	/* read in here */
	if ( argc == 4 ) {
		kdReadART(kd,argv[1],argv[2],&period);
	} else {
		kdReadARTStars(kd,argv[1],argv[2],argv[3],&period);
	}

	fPeriod[0] = period;
	fPeriod[1] = period;
	fPeriod[2] = period;

	printf("bulding tree...\n");
	kdBuildTree(kd);
	printf("initializing density...\n");
	smInit(&smx,kd,nSmooth,fPeriod);
	printf("smoothing density...\n");
	smSmooth(smx,smDensitySym);
	printf("ordering tree...\n");
	kdOrder(kd);

	printf("writing output...\n");

	/* write density file */
	if ( argc == 4 ) {
		fp = fopen(argv[3], "w" );
	} else {
		fp = fopen(argv[4],"wb");	
	}
	assert(fp != NULL);
	smOutDensityBin(smx,fp);
	fclose(fp);

	smFinish(smx);
	kdFinish(kd);
	return 0;
	}
	

