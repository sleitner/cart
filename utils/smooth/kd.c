#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include "kd.h"
#include "tipsydefs.h"


#define MAX_ROOT_ITTR	32


void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket)
{
	KD kd;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	kd->p = NULL;
	kd->kdNodes = NULL;
	*pkd = kd;
	return(1);
	}

void reorder( char *buffer, int size ) {
        int i;
        char tmp;
                                                                                                                                                            
        for ( i = 0; i < (size/2); i++ ) {
                tmp = buffer[i];
                buffer[i] = buffer[size - i - 1];
                buffer[size - i - 1] = tmp;
        }
}

typedef struct {
        float aexpn;
        float aexp0;
        float amplt;
        float astep;
        int   istep;
        float partw;
        float tintg;
        float ekin;
        float ekin1;
        float ekin2;
        float au0;
        float aeu0;
        int   Nrow;
        int   Ngrid;
        int   Nspecies;
        int   Nseed;
        float Om0;
        float Oml0;
        float hubble;
        float Wp5;
        float Ocurv;
        float Omb0;
        float mass[10];
        int   num[10];
        float fill[80];
} particle_header;
        
typedef struct {
        float aexpn;
        float aexp0;
        float amplt;
        float astep;
        int   istep;
        float partw;
        float tintg;
        float ekin;
        float ekin1;
        float ekin2;
        float au0;
        float aeu0;
        int   Nrow;
        int   Ngrid;
        int   Nspecies;
        int   Nseed;
        float Om0;
        float Oml0;
        float hubble;
        float Wp5;
        float Ocurv;
        float mass[10];
        int   num[10];
        float fill[80];
} nbody_particle_header;
                                                                                                                                                    
int kdReadARTStars(KD kd,char *header, char *data, char *stardata, float *period ) 
{
	int i,j, nCnt;

	FILE *header_file;
	FILE *data_file;
	FILE *star_file;
	particle_header gas_header;
	int size;
	int endian;
	int num_grid;
	int num_particles;
	int page_size;
	int num_pages;	
	float *input_page;
	int num_parts_in_page;
	int current_specie;
	char desc[45];
	int num_stars;
	double t, aexp;
	double total_stellar_mass, total_initial_stellar_mass;
	int num_parts;

	endian = 0;

	header_file = fopen(header, "r");
	if ( header_file == NULL ) {
		fprintf(stderr,"Unable to open %s!\n", header );
		exit(1);
	}

	data_file = fopen(data,"rb");
	if ( data_file == NULL ) {
		fprintf(stderr,"Unable to open %s!\n", data );
		exit(1);
	}

	fread( &size, sizeof(int), 1, header_file );

	if ( size != sizeof(particle_header)+45 ) {
		reorder( (char *)&size, sizeof(int) );

		if ( size == sizeof(particle_header)+45 ) {
			endian = 1;
		} else {
			fprintf(stderr,"Error reading from particle header, size = %u!\n", size );
			exit(1);
		}
	}

	fread( desc, sizeof(char), 45, header_file );

	fread( &gas_header, sizeof(particle_header), 1, header_file );

	if ( endian ) {
		reorder( (char *)&gas_header.Ngrid, sizeof(int) );
		reorder( (char *)&gas_header.Nrow, sizeof(int) );
		reorder( (char *)&gas_header.Nspecies, sizeof(int) );

		for ( i = 0; i < gas_header.Nspecies; i++ ) {
			reorder( (char *)&gas_header.num[i], sizeof(int) );
			reorder( (char *)&gas_header.mass[i], sizeof(float) );
		}
	}

	page_size = gas_header.Nrow*gas_header.Nrow;
	num_grid = gas_header.Ngrid;
	num_particles = gas_header.num[ gas_header.Nspecies - 1 ];

	printf("num_particles = %u\n", num_particles );

	fclose( header_file );

	star_file = fopen( stardata, "r" );
        if ( star_file == NULL ) {
                fprintf(stderr,"Unable to open %s!\n", star_file );
                exit(1);
        }

	/* read header information */
	fread( &size, sizeof(int), 1, star_file );
	fread( &t, sizeof(double), 1, star_file );
	fread( &aexp, sizeof(double), 1, star_file );
	fread( &size, sizeof(int), 1, star_file );

	fread( &size, sizeof(int), 1, star_file );
	fread( &num_stars, sizeof(int), 1, star_file );
	
	if ( endian ) {
		reorder( (char *)&num_stars, sizeof(int) );
	}
	fread( &size, sizeof(int), 1, star_file );
	
	fread( &size, sizeof(int), 1, star_file );
	fread( &total_stellar_mass, sizeof(double), 1, star_file );
	fread( &total_initial_stellar_mass, sizeof(double), 1, star_file );
	fread( &size, sizeof(int), 1, star_file );

	printf("num_stars = %u\n", num_stars );

	*period = (float)num_grid;

	num_pages = (num_particles-1) / page_size + 1;

	printf("endian = %u\n", endian );
	printf("num_grid = %u\n", num_grid );
	printf("%u particles\n", num_particles );
	printf("Page size = %u\n", page_size );
	printf("%u pages\n", num_pages );

	/* allocate page */
	input_page = malloc( 6*page_size * sizeof(float) );
	assert( input_page != NULL );

	kd->nParticles = num_stars;
	kd->nDark = num_stars;
	kd->nGas = 0;
	kd->nStar = 0;
	kd->fTime = 0;
	kd->nActive = num_stars;
	kd->bDark = 1;
	kd->bGas = 0;
	kd->bStar = 0;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	num_parts = 0;
	current_specie = 0;
	for (i=0;i<num_pages;++i) {
		if ( i == num_pages - 1 ) {
			num_parts_in_page = num_particles - page_size*(num_pages-1);
		} else {
			num_parts_in_page = page_size;
		}

		/* easiest to just read in entire page even though we're only using positions */
		if ( fread( input_page, sizeof(float), 6*page_size, data_file ) != 6*page_size ) {
			fprintf(stderr,"Error reading page %u from data file!\n", i );
			exit(1);
		}

		if ( endian ) {
			for ( j = 0; j < num_parts_in_page; j++ ) {
				reorder( (char *)&input_page[j], sizeof(float) );
				reorder( (char *)&input_page[j+page_size], sizeof(float) );
				reorder( (char *)&input_page[j+2*page_size], sizeof(float) );
			}
		}
			
		for ( j = 0; j < num_parts_in_page; j++ ) {
			if ( num_parts == gas_header.num[current_specie] ) {
				current_specie++;
			}

			if ( current_specie == gas_header.Nspecies - 1 ) {
				kd->p[nCnt].iOrder = nCnt;
				kd->p[nCnt].iMark = 1;
				kd->p[nCnt].r[0] = input_page[j]-1.0;
				kd->p[nCnt].r[1] = input_page[j+page_size]-1.0;
				kd->p[nCnt].r[2] = input_page[j+2*page_size]-1.0;

				++nCnt;
			}

			num_parts++;
			}
		}
	for ( j = 0; j < num_stars; j++ ) {
		fread( &kd->p[j].fMass, sizeof(float), 1, star_file );
	}
	free( input_page );
	fclose(data_file);
	return(kd->nParticles);
	}


int kdReadART(KD kd,char *header, char *data, float *period ) 
{
	int i,j, nCnt;

	FILE *header_file;
	FILE *data_file;
	particle_header gas_header;
	nbody_particle_header nbody_header;
	int size;
	int nbody_flag;
	int endian;
	int num_grid;
	int num_particles;
	int page_size;
	int num_pages;	
	int num[10];
	float mass[10];
	float *input_page;
	int num_parts_in_page;
	int current_specie;
	char desc[45];

	nbody_flag = 0;
	endian = 0;

	header_file = fopen(header, "r");
	if ( header_file == NULL ) {
		fprintf(stderr,"Unable to open %s!\n", header );
		exit(1);
	}

	data_file = fopen(data,"rb");
	if ( data_file == NULL ) {
		fprintf(stderr,"Unable to open %s!\n", data );
		exit(1);
	}

	fread( &size, sizeof(int), 1, header_file );

	if ( size == sizeof(nbody_particle_header)+45 ) {
		nbody_flag = 1;
	} else if ( size != sizeof(particle_header)+45 ) {
		reorder( (char *)&size, sizeof(int) );

		if ( size == sizeof(particle_header)+45 ) {
			endian = 1;
		} else if ( size == sizeof(nbody_particle_header)+45 ) {
			endian = 1;
			nbody_flag = 1;
		} else {
			fprintf(stderr,"Error reading from particle header, size = %u!\n", size );
			exit(1);
		}
	}

	fread( desc, sizeof(char), 45, header_file );

	if ( nbody_flag ) {
		fread( &nbody_header, sizeof(nbody_particle_header), 1, header_file );

		if ( endian ) {
			reorder( (char *)&nbody_header.Ngrid, sizeof(int) );
			reorder( (char *)&nbody_header.Nrow, sizeof(int) );
			reorder( (char *)&nbody_header.Nspecies, sizeof(int) );

			for ( i = 0; i < nbody_header.Nspecies; i++ ) {
				reorder( (char *)&nbody_header.num[i], sizeof(int) );
				reorder( (char *)&nbody_header.mass[i], sizeof(float) );
			}
		}

		for ( i = 0; i < nbody_header.Nspecies; i++ ) {
			num[i] = nbody_header.num[i];
			mass[i] = nbody_header.mass[i];
		}

		num_particles = nbody_header.num[0];

		printf("Nrow = %u\n", nbody_header.Nrow );
		page_size = nbody_header.Nrow*nbody_header.Nrow;
		num_grid = nbody_header.Ngrid;
	} else {
		fread( &gas_header, sizeof(particle_header), 1, header_file );

		if ( endian ) {
			reorder( (char *)&gas_header.Ngrid, sizeof(int) );
			reorder( (char *)&gas_header.Nrow, sizeof(int) );
			reorder( (char *)&gas_header.Nspecies, sizeof(int) );

			for ( i = 0; i < gas_header.Nspecies; i++ ) {
				reorder( (char *)&gas_header.num[i], sizeof(int) );
				reorder( (char *)&gas_header.mass[i], sizeof(float) );
			}
		}

		/* check for star particles */
		if ( gas_header.mass[ gas_header.Nspecies -1 ] == 0.0 ) {
			gas_header.Nspecies--;
		}

		for ( i = 0; i < gas_header.Nspecies; i++ ) {
			num[i] = gas_header.num[i];
			mass[i] = gas_header.mass[i];
		}

		/* only take first specie */
		num_particles = gas_header.num[0];
		page_size = gas_header.Nrow*gas_header.Nrow;
		num_grid = gas_header.Ngrid;
	}

	fclose( header_file );

	*period = (float)num_grid;

	num_pages = (num_particles-1) / page_size + 1;

	printf("endian = %u\n", endian );
	printf("num_grid = %u\n", num_grid );
	printf("%u particles\n", num_particles );
	printf("Page size = %u\n", page_size );
	printf("%u pages\n", num_pages );

	/* allocate page */
	input_page = malloc( 6*page_size * sizeof(float) );
	assert( input_page != NULL );

	kd->nParticles = num_particles;
	kd->nDark = num_particles;
	kd->nGas = 0;
	kd->nStar = 0;
	kd->fTime = 0;
	kd->nActive = num_particles;
	kd->bDark = 1;
	kd->bGas = 0;
	kd->bStar = 0;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<num_pages;++i) {
		if ( i == num_pages - 1 ) {
			num_parts_in_page = num_particles - page_size*(num_pages-1);
		} else {
			num_parts_in_page = page_size;
		}

		/* easiest to just read in entire page even though we're only using positions */
		if ( fread( input_page, sizeof(float), 6*page_size, data_file ) != 6*page_size ) {
			fprintf(stderr,"Error reading page %u from data file!\n", i );
			exit(1);
		}

		if ( endian ) {
			for ( j = 0; j < num_parts_in_page; j++ ) {
				reorder( (char *)&input_page[j], sizeof(float) );
				reorder( (char *)&input_page[j+page_size], sizeof(float) );
				reorder( (char *)&input_page[j+2*page_size], sizeof(float) );
			}
		}
			
		for ( j = 0; j < num_parts_in_page; j++ ) {
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			kd->p[nCnt].r[0] = input_page[j]-1.0;
			kd->p[nCnt].r[1] = input_page[j+page_size]-1.0;
			kd->p[nCnt].r[2] = input_page[j+2*page_size]-1.0;

			++nCnt;
			}
		}
	current_specie = 0;
	for ( j = 0; j < num_particles; j++ ) {
		if ( j == num[current_specie] ) {
			current_specie++;
		}

		kd->p[j].fMass = mass[current_specie];
	}
	free( input_page );
	fclose(data_file);
	return(kd->nParticles);
	}

void kdInMark(KD kd,char *pszFile)
{
	FILE *fp;
	char ach[80];
	int i,iCnt,iDum;

	fp = fopen(pszFile,"r");
	if (!fp) {
		fprintf(stderr,"Could not open mark array, %s\n",pszFile);
		exit(1);
		}
	fgets(ach,80,fp);	/* ignore the array header! */
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	fclose(fp);
	}


void kdSelect(KD kd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	double v;
	int i,j;

	p = kd->p;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void kdUpPass(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		kdUpPass(kd,l);
		kdUpPass(kd,u);
		kdCombine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->p[u].r[j];
			c[iCell].bnd.fMax[j] = kd->p[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->p[pj].r[j];
				if (kd->p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->p[pj].r[j];
				}
			}
		}
	}


void kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes != NULL) free(kd->kdNodes);
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->p[0].r[j];
		bnd.fMax[j] = kd->p[0].r[j];
		}
	for (i=1;i<kd->nActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->p[i].r[j]) 
				bnd.fMin[j] = kd->p[i].r[j];
			else if (bnd.fMax[j] < kd->p[i].r[j])
				bnd.fMax[j] = kd->p[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelect(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->p[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	kdUpPass(kd,ROOT);
	}


int cmpParticles(const void *v1,const void *v2)
{
	PARTICLE *p1=(PARTICLE *)v1,*p2=(PARTICLE *)v2;
	
	return(p1->iOrder - p2->iOrder);
	}


void kdOrder(KD kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLE),cmpParticles);
	}


void kdFinish(KD kd)
{
	free(kd->p);
	free(kd->kdNodes);
	free(kd);
	}

