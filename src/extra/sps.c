#include "config.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "hydro.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rt.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "halo_tools.h"
#include "sps.h"


#if defined(PARTICLES) && defined(STAR_FORMATION)
void extGalLums(const char *fname, const struct HALO_LIST *halos, float rmax, float tSFR, const char *spsfile)
{
  int ih, i, j, nz, nt; 
  double r, rad, sum;
  char str[999];
  FILE *f;
  float *lfnu, *lz, *lt;
  char *is_star;
#ifdef COSMOLOGY
  double tNow = tphys_from_auni(auni[min_level]);
#else
  double tNow = units->time*tl[min_level]/constants->yr;
#endif
  float ltSFR = log10(tSFR);

  struct GalData
  {
    double Mstar;
    double Zstar;
    double SFR;
    double Lnu;
  };

  struct GalData *gals, *buf;


  if(halos == NULL)
    {
      cart_debug("No halo file is loaded. Skipping computing galaxy luminosities.");
      return;
    }

  f = fopen(spsfile,"r");
  if(f == NULL)
    {
      cart_debug("FSPS data file %s is not found.",spsfile);
      return;
    }

  if(fscanf(f,"%d %d",&nz,&nt) != 2) cart_error("Corrupted FSPS file at header: %d %d.",nz,nt);
 

  lz = cart_alloc(float,nz);
  lt = cart_alloc(float,nt);
  lfnu  = cart_alloc(float,nt*nz);

  if(fscanf(f,"%f",lz)!=1 || lz[0]!=0.0) cart_error("Corrupted FSPS file at loc(0,0).");
  for(i=0; i<nz; i++)
    {
      if(fscanf(f,"%g",lz+i)!=1 || lz[i]<1.0e-30) cart_error("Corrupted FSPS file at loc(0,%d).",i+1);
      /* cart_debug("Z[%d] = %g",i,lz[i]); */
      lz[i] = log10(lz[i]);
    }
  for(j=0; j<nt; j++)
    {
      if(fscanf(f,"%f",lt+j)!=1 || lt[j]<3) cart_error("Corrupted FSPS file at loc(%d,0).",j+1);
      /* cart_debug("lt[%d] = %g",j,lt[j]); */
      for(i=0; i<nz; i++)
	{
	  if(fscanf(f,"%g",lfnu+nz*j+i)!=1 || lfnu[nz*j+i]<1.0e-30) cart_error("Corrupted FSPS file at loc(%d,%d).",j+1,i+1);
	}
    }

  fclose(f);

  gals = cart_alloc(struct GalData,halos->num_halos);
  memset(gals,0,sizeof(struct GalData)*halos->num_halos);

  is_star = cart_alloc(char,num_particles);
#pragma omp parallel for default(none), private(j), shared(is_star,particle_species_indices,particle_id,num_particle_species,particle_level)
  for(j=0; j<num_particles; j++)
    {
      if(particle_level[j] == FREE_PARTICLE_LEVEL)
	{
	  is_star[j] = 0;
	}
      else
	{
	  is_star[j] = particle_is_star(j);
	}
    }

  /*
  // Loop over halos in reverse order, that way we mostly account for satellites
  */
  for(ih=halos->num_halos-1; ih>=0; ih--)
    {
      halo *h = halos->list + ih;
      double sumM = 0.0, sumZ = 0.0, sumS = 0.0, sumL = 0.0;

      rad = rmax*h->rhalo;

      /*
      //  Loop over particles
      */
#pragma omp parallel for default(none), private(j), shared(star_tbirth,star_metallicity_II,star_metallicity_Ia,star_initial_mass,particle_x,h,rad,lt,nt,lz,nz,lfnu,units,constants,tNow,tSFR,ltSFR,is_star), reduction(+:sumM,sumZ,sumS,sumL)
      for(j=0; j<num_particles; j++) if(is_star[j]==1 && (compute_distance_periodic(particle_x[j],(double *)h->pos)<rad))
	{
	  float ltStar, lzStar, fnuStar;
	  int iz, it;
	  float ms = star_initial_mass[j]*units->mass/constants->Msun;

	  sumM += ms;
#ifdef COSMOLOGY
	  ltStar = log10(1.0e-10+fabs(tNow-tphys_from_tcode(star_tbirth[j])));
#else
	  ltStar = log10(1.0e-10+fabs(tNow-units->time*star_tbirth[j]/constants->yr));
#endif /* COSMOLOGY */
	  if(ltStar < ltSFR) sumS += ms/tSFR;

	  if(ltStar < lt[0]) ltStar = lt[0];
	  if(ltStar > lt[nt-1]) ltStar = lt[nt-1];

	  lzStar = 0.0;
#ifdef ENRICHMENT
	  lzStar += star_metallicity_II[j];
#ifdef ENRICHMENT_SNIa
	  lzStar += star_metallicity_Ia[j];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
	  sumZ += star_initial_mass[j]*units->mass/constants->Msun*lzStar;

	  lzStar = log10(1.0e-30+lzStar);
	  if(lzStar < lz[0]) lzStar = lz[0];
	  if(lzStar > lz[nz-1]) lzStar = lz[nz-1];
	  
	  for(it=1; it<nt-1 && ltStar>lt[it]; it++); cart_assert(it>0 && it<nt);
	  for(iz=1; iz<nz-1 && lzStar>lz[iz]; iz++); cart_assert(iz>0 && iz<nz);

	  cart_assert(lt[it-1]<=ltStar && ltStar<=lt[it]);
	  cart_assert(lz[iz-1]<=lzStar && lzStar<=lz[iz]);

	  fnuStar = 
	    (lt[it]-ltStar)/(lt[it]-lt[it-1])*((lz[iz]-lzStar)/(lz[iz]-lz[iz-1])*lfnu[nz*(it-1)+iz-1]+(lzStar-lz[iz-1])/(lz[iz]-lz[iz-1])*lfnu[nz*(it-1)+iz]) +
	    (ltStar-lt[it-1])/(lt[it]-lt[it-1])*((lz[iz]-lzStar)/(lz[iz]-lz[iz-1])*lfnu[nz*it+iz-1]+(lzStar-lz[iz-1])/(lz[iz]-lz[iz-1])*lfnu[nz*it+iz]);

	  sumL += fnuStar*ms;

	  is_star[j] = 0;
	}

      gals[ih].Mstar = sumM;
      gals[ih].Zstar = sumZ/sumM;
      gals[ih].SFR = sumS;
      gals[ih].Lnu = sumL;
    }

  cart_free(lz);
  cart_free(lt);
  cart_free(lfnu);

  sum = 0.0;
#pragma omp parallel for default(none), private(j), shared(star_initial_mass,is_star), reduction(+:sum)
  for(j=0; j<num_particles; j++) if(is_star[j] == 1)
    {
      sum +=  star_initial_mass[j];
    }
  cart_free(is_star);
  cart_debug("Unaccounted stellar mass: %lg",sum*units->mass/constants->Msun);

  buf = cart_alloc(struct GalData,halos->num_halos);
  MPI_Reduce(gals,buf,halos->num_halos*sizeof(struct GalData)/sizeof(double),MPI_DOUBLE,MPI_SUM,MASTER_NODE,mpi.comm.run);

  cart_free(gals);
  gals = buf;

  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(fname,"w");
      if(f == NULL)
        {
          cart_error("Unable to open output file.");
        }
      fprintf(f,"# id.. Mstar.... Zstar.... SFR...... Lnu......\n");
      fprintf(f,"#\n");

      for(ih=0; ih<halos->num_halos; ih++) if(gals[ih].SFR > 0.0)
        {
          fprintf(f,"%6d %9.3e %9.3e %9.3e %9.3e\n",halos->list[ih].id,gals[ih].Mstar,gals[ih].Zstar,gals[ih].SFR,gals[ih].Lnu);
        }
      
      fclose(f);
    }

  cart_free(gals);
}
#endif /* PARTICLES && STAR_FORMATION */

