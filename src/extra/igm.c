#include "config.h"

#ifdef COSMOLOGY

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "density.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rt.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "F/frt_c.h"

#include "los.h"
#include "halo_finder.h"
#include "ifrit.h"
#include "ism.h"




#ifndef ISM_BUFFER_SIZE
#define ISM_BUFFER_SIZE 100
#endif


#ifdef RADIATIVE_TRANSFER

void extDumpRadiationBackground(const char *filename)
{
  int i, n;
  FILE *f;
  float *wlen, *ngxi, uJnu = 1.5808e-17/pow(auni[min_level],3.0);

  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(filename,"w");
      cart_assert(f != NULL);

      fprintf(f,"# wlen[A]  jnu[cgs] at a=%lf\n",auni[min_level]);

      n = 501;
      wlen = cart_alloc(float,n);
      ngxi = cart_alloc(float,n);
      
      for(i=0; i<n; i++) wlen[i] = 10*912/pow(10.0,0.01*i);

      rtGetRadiationField(-1,n,wlen,ngxi);

      for(i=0; i<n; i++) fprintf(f,"%9.3e %9.3e\n",wlen[i],uJnu*ngxi[i]);

      cart_free(wlen);
      cart_free(ngxi);

      fclose(f);
    }
}


void extExtractRadiationField(int n, const float *wlen, float *mean_rf)
{
  int i;
  int vars[rt_num_vars];
  float ngxi[rt_num_vars];
  MESH_RUN_DECLARE(level,cell);
  float uJnu = 1.5808e-17/pow(auni[min_level],3.0);

  cart_assert(n>0 && n<=rt_num_vars);
  
  /*
  //  Find the mean field
  */
  if(mean_rf != NULL)
    {
      rtGetRadiationField(-1,n,wlen,mean_rf);
      if(local_proc_id == MASTER_NODE)
	{
	  for(i=0; i<n; i++)
	    {
	      cart_debug("Radiation background[%e] = %e [cgs] = %e [CU]",wlen[i],uJnu*mean_rf[i],mean_rf[i]);
	    }
	}
    }

  for(i=0; i<n; i++)
    {
      vars[i] = rt_disk_offset + i;
    }

  /*
  //  Extract radiation field and store it in the cell_vars array
  */
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Extracting radiation field on level %d...",level);
    }
#pragma omp parallel for default(none), private(_Index,cell,i,ngxi), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,vars,n,wlen)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
      rtGetRadiationField(cell,n,wlen,ngxi);
      for(i=0; i<n; i++)
	{
	  cell_var(cell,vars[i]) = ngxi[i];
	}

    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  
  update_buffer_level(level,vars,n);

  MESH_RUN_OVER_LEVELS_END;
}  


void extExtractPhotoRates(int n, const int *idx, float *mean_rate)
{
  int i;
  int vars[rt_num_vars];
  float rate[FRT_RATE_DIM], rate0[FRT_RATE_DIM];
  MESH_RUN_DECLARE(level,cell);

  cart_assert(n>0 && n<=rt_num_vars);

  /*
  //  Find the mean rates
  */
  if(mean_rate != NULL)
    {
      rtGetPhotoRates(-1,rate0);
      for(i=0; i<n; i++)
	{
	  mean_rate[i] = rate0[idx[i]];
	}

      if(local_proc_id == MASTER_NODE)
	{
	  for(i=0; i<n; i++)
	    {
	      cart_debug("Mean rate[%d] = %e [cgs]",i,mean_rate[i]);
	    }
	}
    }

  for(i=0; i<n; i++)
    {
      vars[i] = rt_disk_offset + i;
    }

  /*
  //  Extract radiation field and store it in the cell_vars array
  */
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Extracting photo rates on level %d...",level);
    }
#pragma omp parallel for default(none), private(_Index,cell,i,rate), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,vars,n,idx)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
      rtGetPhotoRates(cell,rate);
      for(i=0; i<n; i++)
	{
	  cell_var(cell,vars[i]) = rate[idx[i]];
	}
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  
  update_buffer_level(level,vars,n);

  MESH_RUN_OVER_LEVELS_END;
}  


#if defined(RT_TRANSFER) && defined(RT_EXTERNAL_BACKGROUND)

/*
//  Find proximity zones around galaxies
*/
float *extProximityZone_MeanRF;
#define NBINS 4

typedef struct extProximityZone_DataType
{
  float R[NBINS];
}
extProximityZone_Data;


int extProximityZone_Worker(int id, int cell, double r1, double r2, losBuffer data)
{
  float rf;
  int i, j, ret, nb[num_neighbors];
  extProximityZone_Data *d = (extProximityZone_Data *)data.Data;

  cell_all_neighbors(cell,nb);

  ret = 1;
  for(j=0; j<NBINS; j++) if(d->R[j] < 0.0)
    {
      /*
      //  Average over neighbors to reduce numerical noise
      */
      rf = 2*cell_var(cell,rt_disk_offset+j);
      for(i=0; i<num_neighbors; i++)
	{
	  rf += cell_var(nb[i],rt_disk_offset+j);
	}
      rf /= (2+num_neighbors);
      if(rf < 2*extProximityZone_MeanRF[j])
	{
	  d->R[j] = r2;
	}
      else
	{
	  ret = 0;
	}
      //if(j == 0) cart_debug("%le %le %e %e",r1,r2,rf,extProximityZone_MeanRF[j]);
    }

  return ret;
}


void extProximityZone_Collector(losBuffer *result, int num_segments, const losSegment *segments)
{
  int i, j;
  extProximityZone_Data *dr = (extProximityZone_Data *)result->Data;
  extProximityZone_Data *dn;

  for(i=0; i<num_segments; i++)
    {
      dn = (extProximityZone_Data *)segments[i].Buffer.Data;

      for(j=0; j<NBINS; j++) if(dn->R[j] > 0.0)
	{
	  if(dr->R[j] < 0.0) /* first assignment */
	    {
	      dr->R[j] = dn->R[j];
	    }
	  else
	    {
	      dr->R[j] = min(dr->R[j],dn->R[j]);     
	    }
	}
    }
}


void extProximityZones(const char *fname, int resolution_level, int nside, int halo_id, const halo_list *halos)
{
  //const float wbins[NBINS] = { 1000.0, 911.75, 504.25, 227.83 };
  const int lbins[NBINS] = { 12, 5, 3, 1 };

  int npix = 12*nside*nside;
  int j, ipix, ih;
  float mean_rf[NBINS], uKpc = units->length/constants->kpc;
  FILE *f;

  losBuffer *lines;
  extProximityZone_Data *data;

  int nd[NBINS];
  float Ravg[NBINS], Rmin[NBINS], Rmax[NBINS];

  if(halos == NULL)
    {
      cart_debug("No halo file is loaded. Skipping proximity zone search.");
      return;
    }

  /*
  // Open output file
  */
  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(fname,"w");
      if(f == NULL)
	{
	  cart_error("Unable to open output file.");
	}
      fprintf(f,"# id Mass..... Vmax..... Rvir.....  RUVavg... RUVlow... RUVhigh..  RHIavg... RHIlow... RHIhigh..  RHeIavg.. RHeIlow.. RHeIhigh.  RHeIIavg. RHeIIlow. RHeIIhigh\n");
      fprintf(f,"#    Msun      km/s      kpc        kpc       kpc       kpc        kpc       kpc       kpc        kpc       kpc       kpc        kpc       kpc       kpc      \n");
      fprintf(f,"#\n");
    }

  /*
  //  Extract radiation field and store it in the cell_vars array
  */
  //extExtractRadiationField(NBINS,wbins,mean_rf);
  extExtractPhotoRates(NBINS,lbins,mean_rf);
  extProximityZone_MeanRF = mean_rf;

  /*
  //  Prepare buffers for the LOS traversal
  */
  lines = cart_alloc(losBuffer,npix);
  data = cart_alloc(extProximityZone_Data,npix);

#pragma omp parallel for default(none), private(ipix), shared(npix,lines,data)
  for(ipix=0; ipix<npix; ipix++)
    {
      lines[ipix].Size = sizeof(extProximityZone_Data);
      lines[ipix].Data = data + ipix;
    }

  /*
  //  Loop over all halos resoved to the lowest possible level at the center
  */
  for(ih=0; ih<halos->num_halos; ih++) if((halo_id==0 || halo_id==halos->list[ih].id) && halo_level(&halos->list[ih],mpi.comm.run)>=resolution_level)
    {
      /*
      //  Init LOS buffers
      */
#pragma omp parallel for default(none), private(ipix,j), shared(npix,lines,data)
      for(ipix=0; ipix<npix; ipix++)
	{
	  for(j=0; j<NBINS; j++) data[ipix].R[j] = -1;
	}

      losTraverseSky(nside,halos->list[ih].pos,0.5*num_grid,max_level,lines,extProximityZone_Worker,extProximityZone_Collector);

      for(j=0; j<NBINS; j++)
	{
	  nd[j] = 0;
	  Ravg[j] = 0.0;
	  Rmin[j] = num_grid;
	  Rmax[j] = 0.0;
	}

      for(ipix=0; ipix<npix; ipix++)
	{
	  for(j=0; j<NBINS; j++) if(data[ipix].R[j] > 0.0)
	    {
	      nd[j]++;
	      Ravg[j] += data[ipix].R[j];
	      Rmin[j] = min(Rmin[j],data[ipix].R[j]);
	      Rmax[j] = max(Rmax[j],data[ipix].R[j]);
	    }
	}

      for(j=0; j<NBINS; j++) if(nd[j] > 0) Ravg[j] /= nd[j];

      if(local_proc_id == MASTER_NODE)
	{
	  printf("\rHalo #%d proximity zone: %g %g %g\n",halos->list[ih].id,Ravg[0]*uKpc,Ravg[1]*uKpc,Ravg[2]*uKpc);
	  fprintf(f,"%4d %9.3e %9.3e %9.3e  %9.3e %9.3e %9.3e  %9.3e %9.3e %9.3e  %9.3e %9.3e %9.3e  %9.3e %9.3e %9.3e\n",halos->list[ih].id,units->mass/constants->Msun*halos->list[ih].mvir,units->velocity/constants->kms*halos->list[ih].vmax,units->length/constants->kpc*halos->list[ih].rvir,Ravg[0]*uKpc,Rmin[0]*uKpc,Rmax[0]*uKpc,Ravg[1]*uKpc,Rmin[1]*uKpc,Rmax[1]*uKpc,Ravg[2]*uKpc,Rmin[2]*uKpc,Rmax[2]*uKpc,Ravg[3]*uKpc,Rmin[3]*uKpc,Rmax[3]*uKpc);
	}
    }

  if(local_proc_id == MASTER_NODE)
    {
      fclose(f);
    }

  cart_free(lines);
  cart_free(data);
}


#endif  /* RT_TRANSFER && RT_EXTERNAL_BACKGROUND */
#endif  /* RADIATIVE_TRANSFER */


#ifdef HYDRO
void extMassFractions(const char *fname, halo_list *halos)
{
  const int nmass = 7;
  int j, ih, *hlev;
  float *massl[nmass], *mass[nmass];
  float v, dv;
  FILE *f;
  MESH_RUN_DECLARE(level,cell);

  if(halos == NULL)
    {
      cart_debug("No halo file is loaded. Skipping computing gas fraction.");
      return;
    }

  /*
  //  Map cells
  */
  if(halos->map == NULL)
    {
      map_halos(min_level,halos,1.0);
    }

  for(j=0; j<nmass; j++)
    {
      mass[j] = cart_alloc(float,halos->num_halos);
      massl[j] = cart_alloc(float,halos->num_halos);
      for(ih=0; ih<halos->num_halos; ih++) massl[j][ih] = 0.0;
    }

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  cart_debug("Analysing level %d...",level);

  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
  if(cell_is_leaf(cell))
    {
      ih = halos->map[cell];
      cart_assert(ih>=0 && ih<=halos->num_halos);
      if(ih > 0) 
	{
	  ih--;
	  massl[0][ih] += cell_total_mass(cell)*cell_volume[level];
	  massl[1][ih] += cell_gas_density(cell)*cell_volume[level];
#ifdef RADIATIVE_TRANSFER
	  massl[2][ih] += cell_HI_density(cell)*cell_volume[level];
	  massl[3][ih] += cell_HII_density(cell)*cell_volume[level];
	  massl[4][ih] += 2*cell_H2_density(cell)*cell_volume[level];

	  /*
	  //  HI line width
	  */
	  dv = 0.0;
	  for(j=0; j<nDim; j++)
	    {
	      v = cell_momentum(cell,j)/cell_gas_density(cell) - halos->list[ih].vel[j];
	      dv += v*v;
	    } 
	  massl[6][ih] += cell_HI_density(cell)*cell_volume[level]*dv;
#endif /* RADIATIVE_TRANSFER */
	}
    }

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  /*
  // Measure stellar masses
  */
#ifdef STARFORM
  for(j=0; j<num_particles; j++) if(particle_is_star(j))
    {
      cell = cell_find_position(particle_x[j]);
      if(cell > -1)
	{
	  ih = halos->map[cell];
	  cart_assert(ih>=0 && ih<=halos->num_halos);
	  if(ih > 0) 
	    {
	      ih--;
	      massl[5][ih] += particle_mass[j];
	    }
	}
    }
#endif

  for(j=0; j<nmass; j++)
    {
      MPI_Allreduce(massl[j],mass[j],halos->num_halos,MPI_FLOAT,MPI_SUM,mpi.comm.run);
    }

  for(j=0; j<nmass; j++)
    {
      cart_free(massl[j]);
    }

  hlev = cart_alloc(int,halos->num_halos);
  for(ih=0; ih<halos->num_halos; ih++) 
    {
      hlev[ih] = halo_level(&halos->list[ih],mpi.comm.run);
    }

  /*
  // Open output file
  */
  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(fname,"w");
      if(f == NULL)
	{
	  cart_error("Unable to open output file.");
	}
      fprintf(f,"# id Mvir..... Mtot..... Mgas..... MHI...... MHII..... MH2...... Mstars... Level.... Vmax..... sigHI.... (all masses are in Msun, velocities are in km/s)\n");
      fprintf(f,"#\n");

      for(ih=0; ih<halos->num_halos; ih++) 
	{
	  fprintf(f,"%4d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %2d %9.3e %9.3e\n",halos->list[ih].id,units->mass/constants->Msun*halos->list[ih].mvir,mass[0][ih]*units->mass/constants->Msun,mass[1][ih]*units->mass/constants->Msun,mass[2][ih]*units->mass/constants->Msun,mass[3][ih]*units->mass/constants->Msun,mass[4][ih]*units->mass/constants->Msun,mass[5][ih]*units->mass/constants->Msun,hlev[ih],units->velocity/constants->kms*halos->list[ih].vmax,mass[2][ih]>0.0?sqrt(mass[6][ih]/mass[2][ih])*units->velocity/constants->kms:0.0);
	}
      
      fclose(f);
    }

  cart_free(hlev);
  for(j=0; j<nmass; j++)
    {
      cart_free(mass[j]);
    }
}
#endif /* HYDRO */

#endif /* COSMOLOGY */

