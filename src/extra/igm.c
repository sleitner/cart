#include "config.h"

#ifdef COSMOLOGY

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "iterators.h"
#include "parallel.h"
#include "rt_utilities.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "los.h"
#include "halo_finder.h"
#include "ifrit.h"
#include "ism.h"


#ifdef RADIATIVE_TRANSFER

#include "rt_solver.h"
#include "F/frt_c.h"


void extDumpRadiationBackground(const char *filename)
{
  int i;
  FILE *f;
  float w, nxi, uJnu = 1.5808e-17/pow(auni[min_level],3.0);

  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(filename,"w");
      cart_assert(f != NULL);

      fprintf(f,"# wlen[A]  jnu[cgs] at a=%lf\n",auni[min_level]);
      i = 0;
      do
	{
	  i++;
	  rtGetBinWavelengths(1,&i,&w);
	  rtGetRadiationField(-1,1,&i,&nxi);
	  if(w > 0.0) fprintf(f,"%9.3e %9.3e\n",w,uJnu*nxi);
	}
      while(w > 0.0);

      fclose(f);
    }
}


void extExtractRadiationField(int nbins, const float wbins[], float *mean_rf)
{
  int i;
  int lbins[rt_num_frequencies], vars[rt_num_frequencies];
  float nxi[rt_num_frequencies];
  MESH_RUN_DECLARE(level,cell);
  float uJnu = 1.5808e-17/pow(auni[min_level],3.0);

  cart_assert(nbins>0 && nbins<=rt_num_frequencies);

  /*
  //  Find the frequency bin indecies
  */
  rtGetBinIds(nbins,wbins,lbins);
  
  /*
  //  Find the mean field
  */
  if(mean_rf != NULL)
    {
      rtGetRadiationField(-1,nbins,lbins,mean_rf);
      if(local_proc_id == MASTER_NODE)
	{
	  for(i=0; i<nbins; i++)
	    {
	      cart_debug("Selecting bin %d for wavelength %f; BG = %e [cgs] = %e [CU]",lbins[i],wbins[i],uJnu*mean_rf[i],mean_rf[i]);
	    }
	}
    }
  else
    {
      if(local_proc_id == MASTER_NODE)
	{
	  for(i=0; i<nbins; i++)
	    {
	      cart_debug("Selecting bin %d for wavelength %f ",lbins[i],wbins[i]);
	    }
	}
    }

  for(i=0; i<nbins; i++)
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
#pragma omp parallel for default(none), private(_Index,cell,i,nxi), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,vars,nbins,lbins)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
      rtGetRadiationField(cell,nbins,lbins,nxi);
      for(i=0; i<nbins; i++)
	{
	  cell_var(cell,vars[i]) = nxi[i];
	}

    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  
  update_buffer_level(level,vars,nbins);

  MESH_RUN_OVER_LEVELS_END;
}  


void extExtractPhotoRates(int nbins, const int lbins[], float mean_rate[])
{
  int i;
  int vars[rt_num_frequencies];
  float rate[frtRATE_DIM], rate0[frtRATE_DIM];
  MESH_RUN_DECLARE(level,cell);

  cart_assert(nbins>0 && nbins<=rt_num_frequencies);

  /*
  //  Find the mean rates
  */
  rtGetPhotoRates(-1,rate0);
  for(i=0; i<nbins; i++)
    {
      mean_rate[i] = rate0[lbins[i]];
    }

  if(local_proc_id == MASTER_NODE)
    {
      for(i=0; i<nbins; i++)
	{
	  cart_debug("Mean rate[%d] = %e [cgs]",i,mean_rate[i]);
	}
    }

  for(i=0; i<nbins; i++)
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
#pragma omp parallel for default(none), private(_Index,cell,i,rate), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,vars,nbins,lbins)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
      rtGetPhotoRates(cell,rate);
      for(i=0; i<nbins; i++)
	{
	  cell_var(cell,vars[i]) = rate[lbins[i]];
	}
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  
  update_buffer_level(level,vars,nbins);

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
  for(ih=0; ih<halos->num_halos; ih++) if((halo_id==0 || halo_id==halos->list[ih].id) && halo_level(&halos->list[ih],MPI_COMM_WORLD)>=resolution_level)
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
void extGasFractions(const char *fname, halo_list *halos)
{
  const int nmass = 5;
  int j, ih, *hlev;
  float *massl[nmass], *mass[nmass];
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
  if(halos->map == -1)
    {
      map_halos(VAR_ACCEL,min_level,halos,1.0);
    }

#ifndef RADIATIVE_TRANSFER
  for(level=min_level; level<=max_level; level++)
    {
      cart_debug("assigning density on level %u",level);
      assign_density(level);
    }
#endif /* RADIATIVE_TRANSFER */

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
      ih = (int)(0.5+cell_var(cell,VAR_ACCEL+0));
      cart_assert(ih>=0 && ih<=halos->num_halos);
      if(ih > 0) 
	{
	  ih--;
	  massl[0][ih] += (cell_volume[level]+cell_density(cell));
	  massl[1][ih] += cell_gas_density(cell)*cell_volume[level];
#ifdef RADIATIVE_TRANSFER
	  massl[2][ih] += cell_HI_density(cell)*cell_volume[level];
	  massl[3][ih] += cell_HII_density(cell)*cell_volume[level];
	  massl[4][ih] += 2*cell_H2_density(cell)*cell_volume[level];
#endif /* RADIATIVE_TRANSFER */
	}
    }

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  for(j=0; j<nmass; j++)
    {
      MPI_Allreduce(massl[j],mass[j],halos->num_halos,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    }

  for(j=0; j<nmass; j++)
    {
      cart_free(massl[j]);
    }

  hlev = cart_alloc(int,halos->num_halos);
  for(ih=0; ih<halos->num_halos; ih++) 
    {
      hlev[ih] = halo_level(&halos->list[ih],MPI_COMM_WORLD);
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
      fprintf(f,"# id Mvir..... Mtot..... Mgas..... MHI...... MHII..... MH2...... Level (all masses are in Msun)\n");
      fprintf(f,"#\n");

      for(ih=0; ih<halos->num_halos; ih++) 
	{
	  fprintf(f,"%4d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %2d %9.3e\n",halos->list[ih].id,units->mass/constants->Msun*halos->list[ih].mvir,mass[0][ih]*units->mass/constants->Msun,mass[1][ih]*units->mass/constants->Msun,mass[2][ih]*units->mass/constants->Msun,mass[3][ih]*units->mass/constants->Msun,mass[4][ih]*units->mass/constants->Msun,hlev[ih],units->velocity/constants->kms*halos->list[ih].vmax);
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

