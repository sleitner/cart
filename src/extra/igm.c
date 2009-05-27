#include "defs.h"

#include <stdio.h>
#include <mpi.h>

#include "auxiliary.h"
#include "parallel.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "los.h"
#include "halo_finder.h"

#include "rt_utilities.h"


#ifdef RADIATIVE_TRANSFER

#include "rt_solver.h"


void extDumpRadiationBackground()
{
  int i, n;
  float *wlen, *ngxi;
  FILE *f;
  float uJnu = 1.5808e-17/pow(auni[min_level],3.0);

  if(local_proc_id == MASTER_NODE)
    {
      rtGetRadiationBackground(&n,&wlen,&ngxi);

      f = fopen("jnu.res","w");
      cart_assert(f != NULL);

      fprintf(f,"# wlen[A]  jnu[cgs] at a=%lf\n",auni[min_level]);
      for(i=0; i<n; i++)
	{
	  fprintf(f,"%9.3e %9.3e\n",wlen[i],uJnu*ngxi[i]);
	}

      fclose(f);

      cart_free(wlen);
      cart_free(ngxi);
    }
}


void extExtractRadiationField(int nbins, const float wbins[], float *mean_rf)
{
  int i, j, nbg;
  float *wbg, *fbg;
  int lbins[rt_num_frequencies];
  float nxi[rt_num_frequencies];
  MESH_RUN_DECLARE(level,cell);

  cart_assert(nbins>0 && nbins<=rt_num_disk_vars);

  /*
  // Init internal data
  */
  rtStepBegin();

  /*
  //  Find the frequency bin indecies
  */
  rtGetRadiationBackground(&nbg,&wbg,&fbg);

  for(i=0; i<nbins; i++)
    {
      lbins[i] = -1;
      for(j=0; lbins[i]==-1 && j<nbg; j++)
	{
	  if(wbg[j] < wbins[i]) lbins[i] = j;
	}
      if(mean_rf!=NULL && lbins[i]>-1)
	{
	  mean_rf[i] = fbg[lbins[i]];
	}
      if(local_proc_id == MASTER_NODE)
	{
	  cart_debug("Selecting bin %d for wavelength %f",lbins[i],wbins[i]);
	}
    }
  
  cart_free(wbg);
  cart_free(fbg);

  /*
  //  Extract radiation field and store it in the cell_vars array
  */
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Extracting radiation field on level %d...",level);
    }
#pragma omp parallel for default(none), private(_Index,cell,i,nxi), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,nbins,lbins)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
      rtGetRadiationField(cell,nbins,lbins,nxi);
      for(i=0; i<nbins; i++)
	{
	  cell_var(cell,rt_disk_offset+i) = nxi[i];
	}

    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
}  


#ifdef RT_TRANSFER

#include "rt_tree.h"

/*
//  Find proximity zones around galaxies
*/
float *extProximityZone_MeanRF;

typedef struct extProximityZone_DataType
{
  float R[3];
}
extProximityZone_Data;


int extProximityZone_Worker(int id, int cell, double r1, double r2, losBuffer data)
{
  int j;
  extProximityZone_Data *d = (extProximityZone_Data *)data.Data;

  for(j=0; j<3; j++)
    {
      if(d->R[j]<0.0 && cell_var(cell,rt_disk_offset+j)<2*extProximityZone_MeanRF[j])
	{
	  d->R[j] = 0.5*(r1+r2);
	}
    }

  for(j=0; j<3; j++)
    {
      if(d->R[j] < 0.0) break;
    }

  return (j==3 ? 1 : 0);
}


void extProximityZone_Collector(losBuffer *result, losSegment *next)
{
  int j;
  extProximityZone_Data *dr = (extProximityZone_Data *)result->Data;
  extProximityZone_Data *dn = (extProximityZone_Data *)next->Buffer.Data;

   for(j=0; j<3; j++)
    {
      if(dr->R[j]<0.0 && dn->R[j]>0.0) /* segments are ordered */
	{
	  dr->R[j] = dn->R[j];
	}
    }
}


void extFindProximityZones(const char *fname, int nside, int halo_id)
{
  const float wbins[] = { 911.75, 504.25, 227.83 };
  const int nbins = sizeof(wbins)/sizeof(float);

  int npix = 12*nside*nside;
  int j, ipix, ih, num_halos;
  int floor_level;
  float mean_rf[nbins], uKpc = 1.0e3*r0;
  FILE *f;

  hfHalo *halos;
  losBuffer *lines;
  extProximityZone_Data *data;

  int nd[nbins];
  float Ravg[nbins], Rmin[nbins], Rmax[nbins];

  /*
  //  Find the lowest level on the mesh
  */
  floor_level = max_level_now_global(MPI_COMM_WORLD);

  /*
  //  Read the halo data from the hlist file
  */
  if(hfReadHFINDHalos(fname,100,0.0,0.0,0.0,&halos,&num_halos) != 0)
    {
      cart_error("Failed to read in the halo list from file %s",fname);
    }

  /*
  // Open output file
  */
  if(local_proc_id == MASTER_NODE)
    {
      f = fopen("pz.res","w");
      if(f == NULL)
	{
	  cart_error("Unable to open output file.");
	}
      fprintf(f,"# id Mass..... Vmax..... Rvir.....  RHIavg... RHIlow... RHIhigh..  RHeIavg.. RHeIlow.. RHeIhigh.  RHeIIavg. RHeIIlow. RHeIIhigh\n");
      fprintf(f,"#    Msun      km/s      chik       chik      chik      chik       chik      chik      chik       chik      chik      chik     \n");
      fprintf(f,"#\n");
    }

  /*
  //  Extract radiation field and store it in the cell_vars array
  */
  extExtractRadiationField(nbins,wbins,mean_rf);
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
  for(ih=0; ih<num_halos; ih++) if((halo_id==0 || halo_id==halos[ih].Id) && hfHaloLevel(halos+ih)==floor_level)
    {

      if(local_proc_id == MASTER_NODE)
	{
	  printf("Halo #%d (%g Msun): analysing...",halos[ih].Id,halos[ih].Mass);
	}

      /*
      //  Init LOS buffers
      */
#pragma omp parallel for default(none), private(ipix,j), shared(npix,lines,data,nbins)
      for(ipix=0; ipix<npix; ipix++)
	{
	  for(j=0; j<nbins; j++) data[ipix].R[j] = -1;
	}

      losTraverseSky(nside,halos[ih].Pos,0.5*num_grid,floor_level,lines,extProximityZone_Worker,extProximityZone_Collector);

      for(j=0; j<nbins; j++)
	{
	  nd[j] = 0;
	  Ravg[j] = 0.0;
	  Rmin[j] = num_grid;
	  Rmax[j] = 0.0;
	}

      for(ipix=0; ipix<npix; ipix++)
	{
	  for(j=0; j<nbins; j++) if(data[ipix].R[j] > 0.0)
	    {
	      nd[j]++;
	      Ravg[j] += data[ipix].R[j];
	      Rmin[j] = min(Rmin[j],data[ipix].R[j]);
	      Rmax[j] = max(Rmax[j],data[ipix].R[j]);
	    }
	}

      for(j=0; j<nbins; j++) if(nd[j] > 0) Ravg[j] /= nd[j];

      if(local_proc_id == MASTER_NODE)
	{
	  printf("\rHalo #%d proximity zone: %g %g %g\n",halos[ih].Id,Rmin[0]*uKpc,Ravg[0]*uKpc,Rmax[0]*uKpc);
	  fprintf(f,"%4d %9.3e %9.3e %9.3e  %9.3e %9.3e %9.3e  %9.3e %9.3e %9.3e  %9.3e %9.3e %9.3e\n",halos[ih].Id,halos[ih].Mass,halos[ih].Vmax,halos[ih].Rvir,Ravg[0]*uKpc,Rmin[0]*uKpc,Rmax[0]*uKpc,Ravg[1]*uKpc,Rmin[1]*uKpc,Rmax[1]*uKpc,Ravg[2]*uKpc,Rmin[2]*uKpc,Rmax[2]*uKpc);
	}
    }

  if(local_proc_id == MASTER_NODE)
    {
      fclose(f);
    }

  cart_free(halos);
  cart_free(lines);
  cart_free(data);
}


#endif  /* RT_TRANSFER */
#endif  /* RADIATIVE_TRANSFER */
