#include "config.h"


#include <stdio.h>
#include <string.h>


#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "density.h"
#include "hydro.h"
#include "io.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "plugin.h"
#include "rt_solver.h"
#include "system.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "extra/halo_finder.h"
#include "extra/ifrit.h"
#include "extra/igm.h"
#include "extra/ism.h"
#include "extra/rfsfr.h"
#include "extra/ng_analysis.h"

#ifdef RADIATIVE_TRANSFER
#include "F/frt_c.h"
#endif


halo_list *ngHalos = NULL;


char *ngOutputDir = NULL;

const char *ngOutputFile(const char *filepath)
{
  int ret;
  static char str[999];
  char *name;

  if(ngOutputDir == NULL)
    {
      strcpy(str,filepath);
    }
  else
    {
      strcpy(str,ngOutputDir);
      strcat(str,"/");
      strcat(str,filepath);
    }

  name = strrchr(str,'/');
  if(name != NULL)
    {
      *name = 0;
      /*
      //  Create output directory
      */
      ret = system_mkdir(str);
      if(ret != 0)
	{
	  cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",str,ret);
	  return NULL;
	}
      *name = '/';
    }

  return str;
}


int ngIsDirty[num_vars] = { 0 };

#ifdef RADIATIVE_TRANSFER
int ngVarIds[] = { I_FRACTION+RT_HVAR_OFFSET+0, I_GAS_NUMBER_DENSITY, I_GAS_TEMPERATURE, I_FRACTION+RT_HVAR_OFFSET+5, RT_HVAR_OFFSET+5, I_CELL_LEVEL, I_LOCAL_PROC };
#else
#ifdef HYDRO
int ngVarIds[] = { I_GAS_NUMBER_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC, I_GAS_TOVERMU };
#else
int ngVarIds[] = { VAR_DENSITY, I_CELL_LEVEL, I_LOCAL_PROC };
#endif /* HYDRO */
#endif /* RADIATIVE_TRANSFER */
ng_var_list_t ngVarList = { sizeof(ngVarIds)/sizeof(int), ngVarIds };


/*
//  Initialize what is needed
*/
int ngInitRTSet = 0;
void ngInitRT()
{
#ifdef RADIATIVE_TRANSFER
  if(!ngInitRTSet)
    {
      ngInitRTSet = 1;

      cart_debug("Initializing RT...");
      rtStepBegin();
      rtUpdateTables(min_level,MPI_COMM_WORLD);
    }
#endif
}


int ngInitDensitySet = 0;
void ngInitDensity()
{
  int level;

  if(!ngInitDensitySet)
    {
      ngInitDensitySet = 1;

      for(level=min_level; level<=max_level; level++)
	{
	  cart_debug("assigning density on level %u", level );
	  assign_density( level );
	}
    }
}


#ifdef RADIATIVE_TRANSFER
const char *ngDumpHeader[] = {
  "baryon number density (cm^{-3})",
  "temperature (K)",
  "baryon column density (cm^{-2})",
  "gas metallicity (solar units)",
  "HI fraction",
  "HII fraction",
  "H_2 fraction",
  "radiation field in MW units (U_{MW})",
  "dust-to-gas ratio in MW units (D_{MW}",
  "density in this cell over the average density in neighboring cells",
  "cooling rate (erg/s)",
  "heating rate (erg/s)"
};
#else
const char *ngDumpHeader[] = {
  "baryon number density (cm^{-3})",
  "temperature (K)",
  "baryon column density (cm^{-2})",
  "gas metallicity (solar units)",
  "total density (g/cm^3)",
};
#endif


const int ngDumpSize = sizeof(ngDumpHeader)/sizeof(const char *);


#ifdef RADIATIVE_TRANSFER
const int ngDumpWeight[] = {
  0,
  2,
  2,
  2,
  2,
  2,
  2,
  0,
  2,
  0,
  2,
  2
};
#else
const int ngDumpWeight[] = {
  0,
  2,
  2,
  2,
  0
};
#endif


void ngDumpWorker(int level, int cell, int num, float *ptr)
{
  float s, soblen, sobvel, cr, hr;
  int i, nb[num_neighbors];

#ifdef RADIATIVE_TRANSFER
  cart_assert(num >= 12);
#else
  cart_assert(num >= 5);
#endif

#ifdef HYDRO
  ptr[0] = units->number_density*cell_gas_density(cell);
  ptr[1] = units->temperature*cell_gas_temperature(cell);

  rtGetSobolevFactors(cell,level,&soblen,&sobvel);
  ptr[2] = cell_gas_density(cell)*soblen*units->number_density*units->length;
#else
  soblen = sobvel = 0.0;
  ptr[0] = ptr[1] = ptr[2] = 0.0;
#endif /* HYDRO */

#ifdef ENRICH
  ptr[3] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  ptr[3] = 0.0;
#endif /* ENRICH */

#ifdef RADIATIVE_TRANSFER
  ptr[4] = cell_HI_fraction(cell);
  ptr[5] = cell_HII_fraction(cell);
  ptr[6] = cell_H2_fraction(cell);
  ptr[7] = rtUmw(cell);
  ptr[8] = rtDmw(cell);

  cell_all_neighbors(cell,nb);
  /*
  //  Mean of all neighbors
  */
  s = 0.0;
  for(i=0; i<num_neighbors; i++)
    {
      s += cell_gas_density(nb[i]);
    }
  s /= num_neighbors;
  ptr[9] = cell_gas_density(cell)/s;

  rtGetCoolingRate(cell,&cr,&hr);
  ptr[10] = cr;
  ptr[11] = hr;

#else  /* RADIATIVE_TRANSFER */

  ptr[4] = units->density*(cell_density(cell)*cell_volume_inverse[level]+1.0)/(cgs->g/pow(cgs->cm,3.0));

#endif /* RADIATIVE_TRANSFER */
}


/*
//  Load hlist file for other command to use
*/
void ngLoadHalos(const char *dirname)
{
#ifdef COSMOLOGY
  char str[999];

  if(dirname == NULL)
    {
      if(ngHalos != NULL)
	{
	  destroy_halo_list(ngHalos);
	  ngHalos = NULL;
	}
    }
  else
    {
      sprintf(str,"%s/hlist_%6.4f.dat",dirname,auni[min_level]);
      ngHalos = load_halo_finder_catalog(str,50,0.0,0.0,0.0);
      if(ngHalos == NULL)
	{
	  cart_error("Failed to read in the halo list from file %s",str);
	}
    }
#else
  cart_debug("COSMOLOGY is not set. Skipping ngLoadHalos.");
#endif /* COSMOLOGY */
}


void ngDumpLevels(const char *filename, int level)
{
  ngInitRT();
  ngInitDensity();
  extDumpLevels(ngOutputFile(filename),ngDumpSize,ngDumpWorker,ngDumpHeader,level,max_level,ngHalos);
  ngIsDirty[VAR_ACCEL] = 1;
}


void ngHaloParticles(const char *filename, int id, float rmax)
{
#ifdef COSMOLOGY
  if(ngHalos == NULL)
    {
      cart_debug("Halos are not loaded. Skipping ngDumpHaloParticles.");
      return;
    }

  if(find_halo_by_id(ngHalos,id) == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ngDumpHaloParticles.",id);
    }
  else
    {
      dump_region_around_halo(ngOutputFile(filename),find_halo_by_id(ngHalos,id),rmax);
    }
#else
  cart_debug("COSMOLOGY is not set. Skipping ngDumpHaloParticles.");
#endif /* COSMOLOGY */
}


void ngDumpProfiles(const char *filename, int resolution_level, float rmin, float rmax)
{
#ifdef COSMOLOGY
  ngInitRT();
  ngInitDensity();
  extDumpHaloProfiles(ngOutputFile(filename),ngDumpSize,ngDumpWorker,ngDumpHeader,ngDumpWeight,rmin,rmax,25,ngHalos,resolution_level,1.0);
  ngIsDirty[VAR_ACCEL] = 1;
#else
  cart_debug("COSMOLOGY is not set. Skipping ngDumpProfiles.");
#endif /* COSMOLOGY */
}


void ngGasFractions(const char *filename)
{
#if defined(COSMOLOGY) && defined(HYDRO)
  ngInitDensity();
  extGasFractions(ngOutputFile(filename),ngHalos);
  ngIsDirty[VAR_ACCEL] = 1;
#else
  cart_debug("COSMOLOGY and HYDRO are not set. Skipping ngGasFractions.");
#endif /* COSMOLOGY && HYDRO */
}


void ngRegion2IFrIT(const char *filename, int nbin1, int floor_level, const double pos[], const ng_var_list_t *list)
{
  int i, nbin[] = { nbin1, nbin1, nbin1 };

  if(list == NULL) list = &ngVarList;
  for(i=0; i<list->Num; i++) if(list->Ids[i] == VAR_DENSITY)
    {
      ngInitDensity();
    }

  ifrit.OutputBox(ngOutputFile(filename),floor_level,nbin,pos,list->Num,list->Ids);
}


void ngHalo2IFrIT(const char *filename, int floor_level, int id, const ng_var_list_t *list)
{
  const int nbin1 = 256;
  int i, nbin[] = { nbin1, nbin1, nbin1 };

  if(ngHalos == NULL)
    {
      cart_debug("Halos are not loaded. Skipping ngHalo2IFrIT.");
      return;
    }
  if(find_halo_by_id(ngHalos,id) == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ngHalo2IFrIT.",id);
    }
  else
    {
      if(list == NULL) list = &ngVarList;
      for(i=0; i<list->Num; i++) if(list->Ids[i] == VAR_DENSITY)
	{
	  ngInitDensity();
	}
      ifrit.OutputHalo(ngOutputFile(filename),floor_level,nbin,find_halo_by_id(ngHalos,id),list->Num,list->Ids);
    }
}


void ngProximityZone(const char *filename, int id, int resolution_level)
{
#if defined(COSMOLOGY) && defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER) && defined(RT_EXTERNAL_BACKGROUND)
  extProximityZones(ngOutputFile(filename),resolution_level,2,id,ngHalos);
#else
  cart_debug("COSMOLOGY, RADIATIVE_TRANSFER, RT_TRANSFER, and RT_EXTERNAL_BACKGROUND are not set. Skipping ngProximityZone.");
#endif
}


void ngRFvsSFR(const char *fileroot, int top_level)
{
#if defined(RADIATIVE_TRANSFER) && defined(HYDRO) && defined(STARFORM)
  extRFvsSFR(fileroot,top_level,ngHalos);
#else
  cart_debug("RADIATIVE_TRANSFER, HYDRO, and STARFORM are not set. Skipping ngRFvsSFR.");
#endif
}

 
void ngStarFormationLaw(const char *filename, float spatial_scale, float time_scale, float stellar_age_limit)
{
#if defined(PARTICLES) && defined(STARFORM)
  ngInitRT();
  extStarFormationLaw(ngOutputFile(filename),spatial_scale,time_scale,stellar_age_limit,ngHalos);
#else
  cart_debug("PARTICLES and STARFORM are not set. Skipping ngStarFormationLaw.");
#endif
}


void ngStarFormationLaw2(const char *filename, float spatial_scale)
{
#if defined(HYDRO) && defined(STARFORM)
  ngInitRT();
  extStarFormationLaw2(ngOutputFile(filename),spatial_scale,ngHalos);
#else
  cart_debug("HYDRO and STARFORM are not set. Skipping ngStarFormationLaw.");
#endif
}


void ngHaloStars(const char *filename, int id, float rmax)
{
  if(ngHalos == NULL)
    {
      cart_debug("Halos are not loaded. Skipping ngHaloStars.");
      return;
    }

  if(find_halo_by_id(ngHalos,id) == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ngHaloStars.",id);
    }
  else
    {
#if defined(PARTICLES) && defined(STARFORM)
      extHaloStars(filename,find_halo_by_id(ngHalos,id),rmax);
#endif
    }
}


/* ------------------------------------------------------------------------------- */

void ngAnalysis();


void run_output()
{
  if(local_proc_id == MASTER_NODE)
    {
      cart_error("This is the analysis utility. It must be run in a non-restart mode.");
    }
}


#ifdef USER_PLUGIN
plugin_t ngPlugin = { NULL };
const plugin_t* add_plugin(int id)
{
  if(id == 0)
    {
      ngPlugin.RunBegin = run_output;
      return &ngPlugin;
    }
  else return NULL;
}
#endif


void init_run()
{
  int i, ret, level;
  double aexp;
  char str[999];

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Entering analysis mode...");
    }

  for(i=0; i<num_options; i++)
    {
      aexp = -1.0;
      if(sscanf(options[i],"%lf%c",&aexp,str)!=1 || aexp<1.0e-10 || aexp>1.1)
	{
	  cart_debug("Invalid scale factor value %s; skipping option #%d",options[i],i);
	  continue;
	}
      
      if(buffer_enabled) destroy_cell_buffer();
#ifdef PARTICLES
      init_particles();
#endif /* PARTICLES */

      ngLoadHalos(NULL);

      ngInitRTSet = 0;
      ngInitDensitySet = 0;

      strcpy(str+1,options[i]);
      str[0] = 'a';
      read_restart(str);
      load_balance();

#ifdef COSMOLOGY
      abox[min_level] = abox_from_tcode(tl[min_level]);
      auni[min_level] = auni_from_tcode(tl[min_level]);
#endif /* COSMOLOGY */

      for(level=min_level+1; level<=max_level; level++)
	{
	  tl[level] = tl[min_level];
	  dtl[level] = 0.5*dtl[level-1];
#ifdef COSMOLOGY
	  abox[level] = abox[min_level];
	  auni[level] = auni[min_level];
#endif /* COSMOLOGY */
	}

      if(num_options > 1)
	{
	  strcpy(str+2,options[i]);
	  str[0] = 'a';
	  str[1] = '=';

	  /*
	  //  Create output directory
	  */
	  ret = system_mkdir(str);
	  if(ret != 0)
	    {
	      cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",str,ret);
	      return;
	    }
	  ngOutputDir = str;
	}
	  
      ngAnalysis();

      if(num_options > 1)
	{
	  ngOutputDir = NULL;
	}
    }

  MPI_Finalize();
  exit(0);
}


