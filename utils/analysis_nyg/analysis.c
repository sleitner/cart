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
#include "rt.h"
#include "system.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "extra/start_analysis.h"

#include "extra/halo_finder.h"
#include "extra/ifrit.h"
#include "extra/igm.h"
#include "extra/ism.h"
#include "extra/rfsfr.h"

#ifdef RADIATIVE_TRANSFER
#include "F/frt_c.h"
#endif

#include "analysis.h"
#include "cell_dump.h"


struct 
{
  halo_list *Halos;
  char *OutputDir;
  int IsDirty[num_vars];
}
ngGlobal = { NULL, NULL, { 0 } };

extern struct ngCellDump ngCD_p_;


const struct HALO_LIST *ngHalos()
{
  return ngGlobal.Halos;
}


const char *ngOutputFile(const char *filepath)
{
  int ret;
  static char str[999];
  char *name;

  if(ngGlobal.OutputDir == NULL)
    {
      strcpy(str,filepath);
    }
  else
    {
      strcpy(str,ngGlobal.OutputDir);
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


#ifdef RADIATIVE_TRANSFER
const int ngVI_p_[] = { I_HI_FRACTION, I_GAS_NUMBER_DENSITY, I_GAS_TEMPERATURE, I_FRACTION+RT_HVAR_OFFSET+5, RT_HVAR_OFFSET+5, I_CELL_LEVEL, I_LOCAL_PROC 
#ifdef RT_UV
		   , rt_field_offset+rt_num_freqs-1
#endif
};
#else
#ifdef HYDRO
const int ngVI_p_[] = { I_GAS_NUMBER_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC, I_GAS_TOVERMU };
#else
const int ngVI_p_[] = { VAR_TOTAL_DENSITY, I_CELL_LEVEL, I_LOCAL_PROC };
#endif /* HYDRO */
#endif /* RADIATIVE_TRANSFER */


/*
//  Load hlist file for other commands to use
*/
void ngLH_p_(const char *path)
{
#ifdef COSMOLOGY
  char str[999];

  if(path == NULL)
    {
      if(ngGlobal.Halos != NULL)
	{
	  destroy_halo_list(ngGlobal.Halos);
	  ngGlobal.Halos = NULL;
	}
    }
  else
    {
      sprintf(str,"%s/hlist_%6.4f.dat",path,auni[min_level]);
      ngGlobal.Halos = load_halo_finder_catalog(str,ngLoadHalos.Nmin,ngLoadHalos.Mvir,ngLoadHalos.Vmax,ngLoadHalos.Rvir);
      if(ngGlobal.Halos == NULL)
	{
	  cart_error("Failed to read in the halo list from file %s",str);
	}
    }
#else
  cart_debug("COSMOLOGY is not set. Skipping ngLoadHalos.");
#endif /* COSMOLOGY */
}
struct NG_OBJECT_LH_TYPE ngLoadHalos = { ngLH_p_, 50, 0.0, 0.0, 0.0 };


void ngDL_p_(const char *filename)
{
  if(ngDumpLevels.Dump == NULL) ngDumpLevels.Dump = &ngCD_p_;
  if(ngDumpLevels.Level < min_level) ngDumpLevels.Level = max_level_now_global(mpi.comm.run);
  if(ngDumpLevels.MaxLevel < ngDumpLevels.Level) ngDumpLevels.MaxLevel = ngDumpLevels.Level;

  init->All();

  //  extDumpLevels(ngOutputFile(filename),ngDumpLevels.Dump->Size,ngDumpLevels.Dump->Worker,ngDumpLevels.Dump->Header,ngDumpLevels.Level,ngDumpLevels.MaxLevel,ngGlobal.Halos);
  ngGlobal.IsDirty[VAR_ACCEL] = 1;
}
struct NG_OBJECT_DL_TYPE ngDumpLevels = { ngDL_p_, -1, -1, NULL };


void ngDP_p_(const char *filename)
{
#ifdef COSMOLOGY
  if(ngDumpProfiles.Dump == NULL) ngDumpProfiles.Dump = &ngCD_p_;

  init->All();

  //  extDumpHaloProfiles(ngOutputFile(filename),ngDumpProfiles.Dump->Size,ngDumpProfiles.Dump->Worker,ngDumpProfiles.Dump->Header,ngDumpProfiles.Dump->Weight,ngDumpProfiles.Rmin,ngDumpProfiles.Rmax,ngDumpProfiles.Ndex,ngGlobal.Halos,ngDumpProfiles.Level,ngDumpProfiles.HaloEdge);
  ngGlobal.IsDirty[VAR_ACCEL] = 1;
#else
  cart_debug("COSMOLOGY is not set. Skipping ngDumpProfiles.");
#endif /* COSMOLOGY */
}
struct NG_OBJECT_DP_TYPE ngDumpProfiles = { ngDP_p_, 25, max_level, 0.1, 100.0, 1.0, NULL };


void ngMF_p_(const char *filename)
{
#if defined(COSMOLOGY) && defined(HYDRO)
  init->TotalDensity();

  extMassFractions(ngOutputFile(filename),ngGlobal.Halos);

  ngGlobal.IsDirty[VAR_ACCEL] = 1;
#else
  cart_debug("COSMOLOGY and HYDRO are not set. Skipping ngMassFractions.");
#endif /* COSMOLOGY && HYDRO */
}
struct NG_OBJECT_MF_TYPE ngMassFractions = { ngMF_p_ };


void ngRI_p_(const char *filename)
{
  int i, nbin[] = { ngRegion2IFrIT.Nbin, ngRegion2IFrIT.Nbin, ngRegion2IFrIT.Nbin };

  if(ngRegion2IFrIT.Vars == NULL) ngRegion2IFrIT.Vars = (const int *)ngVI_p_;

  for(i=0; i<ngRegion2IFrIT.NumVars; i++)
    {
      if(ngRegion2IFrIT.Vars[i] == VAR_TOTAL_DENSITY) init->TotalDensity();
      if(ngRegion2IFrIT.Vars[i] == VAR_STELLAR_DENSITY) init->StellarDensity();
    }

  ifrit.OutputBox(ngOutputFile(filename),ngRegion2IFrIT.Level,nbin,ngRegion2IFrIT.Pos,ngRegion2IFrIT.NumVars,ngRegion2IFrIT.Vars);
}
struct NG_OBJECT_RI_TYPE ngRegion2IFrIT = { ngRI_p_, num_grid, min_level, { 0.5*num_grid, 0.5*num_grid, 0.5*num_grid }, sizeof(ngVI_p_)/sizeof(int), NULL };


void ngHI_p_(const char *filename)
{
  int i, nbin[] = { ngHalo2IFrIT.Nbin, ngHalo2IFrIT.Nbin, ngHalo2IFrIT.Nbin };
  halo *h;

  if(ngGlobal.Halos == NULL)
    {
      cart_debug("Halos are not loaded. Skipping ngHalo2IFrIT.");
      return;
    }
  
  if(ngHalo2IFrIT.Vars == NULL) ngHalo2IFrIT.Vars = (const int *)ngVI_p_;

  h = find_halo_by_id(ngGlobal.Halos,ngHalo2IFrIT.Id);
  if(h == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ngHalo2IFrIT.",ngHalo2IFrIT.Id);
    }
  else
    {
      for(i=0; i<ngHalo2IFrIT.NumVars; i++)
	{
	  if(ngHalo2IFrIT.Vars[i] == VAR_TOTAL_DENSITY) init->TotalDensity();
	  if(ngHalo2IFrIT.Vars[i] == VAR_STELLAR_DENSITY) init->StellarDensity();
	}
      ifrit.OutputHalo(ngOutputFile(filename),ngHalo2IFrIT.Level,nbin,h,ngHalo2IFrIT.NumVars,ngHalo2IFrIT.Vars);
    }
}
struct NG_OBJECT_HI_TYPE ngHalo2IFrIT = { ngHI_p_, 1, min(num_grid,256), max_level, sizeof(ngVI_p_)/sizeof(int), NULL };



void ngPZ_p_(const char *filename)
{
#if defined(COSMOLOGY) && defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER) && defined(RT_EXTERNAL_BACKGROUND)
  extProximityZones(ngOutputFile(filename),ngProximityZone.Level,ngProximityZone.Nside,ngProximityZone.Id,ngGlobal.Halos);
#else
  cart_debug("COSMOLOGY, RADIATIVE_TRANSFER, RT_TRANSFER, and RT_EXTERNAL_BACKGROUND are not set. Skipping ngProximityZone.");
#endif
}
struct NG_OBJECT_PZ_TYPE ngProximityZone = { ngPZ_p_, 0, max_level, 2 };


void ngRS_p_(const char *fileroot)
{
#if defined(RADIATIVE_TRANSFER) && defined(HYDRO) && defined(STARFORM)
  extRFvsSFR(fileroot,ngRFvsSFR.Level,ngGlobal.Halos);
#else
  cart_debug("RADIATIVE_TRANSFER, HYDRO, and STARFORM are not set. Skipping ngRFvsSFR.");
#endif
}
struct NG_OBJECT_RS_TYPE ngRFvsSFR = {  ngRS_p_, min_level };

 
void ngOKS_p_(const char *filename)
{
#if defined(PARTICLES) && defined(STARFORM)
  init->RadiativeTransfer();

  extStarFormationLaw(ngOutputFile(filename),ngObservedKSR.LengthScale,ngObservedKSR.TimeScale,ngObservedKSR.StellarAgeLimit,ngGlobal.Halos);
#else
  cart_debug("PARTICLES and STARFORM are not set. Skipping ngObservedKSR.");
#endif
}
struct NG_OBJECT_OKS_TYPE ngObservedKSR = { ngOKS_p_, 20.0, 0.5, 1.0e30 };


void ngIKS_p_(const char *filename)
{
#if defined(HYDRO) && defined(STARFORM)
  init->RadiativeTransfer();

  extStarFormationLaw2(ngOutputFile(filename),ngInstantaneousKSR.LengthScale,ngGlobal.Halos);
#else
  cart_debug("HYDRO and STARFORM are not set. Skipping ngInstantaneousKSR.");
#endif
}
struct NG_OBJECT_IKS_TYPE ngInstantaneousKSR = { ngIKS_p_, 0.5 };


void ngHP_p_(const char *filename)
{
  halo *h;

#ifdef COSMOLOGY
  if(ngGlobal.Halos == NULL)
    {
      cart_debug("Halos are not loaded. Skipping ngDumpHaloParticles.");
      return;
    }

  h = find_halo_by_id(ngGlobal.Halos,ngHaloParticles.Id);
  if(h == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ngDumpHaloParticles.",ngHaloParticles.Id);
    }
  else
    {
      dump_region_around_halo(ngOutputFile(filename),h,ngHaloParticles.Rmax);
    }
#else
  cart_debug("COSMOLOGY is not set. Skipping ngDumpHaloParticles.");
#endif /* COSMOLOGY */
}
struct NG_OBJECT_HP_TYPE ngHaloParticles = { ngHP_p_, 1, 1.0 };


void ngHS_p_(const char *filename)
{
  halo *h;

  if(ngGlobal.Halos == NULL)
    {
      cart_debug("Halos are not loaded. Skipping ngHaloStars.");
      return;
    }

  h = find_halo_by_id(ngGlobal.Halos,ngHaloStars.Id);
  if(h == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ngHaloStars.",ngHaloStars.Id);
    }
  else
    {
#if defined(PARTICLES) && defined(STARFORM)
      extHaloStars(filename,h,ngHaloStars.Rmax);
#endif
    }
}
struct NG_OBJECT_HS_TYPE ngHaloStars = { ngHS_p_, 1, 1.0 };



/* ---------------------------------------------------------------------- */


void AnalyseSnapshot(int argc, char **argv);


int main_analysis(int argc, char **argv)
{
  int ret;
  char str[999];
  int sepdirs = 1;

  ngLH_p_(NULL);

  if(sepdirs)
    {
      if(argc>0 && strcmp(argv[0],"-l")==0)
	{
	  sprintf(str,"a=%6.4f",auni[min_level]);
	  argc--;
	  argv++;
	}
      else
	{
	  sprintf(str,"a=%4.2f",auni[min_level]);
	}

      /*
      //  Create output directory
      */
      ret = system_mkdir(str);
      if(ret != 0)
	{
	  cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",str,ret);
	  return ret;
	}
      ngGlobal.OutputDir = str;
    }
	  
  AnalyseSnapshot(argc,argv);

  if(sepdirs)
    {
      ngGlobal.OutputDir = NULL;
    }

  return 0;
}

