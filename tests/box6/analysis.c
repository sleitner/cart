#include "config.h"

#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "parallel.h"
#include "rt_solver.h"
#include "tree.h"
#include "units.h"

#include "extra/halo_finder.h"
#include "extra/ifrit.h"
#include "extra/igm.h"
#include "extra/ism.h"

#ifdef RADIATIVE_TRANSFER
#include "F/frt_c.h"
#endif


extern int num_options;
extern char **options;


void dump_worker(int level, int cell, int num, float *ptr)
{
  float s, soblen, sobvel;
#ifdef RADIATIVE_TRANSFER
  float rate[frtRATE_DIM];
#else
  float rate;
#endif
  int i, nb[num_neighbors];

  cart_assert(num >= 9);

#ifdef HYDRO
  rtGetSobolevFactors(cell,level,&soblen,&sobvel);
  ptr[0] = units->temperature*cell_gas_temperature(cell);
  ptr[1] = cell_gas_density(cell)*soblen*units->number_density*units->length;
#else
  soblen = sobvel = 0.0;
#endif /* HYDRO */

#ifdef ENRICH
  ptr[2] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#endif

#ifdef RADIATIVE_TRANSFER
  rtGetPhotoRates(cell,rate);
  ptr[3] = cell_HI_fraction(cell);
  ptr[4] = cell_HII_fraction(cell);
  ptr[5] = cell_H2_fraction(cell);
  ptr[6] = rate[frtRATE_CiLW]*1.05e10;  /* UV field at 12.0eV in units of Draine field (1.0e6 phot/cm^2/s/ster/eV) */
  ptr[7] = rtDustToGas(cell);
#endif  /* RADIATIVE_TRANSFER */

#ifdef HYDRO
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
  ptr[8] = cell_gas_density(cell)/s;
#endif /* HYDRO */
}


void analysis()
{
  int i, id, level;
  float s, t, a, rmin, rmax;
  halo_list *halos = NULL;
  int floor_level;
  const char *tmp;
#ifdef RADIATIVE_TRANSFER
  int varid[] = { I_FRACTION+RT_HVAR_OFFSET+0, HVAR_GAS_DENSITY, I_GAS_TEMPERATURE, I_FRACTION+RT_HVAR_OFFSET+5, I_CELL_LEVEL, I_LOCAL_PROC, I_FRACTION+RT_HVAR_OFFSET+1 };
#else
#ifdef HYDRO
  int varid[] = { HVAR_GAS_DENSITY, HVAR_PRESSURE };
#else
  int varid[] = { VAR_DENSITY };
#endif /* HYDRO */
#endif
  int nvars = sizeof(varid)/sizeof(int);
  const int nbin1 = 256;
  int nbin[] = { nbin1, nbin1, nbin1 };
  double bb[6];

  if(num_options < 1)
    {
      if(local_proc_id == MASTER_NODE)
	{
	  printf("An option is required. Available options:\n");
#ifdef HYDRO
	  printf("  -dl[=<level>]:  dump chemical state of the gas at level <level>, or max_level by default\n");
#endif /* HYDRO */
#ifdef COSMOLOGY
	  printf("  -dp=<rmin>,<rmax>:  dump profiles of gas quantities from <rmin> to <rmax>\n");
	  printf("  -hl=<hlist-file>:  load halo file for other options to use\n");
	  printf("  -gf:   dump gas fractions for halos (which must be loaded)\n");
#endif /* COSMOLOGY */
	  printf("  -ifrit[=<halo-id>,<zoom>]:  produce IFrIR uniform scalars file, rebinned at the lowest level\n");
#ifdef COSMOLOGY
	  printf("  -jnu:  dump average cosmic radiation background\n");
	  printf("  -pz[=<halo-id>]:  compute proximity zones (optionally for a specific halo)\n");
#endif /* COSMOLOGY */
	  printf("  -sfl=<scale>,<time>[,stellar_age]:  dump Kennicutt law on scale <scale> kpc averaged over <time> Myr, and optionally limit stellar age in output M* to <stellar_age>\n");
	  printf("  -sfl2=<scale>:        dump instantaneous Kennicutt law on scale <scale> kpc\n");
	  printf("  -vsl[=<level>]:  verify Sobolev-like approximations to true integrated column densities\n");
	}
    }

  /*
  //  Find the lowest level on the mesh
  */
  floor_level = max_level_now_global(MPI_COMM_WORLD);

  for(i=0; i<num_options; i++)
    {
#ifdef HYDRO
      tmp = check_option1(options[i],"dl","-1");
      if(tmp != NULL)
	{
	  if(sscanf(tmp,"%d",&level)==0 || level < -1)
	    {
	      cart_error("-dl[=<level>] requires non-negative integer <level> as an argument.");
	    }
	  if(level<0 || level>floor_level) level = floor_level;
	  extDumpLevels("dl.res",dump_worker,level,level,halos);
	  continue;
	}
#endif /* HYDRO */


#ifdef COSMOLOGY
      tmp = check_option1(options[i],"dp","-1");
      if(tmp != NULL)
	{
	  if(sscanf(tmp,"%e,%e",&rmin,&rmax) != 2)
	    {
	      cart_error("-dp=<rmin>,<rmax> requires two radii as arguments.");
	    }
	  extDumpProfiles("dp.res",dump_worker,floor_level,rmin,rmax,10,halos);
	  continue;
	}


      /*
      //  Load hlist file for other command to use
      */
      tmp = check_option1(options[i],"hl","");
      if(tmp != NULL)
	{
	  if(strcmp(tmp,"") == 0)
	    {
	      if(halos != NULL)
		{
		  destroy_halo_list(halos);
		  halos = NULL;
		}
	    }
	  else
	    {
	      halos = load_halo_finder_catalog(tmp,100,0.0,0.0,0.0);
	      if(halos == NULL)
		{
		  cart_error("Failed to read in the halo list from file %s",tmp);
		}
	      else
		{
		  // compute_halo_properties(".",-1,halos);
		}
	    }
	  continue;
	}


      tmp = check_option0(options[i],"gf");
      if(tmp != NULL)
	{
	  extDumpGasFractions("gf.res",halos);
	  continue;
	}
#endif /* COSMOLOGY */


      tmp = check_option1(options[i],"ifrit","0,1");
      if(tmp != NULL)
	{
	  id = -1;
	  s = -1.0;
	  if(sscanf(tmp,"%d,%e",&id,&s)==0 || id<0 || s<0.1)
	    {
	      cart_error("-ifrit=<halo-id>,<zoom> requires positive integer <halo-id> and a positive float number <zoom> as arguments.");
	    }
	  if(id>0 && halos==NULL)
	    {
	      cart_debug("Halo list must be loaded for dumping IFrIT file for a given halo. Skipping this command.");
	    }
	  else
	    {
	      if(halos == NULL)
		{
		  bb[0] = bb[2] = bb[4] = 0.0;
		  bb[1] = bb[3] = bb[5] = num_grid;
		  iOutputMesh("ifrit-box.bin",floor_level,nbin,bb,nvars,varid);
		}
	      else
		{
		  if(id >= halos->num_halos)
		    {
		      cart_debug("There is no halo with id=%d. Skipping this command.",id);
		    }
		  else iOutputHalo("ifrit",floor_level,s,&halos->list[id],nvars,varid);
		}
	    }
	  continue;
	}


#ifdef COSMOLOGY
      tmp = check_option0(options[i],"jnu");
      if(tmp != NULL)
	{
#ifdef RADIATIVE_TRANSFER
	  extDumpRadiationBackground();
#else
	  cart_error("-jnu requires RADIATIVE_TRANSFER on");
#endif
	  continue;
	}


      tmp = check_option1(options[i],"pz","0");
      if(tmp != NULL)
	{
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER) && defined(RT_EXTERNAL_BACKGROUND)
	  id = -1;
	  if(sscanf(tmp,"%d",&id)==0 || id<0)
	    {
	      cart_error("-pz=<halo-id> requires positive integer <halo-id> as an argument.");
	    }
	  extFindProximityZones("pz.res",floor_level,2,id,halos);
#else
	  cart_error("-pz requires RADIATIVE_TRANSFER, RT_TRANSFER, and RT_EXTERNAL_BACKGROUND on");
#endif
	  continue;
	}
#endif /* COSMOLOGY */


      tmp = check_option1(options[i],"sfl",NULL);
      if(tmp != NULL)
	{
#if defined(PARTICLES) && defined(STARFORM)
	  if(sscanf(tmp,"%e,%e,%e",&s,&t,&a) != 3)
	    {
	      if(sscanf(tmp,"%e,%e",&s,&t) != 2)
		{
		  cart_error("-sfl option format is: -sfl=<scale>,<time>[,stellar_age]");
		}
	      a = 1.0e10;
	    }
	  extStarFormationLaw("sfl.res",s,t,a,halos);
#else
	  cart_error("-sfl requires PARTICLES && STARFORM on");
#endif
	  continue;
	}


      tmp = check_option1(options[i],"sfl2",NULL);
      if(tmp != NULL)
	{
#if defined(HYDRO) && defined(STARFORM)
	  if(sscanf(tmp,"%e",&s) != 1)
	    {
	      cart_error("-sfl2 option format is: -sfl2=<scale>");
	    }
	  extStarFormationLaw2("sfl2.res",s,halos);
#else
	  cart_error("-sfl2 requires HYDRO && STARFORM on");
#endif
	  continue;
	}


      tmp = check_option1(options[i],"vsl","-1");
      if(tmp != NULL)
	{
#if defined (HYDRO) && defined(RADIATIVE_TRANSFER)
	  if(sscanf(tmp,"%d",&level)==0 || level < -1)
	    {
	      cart_error("-sl[=<level>] requires non-negative integer <level> as an argument.");
	    }
	  if(level<0 || level>floor_level) level = floor_level;
	  extCheckSobolevApproximations("sl.res",level,1,10.0);
#else
	  cart_error("-sl requires HYDRO and RADIATIVE_TRANSFER on");
#endif  /* HYDRO && RADIATIVE_TRANSFER */
	  continue;
	}

      cart_error("Invalid option %s",options[i]);
    }

  if(halos != NULL) destroy_halo_list(halos);
}


