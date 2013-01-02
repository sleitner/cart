#include "config.h"

#include <math.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "tree.h"

#include "cd.h"
#include "los.h"


/*
//  Compute column densities of a set of vars to a given distance
*/
struct
{
  int nvars;
  int *vars;
}
cdPars;


int cdWorker(int id, int cell, double r1, double r2, losBuffer data)
{
  int j;
  cdData *d = (cdData *)data.Data;

  for(j=0; j<cdPars.nvars; j++)
    {
      d->val[j] += cell_var(cell,cdPars.vars[j])*(r2-r1);
    }
  d->len += (r2-r1);

  return 0;
}


void cdCollector(const losBuffer *result, int num_segments, const losSegment *segments)
{
  int i, j;
  cdData *dn, *dr = (cdData *)result->Data;

  for(j=0; j<cdPars.nvars; j++)
    {
      dr->val[j] = 0.0;
    }
  dr->len = 0.0;

  for(i=0; i<num_segments; i++)
    {
      dn = (cdData *)segments[i].Buffer.Data;

      for(j=0; j<cdPars.nvars; j++)
	{
	  dr->val[j] += dn->val[j];
	}
      dr->len += dn->len;
    }
}


/*
//  Compute column densities of <nvars> grid variables (with indicies 
//  supplied in <vars[]> array, sampled over the sky using the HealPIX 
//  binning with <nside> bins per section (total number of rays is 
//  12*nside^2) with common origin <pos>, of maximum length <len>, not deeper 
//  than <floor_level>. The final result is returned in <output[]>, whose
//  dimension must be 12*nside^2
*/
void cdTraverseSky(int nvars, int vars[], int nside, double pos[3], double len, int floor_level, cdData output[])
{
  int npix = 12*nside*nside;
  int ipix, j;
  losBuffer *lines;

  cart_assert(0<nvars && nvars<=CD_MAX_NVARS);

  cdPars.nvars = nvars;
  cdPars.vars = vars;

  /*
  //  Prepare buffers for the LOS traversal
  */
  lines = cart_alloc(losBuffer,npix);
  for(ipix=0; ipix<npix; ipix++)
    {
      lines[ipix].Size = sizeof(cdData);
      lines[ipix].Data = output + ipix;
    }

  /*
  //  Init LOS buffers
  */
#pragma omp parallel for default(none), private(ipix,j), shared(npix,lines,nvars,output)
  for(ipix=0; ipix<npix; ipix++)
    {
      for(j=0; j<nvars; j++) output[ipix].val[j] = 0.0;
      output[ipix].len = 0.0;
    }

  losTraverseSky(nside,pos,len,floor_level,lines,cdWorker,cdCollector);

  cart_free(lines);
}



