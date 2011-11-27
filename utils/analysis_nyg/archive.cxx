#include "config.h"
#include "framework.h"


#ifdef NG_ARCHIVE


extern "C"
{
#include "extra/igm.h"
}


#include "archive.h"
#include "celldata.h"


using namespace ng;


//
//  Class ProximityZones
//
ProximityZones::ProximityZones()
{
  p.Id = 0;
  p.MinHaloLevel = min_level;
  p.Nside = 1;
}


ProximityZones::ProximityZones(int id, int level, int nside)
{
  this->SetHaloId(id);
  this->SetMinHaloLevel(level);
  this->SetResolution(nside);
}

void ProximityZones::SetHaloId(int id)
{
  cart_assert(id>=0 && id<=this->MaxHaloId());

  p.Id = id;
}


void ProximityZones::SetMinHaloLevel(int level)
{
  cart_assert(level>=min_level && level<=max_level);

  p.MinHaloLevel = level;
}


void ProximityZones::SetResolution(int nside)
{
  cart_assert(nside>0 && nside<1000);

  p.Nside = nside;
}


void ProximityZones::Exec(const char *path, int id)
{
  this->SetHaloId(id);
  this->Exec(path);
}


void ProximityZones::Exec(const char *path)
{
#if defined(COSMOLOGY) && defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER) && defined(RT_EXTERNAL_BACKGROUND)
  extProximityZones(this->File(path),p.MinHaloLevel,p.Nside,p.Id,this->Halos());
#else
  cart_debug("COSMOLOGY, RADIATIVE_TRANSFER, RT_TRANSFER, and RT_EXTERNAL_BACKGROUND are not set. Skipping ng::ProximityZone...");
#endif
}

#endif // NG_ARCHIVE
