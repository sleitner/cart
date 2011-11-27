#include "config.h"
#include "framework.h"


extern "C"
{
#include "extra/ism.h"
#include "extra/start_analysis.h"
}


#include "celldata.h"

#include "profiles.h"


using namespace ng;


//
//  Class AbstractProbe
//
CellData AbstractProbe::null;
AbstractProbe::AbstractProbe() : cd(null)
{
}


AbstractProbe::AbstractProbe(const CellData& d) : cd(null)
{
  this->SetCellData(d);
}


int AbstractProbe::Size() const
{
  return cd.Size();
}


const CellDataWorker* AbstractProbe::Data() const
{
  return cd.Data();
}


void AbstractProbe::SetCellData(const CellData& d)
{
  cd.Copy(d);
}


void AbstractProbe::Exec(const char *path)
{
  int i;

  for(i=0; i<cd.Size(); i++)
    {
      if(cd.Worker(i) == cd::total_density)
	{
	  init->TotalDensity();
	}

      if(cd.Worker(i) == cd::stellar_density)
	{
	  init->StellarDensity();
	}

      if(cd.Worker(i) == cd::interstellar_radiation_field)
	{
	  init->RadiativeTransfer();
	}
    }

  this->ExecBody(path);
}


//
//  Class DumpLevels
//
DumpLevels::DumpLevels() : AbstractProbe(), lev(p.MinLevel), maxlev(p.MaxLevel)
{
  p.LowMemoryMode = false;
  this->SetLevel(max_level_now_global(mpi.comm.run));
}


DumpLevels::DumpLevels(const CellData& d) : AbstractProbe(d), lev(p.MinLevel), maxlev(p.MaxLevel)
{
  p.LowMemoryMode = false;
  this->SetLevel(max_level_now_global(mpi.comm.run));
}


DumpLevels::DumpLevels(const CellData& d, int minLevel, int maxLevel) : AbstractProbe(d), lev(p.MinLevel), maxlev(p.MaxLevel)
{
  p.LowMemoryMode = false;
  this->SetLevelRange(minLevel,maxLevel);
}


void  DumpLevels::SetLevel(int level)
{
  this->SetLevelRange(level,level);
}


void DumpLevels::SetLevelRange(int minLevel, int maxLevel)
{
  cart_assert(minLevel>=min_level && minLevel<=maxLevel && maxLevel<=max_level);

  p.MinLevel = minLevel;
  p.MaxLevel = maxLevel;
}


void DumpLevels::SetLowMemoryMode(bool s)
{
  p.LowMemoryMode = s;
}


void DumpLevels::ExecBody(const char *path)
{
  if(p.LowMemoryMode)
    {
      extDumpLevelsLowMemory(this->File(path),this->Size(),this->Data(),p.MinLevel,p.MaxLevel,this->Halos());
    }
  else
    {
      extDumpLevels(this->File(path),this->Size(),this->Data(),p.MinLevel,p.MaxLevel,this->Halos());
    }
}


//
//  Class AbstractProfile
//
AbstractProfile::AbstractProfile() : AbstractProbe(), ndex(p.Ndex), rmin(p.Rmin), rmax(p.Rmax)
{
  p.Ndex = 25;
  p.Rmin = 0.1;
  p.Rmax = 300.0;
}


AbstractProfile::AbstractProfile(const CellData& d) : AbstractProbe(d), ndex(p.Ndex), rmin(p.Rmin), rmax(p.Rmax)
{
  p.Ndex = 25;
  p.Rmin = 0.1;
  p.Rmax = 300.0;
}


void AbstractProfile::SetRadialRange(float rmin, float rmax)
{
  cart_assert(rmin>0.0 && rmin<rmax && rmax<1.0e10);

  p.Rmin = rmin;
  p.Rmax = rmax;
}


void AbstractProfile::SetBinsPerDex(int ndex)
{
  cart_assert(ndex>0 && ndex<1000);

  p.Ndex = ndex;
}


//
//  Class HaloProfiles
//
HaloProfiles::HaloProfiles() : AbstractProfile(), lev(p.MinHaloLevel), edge(p.OuterEdge)
{
  p.MinHaloLevel = min_level;
  p.OuterEdge = 1.0;
}


HaloProfiles::HaloProfiles(const CellData& d, int level, float edge) : AbstractProfile(d), lev(p.MinHaloLevel), edge(p.OuterEdge)
{
  this->SetMinHaloLevel(level);
  this->SetOuterEdge(edge);
}


void HaloProfiles::SetMinHaloLevel(int level)
{
  cart_assert(level>=min_level && level<=max_level);

  p.MinHaloLevel = level;
}


void HaloProfiles::SetOuterEdge(float edge)
{
  cart_assert(edge>0 && edge<100);

  p.OuterEdge = edge;
}


void HaloProfiles::ExecBody(const char *path)
{
#ifdef COSMOLOGY
  extDumpHaloProfiles(this->File(path),this->Size(),this->Data(),AbstractProfile::p.Rmin,AbstractProfile::p.Rmax,AbstractProfile::p.Ndex,this->Halos(),p.MinHaloLevel,p.OuterEdge);
#else
  cart_debug("COSMOLOGY is not set. Skipping writing halo profiles...");
#endif /* COSMOLOGY */
}


//
//  Class PointProfile
//
PointProfile::PointProfile() : AbstractProfile()
{
  int i;
  for(i=0; i<nDim; i++) p.Point[i] = 0.0;
}


PointProfile::PointProfile(const CellData& d, const double *point) : AbstractProfile(d)
{
  this->SetPoint(point);
}


void PointProfile::SetPoint(const double *point)
{
  int i;
  for(i=0; i<nDim; i++) p.Point[i] = point[i];
}


void PointProfile::SetPoint(double x, double y, double z)
{
  p.Point[0] = x;
#if (nDim > 1)
  p.Point[1] = y;
#if (nDim > 2)
  p.Point[2] = z;
#endif
#endif
}


void PointProfile::ExecBody(const char *path)
{
  extDumpPointProfile(this->File(path),this->Size(),this->Data(),AbstractProfile::p.Rmin,AbstractProfile::p.Rmax,AbstractProfile::p.Ndex,p.Point);
}


