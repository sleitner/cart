#include "config.h"
#include "framework.h"


extern "C"
{
#include "../extra/ifrit.h"
#include "../drivers/start_analysis.h"
}


#include "ifrit.h"


using namespace ng;
using namespace ng::ifrit;
  

//
//  Class VarSet
//
bool VarSet::IsValid(const int& w) const
{
  return (w >= 0);
}



//
//  Class Base
//
VarSet Base::null;
Base::Base(const VarSet& set) : p(null)
{
  this->SetVars(set);
  p.NumBins[0] = p.NumBins[1] = p.NumBins[2] = num_grid;
  p.PixelLevel = max_level;
}

void Base::SetNumBins(int nbins)
{
  cart_assert(nbins>0 && nbins<10000);

  p.NumBins[0] = p.NumBins[1] = p.NumBins[2] = nbins; 
}


void Base::SetPixelLevel(int level)
{
  cart_assert(min_level<=level && level<=max_level);

  p.PixelLevel = level;
}


void Base::SetVars(const VarSet& set)
{
  p.Set.Copy(set);
}


//
//  Class Region
//
Region::Region(const VarSet& set) : Base(set)
{
  int i;
  for(i=0; i<nDim; i++) p.Point[i] = 0.0;
}


void Region::SetPoint(const double *point)
{
  int i;
  for(i=0; i<nDim; i++) p.Point[i] = point[i];
}


void Region::SetPoint(double x, double y, double z)
{
  p.Point[0] = x;
#if (nDim > 1)
  p.Point[1] = y;
#if (nDim > 2)
  p.Point[2] = z;
#endif
#endif
}


void Region::Exec(const char *path)
{
  int i;
  const VarSet& s(Base::p.Set);

  for(i=0; i<s.Size(); i++)
    {
#ifdef GRAVITY
      if(s.Data()[i] == VAR_TOTAL_DENSITY) init->TotalDensity();
#ifdef STAR_FORMATION
      if(s.Data()[i] == VAR_STELLAR_DENSITY) init->StellarDensity();
#endif
#endif
    }

  ::ifrit.OutputBox(this->File(path),Base::p.PixelLevel,Base::p.NumBins,p.Point,s.Size(),s.Data());
}


//
//  Class Halo
//
Halo::Halo(const VarSet& set, int i) : Base(set), id(p.Id)
{
  this->SetHaloId(i);
}


void Halo::SetHaloId(int id)
{
  cart_assert(id > 0);

  const halo *hptr = Catalog::GetHaloById(id);
  if(hptr == NULL)
    {
      cart_debug("There is no halo with id=%d.",id);
      return;
    }

  p.Id = id;
  this->SetPixelLevel(halo_level(hptr,mpi.comm.run));
}


void Halo::Exec(const char *path)
{
  if(this->Halos() == 0)
    {
      cart_debug("Halos are not loaded. Skipping ng::ifrit::Halo...");
      return;
    }
  
  const halo *hptr = Catalog::GetHaloById(p.Id);
  if(hptr == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ng::ifrit::Halo...",p.Id);
    }
  else
    {
      int i;
      const VarSet& s(Base::p.Set);

      for(i=0; i<s.Size(); i++)
	{
#ifdef GRAVITY
	  if(s.Data()[i] == VAR_TOTAL_DENSITY) init->TotalDensity();
#ifdef STAR_FORMATION
	  if(s.Data()[i] == VAR_STELLAR_DENSITY) init->StellarDensity();
#endif
#endif
	}

      char str[strlen(path)+20];
      sprintf(str,"%s-id=%05d",path,p.Id);

      ::ifrit.OutputHalo(this->File(str),Base::p.PixelLevel,Base::p.NumBins,hptr,s.Size(),s.Data());
    }
}
