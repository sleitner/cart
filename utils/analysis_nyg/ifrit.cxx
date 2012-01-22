#include "config.h"
#include "framework.h"


extern "C"
{
#include "extra/ifrit.h"
#include "extra/start_analysis.h"
}


#include "ifrit.h"


using namespace ng;
using namespace ifrit;
  

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
  p.TopLevel = min_level;
}

void Base::SetNumBins(int nbins)
{
  cart_assert(nbins>0 && nbins<10000);

  p.NumBins[0] = p.NumBins[1] = p.NumBins[2] = nbins; 
}


void Base::SetTopLevel(int level)
{
  cart_assert(min_level<=level && level<=max_level);

  p.TopLevel = level;
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
      if(s.Data()[i] == VAR_TOTAL_DENSITY) init->TotalDensity();
      if(s.Data()[i] == VAR_STELLAR_DENSITY) init->StellarDensity();
    }

  ::ifrit.OutputBox(this->File(path),Base::p.TopLevel,Base::p.NumBins,p.Point,s.Size(),s.Data());
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
  
  p.Id = id;
}


void Halo::Exec(const char *path)
{
  if(this->Halos() == 0)
    {
      cart_debug("Halos are not loaded. Skipping ng::ifrit::Halo...");
      return;
    }
  
  halo *h = find_halo_by_id(this->Halos(),p.Id);
  if(h == NULL)
    {
      cart_debug("There is no halo with id=%d. Skipping ng::ifrit::Halo...",p.Id);
    }
  else
    {
      int i;
      const VarSet& s(Base::p.Set);

      for(i=0; i<s.Size(); i++)
	{
	  if(s.Data()[i] == VAR_TOTAL_DENSITY) init->TotalDensity();
	  if(s.Data()[i] == VAR_STELLAR_DENSITY) init->StellarDensity();
	}

      char str[strlen(path)+20];
      sprintf(str,"%s-id=%05d",path,p.Id);

      ::ifrit.OutputHalo(this->File(str),Base::p.TopLevel,Base::p.NumBins,h,s.Size(),s.Data());
    }
}