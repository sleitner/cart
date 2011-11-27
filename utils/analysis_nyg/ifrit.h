#ifndef __NG_IFRIT_H__
#define __NG_IFRIT_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "framework.h"
#include "list.h"


namespace ng
{
  namespace ifrit
  {
    class VarSet : public List<int>
    {
    protected:

      virtual bool IsValid(const int& w) const;
    };


    class Base : public Algorithm
    {
    public:

      void SetNumBins(int nbins);
      void SetTopLevel(int level);
      void SetVars(const VarSet& set);

    protected:

      Base(const VarSet& set);

      struct Pars
      {
	VarSet& Set;
	int NumBins[3];
	int TopLevel;
	Pars(VarSet& set) : Set(set){}
      }
      p;

      static VarSet null;
    };


    //
    //  Produce IFrIR uniform scalars file of <Nbin>^3 size, centered at 
    //  <Point>, resolved to <Level>, for <VarSet> variables.
    //
    class Region : public Base
    {
    public:

      Region(const VarSet& set);

      void SetPoint(const double *point);
      void SetPoint(double x, double y, double z);

      void pos(double x, double y, double z){ this->SetPoint(x,y,z); }

      virtual void Exec(const char *path);

    protected:

      struct Pars
      {
	double Point[nDim];
      }
      p;
    };


    //
    //  As above, but now centered at a given halo with id=<Id>.
    //
    class Halo : public Base
    {
    public:

      Halo(const VarSet& set, int id = 1);

      void SetHaloId(int id);

      Shorthand<int> id;

      virtual void Exec(const char *path);

    protected:

      struct Pars
      {
	int Id;
      }
      p;
    };
  };
};


#endif // __NG_IFRIT_H__
