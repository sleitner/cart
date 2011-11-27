#ifndef __NG_PROFILES_H__
#define __NG_PROFILES_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "framework.h"


namespace ng
{
  //
  //  Base class for all cell data probes
  //
  class AbstractProbe : public Algorithm
  {
  public:

    void SetCellData(const CellData& dump);

    virtual void Exec(const char *path);

  protected:

    AbstractProbe();
    AbstractProbe(const CellData& dump);

    virtual void ExecBody(const char *path) = 0;

    int Size() const;
    const struct DUMP_WORKER* Data() const;

  private:

    CellData &cd;

    static CellData null;
  };


  //
  //  Write cell data at a range of levels.
  //
  class DumpLevels : public AbstractProbe
  {
  public:

    DumpLevels();
    DumpLevels(const CellData& dump);
    DumpLevels(const CellData& dump, int minLevel, int maxLevel = max_level);

    void SetLevel(int level);
    void SetLevelRange(int minLevel, int maxLevel);
    void SetLowMemoryMode(bool s);

    Shorthand<int> lev;
    Shorthand<int> maxlev;

  protected:

    virtual void ExecBody(const char *path);

    struct
    {
      int MinLevel;
      int MaxLevel;
      bool LowMemoryMode;
    }
    p;
  };


  //
  //  Base class for spherical profile writers
  //
  class AbstractProfile : public AbstractProbe
  {
  public:

    AbstractProfile();
    AbstractProfile(const CellData& d);

    void SetRadialRange(float rmin, float rmax);
    void SetBinsPerDex(int ndex);

    Shorthand<int> ndex;
    Shorthand<float> rmin;
    Shorthand<float> rmax;

  protected:

    struct Pars
    {
      int Ndex;
      float Rmin;
      float Rmax;
    }
    p;
  };


  //
  //  Write spherical profiles around halos
  //
  class HaloProfiles : public AbstractProfile
  {
  public:

    HaloProfiles();
    HaloProfiles(const CellData& d, int level, float edge = 1.0);

    void SetMinHaloLevel(int level);
    void SetOuterEdge(float edge);

    Shorthand<int> lev;
    Shorthand<float> edge;

  protected:

    virtual void ExecBody(const char *path);

    struct
    {
      int MinHaloLevel;
      float OuterEdge;
    }
    p;
  };


  //
  //  Write spherical profile around a single points
  //
  class PointProfile : public AbstractProfile
  {
  public:

    PointProfile();
    PointProfile(const CellData& d, const double *point);

    void SetPoint(const double *point);
    void SetPoint(double x, double y, double z);

    void pos(double x, double y, double z){ this->SetPoint(x,y,z); }

  protected:

    virtual void ExecBody(const char *path);

    struct Pars
    {
      double Point[nDim];
    }
    p;
  };
};

#endif  // __NG_PROFILES_H__
