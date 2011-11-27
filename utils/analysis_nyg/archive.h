#ifndef __NG_ARCHIVE_H__
#define __NG_ARCHIVE_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


namespace ng
{
  //
  //  Compute proximity zones for a specific halo with id=<Id>, 
  //  or for all halos, if <Id>=0. Only consider halos resolved 
  //  to at least <Level> and use <Nside> parameter in Healpix.
  //
  class ProximityZones : public Algorithm
  {
  public:

    ProximityZones();
    ProximityZones(int id, int level, int nside = 1);

    void SetHaloId(int id);
    void SetMinHaloLevel(int level);
    void SetResolution(int nside);

    void Exec(const char *path, int id);
    virtual void Exec(const char *path);

  protected:

    struct Pars
    {
      int Id;
      int MinHaloLevel;
      int Nside;
    }
    p;
  };
};


#endif // __NG_ARCHIVE_H__
