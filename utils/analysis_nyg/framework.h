#ifndef __NG_FRAMEWORK_H__
#define __NG_FRAMEWORK_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


extern "C"
{
#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "tree.h"
#include "units.h"
#include "../extra/halo_finder.h"

  struct DUMP_WORKER;
}


namespace ng
{
  void execute(int argc, char **argv);

  class CellData;

  template<class T>
    class Shorthand
    {
    public:

      Shorthand(T& v) : ref(v){}
      void operator=(T v){ ref = v; }

    private:

      T& ref;
    };

  class MeshTraversal
  {
  public:

    void Traverse();
    void Traverse(int minLevel, int maxLevel);

  protected:

    virtual void OnLevelBegin(int level){}
    virtual void OnLevelEnd(int level){}
    virtual void OnAnyCell(int cell){}
    virtual void OnLeafCell(int cell){}

    MeshTraversal(bool parallel);

    struct Pars
    {
      bool Parallel;
    } 
    p;

  private:
    MeshTraversal(const MeshTraversal&); // Not implemented.
    void operator=(const MeshTraversal&);  // Not implemented.
  };


  class Catalog
  {

    friend class Algorithm;

  public:

    static int NumHalos();
    static const halo* Halo(int i);
    static const halo* GetHaloById(int id);

    //
    //  Load hlist file for other command to use.
    //
    static void LoadHalos(const char *path, int Nmin = 50, float Mvir = 0.0, float Vmax = 0.0, float Rvir = 0.0, int MaxNumHalos = 1000000000);

  private:

    static halo_list *HaloList;
  };


  class Algorithm
  {
  public:

    virtual void Exec(const char *path) = 0;

    static void SetOutputPath(const char *path);

  protected:

    Algorithm(){}

    const char* File(const char *path);

    struct HALO_LIST* Halos() const { return Catalog::HaloList; }
    int MaxHaloId() const { return Catalog::HaloList->list[Catalog::HaloList->num_halos-1].id; }

  private:

    static const char *OutputPath;

  private:
    Algorithm(const Algorithm&); // Not implemented.
    void operator=(const Algorithm&);  // Not implemented.
  };
};

#endif  // __NG_FRAMEWORK_H__
