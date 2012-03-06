#include "config.h"
#include "framework.h"


extern "C"
{
#include "io.h"
#include "system.h"
#include "times.h"
}


using namespace ng;


struct HALO_LIST* Catalog::HaloList = 0;
const char* Algorithm::OutputPath = 0;


// ----------------------------------------------------------------------


extern "C" int main_analysis(int argc, char **argv)
{
  int ret;
  char str[999];
  int sepdirs = 1;

  Catalog::LoadHalos(0);

  if(sepdirs)
    {
      if(argc>0 && strcmp(argv[0],"-l")==0)
	{
	  sprintf(str,"a=%6.4f",auni[min_level]);
	  argc--;
	  argv++;
	}
      else
	{
	  sprintf(str,"a=%4.2f",auni[min_level]);
	}

      /*
      //  Create output directory
      */
      ret = system_mkdir(str);
      if(ret != 0)
	{
	  cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",str,ret);
	  return ret;
	}
      Algorithm::SetOutputPath(str);
    }
	  
  if(argc>0 && strstr(argv[0],"-h")==argv[0])
    {
      if(strcmp(argv[0],"-h") == 0)
	{
	  char s[strlen(output_directory)+9];
	  sprintf(s,"%s/HC",output_directory);
	  Catalog::LoadHalos(s);
	}
      else
	{
	  Catalog::LoadHalos(argv[0]+3);
	}
      argc--;
      argv++;
    }

  execute(argc,argv);

  if(sepdirs)
    {
      Algorithm::SetOutputPath(0);
    }

  return 0;
}


// ----------------------------------------------------------------------


//
//  MeshTraversal class
//
MeshTraversal::MeshTraversal(bool parallel)
{
  p.Parallel = parallel;
}


void MeshTraversal::Traverse()
{
  this->Traverse(min_level,max_level);
}


void MeshTraversal::Traverse(int minLevel, int maxLevel)
{
  cart_assert(minLevel>=min_level && minLevel<=maxLevel && max_level<=max_level);

  MESH_RUN_DECLARE(level,cell);
  int par = p.Parallel ? 1 : 0;

  MESH_RUN_OVER_LEVELS_BEGIN(level,minLevel,maxLevel);

  this->OnLevelBegin(level);

#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct), if(par)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

  this->OnAnyCell(cell);

  if(cell_is_leaf(cell))
    {
      this->OnLeafCell(cell);
    }

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  this->OnLevelEnd(level);
  
  MESH_RUN_OVER_LEVELS_END;
}


//
//  Catalog class
//
int Catalog::NumHalos()
{
  if(HaloList != 0) return HaloList->num_halos; else return 0;
}

const halo* Catalog::Halo(int i)
{
  if(HaloList!=0 && i>=0 && i<HaloList->num_halos) return HaloList->list+i; else return 0;
}


const halo* Catalog::GetHaloById(int id)
{
  return find_halo_by_id(HaloList,id);
}


//
//  Load hlist file for other commands to use
//
void Catalog::LoadHalos(const char *path, int Nmin, float Mvir, float Vmax, float Rvir, int MaxNumHalos)
{
#ifdef COSMOLOGY
  char str[999];

  if(path == 0)
    {
      if(HaloList != 0)
	{
	  destroy_halo_list(HaloList);
	  HaloList = 0;
	}
    }
  else
    {
      sprintf(str,"%s/hlist_%6.4f.dat",path,auni[min_level]);
      HaloList = load_halo_finder_catalog(str,Nmin,Mvir,Vmax,Rvir,MaxNumHalos);
      if(HaloList == 0)
	{
	  cart_error("Failed to read in the halo list from file %s",str);
	}
    }
#else
  cart_debug("COSMOLOGY is not set. Skipping ng::Catalog::LoadHalos()...");
#endif /* COSMOLOGY */
}


//
//  Algorithm class
//
void Algorithm::SetOutputPath(const char *path)
{
  OutputPath = path;
}


const char* Algorithm::File(const char *path)
{
  static char str[999];

  cart_assert(strlen(path) < 500);

  if(OutputPath == 0)
    {
      strcpy(str,path);
    }
  else
    {
      cart_assert(strlen(OutputPath) < 490);

      strcpy(str,OutputPath);
      strcat(str,"/");
      strcat(str,path);
    }

  char *name = strrchr(str,'/');
  if(name != 0)
    {
      *name = 0;
      /*
      //  Create output directory
      */
      int ret = system_mkdir(str);
      if(ret != 0)
	{
	  cart_error("system_mkdir(\"%s\") returned with the error code %d. Abort.",str,ret);
	  return 0;
	}
      *name = '/';
    }

  return str;
}

