#ifndef __ANALYSIS_NYG_H__
#define __ANALYSIS_NYG_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <math.h>
#include <stdio.h>


#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "tree.h"
#include "units.h"


struct HALO_LIST;
struct ngCellDump;


const struct HALO_LIST *ngHalos();
const char *ngOutputFile(const char *filepath);


/*
//  Load hlist file for other command to use.
*/
struct NG_OBJECT_LH_TYPE
{
  void (*Exec)(const char *path);
  int Nmin;
  float Mvir;
  float Vmax;
  float Rvir;
};
extern struct NG_OBJECT_LH_TYPE ngLoadHalos;


/*
//  Dump chemical state of the gas at levels <Level> to <MaxLevel>.
*/
struct NG_OBJECT_DL_TYPE
{
  void (*Exec)(const char *filename);
  int Level;
  int MaxLevel;
  const struct ngCellDump *Dump;
};
extern struct NG_OBJECT_DL_TYPE ngDumpLevels;


/*
//  Dump profiles of gas quantities from <Rmin> to <Rmax> (in kpc) with
// <Ndex> bins per decade for halos resolved to at least <Level>;
// halos end at <HaloEdge>*Rvir.
*/
struct NG_OBJECT_DP_TYPE
{
  void (*Exec)(const char *filename);
  int Ndex;
  int Level;
  float Rmin;
  float Rmax;
  float HaloEdge;
  const struct ngCellDump *Dump;
};
extern struct NG_OBJECT_DP_TYPE ngDumpProfiles;


/*
//  Dump gas and other mass fractions fractions for halos.
*/
struct NG_OBJECT_MF_TYPE
{
  void (*Exec)(const char *filename);
};
extern struct NG_OBJECT_MF_TYPE ngMassFractions;


/*
//  Produce IFrIR uniform scalars file of <Nbin>^3 size, centered at <Pos>,
//  resolved to <Level>, for <NumVars> variables from <Vars>.
*/
struct NG_OBJECT_RI_TYPE
{
  void (*Exec)(const char *filename);
  int Nbin;
  int Level;
  double Pos[3];
  int NumVars;
  const int *Vars;
};
extern struct NG_OBJECT_RI_TYPE ngRegion2IFrIT;


/*
//  As above, but now centered at a given halo with id=<Id>.
*/
struct NG_OBJECT_HI_TYPE
{
  void (*Exec)(const char *filename);
  int Id;
  int Nbin;
  int Level;
  int NumVars;
  const int *Vars;
};
extern struct NG_OBJECT_HI_TYPE ngHalo2IFrIT;


/*
//  Compute proximity zones for a specific halo with id=<Id>, 
//  or for all halos, if <Id>=0. Only consider halos resolved 
//  to at least <Level> and use <Nside> parameter in Healpix.
*/
struct NG_OBJECT_PZ_TYPE
{
  void (*Exec)(const char *filename);
  int Id;
  int Level;
  int Nside;
};
extern struct NG_OBJECT_PZ_TYPE ngProximityZone;


/*
//  Output radiation field as a function of SFR for all levels at or 
//  below <Level>.
*/
struct NG_OBJECT_RS_TYPE
{
  void (*Exec)(const char *fileroot);
  int Level;
};
extern struct NG_OBJECT_RS_TYPE ngRFvsSFR;


/*
//  Dump the Kennicutt-Schmidth relation on scale <LengthScale> kpc,
//  averaged over <TimeScale> Myr, and limit stellar age in output M*
//  to <StellarAgeLimit> Myr.
*/
struct NG_OBJECT_OKS_TYPE
{
  void (*Exec)(const char *filename);
  float TimeScale;
  float LengthScale;
  float StellarAgeLimit;
};
extern struct NG_OBJECT_OKS_TYPE ngObservedKSR;


/*
//  Same as above, but compute instantaneous SFR.
*/
struct NG_OBJECT_IKS_TYPE
{
  void (*Exec)(const char *filename);
  float LengthScale;
};
extern struct NG_OBJECT_IKS_TYPE ngInstantaneousKSR;


/*
//  Dump ids of all particles within <Rmax>*Rvir of a 
//  given halo with id <Id>.
*/
struct NG_OBJECT_HP_TYPE
{
  void (*Exec)(const char *filename);
  int Id;
  float Rmax;
};
extern struct NG_OBJECT_HP_TYPE ngHaloParticles;


/*
//  Dump properties of stellar particles within <Rmax>*Rvir
//  of a given halo with id <Id>.
*/
struct NG_OBJECT_HS_TYPE
{
  void (*Exec)(const char *filename);
  int Id;
  float Rmax;
};
extern struct NG_OBJECT_HS_TYPE ngHaloStars;

#endif  /* __ANALYSIS_NYG_H__ */
