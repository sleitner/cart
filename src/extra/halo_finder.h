#ifndef __HALO_FINDER_H__
#define __HALO_FINDER_H__


typedef struct hfHaloType
{
  int Id;
  float Mass;
  float Vmax;
  float Rvir;
  float Vel[3];
  double Pos[3];
}
hfHalo;


/*
//  Read an hlist file produced by HFIND and return the list of halos.
//  Limit the set of halos returned by specifying min number of members, mass, vmax, and/or rvir.
//  This function must be called by all procs.
*/
int hfReadHFINDHalos(const char* fname, int nmemMin, float massMin, float vmaxMin, float rvirMin, hfHalo **list, int *list_size);


/*
//  Return the level at the center of the halo.
//  This function must be called by all procs.
*/
int hfHaloLevel(const hfHalo *h);


#endif  /* __HALO_FINDER_H__ */
