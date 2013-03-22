#ifndef __EXT_LOS_H__
#define __EXT_LOS_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/*
//  Buffer for keeping information used during the LOS traversal
*/
typedef struct 
{
  int Size;
  void *Data;
} 
losBuffer;


/*
//  Segment for a LOS: a buffer and a range of distances that it represents
*/
typedef struct 
{
  double Range[2];
  losBuffer Buffer;
  int Id, WorkerReturnCode;
}
losSegment;

typedef int (*losWorkerCallback)(int id, int cell, double r1, double r2, losBuffer data);
typedef void (*losCollectorCallback)(const losBuffer *result, int num_segments, const losSegment *segments);


/*
//  Traverse a segment of a LOS located on a single processor,
//  for a LOS starting at pos0, in the direction (theta,phi), 
//  of maximum length len, not deeper than floor_level, and call
//  worker(id,cell,r1,r2,data) at each step inside this processor domain.
//  Worker arguments:
//    id:    integer id of this LOS (to distinguish different LOS
//           inside an OpenMP loop)
//    cell:  current cell crossed by the LOS
//    r1:    starting position of the cell along the LOS
//    r2:    ending position of the cell along the LOS
//    data:  buffer with optional external data
//  worker should return 0 if the traversal to continue; otherwise, 
//  it will be stopped. The worker return value is saved as segment->WorkerReturnCode.
//  The output is the filled *segment argument.
//  This function is thread-safe, can be called from an OpenMP-parallel loop
*/
void losTraverseSegment(int id, double pos0[3], double theta, double phi, double len, int floor_level, losWorkerCallback worker, losSegment *segment);


/*
//  Collect all LOS segments from different processors on the master node
//  and broadcast their buffer data into back into result. A user-supplied 
//  function collector is used to assemble separate segments into a single LOS.
//  Keep in mind that separate segments are not necessarily continues, as a
//  single line-of-sight can cross a given domain more than once.
//  This function calls MPI inside and is manifestly thread-unsafe.
//  The buffer for result needs to be pre-allocated, this function does not allocate memory.
*/
void losCollectSegments(const losBuffer *result, const losSegment *segment, losCollectorCallback collector);


/*
//  A useful wrapper over the two previous functions: samples the sky using 
//  HealPIX binning with nside bins per section (total number of rays is 
//  12*nside^2) with common origin pos0, of maximum length len, not deeper 
//  than floor_level. The user-supplied functions worker and collector are
//  described above. The final result is returned in lines[].
*/
void losTraverseSky(int nside, double pos0[3], double len, int floor_level, losBuffer *lines, losWorkerCallback worker, losCollectorCallback collector);

#endif  /* __EXT_LOS_H__ */
