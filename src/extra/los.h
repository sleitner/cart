#ifndef __EXT_LOS_H__
#define __EXT_LOS_H__


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


/*
//  Traverse a segment of a LOS located on a single processor,
//  for a LOS starting at pos0, in the direction (theta,phi), 
//  of maximum length len, not deeper than floor_level, and call
//  worker(cell,r1,r2,data) at each step inside this processor domain.
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
void losTraverseSegment(int id, double pos0[3], double theta, double phi, double len, int floor_level, int (*worker)(int id, int cell, double r1, double r2, losBuffer data), losSegment *segment);


/*
//  Collect all LOS segments from different processors on the master node
//  and broadcast their buffer data into back into result. A user-supplied 
//  function collector is used to assemble separate segments into a single LOS.
*/
void losCollectSegments(losBuffer *result, losSegment *segment, void (*collector)(losBuffer *result, losSegment *next));


void losTraverseSky(int nside, double pos0[3], double len, int floor_level, losBuffer *lines, int (*worker)(int id, int cell, double r1, double r2, losBuffer data), void (*collector)(losBuffer *result, losSegment *next));

#endif  /* __EXT_LOS_H__ */
