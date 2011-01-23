#ifndef __EXT_CHULL_H__
#define __EXT_CHULL_H__

void chAddPoints(int n, const int *pos);
void chReset();
void chGetLimits(int min[], int max[]);
int  chIsPointInside(int pos[]);

/*
//  Construct the exact hull
*/
void chMakeFullHull();

/*
//  Construct the approximate hull of minimum possible volume that
//  excludes at most <tolnum> fraction of all points, achieves at 
//  least <tolvol> fractional reduction in the hull volume, and
//  makes at most <numits> iterations.
*/
void chMakeHull(float tolnum, float tolvol, int numits, int loud);

#endif  /* __EXT_CHULL_H__ */
