#ifndef __EXT_CHULL_H__
#define __EXT_CHULL_H__

void chAddPoints(int n, const int *pos);
void chReset();
void chConstruct();
void chGetLimits(int min[], int max[]);
int  chIsPointInside(int pos[]);

#endif  /* __EXT_CHULL_H__ */
