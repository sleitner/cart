#ifndef __RT_OTVET_H__
#define __RT_OTVET_H__

#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)


void rtInitRunTransferOtvet();
void rtStepBeginTransferOtvet();
void rtAfterAssignDensityTransferOtvet(int level, int num_level_cells, int *level_cells);
void rtLevelUpdateTransferOtvet(int level, MPI_Comm local_comm);

#endif

#endif  /* __RT_OTVET_H__ */
