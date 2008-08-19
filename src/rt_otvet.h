#ifndef __RT_OTVET_H__
#define __RT_OTVET_H__

#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)


void rtInitRunTransferOtvet();
void rtStepBeginTransferOtvet();
void rtLevelUpdateTransferOtvet(int level, MPI_Comm local_comm);


#endif

#endif  /* __RT_OTVET_H__ */
