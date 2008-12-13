#ifndef __RT_TRANSFER_H__
#define __RT_TRANSFER_H__


#ifdef RT_TRANSFER


/*
// Global quantities 
*/

struct rtArrayAverageData
{
  float Value;
  float LevelSum[num_refinement_levels+1];
};

extern struct rtArrayAverageData rt_glob_Avg[2];

#define rt_source_Avg     rt_glob_Avg[0]
#define rt_num_glob       1

#ifdef RT_VAR_OT_FIELD
#define rt_ot_field_Avg   rt_glob_Avg[1]
#undef rt_num_glob
#define rt_num_glob       2
#endif


void rtInitRunTransfer();
void rtUpdateTablesTransfer();
void rtStepBeginTransfer();
void rtStepEndTransfer();
void rtAfterAssignDensityTransfer(int level, int num_level_cells, int *level_cells);
void rtLevelUpdateTransfer(int level, MPI_Comm local_comm);
void rtComputeAbsLevel(int ncells, int *cells, int ifreq, float **abc);


#ifdef RT_SINGLE_SOURCE
extern float rtSingleSourceVal;
extern double rtSingleSourcePos[nDim];
#endif


#endif

#endif  /* __RT_TRANSFER_H__ */
