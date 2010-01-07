#ifndef __RT_IO_H__
#define __RT_IO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


void rtWriteRadiationFieldData(const char *fileroot, int fortran_style); 
void rtReadRadiationFieldData(const char *fileroot, int fortran_style); 

#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_IO_H__ */
