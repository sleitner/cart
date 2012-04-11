
#ifndef __FRT_H__
#define __FRT_H__


#define __FRT_OFFSET                0
#include "frt_index.h"


/*
//  Types for C - Fortran interface
*/
#define frt_intg        int
#define frt_real        float

#define frtCall(fun)    frt##fun##_


/* 
//  Fortran interface 
*/
void frtCall(initvar)(frt_real *var);
void frtCall(cooloff)(frt_real *var, frt_real *rawrf, frt_real *time, frt_intg *info);
void frtCall(initrun)();
void frtCall(initrun2)(frt_real *Yp, frt_real *Tmin, frt_real *D2Gmin, frt_real *ClumpH2, frt_real *CohLenH2, frt_real *fGal, frt_real *fQSO, frt_intg *IPOP, frt_intg *IREC, frt_intg *IOUNIT);
void frtCall(stepbegin)(frt_real *uDen, frt_real *uLen, frt_real *uTime, frt_real *aExp, frt_real *Hubble, frt_real *daExp);
void frtCall(stepend)(frt_real *dt, frt_real *vol, frt_real *par);
void frtCall(updatetables)(frt_real *rfAvg);

frt_real frtCall(tem)(frt_real *var);
void frtCall(getcoolingrate)(frt_real *var, frt_real *rawrf, frt_real *rateCool, frt_real *rateHeat);
void frtCall(getphotorates)(frt_real *var, frt_real *rawrf, frt_real *pRate);
void frtCall(getbackgroundphotorates)(frt_real *pRate);
void frtCall(getradiationfield)(frt_real *var, frt_real *rawrf, frt_intg *nout, frt_real *wlen, frt_real *ngxi);
void frtCall(getbackgroundradiationfield)(frt_intg *nout, frt_real *wlen, frt_real *ngxi);

void frtCall(packradiationfields)(frt_intg *n, frt_real *data);
void frtCall(unpackradiationfields)(frt_intg *n, frt_real *data);

#ifdef RT_TRANSFER
void frtCall(initruntransfer)(frt_intg *nfreq, frt_intg *ncomp);
void frtCall(stependtransfer)(frt_real *abcAvg);

frt_real frtCall(transferglobalac)(frt_intg *idx, frt_real *abc);

frt_real frtCall(getrfwithunits)(frt_intg *freq, frt_real *rfNear, frt_real *rfFar);

#ifdef RT_ABSORPTION_CALLBACK_FULL
void frtCall(transfercomputecellabs)(frt_intg *idx, frt_real *denB, frt_real *denHI, frt_real *denHeI, frt_real *denHeII, frt_real *denMH, frt_real *dx, frt_real *abc, frt_real *var);
#else
void frtCall(transfercomputecellabs)(frt_intg *idx, frt_real *denB, frt_real *denHI, frt_real *denHeI, frt_real *denHeII, frt_real *denMH, frt_real *dx, frt_real *abc, frt_real *Zsol);
#endif
#endif

#ifdef RT_TEST
void frtCall(testmodifytimestep)(frt_real *dt, frt_intg *ns);
#endif

#ifdef RT_CUSTOM_DUST_TO_GAS
frt_real frtCall(dusttogas)(frt_real *var);
#endif

extern struct
{
  frt_real Yp;
  frt_real Tmin;
  frt_real D2Gmin;
  frt_real ClumpH2;
  frt_real CohLenH2;
  frt_real fGal;
  frt_real fQSO;
  frt_intg IPOP;
  frt_intg IREC;
  frt_intg IOUNIT;
}
frtCall(setrun);

#endif
