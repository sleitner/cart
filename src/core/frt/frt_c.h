
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


#ifdef __cplusplus
extern "C" {
#endif

/* 
//  Fortran interface 
*/
void frtCall(initvar)(frt_real *var);
void frtCall(cooloff)(frt_real *var, frt_real *rawrf, frt_real *time, frt_intg *info);
void frtCall(initrun)();
void frtCall(initrun2)(frt_real *Yp, frt_real *Tmin, frt_real *D2GminH2, frt_real *CluFacH2, frt_real *CohLenH2, frt_real *fGal, frt_real *fQSO, frt_intg *IPOP, frt_intg *IREC, frt_intg *IOUNIT);
void frtCall(stepbegin)(frt_real *uDen, frt_real *uLen, frt_real *uTime, frt_real *aExp, frt_real *Hubble, frt_real *daExp, frt_real *dtCode);
void frtCall(stepend)(frt_real *vol, frt_real *par, frt_real *abcLoc, frt_real *abcUni);
void frtCall(updatetables)(frt_real *rfAvg);

frt_real frtCall(tem)(frt_real *var);
void frtCall(getcoolingrate)(frt_real *var, frt_real *rawrf, frt_real *rateCool, frt_real *rateHeat);
void frtCall(getphotorates)(frt_real *var, frt_real *rawrf, frt_real *pRate);
void frtCall(getphotoratesfs)(frt_real *var, frt_real *rawrf, frt_real *pRate);
void frtCall(getbackgroundphotorates)(frt_real *pRate);
void frtCall(getradiationfield)(frt_real *var, frt_real *rawrf, frt_intg *nout, frt_real *wlen, frt_real *ngxi);
void frtCall(getradiationfieldfs)(frt_real *var, frt_real *rawrf, frt_intg *nout, frt_real *wlen, frt_real *ngxi);
void frtCall(getbackgroundradiationfield)(frt_intg *nout, frt_real *wlen, frt_real *ngxi);

void frtCall(packradiationfields)(frt_intg *n, frt_real *data);
void frtCall(unpackradiationfields)(frt_intg *n, frt_real *data);

#ifdef RT_TRANSFER
void frtCall(initruntransfer)(frt_intg *nfreq, frt_intg *ncomp);

frt_real frtCall(transferglobalqterm)(frt_intg *idx);

#if (RT_CFI > 0)
void frtCall(transfercomputecellabs)(frt_intg *idx, frt_real *denHI, frt_real *denHeI, frt_real *denHeII, frt_real *denH2, frt_real *denDU, frt_real *dx, frt_real *abc, frt_real *var);
#else
void frtCall(transfercomputecellabs)(frt_intg *idx, frt_real *denHI, frt_real *denHeI, frt_real *denHeII, frt_real *denH2, frt_real *denDU, frt_real *dx, frt_real *abc);
#endif
#endif /* RT_TRANSFER */

#ifdef RT_TEST
void frtCall(testmodifytimestep)(frt_real *dt, frt_intg *ns);
#endif /* RT_TEST */

extern struct
{
  frt_real Yp;
  frt_real Tmin;
  frt_real D2GminH2;
  frt_real CluFacH2;
  frt_real CohLenH2;
  frt_real fGal;
  frt_real fQSO;
  frt_intg IPOP;
  frt_intg IREC;
  frt_intg IOUNIT;
}
frtCall(setrun);

#ifdef __cplusplus
}
#endif

#endif
