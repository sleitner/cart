
#ifndef __FRT_H__
#define __FRT_H__


#define frtOFFSET                0
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
void frtCall(cooloff)(frt_real *rPar, frt_real *rRadField0, frt_real *rRadField1, frt_real *rTime, frt_real *rVar, frt_intg *info);
void frtCall(initrun)();
void frtCall(initrun2)(frt_real *Yp, frt_real *Tmin, frt_real *D2Gmin, frt_real *ClumpH2, frt_real *CohLenH2, frt_real *fGal, frt_real *fQSO, frt_intg *IPOP, frt_intg *IREC);
void frtCall(stepbegin)(frt_real *uDen, frt_real *uLen, frt_real *uTime, frt_real *dtStep, frt_real *aExp, frt_real *Hubble, frt_real *rfAvg);
void frtCall(stepend)(frt_real *vol, frt_real *parAvg, frt_real *abcAvg);
void frtCall(updatetables)(frt_real *rfAvg);

frt_real frtCall(ein)(frt_real *tem, frt_real *y);
frt_real frtCall(tem)(frt_real *Ein, frt_real *y);
frt_real frtCall(gamma)(frt_real *Ein, frt_real *y);
frt_real frtCall(dusttogas)(frt_real *par);
frt_real frtCall(logwithunits)(frt_real *val, frt_intg *pDen, frt_intg *pLen, frt_intg *pTime);

void frtCall(setglobalradiationfields)();

void frtCall(getphotorates)(frt_real *par, frt_real *rf, frt_intg *itab, frt_real *y0, frt_real *pRate);
void frtCall(getbackgroundphotorates)(frt_real *pRate);
frt_real frtCall(getradiationfield)(frt_intg *lr, frt_real *rf);
frt_real frtCall(getbackgroundradiationfield)(frt_intg *lr);
frt_intg frtCall(getbinid)(frt_real *wlen);
frt_real frtCall(getbinwavelength)(frt_intg *lr);

void frtCall(packradiationfields)(frt_intg *n, frt_real *data);
void frtCall(unpackradiationfields)(frt_intg *n, frt_real *data);

#ifdef RT_TRANSFER
void frtCall(initruntransfer)(frt_intg *nfreq);
void frtCall(stepbegintransfer)(frt_real *abcAvg);
void frtCall(preparetablestransfer)(frt_real *rfAvg);
#if (RT_CFI == 1)
void frtCall(transfercomputecellabs)(frt_intg *L, frt_real *Zsol, frt_real *denB, frt_real *denH1, frt_real *denG1, frt_real *denG2, frt_real *denMH, frt_real *abc, frt_real *abc1);
#else
void frtCall(transfercomputecellabs)(frt_intg *L, frt_real *Zsol, frt_real *denB, frt_real *denH1, frt_real *denG1, frt_real *denG2, frt_real *denMH, frt_real *abc);
#endif
void frtCall(transferpackradiationfield)(frt_real *par, frt_real *y0, frt_real *rawRF0, frt_real *rawRF1, frt_real *rf);
#endif /* RT_TRANSFER */

#ifdef RT_TEST
void frtCall(testmodifytimestep)(frt_real *dt, frt_intg *ns);
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
}
frtCall(setrun);

#endif
