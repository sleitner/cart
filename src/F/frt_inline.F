#include "defs.h"
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"    
C
C *********************************************************
C *
C *  This file contains functions which need to be inlined
C *  (if possible) for better performance
C *
C *********************************************************
C
C  Functions for frequency dependence fitting: direct functions
C 
      function QFuni(x)

*      if(x .lt. 1.0e-8) then
*         QFuni = 1.0
*      else
*         y = sqrt(x)
*         QFuni = (1.0-exp(-y))/y
*      endif

      QFuni = exp(-x)

      return
      end
C
      function QFloc(x)

      QFloc = exp(-x)

      return
      end
C
C  Inverse functions
C
      function QIuni(x)

*      parameter(A = 3.970)
*      parameter(B = 0.712)
*      parameter(ALP = 0.994351)
*      parameter(BET = 0.740)
*      parameter(GAM = 0.6889336)
*      if(x.le.0.0 .or. x.ge.1.0) then
*         QIuni = 0.0
*      else
*         QIuni = ((1-exp(-(1-x**GAM)/x*(1+A*x**ALP)/
*     .     (1+B*x**BET)))/x)**2
*      endif

      QIuni = -log(1.0e-35+x)

      return
      end
C
      function QIloc(x)

      QIloc = -log(1.0e-35+x)

      return
      end
#ifdef RT_TRANSFER
C
C  Absoprtion coefficient for 1 cell at frequency bin lr.
C
#if (RT_CFI == 1)
      subroutine frtTransferComputeCellAbs(L,Zsol,denB,denH1,denG1,
     .     denG2,denMH,abc,abc1)
#else
      subroutine frtTransferComputeCellAbs(L,Zsol,denB,denH1,denG1,
     .     denG2,denMH,abc)
#endif
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_transfer.inc'
C
      lr = LRSETS(L)
C
C  Adopt an approximation that H2 absorbs ionizing radiation as two HI
C
      denHA = denH1 + 2*denMH
C
#ifdef RT_DUST
      abc = uColumn*Zsol*csDU(lr)*denHA
#else
      abc = 0.0
#endif

#if (RT_CFI == 0)
      abc = abc + uColumn*(csH1(lr)*denHA+csG1(lr)*denG1+
     .     csG2(lr)*denG2)
#endif

#if (RT_CFI == 1)
      abc1 = abc + uColumn*(csH1(lr)*denHA+csG1(lr)*denG1+
     .     csG2(lr)*denG2)
#endif

#if (RT_CFI==1 || RT_CFI==2)
      abc = abc + uColumn*(
     .     csH1(lr)*denH1*CIonH1(denB,denHA/denB,denG1/denB,denG2/denB)
     .     csG1(lr)*denG1*CIonG1(denB,denHA/denB,denG1/denB,denG2/denB)
     .     csG2(lr)*denG2*CIonG2(denB,denHA/denB,denG1/denB,denG2/denB)
#endif

      return
      end
#endif  // RT_TRANSFER
#endif  // RADIATIVE_TRANSFER
      subroutine frtDummyInline(i)
      i = 0
      end
