C
C  The following parameters are chosen after much testing. If you decide
C  to change them, make sure you know what you are doing.
C
	parameter (NRAD = 300)
C
	parameter (NTEM = 201)
	parameter (acTmin = 1.0e-1)
	parameter (acTmax = 1.0e+9)
C
	parameter (NOPT = 82)
	parameter (acOmin = 1.0e-2)
	parameter (acOmax = 1.0e+6)
C
C  Derived parameters
C
	parameter (NOPT1 = NOPT+1)
	parameter (NOPTL = NOPT1, NOPTT = NOPT1)
	parameter (NOPT3 = NOPTL*NOPTT*NOPTT)
	parameter (NTOL = 100*NOPT)
C
C  Radiation mesh
C
	parameter (TXI0 = 1.0e+4)
	parameter (XISTEP = 15.0/NRAD)
C
C  Temperature tables: base
C
	common/RT_TAB_A/ acLT(NTEM), acTEM(NTEM), acSQT(NTEM),
     .	                 acLTmin, acLTmax, acLTstp
C
C  Temperature tables: atomic rates
C
	parameter(irateCIiH1 =  1)
	parameter(irateCIiG1 =  2)
	parameter(irateCIiG2 =  3)
	parameter(irateDRiG2 =  4)
	parameter(irateRAiH2 =  5)
	parameter(irateRAiG2 =  6)
	parameter(irateRAiG3 =  7)
	parameter(irateRBiH2 =  8)
	parameter(irateRBiG2 =  9)
	parameter(irateRBiG3 = 10)
	parameter(irateCIcH1 = 11)
	parameter(irateCIcG1 = 12)
	parameter(irateCIcG2 = 13)
	parameter(irateDRcG2 = 14)
	parameter(irateRAcH2 = 15)
	parameter(irateRAcG2 = 16)
	parameter(irateRAcG3 = 17)
	parameter(irateRBcH2 = 18)
	parameter(irateRBcG2 = 19)
	parameter(irateRBcG3 = 20)
	parameter(irateLCcH1 = 21)
	parameter(irateLCcG2 = 22)
	parameter(irateAtomN = 22)

	common/RT_TAB_B/ RateAtom(irateAtomN,NTEM)
C
C  Temperature tables: chemical rates
C
	parameter(irateMHi01 =  1)
	parameter(irateMHi02 =  2)
	parameter(irateMHi03 =  3)
	parameter(irateMHi04 =  4)
	parameter(irateMHi05 =  5)
	parameter(irateMHi06 =  6)
	parameter(irateMHi07 =  7)
	parameter(irateMHi08 =  8)
	parameter(irateMHi09 =  9)
	parameter(irateMHi10 = 10)
	parameter(irateMHi11 = 11)
	parameter(irateMHi14 = 12)
	parameter(irateMHi15 = 13)
	parameter(irateMHi16 = 14)
	parameter(irateMHi21 = 15)
	parameter(irateMHi22 = 16)
	parameter(irateMHi23 = 17)
	parameter(irateMHi24 = 18)
	parameter(irateMHi25 = 19)
	parameter(irateMHi26 = 20)
	parameter(irateMHi27 = 21)
	parameter(irateMHi28 = 22)
	parameter(irateMHi29 = 23)
	parameter(irateMHiDU = 24)
	parameter(irateDuGro = 25)
	parameter(irateDuSpu = 26)
	parameter(irateChemN = 26)
C
C  Temperature tables: chemical cooling rates
C
	parameter(irateMHcH1o =  1)
	parameter(irateMHcH1p =  2)
	parameter(irateMHcMOo =  3)
	parameter(irateMHcMOp =  4)
	parameter(irateMHcMPo =  5)
	parameter(irateMHcMPp =  6)
	parameter(irateMHcG1o =  7)
	parameter(irateMHcG1p =  8)
	parameter(irateMHcH2o =  9)
	parameter(irateMHcH2p = 10)
	parameter(irateMHcFEo = 11)
	parameter(irateMHcFEp = 12)
	parameter(irateMHc08  = 13)
	parameter(irateMHc09  = 14)
	parameter(irateMHc10  = 15)
	parameter(irateMHc14  = 16)
	parameter(irateMHc15  = 17)
	parameter(irateMHcDU  = 18)
	parameter(irateChemN2 = 18)

	common/RT_TAB_C/ RateChem(irateChemN,NTEM), 
     .			 RateChem2(irateChemN2,NTEM)
	common/RT_TAB_D/ temZLT, d_pPAH, d_sPAH, d_sDUST
C
C  Temperature tables: high density rates
C
	parameter(irateHDi09 =  1)
	parameter(irateHDi10 =  2)
	parameter(irateHDi30 =  3)
	parameter(irateHDi31 =  4)
	parameter(irateHDi32 =  5)
	parameter(irateHighN =  5)
	parameter(irateHDc02  =  1)
	parameter(irateHDcNC  =  2)
	parameter(irateHighN2 =  2)

	common/RT_TAB_E/ RateHigh(irateHighN,NTEM), 
     .			 RateHigh2(irateHighN2,NTEM)
C
C  Radiation field tables: base (order DOES matter here)
C
	parameter(iptabPhG2 =  1)
	parameter(iptabPiG2 =  2)
	parameter(iptabPhG1 =  3)
	parameter(iptabPiG1 =  4)
	parameter(iptabPhH1 =  5)
	parameter(iptabPiH1 =  6)
	parameter(iptabPiC6 =  7)
	parameter(iptabCi27 =  8)
	parameter(iptabCi28 =  9)
	parameter(iptabCi29 = 10)
	parameter(iptabCi30 = 11)
	parameter(iptabCi31 = 12)
	parameter(iptabCi32 = 13)
	parameter(iptabCiLW = 14)
	parameter(iptabIon  =  7)
	parameter(iptabOT   =  8)
	parameter(iptabH1   =  5)
	parameter(iptabG1   =  3)
	parameter(iptabG2   =  1)
	parameter(iptabDim  = 14)

	common/RT_TAB_F/ acLO(NOPT), acOPT(NOPT), acTAU(NOPT),
     .	                 acLOmin, acLOmax, acLOstp
	common/RT_TAB_G/ acT2O(NTOL), acO2T(NTOL), acTOstp
	common/RT_TAB_H/ tauH1m, tauG1m, tauG2m,
     .	                 altH1m, altG1m, altG2m,
     .		         tH1min, tG1min, tG2min,
     .		         tH1max, tG1max, tG2max,
     .                   iptabMax
	common/RT_TAB_I/ itloc2(NOPTT,NOPTT)
C
C  Radiation field tables: global field
C
	common/RT_TAB_J/ pRate0(iptabDim), 
     .			 pTab3(iptabDim,NOPT3)
C
C  Radiation field tables: local field
C
	common/RT_TAB_K/ qRate0(iptabDim), 
     .			 qTab3(iptabDim,NOPT3)
C
C  Integration weights and normalization
C
	common/RT_TAB_L/ wpTab(iptabDim,NRAD), uvFac(iptabDim)
C
C  Radiation mesh
C
	common/RT_TAB_M/ xig(NRAD), Txi(NRAD), guvEff(NRAD),
     .                   tauLoc(NRAD)
	common/RT_TAB_N/ angUni(NRAD), angLoc(NRAD),
     .	                 angAvg(NRAD), angEff(NRAD),
     .	                 sxiUni(NRAD)
C
C  Cross sections: atomic & dust
C
	common/RT_TAB_O/ csH1(NRAD), csH1n(NRAD),
     .			 csG1(NRAD), csG1n(NRAD),
     .			 csG2(NRAD), csG2n(NRAD),
     .	                 csUV(NRAD), csUVn(NRAD),
     .			 csDust(NRAD),
     .			 csH1tr, csG1tr, csG2tr, csUVtr,
     .                   csH1AB, csG1AB, csG2AB
	common/RT_TAB_P/ lrH1tr, lrG1tr, lrG2tr, lrUVtr,
     .                   lrH1rB, lrG1rB, lrG2rB,
     .                   lrmax, lrDUtr, lrSEtr
C
C  Cross sections: chemistry 
C  (LW = LW bands, 31 = 32 para, 32 = 32 ortho)
C
	common/RT_TAB_Q/ csMH27(NRAD), csMH28(NRAD), csMH29(NRAD),
     .                   csMH30(NRAD), csMH31(NRAD), csMH32(NRAD),
     .                   csMHLW(NRAD)
	common/RT_TAB_R/ lr27tr, lr28tr, lr29tr, 
     .                   lr30tr, lr31tr, lr32tr, lrLWtr
C
C  Cross sections: meatl species for CF fits
C
	common/RT_TAB_S/ csC6(NRAD), csC6tr, lrC6tr
