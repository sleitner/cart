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
	parameter (NRADL = NRAD+1)
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
	parameter(irateMHi07 =  1)
	parameter(irateMHi08 =  2)
	parameter(irateMHi09 =  3)
	parameter(irateMHi10 =  4)
	parameter(irateMHi11 =  5)
	parameter(irateMHi12 =  6)
	parameter(irateMHi13 =  7)
	parameter(irateMHi14 =  8)
	parameter(irateMHi15 =  9)
	parameter(irateMHi16 = 10)
	parameter(irateMHi17 = 11)
	parameter(irateMHi18 = 12)
	parameter(irateMHi19 = 13)
	parameter(irateMHi20 = 14)
	parameter(irateMHi21 = 15)
	parameter(irateMHiDU = 16)
	parameter(irateChemC = 16)
	parameter(irateMHcLV = 17)
	parameter(irateMHcLR = 18)
	parameter(irateMHcUE = 19)
	parameter(irateMHcUH = 20)
	parameter(irateMHcUM = 21)
	parameter(irateMHc16 = 22)
	parameter(irateMHc17 = 23)
	parameter(irateMHcDU = 24)
	parameter(irateZLcTX = 25)
	parameter(irateZLcTA = 26)
	parameter(irateZHcTA = 27)
	parameter(irateZHcTB = 28)
	parameter(irateZHcTC = 29)
	parameter(irateChemN = 29)

	common/RT_TAB_C/ RateChem(irateChemN,NTEM)
	common/RT_TAB_D/ Tmin, temZLT, d_clump, d_pPAH, d_sPAH
C
C  Temperature tables: high density rates
C
	parameter(irateHDi11 =  1)
	parameter(irateHDi15 =  2)
	parameter(irateHDi22 =  3)
	parameter(irateHighC =  3)
	parameter(irateHDcLV =  4)
	parameter(irateHDcLR =  5)
	parameter(irateHDc08 =  6)
	parameter(irateHighN =  6)

	common/RT_TAB_E/ RateHigh(irateHighN,NTEM)
C
C  Radiation field tables: base (order DOES matter here)
C
	parameter(iptabPhG2 =  1)
	parameter(iptabPiG2 =  2)
	parameter(iptabPhG1 =  3)
	parameter(iptabPiG1 =  4)
	parameter(iptabPhH1 =  5)
	parameter(iptabPiH1 =  6)
	parameter(iptabCi27 =  7)
	parameter(iptabCi28 =  8)
	parameter(iptabCi29 =  9)
	parameter(iptabCi30 = 10)
	parameter(iptabCi31 = 11)
	parameter(iptabCi32 = 12)
	parameter(iptabCi99 = 13)
	parameter(iptabIon  =  6)
	parameter(iptabOT   =  7)
	parameter(iptabH1   =  5)
	parameter(iptabG1   =  3)
	parameter(iptabG2   =  1)
	parameter(iptabDim  = 13)

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
	common/RT_TAB_J/ pRate0(iptabDim), pTab(iptabDim,NOPT3)
C
C  Radiation field tables: local field
C
	common/RT_TAB_K/ sgpRate0(iptabDim), sgpTab(iptabDim,NOPT3)
C
C  Integration weights and normalization
C
	common/RT_TAB_L/ wpTab(iptabDim,NRAD)
C
C  Radiation mesh
C
	common/RT_TAB_M/ tauhat_xi(NRADL), qfdhat_xi(NRADL),
     .                   xig(NRADL), ang(NRADL), Txi(NRADL)
	common/RT_TAB_N/ guveff(NRADL)
C
C  Cross sections: atomic & dust
C
	common/RT_TAB_O/ csH1(NRADL),  csG1(NRADL),  csG2(NRADL),
     .                   csH1n(NRADL), csG1n(NRADL), csG2n(NRADL),
     .			 csLyS(NRADL), csLySn(NRADL),
     .	                 csDust(NRADL), csDustn(NRADL),
     .			 csH1tr, csG1tr, csG2tr,
     .                   csH1AB, csG1AB, csG2AB
	common/RT_TAB_P/ lrH1tr, lrG1tr, lrG2tr,
     .                   lrH1rB, lrG1rB, lrG2rB,
     .                   lrmin, lrmax, lrLStr, lrDUtr, lrSEtr
C
C  Cross sections: chemistry 
C  (99 = LW bands, 31 = 32 para, 32 = 32 ortho)
C
	common/RT_TAB_Q/ csMH27(NRADL), csMH28(NRADL), csMH29(NRADL),
     .                   csMH30(NRADL), csMH31(NRADL), csMH32(NRADL),
     .                   csMH99(NRADL)
	common/RT_TAB_R/ lr27tr, lr28tr, lr29tr, 
     .                   lr30tr, lr31tr, lr32tr, lr99tr
