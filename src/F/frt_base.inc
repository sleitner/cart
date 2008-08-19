C
C  Fundamental parameters - should not be changed
C
	parameter (TH1 = 157807.0)
	parameter (TG1 = 285335.0)
	parameter (TG2 = 631515.0)
C
	parameter (CFSTAR = 1.381e-16)
	parameter (CSSTAR = 9.606e-24)
	parameter (T1eV = 11604.5)
	parameter (fnuH1 = TH1/T1eV)
C
	parameter (GFF = 1.3)
	parameter (facFF = (1.381e-16/CFSTAR)*1.032e-11*GFF)
C
	parameter (TAUMAX = 23.0)
	parameter (XVLOW = 1.0e-10)
C
	parameter (SORTOL = 0.0)
C
C  Equation indecies
C
        parameter(ieqEin = 1)
        parameter(ieqXH1 = 2)
        parameter(ieqXH2 = 3)
        parameter(ieqXG1 = 4)
        parameter(ieqXG2 = 5)
        parameter(ieqXG3 = 6)
        parameter(ieqXMH = 7)
        parameter(ieqXHp = 8)
        parameter(ieqXHm = 9)
        parameter(ieqDim = 9)
C
        parameter(ioCXH1 =  1)
        parameter(ioCXG1 =  2)
        parameter(ioCXG2 =  3)
        parameter(ioRAH2 =  4)
        parameter(ioRAG2 =  5)
        parameter(ioRAG3 =  6)
        parameter(ioRBH2 =  7)
        parameter(ioRBG2 =  8)
        parameter(ioRBG3 =  9)
        parameter(ioDust = 10)
        parameter(ioBrem = 11)
        parameter(ioTemp = 12)
        parameter(ioMax1 = 11)
        parameter(ioTabH = 13)
        parameter(ioTabM = 14)
        parameter(ioInfo = 15)
        parameter(ioDim  = 15)
C
C  Parameter indecies
C  (yes, iparZsol = iparDGas)
C
        parameter(iparRhoB =  1)
        parameter(iparZsol =  2)
        parameter(iparDGas =  2)
        parameter(iparDTsh =  3)
        parameter(iparDTad =  4)
        parameter(iparBias =  5)
        parameter(iparDFac =  6)
        parameter(iparSobL =  7)
        parameter(iparNumF =  8)
        parameter(iparCell =  9)
        parameter(iparRecB = 10)
        parameter(iparSDp1 = 11)
        parameter(iparSDp2 = 12)
        parameter(iparSDp3 = 13)
        parameter(iparCFH1 = 14)
        parameter(iparCFG1 = 15)
        parameter(iparCFG2 = 16)
        parameter(iparDelN = 17)
        parameter(iparDeb  = 18)
        parameter(iparDim  = 18)
C
C  Diagnostic info
C
	parameter(infoTabH = 1)
	parameter(infoTabM = 2)
	parameter(infoNAvg = 3)
	parameter(infoNMax = 4)
	parameter(infoFail = 5)
	parameter(infoCall = 6)
	parameter(infoDim  = 6)
C
        common/RT_BASE_A/ XH, XG, XFEMAX
        common/RT_BASE_B/ uDen, uLen, uTime, uRate, uColumn
        common/RT_BASE_C/ facPi, facPh, facCC, facH2CMB, facCool,
     .                    facLX, facLya, facRecB
        common/RT_BASE_D/ modeAB
        common/RT_BASE_E/ aStep, hStep, tStep, TemCMB
C
        common/RT_BASE_F/ yLimit(ieqDim), ySmall(ieqDim)

