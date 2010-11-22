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
C  Parameter indecies
C  (yes, iparZsol = iparDGas)
C  (yes, iparCool = iparSDp1)
C  (yes, iparHeat = iparSDp2)
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
        parameter(iparVolB = 10)
        parameter(iparRecB = 11)
        parameter(iparSDp1 = 12)
        parameter(iparCool = 12)
        parameter(iparSDp2 = 13)
        parameter(iparHeat = 13)
        parameter(iparSDp3 = 14)
        parameter(iparCFH1 = 15)
        parameter(iparCFG1 = 16)
        parameter(iparCFG2 = 17)
        parameter(iparDelN = 18)
        parameter(iparDeb  = 19)
        parameter(iparDim  = 19)
C
C  Globally-averaged quantities
C
        parameter(iavgCXH1 =  1)
        parameter(iavgCXG1 =  2)
        parameter(iavgCXG2 =  3)
        parameter(iavgRAH2 =  4)
        parameter(iavgRAG2 =  5)
        parameter(iavgRAG3 =  6)
        parameter(iavgRBH2 =  7)
        parameter(iavgRBG2 =  8)
        parameter(iavgRBG3 =  9)
        parameter(iavgDust = 10)
        parameter(iavgBrem = 11)
        parameter(iavgNorm = 12)
        parameter(iavgTemp = 13)
        parameter(iavgMax1 = 11)
        parameter(iavgMax2 = 12)
        parameter(iavgDim  = 13)
C
C  Diagnostic info
C
	parameter(infoTabH = 1)
	parameter(infoTabM = 2)
	parameter(infoNAvg = 3)
	parameter(infoFail = 4)
	parameter(infoCall = 5)
	parameter(infoMax1 = 5)
	parameter(infoNMax = 6)
	parameter(infoDim  = 6)
C
        common/RT_BASE_A/ XH, XG, XFEMAX
        common/RT_BASE_B/ uRate, uColumn, TemCMB
        common/RT_BASE_C/ facPhot, facCC, facH2CMB, facCool,
     .                    facLX, facLya, facRecB
        common/RT_BASE_D/ yLimit(ieqDim), ySmall(ieqDim)
	common/RT_BASE_E/ uDen, uLen, uTime, aExp, Hubble, dtStep
C
	logical ioNode
	common/RT_BASE_F/ ioNode
C
	include 'frt_f77.inc'


