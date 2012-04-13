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
	parameter (VSMALL = 1.0e-35)
	parameter (VLARGE = 1.0e+35)
C
C  Variable indecies
C
        parameter(ivarEin  =  1)
        parameter(ivarXH1  =  2)
        parameter(ivarXH2  =  3)
        parameter(ivarXG1  =  4)
        parameter(ivarXG2  =  5)
        parameter(ivarXG3  =  6)
        parameter(ivarXMH  =  7)
        parameter(ivarXHp  =  8)
        parameter(ivarXHm  =  9)
        parameter(ivarRho  = 10)
        parameter(ivarZsol = 11)
        parameter(ivarSobL = 12)
        parameter(ivarNumF = 13)
        parameter(ivarRF0l = 14)
        parameter(ivarRF0g = 15)
        parameter(ivarDelV = 16)
        parameter(ivarDelX = 17)
        parameter(ivarBias = 18)
        parameter(ivarCSup = 19)
        parameter(ivarDTex = 20)
        parameter(ivarLTE  = 21)
        parameter(ivarDeb  = 22)
C
C  Diagnostic variables
C  (returned on exit, overwrite internals)
C
        parameter(ivarOut  = ivarDeb+1)
        parameter(ivarDGas = ivarOut+0)
        parameter(ivarGam  = ivarOut+1)
        parameter(ivarCool = ivarOut+2)
        parameter(ivarHeat = ivarOut+3)
C
C  Internal variables
C
        parameter(ivarInt  = ivarDeb+1)
*       parameter(ivarDGas = ivarInt+0)  !! already defined
        parameter(ivarRecB = ivarInt+1)
        parameter(ivarCFH1 = ivarInt+2)
        parameter(ivarCFG1 = ivarInt+3)
        parameter(ivarCFG2 = ivarInt+4)
        parameter(ivarDelN = ivarInt+5)
        parameter(ivarDim  = ivarInt+5)
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
        common/RT_BASE_D/ varMin(ivarDim), varLow(ivarDim)
	common/RT_BASE_E/ uDen, uLen, uTime, aExp, Hubble, daExp, qExp3
C
	logical ioNode
	common/RT_BASE_F/ ioNode
C
	include 'frt_f77.inc'


