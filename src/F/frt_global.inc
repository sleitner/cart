C
C This file contains some of the global variables which need to be
C accumulated during the loop over all baryonic elements.
C
	parameter(iglobBDim = NRAD) 
	parameter(iglobDDim = ioMax1) 
	parameter(iglobIDim = infoDim) 
	dimension gBrem(iglobBDim), gData(iglobDDim), gInfo(iglobIDim)
	common/RT_GLB_A/ gBrem, gData, gInfo
C

