C
C Working set of parameters
C
	parameter (NRLW = 20001)
	parameter (fnuLW1 = 11.3)
	parameter (fnuLW2 = 13.6)
C
C  Test - make LW insert coincide with the main mesh
C
*	parameter (NRLW = 4)
*	parameter (fnuLW1 = 11.032618)
*	parameter (fnuLW2 = 1.22140276*fnuLW1)
C
	parameter (LRLWTS = 51)
C
	common/RT_CHM_A/ angLW(NRLW), xigLW(NRLW), XILWST
C
	common/RT_CHM_B/ guvLW0(NRLW), guvLWe(NRLW)
C
	common/RT_CHM_C/ csLWMH(NRLW), csLWHI(NRLW)
