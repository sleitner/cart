C
C  Secondary electrons, no contribution to HeII
C
	parameter (SE1H1A = 0.3908)
	parameter (SE1H1B = 0.4092)
	parameter (SE1H1C = 1.7592)
C
	parameter (SE1G1A = 0.0554)
	parameter (SE1G1B = 0.4614)
	parameter (SE1G1C = 1.6660)
C
	parameter (SE1heA = 1.0000)
	parameter (SE1heB = 0.2663)
	parameter (SE1heC = 1.3163)
C
	parameter (SE2H1A = -0.110)
	parameter (SE2H1B = 0.2000)
	parameter (SE2H1C = 0.3800)
C
	parameter (SE2G1A = -0.0156)
	parameter (SE2G1B = 0.2000)
	parameter (SE2G1C = 0.3800)
C
	parameter (SE2heA = -0.100)
	parameter (SE2heB = 0.4000)
	parameter (SE2heC = 0.3400)
C
C  Radiation field tables: SE (order DOES matter here)
C  Values starting with i- are external, with j- are internal
C
        parameter (istabPiH1G2   =  1)
        parameter (istabPiG1G2   =  2)
        parameter (istabPiH1G1   =  3)
        parameter (istabPiG1G1   =  4)
        parameter (istabPhG1     =  5)
        parameter (istabPiH1H1   =  6)
        parameter (istabPiG1H1   =  7)
        parameter (istabPhH1     =  8)
        parameter (istabDim      =  8)
        parameter (jstabPiH1G2w1 =  1)
        parameter (jstabPiH1G2w2 =  2)
        parameter (jstabPiG1G2w1 =  3)
        parameter (jstabPiG1G2w2 =  4)
        parameter (jstabPiH1G1w1 =  5)
        parameter (jstabPiH1G1w2 =  6)
        parameter (jstabPiG1G1w1 =  7)
        parameter (jstabPiG1G1w2 =  8)
        parameter (jstabPhG1w1   =  9)
        parameter (jstabPhG1w2   = 10)
        parameter (jstabPiH1H1w1 = 11)
        parameter (jstabPiH1H1w2 = 12)
        parameter (jstabPiG1H1w1 = 13)
        parameter (jstabPiG1H1w2 = 14)
        parameter (jstabPhH1w1   = 15)
        parameter (jstabPhH1w2   = 16)
        parameter (jstabH1       = 11)
        parameter (jstabG1       =  5)
        parameter (jstabG2       =  1)
        parameter (jstabDim      = 16)
        parameter (jstabDim2     =  8)
C
        common/RT_XRS_A/ se1H1i(NRAD), se1H1h(NRAD),
     .			 se1G1i(NRAD), se1G1h(NRAD)
C
	common/RT_XRS_B/ se2H1i(NRAD), se2H1h(NRAD),
     .			 se2G1i(NRAD), se2G1h(NRAD)
C
C  Global rates tables
C
	common/RT_XRS_C/ pTabSE(jstabDim,NOPT3)
C
C  Local rates tables
C
	common/RT_XRS_D/ qTabSE(jstabDim,NOPT3)
C
C  Integration weights and normalization
C
        common/RT_XRS_E/ wpTabSE(jstabDim,NRAD)

